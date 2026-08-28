"""Long-baseline time series: every catalogued source, every interval."""
from __future__ import annotations

import json
import logging
from datetime import datetime
from pathlib import Path
from typing import List

import numpy as np

log = logging.getLogger(__name__)

C_LIGHT = 299792458.0


def _read(path):
    import h5py
    with h5py.File(path, "r") as f:
        cfg = json.loads(np.asarray(f["config"][()]).ravel()[0].decode())
        pos = np.asarray(f["antenna_positions"][()], float)
        bl = np.asarray(f["baselines"][()], int)
        vis = np.asarray(f["vis"][()])
        g = np.asarray(f["gains"][()], float)
        ph = np.asarray(f["phases"][()], float)
        ts = [t.decode() if isinstance(t, bytes) else str(t)
              for t in np.asarray(f["timestamp"][()]).ravel()]
    return cfg, pos, bl, vis, g, ph, ts


def load_span(files: List[Path]):
    from astropy.time import Time

    stamps, rows, cfg0, b0, corr, ok = [], [], None, None, None, None
    for path in sorted(files):
        cfg, pos, bl, vis, g, ph, ts = _read(path)
        if cfg0 is None:
            cfg0 = cfg
            gc = g * np.exp(1j * ph)
            ok = (g[bl[:, 0]] > 0) & (g[bl[:, 1]] > 0)
            b0 = (pos[bl[:, 0]] - pos[bl[:, 1]])[ok] / (C_LIGHT / cfg["frequency"])
            corr = (gc[bl[:, 0]] * np.conj(gc[bl[:, 1]]))[ok]
        with np.errstate(divide="ignore", invalid="ignore"):
            V = vis[:, ok] / corr[None, :]
        for i, iso in enumerate(ts):
            stamps.append(Time(datetime.fromisoformat(
                iso.replace("Z", "+00:00"))).unix)
            rows.append(V[i])
    order = np.argsort(stamps)
    return cfg0, b0, np.asarray(stamps)[order], np.asarray(rows)[order]


def catalogue_at(site, lat, lon, when_iso, min_elevation_deg=10.0):
    from tart_tools import api_handler

    api = api_handler.APIhandler("https://api.elec.ac.nz/tart/" + site)
    url = api.catalog_url(lon=lon, lat=lat, datestr=when_iso)
    url += "&elevation=" + str(min_elevation_deg)
    out = []
    for e in api.get_url(url):
        el, az = float(e.get("el", -99)), float(e.get("az", 0))
        if el < min_elevation_deg:
            continue
        er, ar = np.radians(el), np.radians(az)
        out.append({"name": str(e.get("name", "?")), "el": el, "az": az,
                    "s": np.array([np.cos(er) * np.sin(ar),
                                   np.cos(er) * np.cos(ar),
                                   np.sin(er)])})
    return out


def run(files, site, interval_s=30.0, min_elevation_deg=10.0,
        catalogue_every=10):
    from astropy.time import Time

    cfg, b, stamps, V = load_span(files)
    t0, t1 = float(stamps.min()), float(stamps.max())
    n_int = int(np.ceil((t1 - t0) / interval_s))
    log.info("%d integrations over %.1f min -> %d intervals of %.0f s",
             len(stamps), (t1 - t0) / 60.0, n_int, interval_s)

    cat_cache, series, mids = {}, {}, []
    for k in range(n_int):
        lo, hi = t0 + k * interval_s, t0 + (k + 1) * interval_s
        m = (stamps >= lo) & (stamps < hi)
        if m.sum() < 4:
            continue
        mid = float(np.mean(stamps[m]))
        mids.append(mid)
        key = k // max(catalogue_every, 1)
        if key not in cat_cache:
            iso = Time(mid, format="unix").utc.isot + "+00:00"
            try:
                cat_cache[key] = catalogue_at(site, cfg["lat"], cfg["lon"],
                                              iso, min_elevation_deg)
            except Exception as exc:
                log.warning("catalogue fetch failed, interval %d: %s", k, exc)
                cat_cache[key] = cat_cache.get(key - 1, [])
        Vk = V[m].reshape(-1)
        n_row = int(m.sum())
        for src in cat_cache[key]:
            M = np.tile(np.exp(2j * np.pi * (b @ src["s"])), (n_row, 1)).reshape(-1)
            den = np.vdot(M, M)
            if abs(den) == 0:
                continue
            amp = np.vdot(M, Vk) / den
            resid = Vk - amp * M
            sig = np.sqrt((np.abs(resid) ** 2).sum() /
                          max(resid.size - 1, 1) / resid.size)
            rec = series.setdefault(src["name"],
                                    {"t": [], "amp": [], "snr": [], "el": []})
            rec["t"].append(mid)
            rec["amp"].append(float(abs(amp)))
            rec["snr"].append(float(abs(amp) / sig) if sig > 0 else 0.0)
            rec["el"].append(src["el"])
        if (k + 1) % 20 == 0:
            log.info("  %d/%d intervals", k + 1, n_int)

    return {"config": cfg, "t_start": t0, "t_end": t1,
            "interval_s": interval_s, "n_intervals": len(mids),
            "series": series}


def summarise(res, min_points=6):
    out = []
    for name, rec in res["series"].items():
        a = np.asarray(rec["amp"], float)
        if len(a) < min_points or not np.isfinite(a).any():
            continue
        med = float(np.nanmedian(a))
        mad = float(np.nanmedian(np.abs(a - med))) * 1.4826
        pk = int(np.nanargmax(a))
        out.append({"name": name, "n": int(len(a)), "median": med, "mad": mad,
                    "peak": float(a[pk]), "peak_t": float(rec["t"][pk]),
                    "ratio": float(a[pk] / med) if med > 0 else float("nan"),
                    "excess": float(a[pk] - med),
                    "sigma": float((a[pk] - med) / mad) if mad > 0 else float("nan"),
                    "median_el": float(np.nanmedian(rec["el"]))})
    out.sort(key=lambda r: -(r["sigma"] if np.isfinite(r["sigma"]) else -1))
    return out


def figure(res, ranked, out_path, n_show=6):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from astropy.time import Time

    t0 = res["t_start"]
    fig, ax = plt.subplots(2, 1, figsize=(14, 9))
    for r in ranked[:n_show]:
        rec = res["series"][r["name"]]
        t = (np.asarray(rec["t"]) - t0) / 60.0
        ax[0].plot(t, rec["amp"], "-", lw=1.6,
                   label="%s  (%.0f sigma)" % (r["name"][:26], r["sigma"]))
    ax[0].set_ylabel("|fitted amplitude|")
    ax[0].set_title("Most deviant sources - %s from %s UTC, %.0f s intervals"
                    % (res["config"].get("name", "?"),
                       Time(t0, format="unix").iso[:19], res["interval_s"]),
                    fontsize=12)
    ax[0].legend(fontsize=8, ncol=2)
    ax[0].grid(alpha=.25)

    allc = []
    for name, rec in res["series"].items():
        a = np.asarray(rec["amp"], float)
        med = np.nanmedian(a)
        if med > 0 and len(a) > 5:
            t = (np.asarray(rec["t"]) - t0) / 60.0
            ax[1].plot(t, a / med, "-", lw=0.6, color="#bbb", zorder=1)
            allc.append((t, a / med))
    if allc:
        grid = np.linspace(0, max(c[0].max() for c in allc), 200)
        stack = np.vstack([np.interp(grid, t, y) for t, y in allc])
        ax[1].plot(grid, np.median(stack, axis=0), lw=2.4, color="#d62728",
                   label="median of all sources (common mode)")
        ax[1].legend(fontsize=9)
    ax[1].set_xlabel("minutes from start")
    ax[1].set_ylabel("relative to own median")
    ax[1].set_title("Every source normalised - a common rise means one bright "
                    "source is leaking into all of them", fontsize=11)
    ax[1].grid(alpha=.25)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(out_path, dpi=120)
    plt.close(fig)
    return out_path


def onset_offset(res, name, k_sigma=3.0):
    rec = res["series"].get(name)
    if not rec:
        return None
    t = np.asarray(rec["t"], float)
    a = np.asarray(rec["amp"], float)
    med = float(np.nanmedian(a))
    mad = float(np.nanmedian(np.abs(a - med))) * 1.4826
    if not (mad > 0):
        return None
    lvl = med + k_sigma * mad
    above = a > lvl
    if not above.any():
        return None
    pk = int(np.nanargmax(a))
    i = pk
    while i > 0 and above[i - 1]:
        i -= 1
    j = pk
    while j < len(a) - 1 and above[j + 1]:
        j += 1
    return {"name": name, "median": med, "mad": mad, "threshold": lvl,
            "onset_t": float(t[i]), "peak_t": float(t[pk]),
            "offset_t": float(t[j]), "peak": float(a[pk]),
            "duration_s": float(t[j] - t[i]),
            "rise_s": float(t[pk] - t[i]), "decay_s": float(t[j] - t[pk]),
            "n_above": int(above.sum()), "k_sigma": k_sigma}


def event_figure(res, ev, out_path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from astropy.time import Time

    rec = res["series"][ev["name"]]
    t = np.asarray(rec["t"], float)
    a = np.asarray(rec["amp"], float)
    s = np.asarray(rec["snr"], float)
    t0 = res["t_start"]
    x = (t - t0) / 60.0

    fig, ax = plt.subplots(2, 1, figsize=(13, 8), sharex=True)
    for k, (y, lab) in enumerate(((a, "|amplitude|"), (s, "SNR"))):
        ax[k].plot(x, y, "-", lw=1.5, color="#1f77b4")
        ax[k].axvspan((ev["onset_t"] - t0) / 60.0, (ev["offset_t"] - t0) / 60.0,
                      color="#ffd966", alpha=.35, zorder=0)
        ax[k].set_ylabel(lab)
        ax[k].grid(alpha=.25)
    ax[0].axhline(ev["median"], color="#888", ls=":", lw=1)
    ax[0].axhline(ev["threshold"], color="#a32020", ls="--", lw=1.4)
    ax[1].set_xlabel("minutes from %s UTC" % Time(t0, format="unix").iso[:19])
    ax[0].set_title("%s   onset %s   peak %s   end %s   (%.0f s)"
                    % (ev["name"],
                       Time(ev["onset_t"], format="unix").iso[11:19],
                       Time(ev["peak_t"], format="unix").iso[11:19],
                       Time(ev["offset_t"], format="unix").iso[11:19],
                       ev["duration_s"]), fontsize=12)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout(); fig.savefig(out_path, dpi=120); plt.close(fig)
    return out_path
