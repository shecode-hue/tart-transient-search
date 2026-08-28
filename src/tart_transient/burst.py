"""Solar radio burst detection, direct from HDF visibilities."""
from __future__ import annotations

import logging
from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import List

import numpy as np

log = logging.getLogger(__name__)

C_LIGHT = 299792458.0

CONTROL_AZEL = [(30., 40.), (120., 55.), (210., 45.), (300., 60.),
                (75., 70.), (165., 35.), (255., 65.), (345., 50.)]


def fetch(site: str, out_dir: Path, start_iso: str, minutes: float) -> List[Path]:
    import tart_tools.archive_handler as ah
    from minio import Minio

    t0 = datetime.fromisoformat(start_iso)
    if t0.tzinfo is None:
        t0 = t0.replace(tzinfo=timezone.utc)
    t1 = t0 + timedelta(minutes=float(minutes))
    out_dir.mkdir(parents=True, exist_ok=True)

    client = Minio(endpoint=ah.MINIO_API_HOST, secure=True)
    picked, day = [], t0
    while day.date() <= t1.date():
        prefix = "%s/vis/%d/%d/%d/" % (site, day.year, day.month, day.day)
        for obj in client.list_objects(ah.BUCKET_NAME, prefix=prefix, recursive=True):
            name = obj.object_name.rsplit("/", 1)[-1]
            try:
                when = datetime.strptime(name.replace(".hdf", "").split("_", 1)[1],
                                         "%Y-%m-%d_%H_%M_%S.%f")
            except (ValueError, IndexError):
                continue
            if t0 <= when.replace(tzinfo=timezone.utc) <= t1:
                picked.append(obj.object_name)
        day += timedelta(days=1)

    out = []
    for obj_name in sorted(picked):
        dest = out_dir / obj_name.rsplit("/", 1)[-1]
        if not dest.exists():
            client.fget_object(ah.BUCKET_NAME, obj_name, str(dest))
        out.append(dest)
    log.info("%s: %d file(s)", site, len(out))
    return out


def _read(path: Path):
    import json
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


def light_curve(files: List[Path], peak_iso: str, controls=CONTROL_AZEL) -> dict:
    from astropy.coordinates import AltAz, EarthLocation, get_body
    from astropy.time import Time
    import astropy.units as u

    peak = Time(datetime.fromisoformat(peak_iso)).unix
    rows, ctrl = [], {i: [] for i in range(len(controls))}
    elev = []

    for path in sorted(files):
        cfg, pos, bl, vis, g, ph, ts = _read(path)
        lam = C_LIGHT / cfg["frequency"]
        loc = EarthLocation(lat=cfg["lat"] * u.deg, lon=cfg["lon"] * u.deg,
                            height=cfg["alt"] * u.m)
        gc = g * np.exp(1j * ph)
        ok = (g[bl[:, 0]] > 0) & (g[bl[:, 1]] > 0)
        b = (pos[bl[:, 0]] - pos[bl[:, 1]])[ok] / lam
        with np.errstate(divide="ignore", invalid="ignore"):
            V = (vis / (gc[bl[:, 0]] * np.conj(gc[bl[:, 1]]))[None, :])[:, ok]

        for i, iso in enumerate(ts):
            when = Time(datetime.fromisoformat(iso.replace("Z", "+00:00")))
            aa = get_body("sun", when).transform_to(
                AltAz(obstime=when, location=loc))
            dirs = [(aa.az.deg, aa.alt.deg)] + list(controls)
            amps = []
            for azd, eld in dirs:
                el, az = np.radians(eld), np.radians(azd)
                s = np.array([np.cos(el) * np.sin(az), np.cos(el) * np.cos(az),
                              np.sin(el)])
                M = np.exp(2j * np.pi * (b @ s))
                amps.append(abs(np.vdot(M, V[i]) / np.vdot(M, M)))
            rows.append(((when.unix - peak) / 60.0, amps[0]))
            elev.append(aa.alt.deg)
            for k in range(len(controls)):
                ctrl[k].append(amps[k + 1])

    rows.sort(key=lambda r: r[0])
    t = np.array([r[0] for r in rows])
    a = np.array([r[1] for r in rows])
    order = np.argsort([r[0] for r in rows])
    return {"t_min": t, "amp": a, "elev_deg": float(np.median(elev)),
            "control": {k: np.array(v)[order] for k, v in ctrl.items()},
            "site": _read(sorted(files)[0])[0]["name"]}


def summarise(lc: dict, smooth_s: int = 30) -> dict:
    def sm(y):
        n = max(1, int(smooth_s))
        return np.convolve(y, np.ones(n) / n, mode="same")

    t, y = lc["t_min"], sm(lc["amp"])
    pre = t < -4
    base = float(np.median(y[pre])) if pre.any() else float(np.median(y))
    peak = float(np.nanmax(y))
    ratios = []
    for k, c in lc["control"].items():
        cs = sm(c)
        cb = float(np.median(cs[pre])) if pre.any() else float(np.median(cs))
        ratios.append(float(np.nanmax(cs) / max(cb, 1e-12)))
    ratios = np.array(ratios)
    ratio = peak / max(base, 1e-12)
    sigma = ((ratio - ratios.mean()) / ratios.std()) if ratios.std() > 0 else float("nan")
    return {"baseline": base, "peak": peak, "ratio": ratio,
            "peak_t_min": float(t[int(np.nanargmax(y))]),
            "control_ratio_mean": float(ratios.mean()),
            "control_ratio_max": float(ratios.max()),
            "sigma_over_controls": float(sigma),
            "elev_deg": lc["elev_deg"], "n_integrations": int(len(t))}


def figure(curves: dict, peak_iso: str, out_path: Path) -> Path:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    sites = list(curves)
    fig, ax = plt.subplots(2, 1, figsize=(12, 8.5), sharex=True)
    cols = ["#1f77b4", "#2ca02c", "#d62728", "#9467bd", "#ff7f0e"]

    for i, s in enumerate(sites):
        lc, st = curves[s]["lc"], curves[s]["summary"]
        y = np.convolve(lc["amp"], np.ones(30) / 30, mode="same")
        ax[0].plot(lc["t_min"], y / st["baseline"], lw=1.7,
                   color=cols[i % len(cols)],
                   label="%s  (x%.0f, sun %.0f deg)" % (s, st["ratio"], st["elev_deg"]))
    ax[0].axvline(0, color="k", ls="--", lw=1.4)
    ax[0].set_ylabel("brightness at the Sun\n(relative to pre-burst)")
    ax[0].set_title("Solar radio burst, peak %s UTC" % peak_iso[:19])
    ax[0].legend(fontsize=9)
    ax[0].grid(alpha=.25)

    s0 = sites[0]
    lc = curves[s0]["lc"]
    base = curves[s0]["summary"]["baseline"]
    ax[1].plot(lc["t_min"], np.convolve(lc["amp"], np.ones(30) / 30, mode="same") / base,
               lw=2.2, color=cols[0], label="the Sun", zorder=5)
    for k, c in lc["control"].items():
        cs = np.convolve(c, np.ones(30) / 30, mode="same")
        cb = np.median(cs[lc["t_min"] < -4]) or 1.0
        ax[1].plot(lc["t_min"], cs / cb, lw=0.9, color="#999",
                   label="control directions" if k == 0 else None)
    ax[1].axvline(0, color="k", ls="--", lw=1.4)
    ax[1].set_xlabel("minutes from the reported burst peak")
    ax[1].set_ylabel("brightness\n(relative to pre-burst)")
    ax[1].set_title("%s: the rise appears only toward the Sun" % s0)
    ax[1].legend(fontsize=9)
    ax[1].grid(alpha=.25)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(out_path, dpi=130)
    plt.close(fig)
    return out_path
