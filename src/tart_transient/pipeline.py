"""End-to-end orchestration."""
from __future__ import annotations

import csv
import logging
import shutil
import time
from pathlib import Path

import numpy as np

from . import calibration, catalogue as catmod, download, fitting, imaging
from . import measurement_set as msmod, report, search
from .config import Config
from .significance import empirical_null

log = logging.getLogger(__name__)


def stage_download(cfg: Config):
    cfg.make_dirs()
    d = cfg.section("download")
    files = download.archive(d["target"], cfg.data_dir, str(d.get("start", "-60")),
                             float(d.get("duration_min", 2)), int(d.get("n_files", 1)))
    return files


def stage_build(cfg: Config):
    cfg.make_dirs()
    files = sorted(cfg.data_dir.glob("*.hdf"))
    if not files:
        raise RuntimeError(f"no .hdf in {cfg.data_dir}; run `download` first")
    info = download.describe_all(files)
    log.info("input: %d integrations x %d baselines, gains stored=%s",
             info["n_integrations"], info["n_baselines"], info["gains_stored"])

    ms_cfg = cfg.section("measurement_set")
    have_ms = cfg.ms_path.exists() and any(cfg.ms_dir.glob("model_sources_*.txt"))
    if ms_cfg.get("reuse_existing", False) and have_ms:
        log.info("reusing existing measurement set at %s (reuse_existing: true)",
                 cfg.ms_path)
    else:
        msmod.build(files, cfg.ms_path, cfg.ms_dir / "model_sources_",
                    float(cfg.section("catalogue").get("filter_elevation_deg", 5.0)),
                    ms_cfg.get("timeout_s"))

    cal = cfg.section("calibration")
    column = calibration.decide(info, cal.get("mode", "auto"))
    if column == "CORRECTED_DATA":
        msmod.calibrate(cfg.ms_path, cfg.ms_dir / "caltables",
                        float(cal.get("minsnr", 2.0)), str(cal.get("solint", "inf")),
                        str(cal.get("applymode", "calonly")))
    return info, column


def stage_search(cfg: Config):
    from casacore.tables import table
    from astropy.io import fits as pyfits

    t0 = time.time()
    files = sorted(cfg.data_dir.glob("*.hdf"))
    info = download.describe_all(files)
    column = calibration.decide(info, cfg.section("calibration").get("mode", "auto"))

    work = cfg.ms_dir / "work.ms"
    if work.exists():
        shutil.rmtree(work)
    shutil.copytree(cfg.ms_path, work)

    geom = msmod.observation_geometry(work)
    ra0, dec0, freqs = geom["ra0_deg"], geom["dec0_deg"], geom["freqs_hz"]
    beam = geom["beam_deg"]
    log.info("beam %.2f deg from %.1f lambda max baseline", beam, geom["max_baseline_lambda"])

    with table(str(work), ack=False) as t:
        times, uvw = t.getcol("TIME"), t.getcol("UVW")
        data = t.getcol(column)

    ci = catmod.CatalogueIndex.build(cfg.ms_dir, work)
    raw = ci.at(float(times.mean()), cfg.section("catalogue").get("min_flux", 100.0))
    cat, blends = catmod.merge_within_beam(raw, beam)
    names = list(cat["name"].values)
    tracks = ci.positions_for(times, list(cat["track_name"].values))
    if not np.isfinite(tracks[:, :, 0]).any(axis=1).all():
        missing = [names[i] for i in np.where(~np.isfinite(tracks[:, :, 0]).any(axis=1))[0]]
        raise RuntimeError(f"no catalogue positions for: {missing}")

    if cfg.section("catalogue").get("subtract_unmodelled", True):
        try:
            from astropy.coordinates import EarthLocation
            from astropy.time import Time
            import astropy.units as u
            with table(str(work) + "/ANTENNA", ack=False) as t:
                xyz = t.getcol("POSITION").mean(axis=0)
            here = EarthLocation.from_geocentric(*(xyz * u.m))
            site = cfg.section("download").get("target", cfg.name)
            ut = np.unique(times)
            n_ep = int(min(6, max(2, len(ut))))
            samp = ut[np.linspace(0, len(ut) - 1, n_ep).astype(int)]
            epochs = [Time(t_ / 86400.0, format="mjd").utc.isot + "+00:00"
                      for t_ in samp]
            extra_names, extra_ep = catmod.unmodelled_tracks(
                here.lat.deg, here.lon.deg, site, epochs, names,
                cfg.section("catalogue").get("filter_elevation_deg", 5.0))
            if len(extra_names):
                xs = Time(samp / 86400.0, format="mjd").unix
                xt = Time(times / 86400.0, format="mjd").unix
                extra = np.empty((len(extra_names), len(times), 2))
                for k in range(len(extra_names)):
                    ra_u = np.unwrap(np.radians(extra_ep[k, :, 0]))
                    extra[k, :, 0] = np.degrees(np.interp(xt, xs, ra_u)) % 360.0
                    extra[k, :, 1] = np.interp(xt, xs, extra_ep[k, :, 1])
                tracks = np.concatenate([tracks, extra], axis=0)
                names = names + [f"UNMODELLED:{n}" for n in extra_names]
                log.info("subtracting %d objects the name filter dropped: %s",
                         len(extra_names), ", ".join(extra_names[:6]))
        except Exception as exc:
            log.warning("could not add unmodelled objects (%s); "
                        "they will be vetoed but not subtracted", exc)

    if cfg.section("catalogue").get("fix_solar_system", True):
        try:
            from astropy.coordinates import EarthLocation
            import astropy.units as u
            with table(str(work) + "/ANTENNA", ack=False) as _t:
                _xyz = _t.getcol("POSITION").mean(axis=0)
            _here = EarthLocation.from_geocentric(*(_xyz * u.m))
            catmod.fix_solar_system_positions(names, tracks, times,
                                              _here.lat.deg, _here.lon.deg)
        except Exception as exc:
            log.warning("solar-system position correction failed: %s", exc)

    model = fitting.model_matrix(uvw, freqs, tracks, ra0, dec0)
    cond = fitting.gram_condition(model)
    if cond > 1e6:
        log.warning("Gram condition %.2e -- model columns are not distinguishable; "
                    "the fit will split flux arbitrarily", cond)

    img_cfg = {k.replace("_", "-"): v for k, v in cfg.section("imaging").items()}
    r = imaging.image(work, cfg.fits_dir / "01_before.fits", column, img_cfg)
    if r.returncode != 0 or not r.fits_path.exists():
        raise RuntimeError("DiSkO failed imaging " + column + ":\n"
                           + r.failure_message())
    before, hdr = imaging.load_image(r.fits_path)
    wcs = imaging.load_wcs(hdr)
    npix = before.shape[0]; pix = abs(wcs.wcs.cdelt[0]); half = npix / 2 * pix
    _s, _m, rms_b = imaging.find_peaks(before, beam / pix)
    n_before = 0 if _s is None else len(_s)

    amp, sigma = fitting.fit_amplitudes(data, model)
    snr = fitting.snr_of(amp, sigma)
    sig_cfg = cfg.section("significance")
    thr, samples = empirical_null(uvw, freqs, data, ra0, dec0,
                                  n_trials=sig_cfg.get("null_trials", 400),
                                  percentile=sig_cfg.get("null_percentile", 99.0),
                                  base_model=model)
    log.info("null median %.2f threshold %.2f -- %d/%d sources above",
             float(np.median(samples)), thr, int((snr >= thr).sum()), len(names))

    keep = [str(k).strip().lower()
            for k in cfg.section("catalogue").get("keep_sources", [])]
    amp_sub = amp.copy()
    if keep:
        kept = [i for i, n in enumerate(names)
                if any(k in str(n).lower() for k in keep)]
        if kept:
            amp_sub[kept] = 0.0
            log.info("keeping in the residual (not subtracted): %s",
                     ", ".join(str(names[i]) for i in kept))
    residual = fitting.peel(data, model, amp_sub)
    power = 1.0 - float(np.mean(np.abs(residual) ** 2) / np.mean(np.abs(data) ** 2))
    with table(str(work), readonly=False, ack=False) as t:
        if "RESIDUAL_DATA" not in t.colnames():
            desc = t.getcoldesc(column); desc["name"] = "RESIDUAL_DATA"
            t.addcols(desc)
        t.putcol("RESIDUAL_DATA", residual)

    r2 = imaging.image(work, cfg.fits_dir / "02_after.fits", "RESIDUAL_DATA", img_cfg)
    if r2.returncode != 0 or not r2.fits_path.exists():
        raise RuntimeError("DiSkO failed imaging RESIDUAL_DATA:\n"
                           + r2.failure_message())
    after, _h = imaging.load_image(r2.fits_path)
    removed = before - after
    pyfits.writeto(str(cfg.fits_dir / "03_removed.fits"), removed.astype(np.float32),
                   header=pyfits.getheader(str(r.fits_path)), overwrite=True)
    _s2, _m2, rms_a = imaging.find_peaks(after, beam / pix)
    n_after = 0 if _s2 is None else len(_s2)

    veto = None
    try:
        from astropy.coordinates import EarthLocation
        import astropy.units as u
        with table(str(work) + "/ANTENNA", ack=False) as t:
            xyz = t.getcol("POSITION").mean(axis=0)
        loc = EarthLocation.from_geocentric(*(xyz * u.m))
        site = cfg.section("download").get("target", cfg.name)
        veto = catmod.full_sky(loc.lat.deg, loc.lon.deg,
                               info["t_start"], site)
        log.info("veto catalogue: %d objects above the horizon (%d modelled)",
                 len(veto), len(cat))
    except Exception as exc:
        log.warning("full-sky veto catalogue unavailable (%s); "
                    "falling back to the modelled subset", exc)

    result = search.run(after, wcs, uvw, freqs, residual, ra0, dec0, cat,
                        beam, beam / pix, cfg, veto_catalogue=veto,
                        times=times)

    import h5py
    with h5py.File(files[0], "r") as f:
        vis_abs = np.abs(np.asarray(f["vis"][()])).ravel()
        gains = np.asarray(f["gains"][()]); phases = np.asarray(f["phases"][()])
    report.input_data(info, vis_abs, gains, phases, cfg.plots_dir / "01_input_data.png")

    lam = 299792458.0 / float(freqs[0])
    l, m, _n = fitting.lmn(raw["ra_d"].values, raw["dec_d"].values, ra0, dec0)
    A = model.reshape(-1, model.shape[2])
    report.sky_model(uvw / lam, np.column_stack([l, m]),
                     np.linalg.eigvalsh(A.conj().T @ A)[::-1], beam,
                     len(raw), len(cat), cond, cfg.plots_dir / "02_sky_model.png")
    report.snr_vs_null(names, snr, thr, float(np.median(samples)),
                       cfg.plots_dir / "03_snr_vs_null.png")
    report.before_after(before, after, removed, half, len(names), power,
                        cfg.plots_dir / "04_before_after_removed.png")
    report.effectiveness(before, after, power,
                         100 * (rms_a - rms_b) / rms_b,
                         100 * (n_after - n_before) / max(n_before, 1),
                         cfg.plots_dir / "05_clean_effectiveness.png")
    report.transient_search(after, half, wcs, npix, pix, result,
                            cfg.plots_dir / "06_transient_search.png")

    with open(cfg.results_dir / "sources.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["name", "ra_deg", "dec_deg", "amp", "sigma", "snr", "above_threshold"])
        mid = tracks.shape[1] // 2
        for k, nm in enumerate(names):
            w.writerow([nm, float(tracks[k, mid, 0]), float(tracks[k, mid, 1]),
                        float(np.abs(amp[k])), float(sigma[k]), float(snr[k]),
                        bool(snr[k] >= thr)])
    with open(cfg.results_dir / "candidates.csv", "w", newline="") as f:
        rows = result["candidates"]
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()) if rows else ["none"])
        w.writeheader(); w.writerows(rows)

    record = dict(
        name=cfg.name, column_used=column, elapsed_s=round(time.time() - t0, 1),
        input=info, geometry={k: v for k, v in geom.items() if k != "freqs_hz"},
        model=dict(catalogue_raw=len(raw), catalogue_merged=len(cat),
                   gram_condition=cond, blends={k: v for k, v in blends.items() if len(v) > 1}),
        fit=dict(null_median=float(np.median(samples)), threshold=float(thr),
                 n_above_threshold=int((snr >= thr).sum()),
                 snr_min=float(snr.min()), snr_max=float(snr.max()),
                 snr_median=float(np.median(snr)), power_removed=float(power)),
        images=dict(npix=npix, cdelt_deg=float(pix), rms_before=float(rms_b),
                    rms_after=float(rms_a), rms_removed=float(np.nanstd(removed)),
                    peaks_before=n_before, peaks_after=n_after),
        search=result)
    result.pop("window_map", None)
    report.write_summary(record, cfg.results_dir / "summary.json")
    return record


def stage_compare(quiet_dir, burst_dir, out_dir, target_ra=None, target_dec=None):
    import json
    from pathlib import Path as _P

    import numpy as np
    from astropy.io import fits as pyfits

    quiet_dir, burst_dir, out_dir = _P(quiet_dir), _P(burst_dir), _P(out_dir)
    (out_dir / "fits").mkdir(parents=True, exist_ok=True)
    (out_dir / "plots").mkdir(parents=True, exist_ok=True)
    (out_dir / "results").mkdir(parents=True, exist_ok=True)

    def _load(d):
        img, hdr = pyfits.getdata(d / "fits" / "01_before.fits", header=True)
        rec = json.load(open(d / "results" / "summary.json"))
        return np.squeeze(img).astype(float), hdr, rec

    iq, hq, rq = _load(quiet_dir)
    ib, hb, rb = _load(burst_dir)
    if iq.shape != ib.shape:
        raise RuntimeError("epoch images differ in size: %s vs %s"
                           % (iq.shape, ib.shape))

    diff = ib - iq
    pyfits.writeto(str(out_dir / "fits" / "04_epoch_difference.fits"),
                   diff.astype(np.float32), header=hb, overwrite=True)

    wcs = imaging.load_wcs(hb)
    beam = float(rb["geometry"]["beam_deg"])
    pix = abs(wcs.wcs.cdelt[0])

    srcs, _med, rms = imaging.find_peaks(diff, beam / pix, 5.0)
    changed = []
    if srcs is not None:
        for row in srcs:
            x, y = float(row["x_centroid"]), float(row["y_centroid"])
            ra, dec = wcs.wcs_pix2world(x, y, 0)
            if not (np.isfinite(ra) and np.isfinite(dec)):
                continue
            i, j = int(round(y)), int(round(x))
            changed.append({
                "ra_deg": float(ra), "dec_deg": float(dec),
                "quiet": float(np.nanmax(iq[max(0, i-15):i+16, max(0, j-15):j+16])),
                "burst": float(np.nanmax(ib[max(0, i-15):i+16, max(0, j-15):j+16])),
                "change": float(np.nanmax(diff[max(0, i-15):i+16, max(0, j-15):j+16])),
                "change_rms": float(np.nanmax(
                    diff[max(0, i-15):i+16, max(0, j-15):j+16]) / max(rms, 1e-12)),
            })
    changed.sort(key=lambda c: -c["change"])
    for c in changed:
        c["ratio"] = c["burst"] / max(c["quiet"], 1e-12)

    record = {
        "quiet": {"run": quiet_dir.name, "t_start": rq["input"]["t_start"]},
        "burst": {"run": burst_dir.name, "t_start": rb["input"]["t_start"]},
        "beam_deg": beam, "difference_rms": float(rms),
        "n_changed": len(changed), "changed": changed[:20],
    }
    if target_ra is not None:
        x, y = [float(v) for v in wcs.wcs_world2pix(target_ra, target_dec, 0)]
        i, j = int(round(y)), int(round(x))
        box = lambda im: float(np.nanmax(im[max(0, i-30):i+31, max(0, j-30):j+31]))
        record["target"] = {
            "ra_deg": target_ra, "dec_deg": target_dec, "pixel": [x, y],
            "quiet": box(iq), "burst": box(ib), "change": box(diff),
            "ratio": box(ib) / max(box(iq), 1e-12),
            "change_rms": box(diff) / max(rms, 1e-12),
        }

    report.epoch_comparison(iq, ib, diff, record,
                            out_dir / "plots" / "07_epoch_comparison.png")
    with open(out_dir / "results" / "comparison.json", "w") as f:
        json.dump(record, f, indent=2, default=float)
    log.info("epoch difference: %d changed positions, rms %.5f",
             len(changed), rms)
    return record
