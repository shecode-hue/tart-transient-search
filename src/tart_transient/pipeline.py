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
                             int(d.get("duration_min", 2)), int(d.get("n_files", 1)))
    return files


def stage_build(cfg: Config):
    cfg.make_dirs()
    files = sorted(cfg.data_dir.glob("*.hdf"))
    if not files:
        raise RuntimeError(f"no .hdf in {cfg.data_dir}; run `download` first")
    info = download.describe(files[0])
    log.info("input: %d integrations x %d baselines, gains stored=%s",
             info["n_integrations"], info["n_baselines"], info["gains_stored"])

    msmod.build(files, cfg.ms_path, cfg.ms_dir / "model_sources_",
                float(cfg.section("catalogue").get("filter_elevation_deg", 5.0)))

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
    info = download.describe(files[0])
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

    residual = fitting.peel(data, model, amp)
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

    result = search.run(after, wcs, uvw, freqs, residual, ra0, dec0, cat,
                        beam, beam / pix, cfg)

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
        for k, nm in enumerate(names):
            w.writerow([nm, float(cat.iloc[k]["ra_d"]), float(cat.iloc[k]["dec_d"]),
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
    report.write_summary(record, cfg.results_dir / "summary.json")
    return record
