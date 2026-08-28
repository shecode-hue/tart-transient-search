"""Transient search on the residual visibilities."""
from __future__ import annotations

import logging
from typing import List

import numpy as np

from .fitting import C_LIGHT, PHASE_SIGN, lmn
from .significance import (ZenithNull, empirical_null, null_with_zenith,
                           trials_threshold, zenith_angle_deg)

log = logging.getLogger(__name__)


def snr_at(uvw, freqs_hz, data, ra0_deg, dec0_deg, ra_deg, dec_deg) -> float:
    wl = C_LIGHT / freqs_hz
    u = uvw[:, 0][:, None] / wl[None, :]
    v = uvw[:, 1][:, None] / wl[None, :]
    w = uvw[:, 2][:, None] / wl[None, :]
    l, m, n = lmn(ra_deg, dec_deg, ra0_deg, dec0_deg)
    M = np.exp(PHASE_SIGN * 2.0j * np.pi * (u * l + v * m + w * (n - 1.0)))
    V = data.reshape(data.shape[0], data.shape[1], -1)[:, :, 0].astype(np.complex128)
    denom = np.vdot(M, M)
    if abs(denom) == 0:
        return 0.0
    amp = np.vdot(M, V) / denom
    resid = V - amp * M
    sig = np.sqrt((np.abs(resid) ** 2).sum() / (resid.size - 1) / resid.size)
    return float(np.abs(amp) / sig) if sig > 0 else 0.0


def angular_sep(ra1, dec1, ra2, dec2):
    r1, d1, r2, d2 = map(np.radians, (ra1, dec1, ra2, dec2))
    return np.degrees(np.arccos(np.clip(
        np.sin(d1) * np.sin(d2) + np.cos(d1) * np.cos(d2) * np.cos(r1 - r2), -1, 1)))


def run(residual_img, wcs, uvw, freqs_hz, residual_vis, ra0_deg, dec0_deg,
        catalogue, beam_deg, fwhm_pix, cfg, veto_catalogue=None,
        times=None) -> dict:
    from .imaging import find_peaks

    sig_cfg = cfg.section("significance")
    srcs, _median, rms = find_peaks(residual_img, fwhm_pix,
                                    cfg.section("search").get("source_finder_sigma", 5.0))
    npix = residual_img.shape[0]

    positions = []
    if srcs is not None:
        for s in srcs:
            ra, dec = wcs.wcs_pix2world(float(s["x_centroid"]), float(s["y_centroid"]), 0)
            if np.isfinite(ra) and np.isfinite(dec):
                peak = float(s["peak"]) if "peak" in s.dtype.names else float("nan")
                positions.append((float(ra), float(dec), peak))
    n_image = len(positions)

    window_map = None
    wm_cfg = cfg.section("search").get("window_map", {}) or {}
    if wm_cfg.get("enabled", True) and times is not None:
        from .imaging import peaks_from_map, time_window_peak_map
        amap, wmap, gx, _gy = time_window_peak_map(
            uvw, freqs_hz, residual_vis, times, ra0_deg, dec0_deg,
            n_side=int(wm_cfg.get("n_side", 110)),
            n_windows=int(wm_cfg.get("n_windows", 6)))
        extra = peaks_from_map(amap, wmap, gx, ra0_deg, dec0_deg,
                               n_peaks=int(wm_cfg.get("n_peaks", 40)))
        window_map = (amap, gx)
        added = 0
        for e in extra:
            d = angular_sep(e["ra_deg"], e["dec_deg"],
                            np.array([p[0] for p in positions] or [1e9]),
                            np.array([p[1] for p in positions] or [1e9]))
            if float(np.min(d)) > beam_deg:
                positions.append((e["ra_deg"], e["dec_deg"], e["map_value"]))
                added += 1
        log.info("window-max map added %d position(s) the averaged image missed "
                 "(%d -> %d)", added, n_image, len(positions))

    log.info("residual peaks to test: %d", len(positions))
    if not positions:
        return dict(n_peaks=0, n_passed_single_look=0, n_unexplained=0,
                    n_confirmed=0, candidates=[], image_rms=float(rms))

    splits = list(cfg.section("search").get("time_windows", [1, 2, 4, 8]))
    windowed = bool(times is not None and len(splits) > 1)
    if windowed:
        log.info("time-windowed search, splits %s", splits)

    samples, zeniths = null_with_zenith(
        uvw, freqs_hz, residual_vis, ra0_deg, dec0_deg,
        n_trials=sig_cfg.get("trials_null_samples", 8000), seed=99,
        times=times if windowed else None,
        splits=splits if windowed else None)
    znull = ZenithNull(samples, zeniths)
    alpha = sig_cfg.get("false_alarm_rate", 0.01)
    thr_norm = znull.threshold(len(positions), alpha)
    thr_single_norm = float(np.percentile(
        znull.normalised, sig_cfg.get("null_percentile", 99.0)))
    log.info("zenith-matched null: %s", znull.describe())

    veto = catalogue if veto_catalogue is None else veto_catalogue
    cat_ra = veto["ra_d"].values.astype(float)
    cat_dec = veto["dec_d"].values.astype(float)
    cat_names = list(veto["name"].values) if "name" in veto else\
        [""] * len(cat_ra)

    cands = []
    refine = cfg.section("search").get("refine_positions", True)
    for ra, dec, peak in positions:
        if refine:
            ra, dec, s, n_split, i_win = refine_position(
                uvw, freqs_hz, residual_vis, times, ra0_deg, dec0_deg,
                ra, dec, splits=splits if windowed else None)
        elif windowed:
            s, n_split, i_win = snr_at_windows(
                uvw, freqs_hz, residual_vis, times, ra0_deg, dec0_deg,
                ra, dec, splits)
        else:
            s = snr_at(uvw, freqs_hz, residual_vis, ra0_deg, dec0_deg, ra, dec)
            n_split, i_win = 1, 0
        seps = angular_sep(ra, dec, cat_ra, cat_dec)
        j = int(np.argmin(seps)); sep = float(seps[j])
        z = float(zenith_angle_deg(ra, dec, ra0_deg, dec0_deg))
        s_norm = znull.snr_normalised(s, z)
        local_scale = float(znull.scale(z))
        cands.append(dict(
            ra_deg=ra, dec_deg=dec, image_peak=peak, vis_snr=s,
            best_window_split=int(n_split),
            best_window_index=int(i_win),
            zenith_deg=z, elevation_deg=90.0 - z,
            local_null_scale=local_scale,
            snr_normalised=s_norm,
            threshold_here=thr_norm * local_scale,
            nearest_catalogue_deg=sep, nearest_catalogue=cat_names[j],
            passes_single_look=bool(s_norm >= thr_single_norm),
            not_near_catalogue=bool(sep > beam_deg),
            confirmed=bool(s_norm >= thr_norm and sep > beam_deg)))

    thr_single = thr_single_norm * float(znull.scale(0.0))
    thr_corr = thr_norm * float(znull.scale(0.0))
    big = samples

    n_single = sum(c["passes_single_look"] for c in cands)
    n_unexp = sum(c["passes_single_look"] and c["not_near_catalogue"] for c in cands)
    n_conf = sum(c["confirmed"] for c in cands)
    log.info("%d peaks -> %d passed single-look -> %d unexplained -> %d confirmed",
             len(cands), n_single, n_unexp, n_conf)

    return dict(
        window_map=window_map,
        ra0_deg=float(ra0_deg), dec0_deg=float(dec0_deg),
        n_peaks=len(cands), n_peaks_from_image=int(n_image),
        n_passed_single_look=int(n_single),
        n_unexplained=int(n_unexp), n_confirmed=int(n_conf),
        time_windowed=bool(windowed), time_windows=splits,
        threshold_single_look_at_zenith=float(thr_single),
        threshold_trials_corrected_at_zenith=float(thr_corr),
        threshold_normalised=float(thr_norm),
        zenith_null=znull.describe(),
        null_median=float(np.median(samples)) if len(samples) else None,
        null_samples=int(len(big)),
        expected_false_positives=float(len(cands) * sig_cfg.get("false_alarm_rate", 0.01)),
        image_rms=float(rms), candidates=cands)

def snr_at_windows(uvw, freqs_hz, data, times, ra0_deg, dec0_deg,
                   ra_deg, dec_deg, splits=(1, 2, 4, 8)):
    ut = np.unique(times)
    best = (0.0, 1, 0)
    for k in splits:
        if k > 1 and len(ut) < 4 * k:
            continue
        for j, chunk in enumerate(np.array_split(ut, k)):
            m = np.isin(times, chunk)
            if m.sum() < 8:
                continue
            v = snr_at(uvw[m], freqs_hz, data[m], ra0_deg, dec0_deg,
                       ra_deg, dec_deg)
            if v > best[0]:
                best = (float(v), int(k), int(j))
    return best

def refine_position(uvw, freqs_hz, data, times, ra0_deg, dec0_deg,
                    ra_deg, dec_deg, radius_deg=2.0, steps=5, rounds=3,
                    splits=None):
    best = (ra_deg, dec_deg, -1.0, 1, 0)
    span = float(radius_deg)
    for _ in range(int(rounds)):
        ra_c, dec_c = best[0], best[1]
        cosd = max(np.cos(np.radians(dec_c)), 1e-3)
        for dra in np.linspace(-span, span, steps):
            for ddec in np.linspace(-span, span, steps):
                ra = ra_c + dra / cosd
                dec = float(np.clip(dec_c + ddec, -89.9, 89.9))
                if splits and times is not None:
                    v, k, i = snr_at_windows(uvw, freqs_hz, data, times,
                                             ra0_deg, dec0_deg, ra, dec, splits)
                else:
                    v = snr_at(uvw, freqs_hz, data, ra0_deg, dec0_deg, ra, dec)
                    k, i = 1, 0
                if v > best[2]:
                    best = (ra, dec, v, k, i)
        span /= float(steps - 1)
    return best
