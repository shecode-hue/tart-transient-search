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
    """Coherent fit SNR at one sky position."""
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
        catalogue, beam_deg, fwhm_pix, cfg) -> dict:
    """Find residual peaks and decide which, if any, are real."""
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

    log.info("residual peaks to test: %d", len(positions))
    if not positions:
        return dict(n_peaks=0, n_passed_single_look=0, n_unexplained=0,
                    n_confirmed=0, candidates=[], image_rms=float(rms))

    samples, zeniths = null_with_zenith(
        uvw, freqs_hz, residual_vis, ra0_deg, dec0_deg,
        n_trials=sig_cfg.get("trials_null_samples", 8000), seed=99)
    znull = ZenithNull(samples, zeniths)
    alpha = sig_cfg.get("false_alarm_rate", 0.01)
    thr_norm = znull.threshold(len(positions), alpha)
    thr_single_norm = float(np.percentile(
        znull.normalised, sig_cfg.get("null_percentile", 99.0)))
    log.info("zenith-matched null: %s", znull.describe())

    cat_ra = catalogue["ra_d"].values.astype(float)
    cat_dec = catalogue["dec_d"].values.astype(float)

    cands = []
    for ra, dec, peak in positions:
        s = snr_at(uvw, freqs_hz, residual_vis, ra0_deg, dec0_deg, ra, dec)
        sep = float(angular_sep(ra, dec, cat_ra, cat_dec).min())
        z = float(zenith_angle_deg(ra, dec, ra0_deg, dec0_deg))
        s_norm = znull.snr_normalised(s, z)
        local_scale = float(znull.scale(z))
        cands.append(dict(
            ra_deg=ra, dec_deg=dec, image_peak=peak, vis_snr=s,
            zenith_deg=z, elevation_deg=90.0 - z,
            local_null_scale=local_scale,
            snr_normalised=s_norm,
            threshold_here=thr_norm * local_scale,
            nearest_catalogue_deg=sep,
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
        n_peaks=len(cands), n_passed_single_look=int(n_single),
        n_unexplained=int(n_unexp), n_confirmed=int(n_conf),
        threshold_single_look_at_zenith=float(thr_single),
        threshold_trials_corrected_at_zenith=float(thr_corr),
        threshold_normalised=float(thr_norm),
        zenith_null=znull.describe(),
        null_median=float(np.median(samples)) if len(samples) else None,
        null_samples=int(len(big)),
        expected_false_positives=float(len(cands) * sig_cfg.get("false_alarm_rate", 0.01)),
        image_rms=float(rms), candidates=cands)
