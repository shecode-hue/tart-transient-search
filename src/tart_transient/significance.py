"""Detection thresholds measured from the data"""
from __future__ import annotations

import logging
from typing import Optional, Tuple

import numpy as np

from .fitting import C_LIGHT, PHASE_SIGN

log = logging.getLogger(__name__)


def empirical_null(uvw: np.ndarray, freqs_hz: np.ndarray, data: np.ndarray,
                   ra0_deg: float, dec0_deg: float, n_trials: int = 400,
                   percentile: float = 99.0, seed: int = 12345,
                   base_model: Optional[np.ndarray] = None) -> Tuple[float, np.ndarray]:
    """SNR distribution of the coherent fit at positions with NO source."""
    rng = np.random.default_rng(seed)
    wl = C_LIGHT / freqs_hz
    u = uvw[:, 0][:, None] / wl[None, :]
    v = uvw[:, 1][:, None] / wl[None, :]
    w = uvw[:, 2][:, None] / wl[None, :]
    V = data.reshape(data.shape[0], data.shape[1], -1)[:, :, 0].astype(np.complex128)

    ra0, dec0 = np.radians(ra0_deg), np.radians(dec0_deg)
    snrs, guard = [], 0
    while len(snrs) < n_trials and guard < n_trials * 20:
        guard += 1
        ra = rng.uniform(0.0, 360.0)
        dec = np.degrees(np.arcsin(rng.uniform(-1.0, 1.0)))
        rr, dd = np.radians(ra), np.radians(dec)
        l = np.cos(dd) * np.sin(rr - ra0)
        m = np.sin(dd) * np.cos(dec0) - np.cos(dd) * np.sin(dec0) * np.cos(rr - ra0)
        if l * l + m * m >= 1.0:
            continue                      # below the horizon of this projection
        n = np.sqrt(max(1.0 - l * l - m * m, 0.0))
        M = np.exp(PHASE_SIGN * 2.0j * np.pi * (u * l + v * m + w * (n - 1.0)))

        if base_model is None:
            denom = np.vdot(M, M)
            if abs(denom) == 0:
                continue
            amp = np.vdot(M, V) / denom
            resid = V - amp * M
            sig = np.sqrt((np.abs(resid) ** 2).sum() / (resid.size - 1) / resid.size)
            if sig > 0 and np.isfinite(sig):
                snrs.append(float(np.abs(amp) / sig))
        else:
            B = base_model.reshape(-1, base_model.shape[2])
            A = np.concatenate([B, M.reshape(-1, 1)], axis=1)
            Vf = V.reshape(-1)
            try:
                coef, *_ = np.linalg.lstsq(A, Vf, rcond=None)
                cov = np.linalg.pinv(A.conj().T @ A)
            except np.linalg.LinAlgError:
                continue
            resid = Vf - A @ coef
            s2 = (np.abs(resid) ** 2).sum() / max(len(Vf) - A.shape[1], 1)
            var = float(np.real(cov[-1, -1])) * s2
            if var > 0 and np.isfinite(var):
                snrs.append(float(np.abs(coef[-1]) / np.sqrt(var)))

    if not snrs:
        return float("inf"), np.array([])
    samples = np.array(snrs)
    return float(np.percentile(samples, percentile)), samples


def trials_threshold(samples: np.ndarray, n_looks: int,
                     false_alarm_rate: float = 0.01) -> float:
    """Threshold corrected for having tested ``n_looks`` positions."""
    if not len(samples):
        return float("inf")
    q = 100.0 * (1.0 - false_alarm_rate / max(n_looks, 1))
    return float(np.percentile(samples, min(q, 100.0)))

def zenith_angle_deg(ra_deg, dec_deg, ra0_deg: float, dec0_deg: float):
    """Angle from the phase centre, which for a zenith-pointing TART is the"""
    r1, d1 = np.radians(ra_deg), np.radians(dec_deg)
    r0, d0 = np.radians(ra0_deg), np.radians(dec0_deg)
    cosang = (np.sin(d1) * np.sin(d0)
              + np.cos(d1) * np.cos(d0) * np.cos(r1 - r0))
    return np.degrees(np.arccos(np.clip(cosang, -1.0, 1.0)))


def null_with_zenith(uvw, freqs_hz, data, ra0_deg: float, dec0_deg: float,
                     n_trials: int = 20000, seed: int = 4242,
                     max_zenith_deg: float = 85.0):
    """Null SNR draws, each tagged with the zenith angle it was drawn at."""
    rng = np.random.default_rng(seed)
    wl = C_LIGHT / freqs_hz
    u = uvw[:, 0][:, None] / wl[None, :]
    v = uvw[:, 1][:, None] / wl[None, :]
    w = uvw[:, 2][:, None] / wl[None, :]
    V = data.reshape(data.shape[0], data.shape[1], -1)[:, :, 0].astype(np.complex128)
    ra0, dec0 = np.radians(ra0_deg), np.radians(dec0_deg)

    snrs, zen, guard = [], [], 0
    while len(snrs) < n_trials and guard < n_trials * 20:
        guard += 1
        ra = rng.uniform(0.0, 360.0)
        dec = np.degrees(np.arcsin(rng.uniform(-1.0, 1.0)))
        z = float(zenith_angle_deg(ra, dec, ra0_deg, dec0_deg))
        if z > max_zenith_deg:
            continue
        rr, dd = np.radians(ra), np.radians(dec)
        l = np.cos(dd) * np.sin(rr - ra0)
        m = np.sin(dd) * np.cos(dec0) - np.cos(dd) * np.sin(dec0) * np.cos(rr - ra0)
        if l * l + m * m >= 1.0:
            continue
        n = np.sqrt(max(1.0 - l * l - m * m, 0.0))
        M = np.exp(PHASE_SIGN * 2.0j * np.pi * (u * l + v * m + w * (n - 1.0)))
        denom = np.vdot(M, M)
        if abs(denom) == 0:
            continue
        amp = np.vdot(M, V) / denom
        resid = V - amp * M
        sig = np.sqrt((np.abs(resid) ** 2).sum() / (resid.size - 1) / resid.size)
        if sig > 0 and np.isfinite(sig):
            snrs.append(float(np.abs(amp) / sig))
            zen.append(z)
    return np.array(snrs), np.array(zen)


class ZenithNull:
    """Position-matched detection threshold."""

    # scale_percentile must be a tail quantile: 90 carries none of the effect
    def __init__(self, samples, zeniths, n_bins: int = 6,
                 scale_percentile: float = 99.0, max_zenith_deg: float = 85.0,
                 min_per_bin: int = 500):
        self.samples, self.zeniths = samples, zeniths
        self.scale_percentile = scale_percentile
        n_bins = max(1, min(n_bins, len(samples) // max(min_per_bin, 1)))
        qs = np.linspace(0.0, 100.0, n_bins + 1)
        edges = np.percentile(zeniths, qs) if len(zeniths) else np.array([0.0, max_zenith_deg])
        edges[0], edges[-1] = 0.0, max(max_zenith_deg, float(edges[-1]))
        centres, scales = [], []
        for lo, hi in zip(edges[:-1], edges[1:]):
            m = (zeniths >= lo) & (zeniths < hi)
            if m.sum() >= 30:
                centres.append(float(np.median(zeniths[m])))
                scales.append(float(np.percentile(samples[m], scale_percentile)))
        if not centres:                       # degenerate: fall back to global
            centres = [0.0, max_zenith_deg]
            g = float(np.percentile(samples, scale_percentile)) if len(samples) else 1.0
            scales = [g, g]
        self._centres = np.array(centres, dtype=float)
        self._scales = np.array(scales, dtype=float)
        self._scales[self._scales <= 0] = np.median(self._scales[self._scales > 0]) \
            if np.any(self._scales > 0) else 1.0
        self.normalised = samples / self.scale(zeniths)

    def scale(self, zenith_deg):
        """Local noise scale, linearly interpolated between bin centres."""
        return np.interp(zenith_deg, self._centres, self._scales)

    def threshold(self, n_looks: int, false_alarm_rate: float = 0.01) -> float:
        """Trials-corrected threshold in NORMALISED units."""
        return trials_threshold(self.normalised, n_looks, false_alarm_rate)

    def snr_normalised(self, snr: float, zenith_deg: float) -> float:
        return float(snr / float(self.scale(zenith_deg)))

    def describe(self) -> dict:
        return {"n_samples": int(len(self.samples)),
                "scale_percentile": self.scale_percentile,
                "bin_centres_deg": [round(float(c), 1) for c in self._centres],
                "bin_scales": [round(float(s), 2) for s in self._scales],
                "scale_ratio_rim_to_zenith":
                    round(float(self._scales[-1] / self._scales[0]), 3)
                    if len(self._scales) > 1 and self._scales[0] > 0 else None}
