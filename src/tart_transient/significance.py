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
            continue
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
    if not len(samples):
        return float("inf")
    q = 100.0 * (1.0 - false_alarm_rate / max(n_looks, 1))
    return float(np.percentile(samples, min(q, 100.0)))

def zenith_angle_deg(ra_deg, dec_deg, ra0_deg: float, dec0_deg: float):
    r1, d1 = np.radians(ra_deg), np.radians(dec_deg)
    r0, d0 = np.radians(ra0_deg), np.radians(dec0_deg)
    cosang = (np.sin(d1) * np.sin(d0)
              + np.cos(d1) * np.cos(d0) * np.cos(r1 - r0))
    return np.degrees(np.arccos(np.clip(cosang, -1.0, 1.0)))


def null_with_zenith(uvw, freqs_hz, data, ra0_deg: float, dec0_deg: float,
                     n_trials: int = 20000, seed: int = 4242,
                     max_zenith_deg: float = 85.0, times=None, splits=None):
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
        if splits and times is not None:
            from .search import snr_at_windows
            val = snr_at_windows(uvw, freqs_hz, data, times, ra0_deg, dec0_deg,
                                 ra, dec, splits)[0]
            if val > 0:
                snrs.append(float(val))
                zen.append(z)
            continue
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
        if not centres:
            centres = [0.0, max_zenith_deg]
            g = float(np.percentile(samples, scale_percentile)) if len(samples) else 1.0
            scales = [g, g]
        self._centres = np.array(centres, dtype=float)
        self._scales = np.array(scales, dtype=float)
        self._scales[self._scales <= 0] = np.median(self._scales[self._scales > 0]) \
            if np.any(self._scales > 0) else 1.0
        self.normalised = samples / self.scale(zeniths)

    def scale(self, zenith_deg):
        return np.interp(zenith_deg, self._centres, self._scales)

    def threshold(self, n_looks: int, false_alarm_rate: float = 0.01,
                  method: str = "tail") -> float:
        if method == "tail":
            t = tail_threshold(self.normalised, n_looks, false_alarm_rate)
            if np.isfinite(t):
                return t
            log.warning("tail fit failed; falling back to empirical quantile")
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

def tail_threshold(samples, n_looks: int, false_alarm_rate: float = 0.01,
                   n_tail: int = 500, return_fit: bool = False):
    samples = np.asarray(samples, dtype=float)
    samples = samples[np.isfinite(samples)]
    if len(samples) < 20:
        return (float("inf"), {}) if return_fit else float("inf")

    ordered = np.sort(samples)[::-1]
    n_tail = int(max(20, min(n_tail, len(ordered) // 2)))
    x = ordered[:n_tail]
    y = np.arange(1, n_tail + 1, dtype=float)
    keep = x > 0
    x, y = x[keep], y[keep]
    if len(x) < 20 or np.ptp(x) <= 0:
        return (float("inf"), {}) if return_fit else float("inf")

    A = np.column_stack([np.ones_like(x), -x])
    coef, *_ = np.linalg.lstsq(A, np.log(y), rcond=None)
    log_nhat, inv_rho = float(coef[0]), float(coef[1])
    if not np.isfinite(inv_rho) or inv_rho <= 0:
        return (float("inf"), {}) if return_fit else float("inf")
    rhohat = 1.0 / inv_rho

    n_target = len(samples) * false_alarm_rate / max(n_looks, 1)
    thr = rhohat * (log_nhat - np.log(max(n_target, 1e-300)))

    if return_fit:
        pred = np.exp(log_nhat - x / rhohat)
        ss = 1.0 - np.sum((np.log(y) - np.log(pred)) ** 2) / \
            max(np.sum((np.log(y) - np.log(y).mean()) ** 2), 1e-12)
        return float(thr), {"rhohat": rhohat, "log_nhat": log_nhat,
                            "n_tail": int(len(x)), "r2_log": float(ss),
                            "n_target": float(n_target)}
    return float(thr)
