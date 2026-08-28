"""Coherent DFT model, joint amplitude fit, and peeling."""
from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import List, Optional, Tuple

import numpy as np

log = logging.getLogger(__name__)

C_LIGHT = 299792458.0


@dataclass
class PointSource:
    name: str
    ra_deg: float
    dec_deg: float
    flux_jy: float = 1.0


def lmn(ra_deg, dec_deg, ra0_deg, dec0_deg):
    ra, dec = np.radians(ra_deg), np.radians(dec_deg)
    ra0, dec0 = np.radians(ra0_deg), np.radians(dec0_deg)
    l = np.cos(dec) * np.sin(ra - ra0)
    m = np.sin(dec) * np.cos(dec0) - np.cos(dec) * np.sin(dec0) * np.cos(ra - ra0)
    n = np.sqrt(np.clip(1.0 - l * l - m * m, 0.0, None))
    return l, m, n


PHASE_SIGN = +1.0


def model_matrix(uvw: np.ndarray, freqs_hz: np.ndarray,
                 tracks: np.ndarray, ra0_deg: float, dec0_deg: float) -> np.ndarray:
    n_row, n_chan, n_src = uvw.shape[0], freqs_hz.shape[0], tracks.shape[0]
    mat = np.zeros((n_row, n_chan, n_src), dtype=np.complex128)
    wl = C_LIGHT / freqs_hz
    u = uvw[:, 0][:, None] / wl[None, :]
    v = uvw[:, 1][:, None] / wl[None, :]
    w = uvw[:, 2][:, None] / wl[None, :]

    for i in range(n_src):
        ra, dec = tracks[i, :, 0], tracks[i, :, 1]
        ok = np.isfinite(ra) & np.isfinite(dec)
        if not ok.any():
            continue
        l, m, n = lmn(ra, dec, ra0_deg, dec0_deg)
        block = np.exp(PHASE_SIGN * 2.0j * np.pi *
                       (u * l[:, None] + v * m[:, None] + w * (n[:, None] - 1.0)))
        block[~ok, :] = 0.0
        mat[:, :, i] = block
    return mat


def fit_amplitudes(data: np.ndarray, model: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    n_row, n_chan, _ = data.shape
    n_src = model.shape[2]
    A = model.reshape(n_row * n_chan, n_src)
    b = data.mean(axis=2).reshape(n_row * n_chan)

    good = np.isfinite(b) & np.all(np.isfinite(A), axis=1)
    A, b = A[good], b[good]
    if A.shape[0] <= n_src:
        return np.zeros(n_src, dtype=np.complex128), np.full(n_src, np.inf)

    amp, *_ = np.linalg.lstsq(A, b, rcond=None)
    resid = b - A @ amp
    dof = max(A.shape[0] - n_src, 1)
    noise_var = float(np.sum(np.abs(resid) ** 2) / dof)
    try:
        cov = noise_var * np.linalg.inv(A.conj().T @ A)
        sigma = np.sqrt(np.abs(np.diag(cov)))
    except np.linalg.LinAlgError:
        sigma = np.full(n_src, np.inf)
    return amp, sigma


def snr_of(amp: np.ndarray, sigma: np.ndarray) -> np.ndarray:
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(sigma > 0, np.abs(amp) / sigma, 0.0)


def peel(data: np.ndarray, model: np.ndarray, amp: np.ndarray) -> np.ndarray:
    n_row, n_chan, n_corr = data.shape
    stokes_i = (model.reshape(n_row * n_chan, -1) @ amp).reshape(n_row, n_chan)
    return data - np.repeat(stokes_i[:, :, None], n_corr, axis=2)


def gram_condition(model: np.ndarray) -> float:
    A = model.reshape(-1, model.shape[2])
    return float(np.linalg.cond(A.conj().T @ A))
