"""DiSkO imaging, and the FITS pixel-scale correction it needs."""
from __future__ import annotations

import logging
import shutil
import subprocess
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

log = logging.getLogger(__name__)

DEG_PER_RADIAN = 57.29577951308232


def disko_cdelt(npix: int) -> float:
    return DEG_PER_RADIAN / (npix / 2.0)


def load_wcs(header) -> WCS:
    origin = str(header.get("ORIGIN", "") or "")
    npix = header.get("NAXIS1")
    if origin.strip().startswith("DiSkO") and npix:
        header = header.copy()
        scale = disko_cdelt(int(npix))
        header["CDELT1"] = -abs(scale)
        header["CDELT2"] = abs(scale)
    return WCS(header).celestial


def load_image(path: Path):
    with fits.open(str(path)) as hdul:
        data = np.squeeze(hdul[0].data).astype(float)
        return data, hdul[0].header


@dataclass
class DiskoResult:
    fits_path: Path
    returncode: int
    stdout: str
    stderr: str
    cmd: List[str] = field(default_factory=list)
    siblings: List[str] = field(default_factory=list)

    def failure_message(self) -> str:
        if self.returncode in (137, -9):
            head = ("DiSkO was killed by SIGKILL (exit 137) -- this is the "
                    "out-of-memory killer, not a DiSkO error. It writes no "
                    "message when this happens. Reduce 'res' or 'nvis'.")
        elif self.returncode == 0:
            head = ("DiSkO exited 0 but wrote no FITS matching "
                    f"{self.fits_path.name}")
        else:
            head = f"DiSkO exited {self.returncode}"
        parts = [head, "command: " + " ".join(self.cmd)]
        found = ", ".join(self.siblings) if self.siblings else "(nothing)"
        parts.append(f"files in {self.fits_path.parent}: {found}")
        for name, text in (("stderr", self.stderr), ("stdout", self.stdout)):
            text = (text or "").strip()
            if text:
                parts.append(f"--- {name} (tail) ---\n" + text[-2000:])
        return "\n".join(parts)


def healpix_npix_in_fov(fov_deg: float, res_deg: float) -> int:
    import healpy as hp
    res_arcmin = res_deg * 60.0
    nside = 1
    while hp.nside2resol(nside, arcmin=True) > res_arcmin:
        nside *= 2
    cap = (1.0 - np.cos(np.radians(fov_deg / 2.0))) / 2.0
    return int(hp.nside2npix(nside) * cap)


def estimate_memory_gb(fov_deg: float, res_deg: float, n_vis: int) -> float:
    return n_vis * healpix_npix_in_fov(fov_deg, res_deg) * 48 / 1e9


def _available_gb() -> Optional[float]:
    try:
        for line in Path("/proc/meminfo").read_text().splitlines():
            if line.startswith("MemAvailable:"):
                return float(line.split()[1]) * 1024 / 1e9
    except OSError:
        pass
    return None


def _parse_deg(v) -> Optional[float]:
    t = str(v).strip().lower().replace("deg", "").strip()
    try:
        return float(t)
    except ValueError:
        return None


def check_memory(params: dict, safety: float = 2.0) -> None:
    fov = _parse_deg(params.get("fov", params.get("fov-deg", "")))
    res = _parse_deg(params.get("res", params.get("res-deg", "")))
    nvis = params.get("nvis")
    if not (fov and res and nvis) or not params.get("healpix"):
        return
    need = estimate_memory_gb(fov, res, int(nvis)) * safety
    have = _available_gb()
    npix = healpix_npix_in_fov(fov, res)
    log.info("DiSkO solve grid: %d cells x %s visibilities -> ~%.1f GB needed"
             " (%.1f GB available)", npix, nvis, need,
             have if have else float("nan"))
    if have is not None and need > have:
        raise MemoryError(
            f"DiSkO would need ~{need:.1f} GB but only {have:.1f} GB is "
            f"available.\n"
            f"  fov={fov}deg res={res}deg -> {npix} HEALPix cells\n"
            f"  nvis={nvis} visibilities\n"
            f"  peak = n_vis * n_pix * 48 bytes (three operator copies)\n"
            f"Reduce 'res' (coarser grid) or 'nvis' in the config, or give "
            f"Docker more memory. Left alone this is killed by the OOM killer "
            f"with an empty error message and exit code 137.")

def image(ms_path: Path, out_fits: Path, column: str = "DATA",
          params: Optional[dict] = None, timeout_s: int = 1800) -> DiskoResult:
    if shutil.which("disko") is None:
        raise RuntimeError("disko not found. pip install disko")

    out_fits = Path(out_fits)
    out_fits.parent.mkdir(parents=True, exist_ok=True)
    if out_fits.exists():
        out_fits.unlink()
    for stale in out_fits.parent.glob(f"{out_fits.stem}_*.fits"):
        stale.unlink()

    args = dict(params or {})
    cmd = ["disko", "--ms", str(ms_path), "--column", column,
           "--dir", str(out_fits.parent), "--title", out_fits.stem, "--FITS"]
    for k, v in args.items():
        flag = "--" + k.replace("_", "-")
        if isinstance(v, bool):
            if v:
                cmd.append(flag)
        else:
            cmd += [flag, str(v)]

    check_memory(args)

    started = time.time()
    proc = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout_s)

    if not out_fits.exists():
        fresh = [p for p in out_fits.parent.glob(f"{out_fits.stem}_*.fits")
                 if p.stat().st_mtime >= started]
        fresh.sort(key=lambda p: p.stat().st_mtime, reverse=True)
        if fresh:
            fresh[0].rename(out_fits)

    siblings = sorted(q.name for q in out_fits.parent.iterdir()) \
        if out_fits.parent.is_dir() else []
    rc = proc.returncode
    if rc == 0 and not out_fits.exists():
        rc = 0
        return DiskoResult(out_fits, rc, proc.stdout, proc.stderr, cmd, siblings)
    return DiskoResult(out_fits, rc, proc.stdout, proc.stderr, cmd, siblings)


def find_peaks(img: np.ndarray, fwhm_pix: float, sigma: float = 5.0):
    from astropy.stats import sigma_clipped_stats
    from photutils.detection import DAOStarFinder

    mask = np.isfinite(img)
    if not mask.any():
        return None, np.nan, np.nan
    _mean, median, std = sigma_clipped_stats(img[mask], sigma=3.0)
    if not np.isfinite(std) or std <= 0:
        return None, median, std
    finder = DAOStarFinder(fwhm=float(fwhm_pix), threshold=sigma * std)
    return finder(np.where(mask, img, median) - median), median, std


def lm_to_radec(l, m, ra0_deg, dec0_deg):
    ra0, dec0 = np.radians(ra0_deg), np.radians(dec0_deg)
    n = np.sqrt(np.clip(1.0 - l*l - m*m, 0.0, 1.0))
    dec = np.arcsin(np.clip(m*np.cos(dec0) + n*np.sin(dec0), -1.0, 1.0))
    ra = ra0 + np.arctan2(l, n*np.cos(dec0) - m*np.sin(dec0))
    return np.degrees(ra) % 360.0, np.degrees(dec)


def time_window_peak_map(uvw, freqs_hz, data, times, ra0_deg, dec0_deg,
                         n_side=110, n_windows=6, max_zenith_deg=85.0,
                         chunk=400):
    from .fitting import C_LIGHT, PHASE_SIGN

    ut = np.unique(times)
    n_windows = int(max(1, min(n_windows, len(ut))))
    gx = np.linspace(-1.0, 1.0, int(n_side))
    L, M = np.meshgrid(gx, gx)
    lim = np.sin(np.radians(max_zenith_deg))
    good = (L*L + M*M) <= lim*lim
    N = np.sqrt(np.clip(1.0 - L*L - M*M, 0.0, 1.0))
    S = np.stack([L[good], M[good], N[good] - 1.0], axis=1)

    wl = C_LIGHT / np.asarray(freqs_hz, dtype=float)
    best = np.full(S.shape[0], -np.inf)
    which = np.zeros(S.shape[0], dtype=int)

    for w, ch in enumerate(np.array_split(ut, n_windows)):
        sel = np.isin(times, ch)
        if sel.sum() < 4:
            continue
        u = uvw[sel][:, 0][:, None] / wl[None, :]
        v = uvw[sel][:, 1][:, None] / wl[None, :]
        ww = uvw[sel][:, 2][:, None] / wl[None, :]
        V = data[sel].reshape(sel.sum(), data.shape[1], -1)[:, :, 0]
        V = V.astype(np.complex128).reshape(-1)
        base = np.column_stack([u.reshape(-1), v.reshape(-1), ww.reshape(-1)])
        for i in range(0, S.shape[0], chunk):
            Mm = np.exp(PHASE_SIGN * 2.0j * np.pi * (base @ S[i:i+chunk].T))
            num = np.abs((np.conj(Mm) * V[:, None]).sum(axis=0))
            den = (np.abs(Mm) ** 2).sum(axis=0)
            val = num / np.maximum(den, 1e-30)
            upd = val > best[i:i+chunk]
            best[i:i+chunk][upd] = val[upd]
            which[i:i+chunk][upd] = w
        log.info("  window %d/%d imaged", w + 1, n_windows)

    out = np.full(L.shape, np.nan)
    win = np.full(L.shape, -1, dtype=int)
    out[good] = best
    win[good] = which
    return out, win, gx, gx


def peaks_from_map(amap, wmap, gx, ra0_deg, dec0_deg, n_peaks=60, min_sep_pix=3):
    finite = np.isfinite(amap)
    if not finite.any():
        return []
    flat = np.where(finite, amap, -np.inf)
    order = np.argsort(flat, axis=None)[::-1]
    taken, out = [], []
    for idx in order:
        if len(out) >= n_peaks:
            break
        i, j = np.unravel_index(idx, flat.shape)
        if not np.isfinite(flat[i, j]):
            continue
        if any((i-a)**2 + (j-b)**2 < min_sep_pix**2 for a, b in taken):
            continue
        taken.append((i, j))
        ra, dec = lm_to_radec(gx[j], gx[i], ra0_deg, dec0_deg)
        out.append({"ra_deg": float(ra), "dec_deg": float(dec),
                    "map_value": float(flat[i, j]),
                    "window": int(wmap[i, j])})
    return out
