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
    """Degrees per pixel that actually describe a DiSkO image.

    The header's own CDELT overstates this by 1.485x. See README.
    """
    return DEG_PER_RADIAN / (npix / 2.0)


def load_wcs(header) -> WCS:
    """Pixel<->sky mapping, with DiSkO's header scale corrected."""
    origin = str(header.get("ORIGIN", "") or "")
    npix = header.get("NAXIS1")
    if origin.strip().startswith("DiSkO") and npix:
        header = header.copy()
        scale = disko_cdelt(int(npix))
        header["CDELT1"] = -abs(scale)   # RA increases toward smaller x
        header["CDELT2"] = abs(scale)
    return WCS(header).celestial


def load_image(path: Path):
    """Return (image, header). NaN outside DiSkO's circular field is preserved."""
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
        """A message that is useful even when DiSkO says nothing at all."""
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
    """Number of HEALPix cells DiSkO will actually solve for."""
    import healpy as hp
    res_arcmin = res_deg * 60.0
    nside = 1
    while hp.nside2resol(nside, arcmin=True) > res_arcmin:
        nside *= 2
    cap = (1.0 - np.cos(np.radians(fov_deg / 2.0))) / 2.0
    return int(hp.nside2npix(nside) * cap)


def estimate_memory_gb(fov_deg: float, res_deg: float, n_vis: int) -> float:
    """Peak DiSkO memory in GB for these parameters."""
    # 48 = three operator copies; calibrated, not derived
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
    """Refuse to start a solve that cannot fit, with the arithmetic shown."""
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
    """Run DiSkO and return the FITS it wrote during THIS call."""
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
        rc = 0          # keep 0: failure_message distinguishes this mode
        return DiskoResult(out_fits, rc, proc.stdout, proc.stderr, cmd, siblings)
    return DiskoResult(out_fits, rc, proc.stdout, proc.stderr, cmd, siblings)


def find_peaks(img: np.ndarray, fwhm_pix: float, sigma: float = 5.0):
    """DAOStarFinder, NaN-aware."""
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
