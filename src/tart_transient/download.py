"""Fetch visibilities from the public TART archive."""
from __future__ import annotations

import logging
import shutil
import subprocess
from pathlib import Path
from typing import List

log = logging.getLogger(__name__)


def archive(target: str, out_dir: Path, start: str = "-60",
            duration_min: int = 2, n_files: int = 1) -> List[Path]:
    """Download HDF visibility files from the TART S3 archive."""
    if shutil.which("tart_get_archive_data") is None:
        raise RuntimeError(
            "tart_get_archive_data not found. pip install tart_tools")

    out_dir.mkdir(parents=True, exist_ok=True)
    cmd = ["tart_get_archive_data", "--target", target, "--start", str(start),
           "--duration", str(duration_min), "--n", str(n_files),
           "--dir", str(out_dir)]
    log.info("downloading: %s", " ".join(cmd))
    proc = subprocess.run(cmd, capture_output=True, text=True, timeout=900)
    if proc.returncode != 0:
        raise RuntimeError("archive download failed:\n" + proc.stderr[-2000:])

    files = sorted(out_dir.glob("*.hdf"))
    if not files:
        raise RuntimeError(
            "download reported success but produced no .hdf files. "
            "Check the telescope name and that data exists for that window.")
    log.info("downloaded %d file(s) to %s", len(files), out_dir)
    return files


def describe(hdf_path: Path) -> dict:
    """Summarise one HDF file, including whether calibration was recorded."""
    import h5py
    import numpy as np

    with h5py.File(hdf_path, "r") as f:
        vis = np.asarray(f["vis"][()])
        gains = np.asarray(f["gains"][()])
        phases = np.asarray(f["phases"][()])
        stamps = [t.decode() for t in f["timestamp"][()]]

    has_gains = not (np.allclose(gains, 1.0) and np.allclose(phases, 0.0))
    return dict(
        path=str(hdf_path),
        n_integrations=int(vis.shape[0]),
        n_baselines=int(vis.shape[1]),
        t_start=stamps[0],
        t_end=stamps[-1],
        gains_stored=bool(has_gains),
        gain_min=float(gains.min()),
        gain_max=float(gains.max()),
        phase_min_rad=float(phases.min()),
        phase_max_rad=float(phases.max()),
    )
