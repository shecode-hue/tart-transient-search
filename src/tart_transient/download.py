"""Fetch visibilities from the public TART archive."""
from __future__ import annotations

import logging
import shutil
import subprocess
from pathlib import Path
from typing import List

log = logging.getLogger(__name__)


def _archive_by_filename(target: str, out_dir: Path, start_iso: str,
                         duration_min: float, n_files: int) -> List[Path]:
    from datetime import datetime, timedelta, timezone
    import tart_tools.archive_handler as ah
    from minio import Minio

    t0 = datetime.fromisoformat(start_iso)
    if t0.tzinfo is None:
        t0 = t0.replace(tzinfo=timezone.utc)
    t1 = t0 + timedelta(minutes=float(duration_min))

    client = Minio(endpoint=ah.MINIO_API_HOST, secure=True)
    picked = []
    day = t0
    while day.date() <= t1.date():
        prefix = "%s/vis/%d/%d/%d/" % (target, day.year, day.month, day.day)
        for obj in client.list_objects(ah.BUCKET_NAME, prefix=prefix,
                                       recursive=True):
            name = obj.object_name.rsplit("/", 1)[-1]
            try:
                stamp = name.replace(".hdf", "").split("_", 1)[1]
                when = datetime.strptime(stamp, "%Y-%m-%d_%H_%M_%S.%f")
            except (ValueError, IndexError):
                continue
            when = when.replace(tzinfo=timezone.utc)
            if t0 <= when <= t1:
                picked.append((when, obj.object_name))
        day = day + timedelta(days=1)

    picked.sort()
    if n_files > 0:
        picked = picked[:n_files]
    out = []
    for when, obj_name in picked:
        dest = out_dir / obj_name.rsplit("/", 1)[-1]
        if not dest.exists():
            client.fget_object(ah.BUCKET_NAME, obj_name, str(dest))
        out.append(dest)
    log.info("archive: selected %d file(s) by observation time, %s .. %s",
             len(out), t0.isoformat(), t1.isoformat())
    return out


def archive(target: str, out_dir: Path, start: str = "-60",
            duration_min: int = 2, n_files: int = 1) -> List[Path]:
    if shutil.which("tart_get_archive_data") is None:
        raise RuntimeError(
            "tart_get_archive_data not found. pip install tart_tools")

    out_dir.mkdir(parents=True, exist_ok=True)

    if not str(start).lstrip().startswith("-"):
        files = _archive_by_filename(target, out_dir, str(start),
                                     duration_min, n_files)
        if files:
            log.info("downloaded %d file(s) to %s", len(files), out_dir)
            return files
        log.warning("no objects matched by observation time; "
                    "falling back to tart_get_archive_data")

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


def describe_all(paths) -> dict:
    parts = [describe(p) for p in paths]
    first, last = parts[0], parts[-1]
    return dict(
        paths=[q["path"] for q in parts],
        n_files=len(parts),
        n_integrations=sum(q["n_integrations"] for q in parts),
        n_baselines=first["n_baselines"],
        t_start=first["t_start"],
        t_end=last["t_end"],
        gains_stored=any(q["gains_stored"] for q in parts),
        gain_min=min(q["gain_min"] for q in parts),
        gain_max=max(q["gain_max"] for q in parts),
        phase_min_rad=min(q["phase_min_rad"] for q in parts),
        phase_max_rad=max(q["phase_max_rad"] for q in parts),
    )
