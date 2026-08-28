"""Build a Measurement Set from HDF, and calibrate it if needed."""
from __future__ import annotations

import logging
import shutil
import subprocess
from pathlib import Path
from typing import List, Optional

log = logging.getLogger(__name__)


def build(hdf_files: List[Path], ms_path: Path, cat_prefix: Path,
          filter_elevation_deg: float = 5.0,
          timeout_s: Optional[int] = None) -> Path:
    if shutil.which("tart2ms") is None:
        raise RuntimeError("tart2ms not found. pip install tart2ms")

    ms_path = Path(ms_path)
    if ms_path.exists():
        shutil.rmtree(ms_path)
    ms_path.parent.mkdir(parents=True, exist_ok=True)
    Path(cat_prefix).parent.mkdir(parents=True, exist_ok=True)

    cmd = ["tart2ms", "--ms", str(ms_path), "--rephase", "obs-midpoint",
           "--add-model", "--clobber", "--write-model-catalog",
           "--model-catalog-name-prefix", str(cat_prefix),
           "--filter-elevation", str(filter_elevation_deg),
           "--hdf"] + [str(p) for p in hdf_files]
    log.info("tart2ms: %d hdf file(s) -> %s", len(hdf_files), ms_path)
    if timeout_s is None:
        timeout_s = 1800 + 1800 * len(hdf_files)
    log.info("tart2ms timeout %d s for %d file(s)", timeout_s, len(hdf_files))
    proc = subprocess.run(cmd, capture_output=True, text=True,
                          timeout=timeout_s)
    if proc.returncode != 0:
        if "--filter-elevation" in proc.stderr or "unrecognized" in proc.stderr:
            raise RuntimeError(
                "tart2ms does not accept --filter-elevation. Run "
                "scripts/apply_patches.py first (see PATCHES.md).")
        raise RuntimeError("tart2ms failed:\n" + proc.stderr[-3000:])
    if not ms_path.exists():
        raise RuntimeError("tart2ms reported success but wrote no MS")
    return ms_path


def calibrate(ms_path: Path, caltable_dir: Path, minsnr: float = 2.0,
              solint: str = "inf", applymode: str = "calonly") -> None:
    from casatasks import applycal, gaincal

    caltable_dir = Path(caltable_dir)
    if caltable_dir.exists():
        shutil.rmtree(caltable_dir)
    caltable_dir.mkdir(parents=True, exist_ok=True)
    amp, phase = caltable_dir / "G0a", caltable_dir / "G0p"

    log.info("gaincal: amplitude")
    gaincal(vis=str(ms_path), caltable=str(amp), gaintype="G", calmode="a",
            solint=solint, minsnr=minsnr, solnorm=True, refant="0")
    log.info("gaincal: phase")
    gaincal(vis=str(ms_path), caltable=str(phase), gaintype="G", calmode="p",
            solint=solint, minsnr=minsnr, refant="0", gaintable=[str(amp)])
    log.info("applycal")
    applycal(vis=str(ms_path), gaintable=[str(amp), str(phase)],
             applymode=applymode, flagbackup=False)


def observation_geometry(ms_path: Path):
    import numpy as np
    from casacore.tables import table

    with table(str(ms_path) + "/FIELD", ack=False) as t:
        pd = t.getcol("PHASE_DIR")
    ra0, dec0 = np.degrees(pd.reshape(-1, 2)[0])
    with table(str(ms_path) + "/SPECTRAL_WINDOW", ack=False) as t:
        freqs = np.asarray(t.getcol("CHAN_FREQ")[0], dtype=float)
    with table(str(ms_path), ack=False) as t:
        uvw = t.getcol("UVW")
        times = t.getcol("TIME")
        n_rows = t.nrows()

    lam = 299792458.0 / float(freqs[0])
    bmax = float(np.max(np.sqrt((uvw ** 2).sum(axis=1))))
    return dict(ra0_deg=float(ra0), dec0_deg=float(dec0), freqs_hz=freqs,
                n_rows=int(n_rows), n_integrations=int(len(np.unique(times))),
                max_baseline_m=bmax, max_baseline_lambda=bmax / lam,
                beam_deg=float(np.degrees(lam / bmax)))
