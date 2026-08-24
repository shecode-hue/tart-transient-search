"""Whether to calibrate, and how."""
from __future__ import annotations

import logging
from pathlib import Path

log = logging.getLogger(__name__)


def decide(hdf_info: dict, mode: str = "auto") -> str:
    """Return the MS column to use: ``DATA`` or ``CORRECTED_DATA``."""
    if mode == "never":
        return "DATA"
    if mode == "always":
        return "CORRECTED_DATA"
    if mode != "auto":
        raise ValueError(f"unknown calibration mode: {mode!r}")

    if hdf_info["gains_stored"]:
        log.info("HDF carries real gains (%.2f-%.2f) -> using DATA; "
                 "re-calibrating would suppress sources",
                 hdf_info["gain_min"], hdf_info["gain_max"])
        return "DATA"
    log.info("HDF has no stored calibration -> solving gains, using CORRECTED_DATA")
    return "CORRECTED_DATA"
