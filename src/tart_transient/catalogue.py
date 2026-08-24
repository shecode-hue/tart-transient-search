"""Satellite sky model: epoch mapping, elevation, and beam merging."""
from __future__ import annotations

import logging
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

log = logging.getLogger(__name__)

CAT_COLUMNS = ["name", "ra_d", "dec_d", "flux", "spi", "freq0"]
CAT_FNAME_RE = re.compile(r"model_sources_(\d+)\.txt$")

_cache: Dict[str, pd.DataFrame] = {}


def load(path: Path) -> pd.DataFrame:
    """Read one ``model_sources_<N>.txt``, with caching."""
    key = str(path)
    if key in _cache:
        return _cache[key]
    try:
        df = pd.read_csv(path, sep=r"\s+", comment="#", names=CAT_COLUMNS, skiprows=1)
        df = df.dropna(subset=["ra_d", "dec_d", "flux"]).reset_index(drop=True)
    except Exception:
        df = pd.DataFrame(columns=CAT_COLUMNS)
    _cache[key] = df
    return df


@dataclass
class CatalogueIndex:
    """Maps observation time to the right catalogue snapshot."""

    indices: List[int]
    paths: List[Path]
    unique_times: Optional[np.ndarray] = field(default=None, repr=False)

    @classmethod
    def build(cls, msdir: Path, ms_path: Optional[Path] = None) -> "CatalogueIndex":
        found = []
        for p in sorted(Path(msdir).glob("model_sources_*.txt")):
            m = CAT_FNAME_RE.search(p.name)
            if m:
                found.append((int(m.group(1)), p))
        if not found:
            raise FileNotFoundError(f"no model_sources_*.txt under {msdir}")
        found.sort(key=lambda t: t[0])

        times = None
        if ms_path is not None:
            from casacore.tables import table
            with table(str(ms_path), ack=False) as t:
                times = np.unique(t.getcol("TIME"))
            if len(times) != len(found):
                log.warning(
                    "catalogue/MS epoch mismatch: %d files vs %d distinct TIMEs. "
                    "The 1:1 correspondence is the basis of the exact mapping.",
                    len(found), len(times))
        return cls(indices=[i for i, _ in found],
                   paths=[p for _, p in found], unique_times=times)

    def index_at(self, t: float) -> int:
        if self.unique_times is None or not len(self.unique_times):
            raise RuntimeError("build the index with ms_path= for exact mapping")
        pos = int(np.searchsorted(self.unique_times, t))
        if pos > 0 and (pos >= len(self.unique_times) or
                        abs(self.unique_times[pos - 1] - t) <=
                        abs(self.unique_times[pos] - t)):
            pos -= 1
        return self.indices[int(np.clip(pos, 0, len(self.indices) - 1))]

    def path_at(self, t: float) -> Path:
        return self.paths[self.indices.index(self.index_at(t))]

    def at(self, t: float, min_flux: float = 0.0) -> pd.DataFrame:
        df = load(self.path_at(t))
        if min_flux > 0 and len(df):
            df = df[df["flux"] >= min_flux].reset_index(drop=True)
        return df

    def positions_for(self, times: np.ndarray, names: List[str]) -> np.ndarray:
        """Per-row (RA, Dec) for each named source. Shape (n_src, n_row, 2)."""
        if self.unique_times is None:
            raise RuntimeError("build the index with ms_path= for exact mapping")
        ep = np.clip(np.searchsorted(self.unique_times, times),
                     0, len(self.unique_times) - 1)
        lo = np.clip(ep - 1, 0, len(self.unique_times) - 1)
        take_lo = (np.abs(self.unique_times[lo] - times) <=
                   np.abs(self.unique_times[ep] - times))
        ep = np.where(take_lo, lo, ep)

        idx = np.asarray(self.indices)
        out = np.full((len(names), len(times), 2), np.nan)
        for e in np.unique(ep):
            sel = np.where(ep == e)[0]
            df = load(self.path_at(float(self.unique_times[e])))
            if not len(df):
                continue
            lut = {n: (r, d) for n, r, d in
                   zip(df["name"].values, df["ra_d"].values, df["dec_d"].values)}
            for si, nm in enumerate(names):
                hit = lut.get(nm)
                if hit is not None:
                    out[si, sel, 0] = hit[0]
                    out[si, sel, 1] = hit[1]
        return out


def merge_within_beam(cat: pd.DataFrame, beam_deg: float) -> Tuple[pd.DataFrame, dict]:
    """Collapse sources the array cannot separate."""
    if not len(cat):
        return cat, {}

    ra = np.radians(cat["ra_d"].values.astype(float))
    dec = np.radians(cat["dec_d"].values.astype(float))
    n = len(cat)
    radius = beam_deg / 2.0
    assigned = -np.ones(n, dtype=int)
    reps: List[int] = []

    for i in np.argsort(-cat["flux"].values.astype(float)):
        if assigned[i] >= 0:
            continue
        r = len(reps)
        reps.append(int(i))
        assigned[i] = r
        cosd = (np.sin(dec[i]) * np.sin(dec) +
                np.cos(dec[i]) * np.cos(dec) * np.cos(ra[i] - ra))
        sep = np.degrees(np.arccos(np.clip(cosd, -1.0, 1.0)))
        assigned[(sep < radius) & (assigned < 0)] = r

    groups, rows = {}, []
    for r, i in enumerate(reps):
        members = [str(cat.iloc[j]["name"]) for j in np.where(assigned == r)[0]]
        rep = cat.iloc[i].copy()
        rep["track_name"] = str(cat.iloc[i]["name"])
        if len(members) > 1:
            k = np.where(assigned == r)[0]
            w = cat["flux"].values[k].astype(float)
            w = w / w.sum() if w.sum() > 0 else np.ones(len(k)) / len(k)
            rep["ra_d"] = float((cat["ra_d"].values[k].astype(float) * w).sum())
            rep["dec_d"] = float((cat["dec_d"].values[k].astype(float) * w).sum())
            rep["name"] = f"{rep['name']} [+{len(members) - 1} blended]"
        groups[str(rep["name"])] = members
        rows.append(rep)

    merged = pd.DataFrame(rows).reset_index(drop=True)
    log.info("catalogue %d -> %d after beam merge (%.2f deg)",
             len(cat), len(merged), beam_deg)
    return merged, groups
