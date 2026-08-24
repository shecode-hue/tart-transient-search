"""Run configuration."""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict

import yaml


@dataclass
class Config:
    """Parsed run configuration, plus the derived output paths."""

    name: str
    raw: Dict[str, Any] = field(repr=False)
    root: Path

    @classmethod
    def load(cls, path: str | Path, runs_dir: str | Path = "runs") -> "Config":
        path = Path(path)
        with open(path) as f:
            raw = yaml.safe_load(f)
        name = raw.get("name") or path.stem
        return cls(name=name, raw=raw, root=Path(runs_dir) / name)

    def section(self, key: str) -> Dict[str, Any]:
        return dict(self.raw.get(key, {}) or {})

    @property
    def data_dir(self) -> Path:
        return self.root / "data"

    @property
    def ms_dir(self) -> Path:
        return self.root / "ms"

    @property
    def fits_dir(self) -> Path:
        return self.root / "fits"

    @property
    def plots_dir(self) -> Path:
        return self.root / "plots"

    @property
    def results_dir(self) -> Path:
        return self.root / "results"

    @property
    def ms_path(self) -> Path:
        return self.ms_dir / f"{self.name}.ms"

    def make_dirs(self) -> None:
        for d in (self.data_dir, self.ms_dir, self.fits_dir,
                  self.plots_dir, self.results_dir):
            d.mkdir(parents=True, exist_ok=True)
