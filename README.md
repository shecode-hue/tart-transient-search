# TART transient search

A transient-detection pipeline for the [TART](https://tart.elec.ac.nz) radio
telescope array, and the detection of a solar radio burst on 2025-11-11.

TART is a 24-antenna all-sky array observing at 1575.42 MHz. It records
correlations between antenna pairs; images are computed afterwards and are used
only to generate candidate positions. Every detection decision is made on the
raw visibilities.

The 2025-11-11 burst (onset 10:01:11 UTC, peak 10:03:11 UTC, 13.3 sigma at
na-unam, seen at four sites, recorded independently by NOAA/RSTN) is in
`results/unam_event_20251111/`, with figures, FITS and tables per epoch.

## Quick start

```bash
python3 scripts/apply_patches.py
PYTHONPATH=src python3 -m pytest tests/ -q

PYTHONPATH=src python3 -m tart_transient timeseries \
  --sites na-unam --peak "2025-11-11T10:00:00+00:00" \
  --window 180 --interval 30 --out runs/ts

PYTHONPATH=src python3 -m tart_transient run --config config/unam-ev-peak.yaml
```

## Stages

| stage | what it does |
|---|---|
| download | fetch HDF, selected by observation time |
| build | measurement set, per-second satellite positions, calibration decision |
| search | model, subtract, image, find peaks, test each |
| run | all three |
| timeseries | fit every catalogued source per interval across hours |
| transient | light curve at one direction, plus control directions |
| compare | difference two completed runs |

## Layout

```
src/tart_transient/   the pipeline
tests/                convention tests
config/               epoch configs, plus example.yaml
scripts/              apply_patches.py
results/              figures, FITS and tables per run
```

## Requirements

Python 3.12, tart2ms, disko, casacore, casatasks, astropy, h5py, minio.
Two upstream patches are required; see PATCHES.md, applied by
`scripts/apply_patches.py`.
