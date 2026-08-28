# TART transient search

A transient-detection pipeline for the TART radio telescope array, and the
detection of a solar radio burst on 2025-11-11.

TART is a 24-antenna all-sky array observing at 1575.42 MHz. It records
correlations between antenna pairs, not images — 276 numbers per second. Images
are computed afterwards and are used here only to generate candidate positions;
every detection decision is made on the raw visibilities.

---

## The result

A solar radio burst, found with no prior knowledge of its position or time, and
independently recorded by NOAA/RSTN at 53,000 sfu, 1415 MHz.

**Time series** — 3 hours at na-unam, 176 files, 361 intervals of 30 s:

| | |
|---|---|
| onset | 10:01:11 UTC |
| peak | 10:03:11 UTC |
| end | 10:03:41 UTC |
| duration | 150 s (120 s rise, 30 s decay) |
| significance | 13.3 sigma; next source 5.7 |

**Imaging** — the full pipeline run separately at five moments:

| epoch | time UTC | separation | SNR | threshold | rank |
|---|---|---|---|---|---|
| pre | 09:54:30 | — | — | — | nothing within a beam |
| onset | 10:00:39 | 0.77 deg | 29.6 | 28.1 | 2 of 63 |
| peak | 10:02:42 | 0.20 deg | 55.5 | 38.5 | 1 of 59 |
| decay | 10:03:44 | 0.40 deg | 26.7 | 29.2 | 2 of 68 |
| post | 10:09:53 | — | — | — | nothing within a beam |

`results/unam_event_20251111/figures/epochs.png`

**Multi-site** — the same burst at na-unam, za-rhodes, ghana and mu-udm, peak
times agreeing to 0.01 min, cross-correlation r = +0.978, and a 14.5 sigma
excess toward the Sun over eight control directions.

---

## Quick start

```bash
python3 scripts/apply_patches.py
PYTHONPATH=src python3 -m pytest tests/ -q

PYTHONPATH=src python3 -m tart_transient timeseries \
  --sites na-unam --peak "2025-11-11T10:00:00+00:00" \
  --window 180 --interval 30 --out runs/ts

PYTHONPATH=src python3 -m tart_transient run --config config/unam-ev-peak.yaml
```

`docs/COMMANDS.md` has the full sequence from a cold machine.

---

## Stages

| stage | what it does |
|---|---|
| `download` | fetch HDF, selected by observation time |
| `build` | measurement set, per-second satellite positions, calibration decision |
| `search` | model, subtract, image, find peaks, test each |
| `run` | all three |
| `timeseries` | fit every catalogued source per interval across hours |
| `transient` | light curve at one direction, plus control directions |
| `compare` | difference two completed runs |

---

## How a detection is decided

Images generate candidate positions. Nothing else. Every significance test runs
on the raw visibilities.

A candidate must pass four gates:

1. **Above the measured noise.** 20,000 random empty directions are fitted in
   the same data. A fixed 5 sigma cut is invalid here — 64-86% of empty sky
   passes it, because the same 276 baselines are re-measured every second.
2. **Above the noise at its own elevation.** The null tail is ~35% heavier below
   20 degrees.
3. **Above a trials-corrected threshold**, obtained by fitting the tail of the
   null rather than reading an order statistic.
4. **Not within one beam of a catalogued object**, checked against the full
   catalogue rather than the modelled subset.

A source that brightens partway through an observation is invisible to a single
image of the whole span: one amplitude is fitted, so subtraction removes too
much early and too little at peak, and the average cancels. Each time window is
therefore imaged separately and the pixel-wise maximum kept.

---

## Layout

```
README.md             this file
src/tart_transient/   the pipeline, 16 modules
tests/                convention tests
config/               the five epoch configs, plus example.yaml
scripts/              apply_patches.py
docs/                 RESEARCH.md TRANSIENTS.md COMMANDS.md FLOW.md RUNBOOK.md
results/              figures, FITS and tables
```

`results/unam_event_20251111/` is organised by epoch:

```
fits/pre/  fits/onset/  fits/peak/  fits/decay/  fits/post/
figures/pre/ ...        figures/epochs.png  lightcurve.png  timeseries.png
tables/pre/ ...         events.json  ranked_variability.json
```

---

## Requirements

Python 3.12, `tart2ms`, `disko`, `casacore`, `casatasks`, `astropy`, `h5py`,
`minio`. Two upstream patches are required — see `PATCHES.md`, applied by
`scripts/apply_patches.py`.

---

## Instrument limits

TART's detection floor is of order 10^6 Jy in a 1-second sample. Cassiopeia A,
the brightest steady radio source in the sky, is ~500x below it. Of the known
radio transient classes only solar bursts are reachable; `docs/TRANSIENTS.md`
gives the full comparison. Fast transients are additionally diluted by the
1-second integration.
