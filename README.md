# tart-transient-search

A transient search pipeline for the [TART](https://tart.elec.ac.nz) radio
telescope array. Downloads visibilities from the public archive, builds a
Measurement Set, models and subtracts every catalogued GNSS satellite, then
tests whatever is left for real astrophysical signal.

Built on the standard TART toolchain: `tart_tools`, `tart2ms`, `DiSkO`, CASA.

---

## Quick start

```bash
pip install -r requirements.txt

# fetch data, build the MS, calibrate, search — one command
tart-transient run --config config/za-hammanskraal.yaml
```

Or step by step:

```bash
tart-transient download --config config/za-hammanskraal.yaml
tart-transient build    --config config/za-hammanskraal.yaml
tart-transient search   --config config/za-hammanskraal.yaml
tart-transient report   --config config/za-hammanskraal.yaml
```

Outputs land in `runs/<name>/` — Measurement Set, FITS for CARTA, a figure per
stage, and a JSON record of every number.

---

## What it does

| stage | what happens |
|---|---|
| **download** | Pull HDF visibilities from the TART S3 archive |
| **build** | `tart2ms` → Measurement Set + per-epoch source catalogues; decide whether to calibrate |
| **model** | Catalogue at the horizon, merge sources closer than the beam |
| **fit** | Coherent DFT fit of every source, positions tracked per visibility row |
| **peel** | Subtract the whole known sky model |
| **image** | DiSkO before and after |
| **search** | Test every residual peak in the visibility domain, with a trials-corrected threshold |

---

## Seven things this pipeline does differently

These are not stylistic choices. Each was established by measurement on real
TART data; the numbers below are the record, and the code is kept terse.

**1. The image is not the detector.** DiSkO's regularised reconstruction
correlates only r ≈ 0.55 with a direct transform of the same visibilities and
shares just 4 of its 20 brightest positions — the lasso prior redistributes
flux. Image peaks are treated as *candidate positions*; every one is then tested
by a coherent fit against the visibilities. See `search.py`.

**2. The detection threshold is measured, not assumed.** A fixed 5σ cut is
invalid here: 64–86% of random empty-sky positions pass it, because the
residuals are not independent (276 baselines re-measured every integration
understate σ by roughly √60). The threshold comes from fitting ~400 random sky
directions in the same data. See `significance.py`.

**3. Known satellites are subtracted without a detection test.** They have
published ephemerides — their existence is not in question, only their
brightness. Gating subtraction on significance leaves known sources in the map,
each spraying 10–30% sidelobes. See `fitting.py`.

**4. The satellite catalogue reaches the horizon.** `tart2ms` defaults to a
hardcoded 45° elevation cut; a TART image spans 170° FOV, down to ~5°. On one
test file that meant 14 catalogued satellites out of 67 actually above the
horizon. See `catalogue.py` and `PATCHES.md`.

**5. The detection threshold rises toward the horizon.** The null SNR
distribution is not uniform across a 170 deg field. From 24,000 null draws on
real residual visibilities, the median and 90th percentile are flat to within a
few percent, but the 99th percentile jumps ~35% below 20 deg elevation:

| elevation | median | 90th | 99th |
|-----------|--------|------|------|
| 70-90 deg | 5.68   | 11.18| 18.95|
| 40-55 deg | 5.62   | 11.70| 18.09|
| 20-30 deg | 5.58   | 11.05| 17.34|
| 12-20 deg | 5.80   | 11.72| 25.28|
| 5-12 deg  | 5.89   | 11.38| 24.13|

Near the horizon the same 276 baselines constrain a direction weakly, so the fit
occasionally finds a large spurious amplitude. A single global cut judges rim
peaks against a bar set by well-behaved sky. `ZenithNull` in `significance.py`
normalises by the local *tail* quantile in equal-count elevation bins --
normalising by the 90th percentile corrects nothing, because the 90th percentile
carries none of the effect. Beware small samples: a 99th percentile from 200
draws is 19.0 +/- 2.5 on this data, wide enough to invent the whole effect.

**6. DiSkO's memory is checked before the solve starts.** Peak usage is
`n_vis * n_pix * 48` bytes -- it holds three copies of the operator (complex
`gamma`, the real augmented `concatenate(real, imag)`, and sklearn's lasso copy).
At `res=0.5deg, nvis=10000` over a 170 deg field that is 43 GB; the OOM killer
takes it with **no error message, no output, and exit code 137**, which reads as
an unexplained failure after several minutes of work. `imaging.check_memory()`
refuses to start and shows the arithmetic. The shipped `res=1.0deg` is also the
physically honest choice: the array resolves lambda/B_max = 3.21 deg, so 0.5 deg
cells oversample the beam 6x and are constrained by the lasso prior, not by data.

**7. DiSkO's FITS pixel scale is corrected on read.** The header declares
`CDELT = fov/npix` with `CTYPE = RA---SIN`, but the image is a linear
direction-cosine grid, so the correct value is `(180/π)/(npix/2)` — the header
overstates the scale by 1.485×. See `imaging.py`, and `tests/test_conventions.py`
which locks this down.

---

## Upstream patches

Two changes to installed packages are required. `PATCHES.md` documents them and
`scripts/apply_patches.py` applies them idempotently.

---

## Layout

```
config/            run configuration, one YAML per site
recipes/           Stimela recipe for tart2ms + CASA calibration
src/tart_transient/
    download.py    TART archive
    catalogue.py   sky model, epoch mapping, beam merging
    imaging.py     DiSkO wrapper, WCS correction
    fitting.py     DFT model, joint fit, peeling
    significance.py empirical null, trials correction
    search.py      transient search
    report.py      figures and summary
    cli.py         command line entry point
tests/             convention regression tests
runs/              outputs (gitignored)
```

## Requirements

Python 3.10+, and the TART toolchain (`tart_tools`, `tart2ms`, `disko`,
`python-casacore`, `casatasks`). See `requirements.txt`.

---

## License

MIT -- see `LICENSE`.
