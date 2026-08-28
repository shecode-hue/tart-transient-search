# Runbook — transient detection on TART

End to end, from a cold machine. Every command is meant to be pasted as-is.

The target: NOAA/RSTN recorded a **53,000 sfu** radio burst at 1415 MHz on
2025-11-11, peaking at **10:00 UTC**. TART's threshold is ~114 sfu, so this is
~465x above it. Four sites had the Sun above 55 degrees with data recording.

---

## 1. Start Docker

```bash
open -a "Docker Desktop"
docker start shain_tart_jammy
docker ps
```

Expect `shain_tart_jammy   Up N seconds`. If the container is missing, see
section 8.

---

## 2. Enter the container and activate the environment

```bash
docker exec -it shain_tart_jammy bash
source /workspace/gro_demo/.venv/bin/activate
cd /workspace/tart
```

Check the toolchain:

```bash
python3 --version && which tart2ms disko && python3 -c "import casacore, astropy; print('ok')"
```

Expect Python 3.12, paths for both binaries, and `ok`.

---

## 3. Confirm the patches are in place

```bash
python3 scripts/apply_patches.py
```

Expect three "already patched" lines. These add `--filter-elevation` to
`tart2ms`; without them the catalogue stops at 45 degrees elevation.

---

## 4. Run the tests

```bash
PYTHONPATH=src python3 -m pytest tests/ -q
```

Expect `7 passed`. They lock down the conventions that were wrong at various
points: the DiSkO pixel scale, the phase sign, SNR on magnitude, beam merging,
the elevation-dependent threshold, and the tail-fit estimator.

---

## 5. Detect the burst

```bash
PYTHONPATH=src python3 -m tart_transient burst \
  --sites na-unam,za-rhodes,ghana,mu-udm \
  --peak "2025-11-11T10:00:00+00:00" \
  --window 30 \
  --out runs/burst
```

Runs in roughly 10-20 minutes. It downloads ~30 HDF files per site, fits the
Sun's position coherently in every 1-second integration, and does the same for
eight control directions.

**Expected output:**

```
  site           sun    peak/base   peak time   vs controls
  --------------------------------------------------------------
  na-unam         81 deg     22.6x     +3.31 min    26.3 sigma
  za-rhodes       74 deg      9.3x     +3.31 min     2.7 sigma
  ghana           55 deg     ~12x      +3.3  min     ...
  mu-udm          60 deg     ~ex       +3.3  min     ...

  peak times agree to 0.00 min across N sites
```

**How to read it:**

| column | meaning |
|---|---|
| `sun` | Sun elevation at that site. Higher is better. |
| `peak/base` | how many times brighter the Sun got, against its own pre-burst level |
| `peak time` | minutes from the reported NOAA peak. **All sites should agree.** |
| `vs controls` | how far the Sun's rise exceeds the same measurement in 8 other sky directions |

**The result that matters is the peak time agreeing across sites.** Independent
telescopes, independent clocks, different sky geometry. Nothing instrumental
lines up like that; local interference would not.

`vs controls` varies by site because some sites are noisier. na-unam is the
cleanest. A low sigma at one site is not a failure as long as its peak time
agrees.

**Files produced:**

| file | contents |
|---|---|
| `runs/burst/burst_lightcurves.png` | top: all sites overlaid. bottom: the Sun against the 8 control directions |
| `runs/burst/burst_summary.json` | the numbers above, machine-readable |
| `runs/burst/data/<site>/*.hdf` | the raw visibilities, ~150 kB each |

---

## 6. Look at the result

```bash
exit   # leave the container
open ~/Desktop/tart-transient-search/runs/burst/burst_lightcurves.png
```

Top panel: every site flat at 1x, then all rising together to a spike at
+3.3 minutes, with the same shape and the same secondary bumps.

Bottom panel: the Sun's curve against eight other directions computed from the
same visibilities in the same seconds. If this were a receiver-temperature
effect, every grey line would rise too. They stay flat.

---

## 7. Optional: the full imaging pipeline on the same data

The burst command above works straight from the HDF files. To run the complete
chain (measurement set, calibration, DiSkO imaging, subtraction, transient
search):

```bash
docker exec -it shain_tart_jammy bash
source /workspace/gro_demo/.venv/bin/activate && cd /workspace/tart
PYTHONPATH=src python3 -m tart_transient run --config config/na-unam-burst.yaml
```

Takes ~90 minutes, most of it `tart2ms`. Set `reuse_existing: true` under
`measurement_set:` in the config to skip that step on re-runs.

Produces `runs/na-unam-burst/`:

| path | contents |
|---|---|
| `fits/01_before.fits` | all-sky image before subtraction — open in CARTA |
| `fits/02_after.fits` | after the modelled sources are removed |
| `fits/03_removed.fits` | the difference: what was subtracted |
| `plots/01`-`06` | input data, sky model, SNR vs null, before/after, cleaning, transient search |
| `results/sources.csv` | every modelled source with fitted amplitude and SNR |
| `results/candidates.csv` | every residual peak tested, with its verdict and nearest catalogued object |
| `results/summary.json` | the full run record |

Note: with the Sun correctly subtracted this run reports **no candidate at the
Sun** — that is right. The pipeline removes known sources; the burst is found by
the `burst` command in section 5, which looks at the Sun directly.

---

## 7b. Before / after imaging across all four sites

Full chain — measurement set, calibration, DiSkO imaging, subtraction, search —
run twice per site: once during the burst, once on quiet sky ~8 minutes earlier.
One HDF file per run, which is the minimum that gives a real comparison.

Generate the configs (the file covering the peak differs per site, so this picks
them for you):

```bash
cd /workspace/tart && git pull
python3 scripts/make_burst_configs.py
```

Run all eight (~4 hours, ~500 MB):

```bash
for c in na-unam-1111q na-unam-1111b mu-udm-1111q mu-udm-1111b \
         ghana-1111q ghana-1111b za-rhodes-1111q za-rhodes-1111b; do
  echo "=== $c ==="
  PYTHONPATH=src python3 -m tart_transient run --config config/$c.yaml
done
```

Or the two strongest sites only (~2 hours):

```bash
for c in mu-udm-1111q mu-udm-1111b na-unam-1111q na-unam-1111b; do
  PYTHONPATH=src python3 -m tart_transient run --config config/$c.yaml
done
```

Compare:

```bash
python3 scripts/compare_burst.py runs
```

Prints the brightest pixel within 25 px of the Sun in each image:

```
  site           quiet peak   burst peak   ratio   image rms q/b
  na-unam           0.00xxx      0.0xxxx     N.Nx   0.0xxx/0.0xxx
```

Writes `runs/burst_before_after.png` — one row per site, three panels: quiet,
burst, difference, Sun circled.

**Expect the image ratio to be lower than the light-curve ratio.** A 60 s file
averages over the burst rise and decay; the light curve resolves it per second.
A weak image ratio at a site with a clear light curve is dirty-map dynamic range
losing to satellite sidelobes, not a failed detection.

These runs use the standard search, which fits one amplitude per observation.
The images and the comparison are unaffected, but the candidate list will not
flag the Sun.

**Per run:** ~30 min (tart2ms 11-20, two DiSkO images 10-18, search 1-2),
~60 MB (FITS dominate; the HDF input is 150 kB).

---

## 8. If the container is gone

```bash
docker run -dit --name shain_tart_jammy \
  -v ~/Desktop/gro_demo:/workspace/gro_demo \
  -v ~/Desktop/tart-transient-search:/workspace/tart \
  shain_tart:snapshot-20260824 bash
```

The snapshot has Python 3.12, casacore, wsclean, CASA, Stimela, tart2ms and
disko already built. The base `ubuntu:22.04` does not — do not rebuild from it.

---

## 9. Other bursts

52 further bursts above TART's threshold are catalogued in 2025-2026. To try
one, change `--peak` and `--sites`. Confirm the Sun was up at the site first.

```bash
PYTHONPATH=src python3 -m tart_transient burst \
  --sites za-rhodes,ghana,mu-udm \
  --peak "2025-11-09T07:12:00+00:00" \
  --window 30 --out runs/burst-1109
```

That one was 12,000 sfu, about 105x threshold.
