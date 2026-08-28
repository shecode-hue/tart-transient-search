# What actually happens, start to finish

## 0. What TART measures

TART does not take pictures. Each of its 24 antennas records a signal, and the
correlator multiplies every antenna against every other:

    24 antennas  ->  24 x 23 / 2  =  276 antenna pairs (baselines)

For each pair it produces one complex number per second — a **visibility**. Its
amplitude and phase encode how much sky brightness there is at one particular
spatial frequency, set by how far apart those two antennas are and in which
direction.

    one HDF file = 60 seconds = 60 x 276 = 16,560 visibilities

That is the raw data. An image is something you *compute* from it afterwards,
and computing it throws information away. **Everything that matters for
detection is done on the visibilities directly.**

A visibility from a 1-bit correlator is a normalised correlation coefficient in
[-1, 1], roughly `S_source / SEFD`. It carries no absolute flux on its own — we
calibrate it against the GNSS satellites, whose transmit power is specified.

---

## 1. Download

```
TART S3 archive  ->  runs/<name>/data/*.hdf
```

Files are named by the time they *end*. Each holds `vis` (60 x 276 complex),
`antenna_positions`, `baselines`, `gains`, `phases`, `timestamp`, and a `config`
blob with the site's latitude, longitude, frequency and bandwidth.

Selection is by the timestamp in the filename, not by S3 upload time — the
upstream tool uses upload time, which silently returns nothing for archival
data.

---

## 2. Two paths from here

The code does two different things with this data. They answer different
questions.

```
                        runs/<name>/data/*.hdf
                          |                |
        PATH A            |                |          PATH B
   "is a known source     |                |     "is there anything
    doing something?"     |                |      unexpected up there?"
                          v                v
                  tart_transient burst   tart_transient run
                  (minutes)              (~30 min per file)
```

---

## PATH A — the burst command

Works straight on the raw visibilities. No measurement set, no imaging.

### A1. Pick a direction

For a chosen moment, compute where the Sun is from that site (astropy
ephemeris), as a unit vector in the local East/North/Up frame.

### A2. Build the expected signal

For a point source in direction `s`, the visibility on a baseline `b` should be

    M = exp( 2*pi*i * (b . s) / lambda )

This is a prediction, one complex number per baseline. It is what the array
*would* see if there were a source exactly there and nothing else.

### A3. Fit

Project the measured visibilities onto that prediction:

    amplitude = <M, V> / <M, M>

This is a matched filter. It asks: how much of what we measured looks like a
source in that direction? Every baseline votes; noise cancels, signal adds.

### A4. Repeat per second, and for control directions

Do that for every 1-second integration, giving a light curve. Do it in parallel
for eight other sky directions that have no bright source.

### A5. Decide

- A **source** brightens in its own direction only. The controls stay flat.
- A **receiver effect** brightens every direction at once.
- **Local interference** appears at one site, not at four.

Output: `burst_lightcurves.png`, `burst_summary.json`.

The 2025-11-11 test: four sites, all peaking at +3.31 min, controls flat.

---

## PATH B — the full pipeline

Six stages. Only one of them involves an image, and it is not the one that
decides anything.

### B1. HDF -> Measurement Set  (`tart2ms`, the slow step)

Converts to the standard radio astronomy format, and at the same time queries
the TART catalogue for **where every satellite was, at every epoch**, writing
`model_sources_0.txt ... model_sources_N.txt` — one per integration, because
satellites move several degrees in a few minutes.

    runs/<name>/ms/<name>.ms
    runs/<name>/ms/model_sources_*.txt

### B2. Calibration decision

Looks at the gains stored in the HDF. If they are real, use `DATA` as is —
re-solving on a field with no strong isolated calibrator suppresses the sources
you are trying to find. If gains are absent, run CASA `gaincal`/`applycal`.

This is a decision, not a step: on most TART files it correctly does nothing.

### B3. Build the sky model

For each catalogued object, compute its direction per integration and build the
predicted visibility column — the same `exp(2*pi*i (b.s)/lambda)` as A2, one
column per source:

    model matrix:  (n_visibilities, n_sources)

Three corrections happen here:
- objects within one beam of each other are merged, or the fit cannot tell them
  apart and splits flux arbitrarily
- objects the upstream name filter drops (EGNOS, INMARSAT, some BeiDou) are
  added back
- solar-system bodies are re-computed with apparent geocentric positions,
  because upstream writes their barycentric direction — 9.3 degrees wrong for
  the Sun

### B4. Fit and subtract

Solve all source amplitudes at once by least squares:

    V_measured  =  model @ amplitudes  +  residual

Subtract the fitted model. What remains is `RESIDUAL_DATA`: everything the known
sky does not explain.

**This is the step that makes it a transient search.** Known objects are removed
regardless of significance — they have published ephemerides, their existence is
not in question.

### B5. Image  (`DiSkO`)

*Now* an image is made, twice: once from the data, once from the residual.

    fits/01_before.fits   the sky as measured
    fits/02_after.fits    after known sources are removed
    fits/03_removed.fits  the difference

The image is used for **one thing only: generating a list of positions to test.**
It is not used to decide anything. DiSkO's regularised reconstruction correlates
only r ~ 0.55 with a direct transform of the same visibilities and shares just 4
of its 20 brightest positions — the lasso prior moves flux around. An image is a
reconstruction, not a measurement.

### B6. Test each position — back on the visibilities

For every peak the source-finder reports, go back to the raw visibilities and do
the same matched fit as A3. Then four gates:

1. **Is it above the noise?** The threshold is measured, not assumed: fit ~20,000
   random empty directions in this same data and look at the distribution. A
   fixed 5-sigma cut is wrong here — 64-86% of empty sky passes it, because the
   same 276 baselines are re-measured every second and the formal sigma is
   understated by about sqrt(60).

2. **Is it above the noise *for its elevation*?** The null tail is ~35% heavier
   below 20 degrees elevation, so a rim peak is compared against rim noise.

3. **Have we looked in many places?** Testing 40 positions means 40 chances to
   be fooled. The threshold is corrected to the (1 - alpha/N) quantile, obtained
   by fitting the tail rather than reading an order statistic.

4. **Is it a known object?** Compare against the **full** catalogue, not the
   modelled subset. A peak within one beam of a known satellite is that
   satellite.

Output: `results/candidates.csv`, one row per peak with its verdict;
`results/sources.csv`, every modelled source and its fitted brightness;
`results/summary.json`, the whole run.

---

## Where the raw data is used, and where it is not

| stage | works on |
|---|---|
| A2-A5 burst fit | raw visibilities |
| B3 sky model | raw visibilities |
| B4 fit and subtract | raw visibilities |
| B5 imaging | visibilities in, image out |
| B6 candidate positions | **image** |
| B6 significance test | **raw visibilities** |

The image is a pointer, never evidence.

---

## What the answer means

The pipeline reports a candidate count alongside the number expected by chance:

    TRANSIENT CANDIDATES 2
    expected false alarms 0.41  (41 looks at alpha=0.01)

Two candidates against 0.41 expected is not a detection — it is roughly what
chance produces. A candidate is only interesting if it also survives:

- a second site at the same moment (a satellite shifts by parallax, a sky source
  does not)
- a check against external catalogues
- a plausible brightness — TART's floor is ~10^6 Jy, so anything real here is
  either solar or artificial

---

## Known limits

- **The search fits one amplitude per observation.** That is optimal for a steady
  source and worst for a brief one: on the 2025-11-11 burst it diluted a real
  SNR 34 signal to 0.01. `snr_at_windows()` searches time windows as well as
  positions and recovers it; it is available but not yet the default path.
- **TART's floor is ~10^6 Jy.** No known astrophysical transient class except
  solar bursts is that bright. The pipeline is correct; the instrument is small.
- **Images have limited dynamic range.** With ~50 bright satellites and PSF
  sidelobes up to 56%, a real source can be hard to see in an image even when
  the visibility fit finds it clearly.
