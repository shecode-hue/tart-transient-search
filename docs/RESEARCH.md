# Research notes: where this pipeline stands, and how to improve it

Literature review and measurements, 2026-08-25. Every number attributed to
"measured" comes from our own two runs (za-hammanskraal and za-rhodes,
2026-08-24); everything else is cited.

---

## 1. The headline result: what TART can actually detect

This is the most important thing I found, and it reframes the project.

**Measured detection threshold: ~1.1 x 10^6 Jy (about 114 sfu).**

> **Correction, 2026-08-25.** An earlier version of this section gave 11,000 Jy
> and concluded the instrument was mysteriously insensitive. That was my
> arithmetic error, not the telescope's. I divided the GNSS specified received
> *power* (-158.5 dBW) by bandwidth and called the result a flux density,
> omitting the antenna's effective area. The correct conversion is
> S = P / (A_e x B) with A_e = 57.5 cm^2 for a 3 dBi antenna at 19 cm, giving
> **1.2 x 10^6 Jy** per GNSS satellite -- 174x larger. Every threshold below is
> revised accordingly, and the apparent "1400x suppression" of the Sun
> disappears entirely.

Derived two independent ways, which agree:

*Radiometer equation.* For a GNSS patch antenna (gain 3–5 dBi, so effective
area 58–91 cm² at 19 cm) and T_sys 200–300 K, SEFD = 2kT/A_e = 6×10⁷ to
1.4×10⁸ Jy per antenna. With N = 24, Δν = 2.5 MHz, τ = 60 s:

    sigma = SEFD / sqrt(N(N-1) · Δν · τ)  =  211 to 501 Jy

*Empirical, from the data.* GPS L1 C/A arrives at Earth at ≥ −158.5 dBW/m²
spread over 2.046 MHz, i.e. **6.9 × 10³ Jy** in TART's band. Scaling our fitted
satellite amplitudes to that gives a trials-corrected threshold of:

| site | threshold |
|---|---|
| za-hammanskraal | ~1.9 x 10^6 Jy |
| za-rhodes | ~2.1 x 10^6 Jy |
| nz-elec (300 integrations) | **1.1 x 10^6 Jy = 114 sfu** |

The instrument is a **1-bit (two-level) correlator** -- confirmed in the source:
`corr_b()` XORs bit streams and `van_vleck_correction(R) = sin(pi/2 R)` is the
two-level correction. Visibilities are therefore normalised correlation
coefficients, and a coefficient is S/SEFD. The radiometer equation predicts
0.013-0.020 for a GNSS satellite; the measured amplitude is 0.056. Same order.
**The instrument is behaving normally.**

Two different arrays, two different amplitude scales, thresholds agreeing to
10%. That is a real flux calibration, and it is the number the project has been
missing.

**What that means, at 1.575 GHz:**

| source | flux | TART? |
|---|---|---|
| large solar burst (~1000 sfu) | 10⁷ Jy | **yes — 9x above threshold** |
| 2025 X5.1 burst (53,000 sfu at 1415 MHz) | 5.3×10⁸ Jy | **yes — 463x above threshold** |
| quiet Sun (~50 sfu, half rejected by RHCP) | 2.5×10⁵ Jy | marginal — measured 6.2×10⁴ Jy |
| GNSS L1 satellite | 6.9×10³ Jy | marginal (single-look yes, trials-corrected no) |
| **Cassiopeia A — brightest steady radio source in the sky** | 2×10³ Jy | **no** |
| bright AGN flare | ~100 Jy | no |
| X-ray binary outburst | ~1 Jy | no |
| flare star | ~0.1 Jy | no |

1 sfu = 10⁴ Jy ([NOAA/standard definition](https://iopscience.iop.org/article/10.3847/1538-4357/ac34ed)).

**TART cannot detect any known class of astrophysical radio transient except
solar bursts.** Cas A, the brightest steady radio source in the entire sky, is
five times *below* our detection threshold. This is not a defect in our
pipeline — it is the instrument. Twenty-four patch antennas with 2.5 MHz of
bandwidth cannot reach the mJy–Jy regime where transients live.

This also resolves a puzzle from the runs: "0/46 sources above threshold" is
correct behaviour, not a bug. The satellites really are at 6.9 kJy and the
trials-corrected threshold really is at 11 kJy.


### 1.1 The solar burst test — a measured negative

The paper-based reasoning above says solar bursts should be TART's one
accessible astrophysical target. **That prediction was tested against a real
flare and it failed.**

The X8.1 flare of 2026-02-01 (Active Region 4366) is the strongest solar event
since 2024. SWPC puts the peak at **23:57 UTC**. At that moment the Sun was
59 deg above the horizon at nz-elec (Dunedin) and the archive holds continuous
1 s data through the whole event. This is as favourable a test as one could ask
for.

Coherent fit at the Sun's position, per 1 s integration, 24 files spanning
-18 to +18 minutes around the peak (1440 integrations):

| minutes from peak | median \|amp\| | median SNR | sun elev |
|---|---|---|---|
| -16 to -12 | 0.01685 | 1.50 | 57.5 |
| -4 to 0 | 0.00620 | 0.54 | 58.7 |
| **0 to +2** | **0.00484** | **0.44** | 59.0 |
| +4 to +8 | 0.00673 | 0.60 | 59.3 |
| +12 to +20 | 0.00193 | 0.17 | 59.9 |

**No enhancement at the flare peak.**

The reason is now clear, and it is not instrumental blindness. TART's threshold
is ~114 sfu at 1.575 GHz. SWPC reports that this flare's radio noise "affected
mainly radio frequencies above 2 GHz" -- so at 1.575 GHz the burst simply did
not deliver 114 sfu. A non-detection is the correct outcome.

An independent check supports this. Because a 1-bit correlator measures
S/T_sys, a broadband power rise would dim *every* satellite at once. Measuring
the common-mode amplitude across all 44 modelled satellites gives a change of
**+0.81% (0.3 sigma)** at the flare peak -- no system-temperature excursion,
limiting the broadband solar flux to roughly 450 sfu by a route that never
touches an image. See `scripts/16_common_mode_radiometer.py`.

Three reasons this is not surprising in hindsight:

1. **The burst was in the wrong band.** SWPC reports the radio noise from this
   flare "affected mainly radio frequencies above 2 GHz". TART observes at
   1.575 GHz.
2. **The front-ends are GNSS receivers.** They are right-hand circularly
   polarised, so roughly half of the (unpolarised) solar emission is rejected
   before anything else happens, and only 2.5 MHz of bandwidth is sampled.
3. **Automatic gain control and coarse quantisation.** A GNSS front-end holds
   its output level constant. A broadband rise in system temperature — which is
   exactly what a solar burst is — gets compensated by the AGC rather than
   recorded as signal, and few-bit quantisation compresses strong inputs
   further. **This is a hypothesis, not something I have verified**, but it
   would explain why the Sun is ~350x weaker than its flux density implies.

**Consequence:** the sfu-to-TART conversion cannot be assumed. The honest
statement is that TART's response to broadband unpolarised emission is
*unmeasured*, and the one measurement available (this flare) is a
non-detection. Establishing that response — by observing the Sun across many
elevations and comparing to RSTN's 1415 MHz flux — would be a genuinely useful
piece of work, and is a prerequisite for any claim about solar sensitivity.

Data: `tart-candidate-analysis/scripts/14_solar_flare_lightcurve.py`,
`plots/solar_x81_lightcurve.png`.
### What this implies for the thesis

The honest framings, in order of strength:

1. **A methods and verification study.** The pipeline demonstrably finds nothing
   and can prove *why* nothing it finds is real — including a two-site parallax
   test that discriminates sky sources from satellites geometrically. That is a
   genuine contribution, and the toolchain bugs found along the way are real
   results.
2. **Solar radio burst monitoring — viable, and the strongest science case.**
   The threshold is ~114 sfu. The 2025 X5.1 burst reached 53,000 sfu at
   1415 MHz, 463x above it. Bursts above ~1000 sfu are routine near solar
   maximum. TART sites span many longitudes, so between them they watch the Sun
   almost continuously. Two independent detection routes now exist: imaging,
   and the common-mode radiometer of §1.1 which needs no imaging at all.
3. **Satellite and RFI characterisation.** TART sees the full GNSS constellation
   continuously in the L1 band. Characterising satellite contamination — which
   we have already had to do — is publishable in its own right and directly
   relevant to SKA-era RFI concerns ([Dawes Review 13](https://www.cambridge.org/core/journals/publications-of-the-astronomical-society-of-australia/article/dawes-review-13-a-new-look-at-the-dynamic-radio-sky/EF458D9211F84397D327610DDCDEEF8E)
   flags satellite reflections as an emerging challenge).
4. **A transient surface-density limit.** Weak, but real and computable — see §3.

What is *not* honest is presenting this as a search likely to find an
astrophysical transient.


---

## 1.2 FIRST DETECTION — a solar radio burst, confirmed at three sites

**2026-08-25. TART detected a solar radio burst.** This is the project's first
genuine astrophysical detection, and it validates the sensitivity estimate of
§1 end to end.

### How the target was chosen

Selecting on X-ray class failed (§1.1) because soft X-ray output is only loosely
coupled to microwave emission. So instead, select on the *measured radio flux at
almost exactly TART's frequency*.

NOAA/SWPC daily event reports (`YYYYMMDDevents.txt`) list `RBR` entries -- radio
bursts at fixed frequency -- measured by the USAF Radio Solar Telescope Network
at 245, 410, 610, **1415**, 2695, 4995, 8800 and 15400 MHz. 1415 MHz is 160 MHz
from TART's 1575.42 MHz.

Parsing 322 daily reports from 2025-2026:

| | count |
|---|---|
| bursts recorded at 1415 MHz | 65 |
| above TART's ~114 sfu threshold | **52** |
| above 1000 sfu | 7 |

Cross-matching the strongest against archive coverage and Sun elevation gave the
target: **2025-11-11 10:00 UTC, 53,000 sfu — 465x above threshold — with the Sun
above 55 deg at four TART sites simultaneously.**

### The detection

Coherent fit at the Sun's position, per 1 s integration, 2040 integrations per
site across a 30-minute window, straight from the HDF files (no measurement set
required):

| site | sun elev | baseline \|amp\| | peak \|amp\| | ratio | peak time |
|---|---|---|---|---|---|
| na-unam | 81 deg | 0.00292 | 0.07400 | **25.3x** | **+3.30 min** |
| za-rhodes | 74 deg | 0.00803 | 0.07484 | **9.3x** | **+3.30 min** |
| ghana | 55 deg | 0.00515 | 0.06402 | **12.4x** | **+3.30 min** |

(30 s smoothed; unsmoothed single-integration peaks reach 46x at na-unam.)

**Three telescopes on three different parts of the continent, with independent
hardware, independent clocks and completely different sky geometry, peak at the
same instant to within 0.01 minutes.**

Cross-correlation of the light curves:

| pair | r | lag |
|---|---|---|
| na-unam vs ghana | **+0.978** | 0.00 min |
| za-rhodes vs ghana | +0.840 | +0.02 min |
| na-unam vs za-rhodes | +0.819 | -0.02 min |

Nothing instrumental produces that. Not RFI -- it would not correlate across
three continents. Not a satellite -- one cannot sit at the Sun's position from
three sites at once (the parallax argument of the candidate analysis, run in
reverse). It is the Sun.

### Flux, and an important caveat

Using the GNSS flux scale (a satellite is 1.2e6 Jy at amplitude 0.056):

| site | peak flux |
|---|---|
| na-unam | 1.59e6 Jy = **159 sfu** |
| za-rhodes | 1.60e6 Jy = **160 sfu** |
| ghana | 1.37e6 Jy = **137 sfu** |

Three independent telescopes agreeing to within 15% is a strong cross-check on
the flux scale derived in §1.

They do *not* agree with RSTN's 53,000 sfu at 1415 MHz, and that gap needs
stating honestly rather than explaining away. Contributing factors, in likely
order:

1. **1-bit compression.** A two-level correlator measures S/(SEFD + S), which
   saturates toward 1 for a very bright source. The numbers above are therefore
   **lower bounds**, not linear flux measurements. Recovering true flux needs
   the van Vleck correction (`sin(pi/2 R)`, already in `tart.imaging.correlator`)
   applied with a proper system-temperature model.
2. **Polarisation.** The antennas are RHCP; solar burst emission is often
   strongly circularly polarised, so the factor is somewhere between 1 and 2 and
   depends on the burst's handedness.
3. **Spectrum and timing.** 1415 and 1575 MHz are not the same frequency, and
   our peak is 3.3 minutes after the RSTN maximum -- so these are not
   measurements of the same instant.

### What this establishes

- TART **can** detect astrophysical transients of the one class within its
  sensitivity: solar radio bursts.
- The §1 threshold of ~114 sfu is validated -- a 465x event is detected at
  9-25x above its own local baseline, which is the right order.
- 52 further bursts above threshold are already catalogued in 2025-2026 alone
  (`results/rbr_1415_catalogue.json`), so this is a repeatable programme, not a
  one-off.
- Multi-site coincidence works as a discriminator in both directions: it refuted
  a satellite masquerading as a transient, and here it confirms a real source.

Scripts: `17_solar_burst_hunt.py`, `18_verify_burst_coincidence.py`,
`19_rbr_catalogue.py`. Figure: `plots/solar_burst_20251111_detection.png`.


### 1.3 End-to-end validation: the pipeline recovers the burst unaided

§1.2 detected the burst with a direct coherent fit read from the HDF files. That
tests the *estimator*. This tests the *pipeline*: the full chain, run on
na-unam for 2025-11-11 10:00:40 - 10:05:46 UTC (5 files, 300 integrations), with
the Sun absent from the model as usual.

    THE SUN WAS FOUND: 0.55 deg from ephemeris, SNR 58.16, rank 2 of 41 peaks
      veto identified it as: 'sun' at 0.55 deg
      passed single-look: True     confirmed as transient: False

Every stage worked without being told the Sun was there:

| requirement | outcome |
|---|---|
| residual peak at the Sun | **0.55 deg** from ephemeris (beam 3.38 deg, so ~beam/6) |
| strong enough to matter | SNR 58.16, **rank 2 of 41** |
| survives the null | passed single-look |
| correctly classified | veto named it **`sun`**, so NOT reported as a transient |

**This is the natural injection test the project needed.** A real astrophysical
transient, of flux measured independently by RSTN, at a position and time fixed
by ephemeris, recovered end-to-end by a pipeline that had no knowledge of it.
Synthetic injection can only test the code against its own assumptions; this
tests it against the sky.

The astrometry deserves note: 0.55 deg on a 3.38 deg beam, from a solve grid
whose cells are ~0.9 deg. That is better than the grid spacing, which is what a
coherent fit should achieve and a useful independent check that the coordinate
conventions of §2 (phase sign, pixel scale) are right.

**What did not go to plan, recorded honestly:**

- I predicted the Sun would be the *strongest* peak. It is rank 2. A peak at
  SNR 61.79, 59.5 deg away and 4.59 deg from BEIDOU-3 M4, is stronger.
- **Sidelobes behaved as predicted.** Three further peaks at 13.5, 14.5 and
  16.8 deg from the Sun are attributed by the veto to `sun` -- they are its
  sidelobes. One of them (SNR 43.39 at 14.51 deg) passes single-look and, being
  beyond one beam, is counted as "unexplained". An unmodelled bright source
  manufactures its own false candidates, which is a direct argument for
  subtracting the Sun rather than merely vetoing it.
- 2 candidates were confirmed against 0.41 expected. Both sit 4.6-5.6 deg from
  catalogued satellites -- just outside the 3.38 deg beam -- so the most likely
  explanation is mis-modelled satellite flux rather than anything new.

**Consequences for the pipeline:**

1. Model the Sun and the other name-filtered objects, do not merely veto them.
   The veto prevents a false *claim*; it does not prevent the false *peaks* a
   bright unmodelled source scatters across the field.
2. Consider widening the veto radius beyond one beam for very bright objects,
   since sidelobe peaks land well outside it.


### 1.4 The estimator was the problem — a transient search that integrates transients away

Correcting the Sun's position (section 1.3) removed its false peak: after the fix
**no peak lies within one beam of the Sun**, power removed rose 20.8% -> 27.7%,
and image RMS improved 0.0005 -> 0.0004. Good -- but the burst then vanished
from the candidate list entirely, which exposed something more fundamental.

**The search fits ONE amplitude across the whole observation.** Measured at the
Sun's position in the residual, na-unam, 300 integrations:

| window | integrations | residual SNR |
|---|---|---|
| **whole observation (what the pipeline used)** | 300 | **0.01** |
| +/- 90 s around the burst | 175 | 12.67 |
| +/- 30 s | 59 | 28.11 |
| **+/- 15 s** | 29 | **34.43** |
| before the burst (control) | 38 | 10.04 |

A constant amplitude is the *optimal* estimator for a steady source and the
*worst* one for a brief event: it absorbs the average and integrates the excess
to nothing. A transient search built this way is structurally unable to find
short transients. This is not a tuning problem.

**Fix 1, implemented: `search.snr_at_windows()`.** Dyadic top-hat windows
(1, 2, 4, 8 splits), taking the best -- the template-bank idea of Feng et al.
(2017) in its simplest form. At the Sun this recovers SNR 0.01 -> **24.98**
(window 4 of 8, exactly the burst), a factor of 3665.

`significance.null_with_zenith(times=..., splits=...)` draws the null the same
way, because maximising over windows inflates the trials factor and a threshold
built for a different test is not a threshold.

**Honest limit of fix 1.** The windowed SNR at the Sun (24.98) does NOT exceed
the windowed control maximum (33.30) in this data. Windowing recovers the
signal but does not, by itself, make this burst stand out in the residual --
because once the Sun is correctly subtracted, only its *deviation* from the
fitted mean remains, competing with imperfectly subtracted satellites.

**Fix 2, the right architecture (not yet implemented).** Searching the residual
is the wrong place to look for a *known* source that varies. Fit per-window
amplitudes for every modelled source and test each for variability. The Sun's
amplitude in ten 30 s windows:

    0.027  0.055  0.047  0.050  0.119  0.082  0.048  0.038  0.025  0.028
                                 ^ the burst

Unambiguous. But ranked by raw max/median the Sun places 3rd, behind sources
whose excursions are tiny in absolute terms (0.017 -> 0.003). The metric must be
significance-weighted -- (max - median) / sigma -- not a bare ratio.

**Summary of what a transient-capable pipeline needs here:**

1. Time-resolved search, not one amplitude per observation (fix 1, done).
2. A null drawn with the same estimator (done).
3. Variability testing of modelled sources, significance-weighted (to do).
4. Correct source positions, or the subtraction manufactures its own
   candidates (section 1.3, done).

Scripts: `22_windowed_search_test.py`, and the variability probe in this
section.

---

## 2. Where our methods sit relative to the literature

Some good news: several of our hard-won choices are what the field does.

**We were right to detect in the visibility domain.** We measured DiSkO's lasso
image correlating only r ≈ 0.55 with a direct transform and sharing 4 of its 20
brightest positions. The compressed-sensing literature explains why: L1
minimisation is a biased surrogate for L0, and the bias is well documented —
recent work develops non-convex (SCAD) penalties specifically to fix it
([Non-convex sparse regularisation, A&A 2025](https://www.aanda.org/articles/aa/full_html/2025/12/aa55737-25/aa55737-25.html)).
A regularised image is a *reconstruction*, not a measurement. Using it only to
generate candidate positions and testing each by coherent fit is correct.

**Our empirical null is standard practice.** Feng et al. use a "playground
region" of ~10% of the image assumed to contain no transients, to characterise
the background distribution
([Matched filter technique, AJ 153, 98](https://arxiv.org/abs/1701.03557)).
Our 20,000 random-position null is the same idea.

**Our trials correction is standard.** The field computes false-alarm rates
from the number of independent *beams*, not pixels, because interferometer
noise is correlated on the beam scale.

**Our wide-field geometry is right.** At 170° FOV the 2D approximation fails
badly; the standard fixes are faceting, 3D FFT, w-projection or w-stacking
([WSClean, Offringa et al.](https://arxiv.org/pdf/1407.1943)). Both DiSkO and
our own `fitting.py` evaluate exp(2πi(ul + vm + w(n−1))) directly, which is the
exact expression — no w-term approximation is being made anywhere in our chain.

---

## 3. Seven concrete improvements, in priority order

### 3.1 Fit the tail, don't take a percentile (highest value, easiest)

**Problem we hit:** the trials-corrected threshold needs the (1 − α/N) quantile.
With 48 looks at α = 0.01 that is the 99.98th percentile, which needs ≳ 4,400
samples *just to have one draw beyond it*. We threw 20,000 draws at it and the
estimate is still noisy — a 99th percentile from 200 draws is 19.0 ± 2.5 on our
data.

**The fix**, from Feng et al. (2017): fit the *tail* of the background
distribution to an exponential and extrapolate:

    N(≥ρ) = N̂ exp(−ρ/ρ̂)        →        ρ* = ρ̂ (log N̂ − log P_FA)

They fit 500 tail points to reach P_FA = 10⁻³. This converts an
order-statistic problem into a two-parameter fit, so the threshold stops being
limited by how many draws we can afford.

**Effort:** ~30 lines in `significance.py`. Should be done first.

### 3.2 Calibrate in janskys, using the GNSS satellites as flux references

We just did this by hand and it worked — two sites agreeing to 10%. It should
be in the pipeline.

GNSS satellites are *ideal* flux calibrators for TART: their transmit power is
specified by treaty, they are always present, there are dozens above the horizon
at any time, and their positions are known to metres. Every TART observation
carries its own flux calibration.

**What it enables:** results in Jy instead of dimensionless SNR; comparison with
other surveys; injection–recovery in physical units; a real detection limit.

**Effort:** moderate. Needs a per-satellite EIRP table and an antenna gain model
versus elevation.

### 3.3 Injection–recovery, in janskys

Still the biggest gap, and now it can be done properly because §3.2 gives units.

Method, following standard practice
([DSA-110](https://arxiv.org/abs/2510.18136), and the completeness definition in
[Carbone et al. 2016](https://academic.oup.com/mnras/article/459/3/3161/2595112)):
inject synthetic point sources of known flux into the *visibilities* at random
positions, run the full chain, and count recoveries in bins of flux and
elevation. Output is a completeness curve C(S, elevation), and the 90%
completeness flux is the survey's quotable limit.

Note the elevation dependence matters for us specifically — we measured the null
tail running ~35% hotter below 20° elevation, so completeness will be strongly
elevation-dependent. OVRO–LWA model their beam as gain ∝ sin^1.6(elevation)
([Anderson et al. 2019](https://iopscience.iop.org/article/10.3847/1538-4357/ab4f87)).

### 3.4 Publishable rate limits

With zero detections, Carbone et al. (2016) give the 95% confidence limit:

    rho < −ln(0.05) / Omega_tot

TART's instantaneous field of view is huge — that is its one real advantage:

    radius 85° → 5.736 sr = 18,829 deg² = 46% of the entire sky in one snapshot

| exposure | Omega_tot | limit |
|---|---|---|
| one 60 s snapshot | 18,829 deg² | 1.6×10⁻⁴ deg⁻² |
| our two runs | 37,658 deg² | 8.0×10⁻⁵ deg⁻² |
| one day at 60 s cadence | 2.7×10⁷ deg² | 1.1×10⁻⁷ deg⁻² |
| one month | 8.1×10⁸ deg² | 3.7×10⁻⁹ deg⁻² |

For context: OVRO–LWA reach 2.5×10⁻⁸ deg⁻² at 10.5 Jy; LOFAR 1.3×10⁻³ deg⁻² at
0.3 Jy. Our limits are set at ~11 kJy, so they are weak in flux but competitive
in *area* — and there is very little all-sky transient work at 1.575 GHz.

Carbone et al. also give a power-law method that uses every image rather than
discarding the noisy ones, worth a factor ~4:

    N* < −ln(0.05)/Omega · (S*/D)^(−gamma) · 1/Sum_i sigma_i^(−gamma)

**Effort:** small — this is arithmetic over metadata we already store.

### 3.5 Multi-epoch matched filter

Right now each file is searched independently. Feng et al. build a time series
per pixel across epochs and correlate against light-curve templates:

    rho = Sum_i (b_i²/sigma_i²) x_i (f_i − <f>),    rho~ = max over templates of rho/sigma_rho

Maximising over templates *automatically incorporates the trials factor* into
the background distribution. They report this beats source-finding on the same
images and detects below the classical confusion limit, because confusion
sources are time-independent.

For TART this is a natural fit: 60 integrations per file and hundreds of files
per day, and the dominant contaminants (satellites) move predictably while a
transient would not.

### 3.6 Closure quantities as a calibration-independent cross-check

Given how much time we have lost to calibration, this is attractive.

The bispectrum b_ijk = a_ij·a_jk·a_ki·exp(i(φ_ij + φ_jk + φ_ki)) is **independent
of per-antenna phase errors** — they cancel around the loop
([Law & Bower 2012](https://iopscience.iop.org/article/10.1088/0004-637X/749/2/143)).
With 24 antennas TART has n_tr = 24·23·22/6 = **2,024 triples**.

Sensitivity scaling: bispectrum S/N ∝ s³·n_tr versus coherent S/N ∝ s·sqrt(n_bl).
For the VLA they measure a 2.2× sensitivity penalty (38 mJy vs 17 mJy). The
penalty is real but the payoff is a detection statistic that *cannot* be faked by
a gain error — exactly the failure mode we have been chasing.

**Use it as a veto**, not the primary detector: any candidate that survives the
coherent fit should also appear in the bispectrum.

### 3.7 Model the Sun, and everything else the name filter drops

`tart2ms` keeps only `^GPS|^QZS|^BEIDOU|^GSAT`. At our epoch that dropped 6 of
54 objects above the horizon, including the Sun and four EGNOS/SDCM payloads
that transmit *in TART's band*.

I measured what the Sun is actually doing, rather than assuming: SNR 3.96
(Hammanskraal) and 8.83 (Rhodes) in the uncalibrated DATA — much weaker than its
~5×10⁵ Jy would suggest, presumably because the patch antennas are RHCP and
strongly attenuated toward a 21°-elevation source. So it is *not* currently the
dominant contaminant. But it should be modelled anyway, because at higher
elevation or during a burst it would swamp everything.

Our `catalogue.full_sky()` fix already vetoes against the complete list. The
remaining work is to *subtract* them, not just veto.

---

## 4. Two things I could not settle

**TART's SEFD is not published anywhere I could find.** The Otago group's
papers (Scheel, Molteno & Brown, ENZCon 2016; the continuous-calibration paper,
doi:10.1109/ICEAA.2019.8879242) cover calibration and aperture synthesis, and
there is a Stellenbosch MSc characterising TART (N. Mtetho) whose full text I
could not retrieve. Our radiometer estimate and the GNSS-based empirical scale
agree, which is reassuring, but a published SEFD would settle it. **Worth asking
Dr. Hugo whether SARAO has measured this, or emailing Tim Molteno directly.**

**Whether DiSkO's lasso hyperparameters are appropriate for detection.** We
inherited alpha = 0.006, l1_ratio = 0.02 from the reference recipe. Nobody has
checked what those do to point-source photometry at TART's resolution. Since we
no longer use the image for detection this is not urgent, but if the image is
ever used quantitatively it needs testing.

---

## 5. Reading list

Ordered by how directly useful each is to this project.

1. Feng et al. 2017, *A Matched Filter Technique for Slow Radio Transient
   Detection and First Demonstration with the MWA*, AJ 153, 98 —
   [arXiv:1701.03557](https://arxiv.org/abs/1701.03557). **Read first.**
   Detection statistic, exponential tail fit, playground region.
2. Carbone et al. 2016, *New methods to constrain the radio transient rate*,
   MNRAS 459, 3161 — [arXiv:1411.7928](https://arxiv.org/abs/1411.7928).
   How to turn a null result into a limit.
3. Anderson et al. 2019, *New Limits on the Low-frequency Radio Transient Sky
   Using 31 hr of All-sky Data with the OVRO–LWA*, ApJ 886, 123 —
   [IOP](https://iopscience.iop.org/article/10.3847/1538-4357/ab4f87).
   The closest analogue to TART: all-sky, zenith-pointing, satellite/RFI-limited.
4. Law & Bower 2012, *All Transients, All the Time*, ApJ 749, 143 —
   [IOP](https://iopscience.iop.org/article/10.1088/0004-637X/749/2/143).
   Closure quantities, calibration-independent detection.
5. *The Dawes Review 13: A new look at the dynamic radio sky*, PASA —
   [Cambridge](https://www.cambridge.org/core/journals/publications-of-the-astronomical-society-of-australia/article/dawes-review-13-a-new-look-at-the-dynamic-radio-sky/EF458D9211F84397D327610DDCDEEF8E).
   Current state of the field; transient classes and rates.
6. Pietka, Fender & Keane 2015, *The variability timescales and brightness
   temperatures of radio flares*, MNRAS 446, 3687 —
   [arXiv:1411.1067](https://arxiv.org/abs/1411.1067). The transient phase
   space; what is bright enough to matter.
7. Offringa et al. 2014, *WSClean* — [arXiv:1407.1943](https://arxiv.org/abs/1407.1943).
   Wide-field imaging and w-term handling.
8. Scheel, Molteno & Brown 2016, *Transient array radio telescope: Calibration
   and aperture synthesis*, ENZCon. TART's own calibration description.

---

## 6. What I would do next, concretely

In order:

1. **Exponential tail fit** for thresholds (§3.1) — half a day, removes a real
   statistical weakness.
2. **Jy calibration** from GNSS satellites (§3.2) — gives every future result
   physical units.
3. **Injection–recovery** in Jy, binned by elevation (§3.3) — produces the
   completeness curve, which is the thesis's missing figure.
4. **Rate limits** via Carbone (§3.4) — converts "we found nothing" into a
   number that can go in a paper.
5. Then reconsider the science framing with supervisors, armed with the fact
   that the instrument's floor is ~11 kJy.
