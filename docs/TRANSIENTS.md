# Radio transient classes, and whether TART can see them

TART's measured detection floor is **~10⁶ Jy** in a 1-second integration at
1575.42 MHz (see `RESEARCH.md` §1). The table below asks, for each known class
of radio transient, whether anything it produces reaches that.

**The critical arithmetic.** TART integrates for **1 second**. A burst shorter
than that is diluted:

    effective flux in a 1 s sample  =  peak flux × (duration / 1 s)

A 1-millisecond burst therefore loses a factor of 1000 before it is even
compared to the threshold. This single fact rules out most of the transient sky
for TART regardless of how bright the events are intrinsically.

---

## The table

| class | typical peak at ~1.4 GHz | duration | effective flux in 1 s | TART? |
|---|---|---|---|---|
| **Solar radio burst** | 10⁴ – 10⁸ Jy | seconds to hours | **same** (resolved) | **YES** |
| Quiet Sun | ~5×10⁵ Jy | steady | 5×10⁵ Jy | marginal |
| GNSS satellite (artificial) | 1.2×10⁶ Jy | steady | 1.2×10⁶ Jy | **YES** |
| Magnetar giant radio burst | 1.5×10⁶ Jy | ~1 ms | 1.5×10³ Jy | no |
| Fast radio burst, brightest ever | 1.2×10³ Jy | ~1 ms | ~1.7 Jy | no |
| Fast radio burst, typical | 0.1 – 10 Jy | ~1 ms | 10⁻⁴ – 10⁻² Jy | no |
| Long-period transient (LPT) | up to ~50 Jy | seconds–minutes | ~50 Jy | no |
| Flare star / UV Ceti burst | 0.01 – 1 Jy | minutes–hours | ≤1 Jy | no |
| X-ray binary outburst | 0.1 – 1 Jy | days | ≤1 Jy | no |
| AGN / blazar flare | 1 – 100 Jy | weeks–months | ≤100 Jy | no |
| GRB afterglow | mJy | days–weeks | ≪1 Jy | no |
| Radio supernova | mJy | weeks–years | ≪1 Jy | no |
| Tidal disruption event | mJy | months | ≪1 Jy | no |
| Cassiopeia A *(steady, for scale)* | 2×10³ Jy | — | 2×10³ Jy | no |

**Only one astrophysical class clears the bar: solar radio bursts.** Everything
else is between three and nine orders of magnitude too faint.

For scale: Cassiopeia A, the brightest steady radio source in the entire sky, is
**500× below** TART's detection floor.

---

## Terrestrial and near-Earth signals

These are not astrophysical, but they are real, detectable, and worth listing —
several are already showing up in our data.

| source | at L-band? | TART? | notes |
|---|---|---|---|
| GNSS satellites | yes, by design | **YES** | dozens visible at all times; the dominant signal |
| EGNOS/SDCM/GAGAN augmentation | yes | **YES** | 13 found unmodelled in one file |
| Geostationary comms (INMARSAT etc.) | yes | **YES** | SES-5 caused our first false candidate |
| Aircraft transponders / reflections | partially | possible | 1030/1090 MHz, off-band, but reflections occur |
| **Lightning** | **mostly no** | **unlikely** | see below |
| Meteor reflections | no | no | VHF phenomenon (54–88 MHz), off-band by 20× |
| Ground RFI, power lines | varies | possible | shows up near the horizon |

### On lightning specifically

Lightning is a broadband radio emitter, but its power is overwhelmingly at
**low frequencies** — the sferic spectrum peaks near 5–10 kHz and falls steeply
with frequency. By 1.5 GHz the emission is many orders of magnitude weaker than
at VLF, and what remains is impulsive on microsecond scales.

Two things work against detection with TART: the frequency is far off the
emission peak, and the **1-second integration** dilutes a microsecond impulse by
a factor of ~10⁶.

That said, this is worth an actual test rather than a calculation, because:

- lightning would appear near the **horizon**, where we already know the noise
  tail is ~35% heavier
- it would be **strongly time-localised**, which the window-max search is now
  built to find
- storm times and locations are publicly catalogued, so it can be tested the
  same way the solar burst was: pick a known event, look for it

**Proposed test:** take a TART site during a documented local thunderstorm
(lightning-detection networks publish strike times and coordinates), run the
window-max search restricted to low elevations, and see whether anything
time-localised appears. A null result is publishable as a limit; a detection
would be genuinely new.

---

## Notable events in 2025–2026

| date | event | flux | TART? | why |
|---|---|---|---|---|
| **2025-11-11** | **solar radio burst** | **53,000 sfu = 5.3×10⁸ Jy at 1415 MHz** | **DETECTED** | our result — four sites |
| 2025-11-09 | solar radio burst | 12,000 sfu | detectable | not yet analysed |
| 2025-08-04 | solar radio burst | 6,000 sfu | detectable | no site had Sun up |
| 2025-12-06 | solar radio burst | 1,700 sfu | detectable | nz-elec had Sun up |
| 2026-01-18 | solar radio burst | 1,300 sfu | detectable | not yet analysed |
| 2026-02-01 | X8.1 flare (AR 4366) | radio mainly >2 GHz | **not detected** | tested; off-band |
| 2025-03-16 | **FRB 20250316A "RBFLOAT"** — brightest FRB ever recorded | 1.2 kJy peak, 1.7 kJy·ms | **no** | ~1.7 Jy once diluted to 1 s; 10⁶× too faint |
| 2025–2026 | repeating FRB 20201124A campaign (uGMRT) | Jy-level, ms | no | same dilution |
| 2026-01 | FRB in a binary system, plasma flare | Jy-level | no | same |

**52 further solar bursts above TART's threshold** are already catalogued in
2025–2026 from the NOAA/RSTN 1415 MHz records — see
`tart-candidate-analysis/results/rbr_1415_catalogue.json`. Seven exceed
1,000 sfu.

---

## What this means for the project

TART cannot contribute to the astrophysical transient field: the instrument is
three to nine orders of magnitude short for every class except the Sun. This is
not a defect of the pipeline; it is 24 small antennas with 2.5 MHz of bandwidth.

Where it **can** contribute:

1. **Solar radio burst monitoring.** 52 detectable events already catalogued in
   two years, and the sites span enough longitude to watch the Sun almost
   continuously. This is a real, repeatable programme.
2. **Satellite and RFI characterisation at L-band.** Directly relevant to
   SKA-era interference concerns, and TART sees the full GNSS constellation
   continuously.
3. **Method development.** Everything the pipeline now does — elevation-matched
   nulls, tail-fitted thresholds, time-windowed search, multi-site parallax
   discrimination — is transferable to instruments that *are* sensitive enough.
4. **Possibly lightning**, untested, see above.

## Sources

- [Anderson et al. 2019, OVRO-LWA all-sky transient limits](https://iopscience.iop.org/article/10.3847/1538-4357/ab4f87)
- [Pietka, Fender & Keane 2015, variability timescales and brightness temperatures](https://arxiv.org/abs/1411.1067)
- [Dawes Review 13: the dynamic radio sky](https://www.cambridge.org/core/journals/publications-of-the-astronomical-society-of-australia/article/dawes-review-13-a-new-look-at-the-dynamic-radio-sky/EF458D9211F84397D327610DDCDEEF8E)
- [Bright millisecond burst from Galactic magnetar SGR 1935+2154](https://www.nature.com/articles/s41586-020-2863-y)
- [FRB 20250316A, brightest FRB recorded](https://www.space.com/astronomy/brightest-ever-fast-radio-burst-challenges-assumptions-about-mysterious-blasts-of-energy-this-marks-the-beginning-of-a-new-era)
- [Long Period Transients review](https://arxiv.org/html/2601.10393v2)
- NOAA/SWPC daily solar event reports, USAF RSTN 1415 MHz
