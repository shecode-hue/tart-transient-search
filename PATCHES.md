# Upstream patches

Two changes to installed packages are required. `scripts/apply_patches.py`
applies both idempotently and is safe to re-run.

Both are defects in the standard toolchain, not in this pipeline. They are
worth reporting to the TART maintainers.

---

## 1. `tart2ms` — expose the satellite elevation cut

`tart2ms` filters the GNSS catalogue at a hardcoded **45° elevation**, with no
CLI or recipe parameter. A TART image spans **170° FOV** — down to ~5°.

Measured at za-hammanskraal, lat −25.24 lon 28.32, 2025-10-09T23:59Z:

| elevation cut | satellites |
|---|---|
| ≥ 45° | **14** ← what the catalogue contained |
| ≥ 10° | 60 |
| ≥ 0° | **67** ← what is actually in the image |

So ~80% of real emitters had no catalogue entry: they cannot be subtracted,
they inflate the empirical noise floor, and they appear as unexplained residual
peaks. Dropping the cut to 5° took visibility power removed from 6.95% to 30.7%.

The catalogue API **ignores its own `elevation` query parameter** and returns
all 67 objects; the cut is applied client-side at `tart2ms.py:1422`.

The patch adds `--filter-elevation` (default 45.0, so nothing changes unless
asked) and threads it through
`ms_from_hdf5`/`ms_from_json` → `ms_create` → `predict_model`, plus
`__fetch_sources`.

### Also worth knowing (not patched)

`tart2ms` additionally filters **by name**:

```python
filter_name = r"(?:^GPS.*)|(?:^QZS.*)|(?:^BEIDOU.*)|(?:^GSAT.*)"
```

Six objects above the horizon in the test field are excluded by this: INMARSAT
4-F2, SES-5, LUCH 5B, LUCH 5V, ASTRA 5B, EUTELSAT 5 West B. They are
geostationary communications satellites and are neither modelled nor subtracted.

---

## 2. `tartcargo` — accept the new parameter in the cab schema

Stimela validates recipe parameters against the cab schema, so
`filter-elevation` must be declared in `tartcargo/tart2ms.yml` for a recipe to
pass it.

---

## Not patched — worked around instead

**DiSkO FITS pixel scale.** The header declares `CTYPE = RA---SIN` with
`CDELT = fov/npix`, but the image is a linear direction-cosine grid, so the
correct value is `(180/π)/(npix/2)`. The header overstates the scale by
**1.485×**. Corrected on read in `imaging.load_wcs()` rather than patching
DiSkO, since the official pipeline never writes FITS (it keeps the `.sphere`
and renders with `disko-draw`) — this sits on a code path upstream does not
exercise. `tests/test_conventions.py` locks the correction down.

**CASA `plotms` in the reference recipe.** Four cosmetic plotting steps sit
before calibration and fail without a display, aborting the recipe before it
ever calibrates. This pipeline calls `tart2ms` and `casatasks` directly instead
(see `measurement_set.py`).
