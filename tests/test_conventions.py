"""Regression tests for the two conventions that cost the most to find."""
import numpy as np
import pytest

from tart_transient import fitting
from tart_transient.imaging import disko_cdelt


def test_disko_pixel_scale():
    assert disko_cdelt(2000) == pytest.approx(0.0572958, rel=1e-4)
    header_claims = 170.0 / 2000
    assert header_claims / disko_cdelt(2000) == pytest.approx(1.485, rel=1e-2)


def test_phase_sign_is_positive():
    assert fitting.PHASE_SIGN == +1.0


def test_injected_source_recovers_at_its_own_position():
    rng = np.random.default_rng(0)
    n_row = 2000
    uvw = rng.normal(0, 1.5, size=(n_row, 3))
    uvw[:, 2] *= 0.1
    freqs = np.array([1.57542e9])
    ra0, dec0 = 30.0, -25.0
    true_ra, true_dec = 42.0, -14.0

    tracks = np.zeros((1, n_row, 2))
    tracks[0, :, 0] = true_ra
    tracks[0, :, 1] = true_dec
    model = fitting.model_matrix(uvw, freqs, tracks, ra0, dec0)
    data = model.copy().reshape(n_row, 1, 1) * 1.0

    amp, sigma = fitting.fit_amplitudes(data, model)
    assert np.abs(amp[0]) == pytest.approx(1.0, rel=1e-6)

    off = np.zeros((1, n_row, 2))
    off[0, :, 0] = true_ra + 25.0
    off[0, :, 1] = true_dec + 25.0
    wrong = fitting.model_matrix(uvw, freqs, off, ra0, dec0)
    amp_wrong, _ = fitting.fit_amplitudes(data, wrong)
    assert np.abs(amp_wrong[0]) < 0.5 * np.abs(amp[0])


def test_snr_uses_magnitude_not_real_part():
    amp = np.array([1.0 * np.exp(1j * np.pi)])
    sigma = np.array([0.1])
    assert fitting.snr_of(amp, sigma)[0] == pytest.approx(10.0)


def test_beam_merge_collapses_duplicates():
    import pandas as pd
    from tart_transient.catalogue import merge_within_beam

    cat = pd.DataFrame({
        "name": ["SAT_A_(PRN_1)", "SAT_A_(C60)", "SAT_B"],
        "ra_d": [10.0, 10.0, 60.0],
        "dec_d": [0.0, 0.0, 0.0],
        "flux": [1e5, 1e5, 1e5],
    })
    merged, groups = merge_within_beam(cat, beam_deg=3.4)
    assert len(merged) == 2
    assert "track_name" in merged.columns
    assert any(len(v) == 2 for v in groups.values())


def test_zenith_null_demands_more_of_low_elevation_peaks():
    import numpy as np
    from tart_transient.significance import ZenithNull

    rng = np.random.default_rng(0)
    n = 12000
    zen = np.degrees(np.arccos(1.0 - rng.random(n) * 0.92))
    snr = rng.gamma(2.0, 3.0, n)
    rim = zen > 70.0
    snr[rim] *= 1.35

    zn = ZenithNull(snr, zen)
    scales = zn.describe()["bin_scales"]
    assert scales[-1] > scales[0], "rim scale must exceed zenith scale"

    thr = zn.threshold(44)
    at_zenith = thr * float(zn.scale(10.0))
    at_rim = thr * float(zn.scale(80.0))
    assert at_rim > at_zenith * 1.1, (
        f"low-elevation threshold {at_rim:.1f} should clearly exceed "
        f"the zenith threshold {at_zenith:.1f}")


def test_tail_fit_beats_quantile_on_a_known_exponential_tail():
    import numpy as np
    from tart_transient.significance import tail_threshold, trials_threshold

    rng = np.random.default_rng(0)
    scale, n_looks, alpha = 4.0, 48, 0.01
    n = 3000
    samples = rng.exponential(scale, n)
    truth = -scale * np.log(alpha / n_looks)

    fit = np.mean([tail_threshold(rng.exponential(scale, n), n_looks, alpha,
                                  n_tail=300) for _ in range(30)])
    quant = np.mean([trials_threshold(rng.exponential(scale, n), n_looks, alpha)
                     for _ in range(30)])

    assert abs(fit - truth) < abs(quant - truth), (
        f"tail fit {fit:.2f} should beat quantile {quant:.2f} "
        f"against truth {truth:.2f}")
    assert quant < truth, "the empirical quantile is expected to be biased low"
    assert abs(fit - truth) / truth < 0.15
