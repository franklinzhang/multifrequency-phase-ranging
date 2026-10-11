import math
import sys
from pathlib import Path

import numpy as np
import pytest
from scipy.optimize import brentq

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from mfpr import bounds, core, estimators, lattice, sim  # noqa: E402

F = core.equal_hz()
M = core.moments(F)


def test_equal_grid_numbers():
    assert M["R_db"] == pytest.approx(40.5288, abs=1e-4)
    assert M["sigma_jump"] == pytest.approx(0.083613, abs=1e-6)
    assert M["S"] == pytest.approx(4200e12)
    n = M["n"]
    fbar, B = M["fbar_hz"], 70e6
    assert M["R"] == pytest.approx(1 + 12 * (n - 1) / (n + 1) * (fbar / B) ** 2, rel=1e-12)


def test_matched_grid_shares_moments():
    m = core.moments(core.matched_hz())
    assert m["fbar_hz"] == pytest.approx(M["fbar_hz"])
    assert m["S"] == pytest.approx(M["S"])


def test_one_wavelength_cost_is_four_sigma_jump_squared():
    assert core.phase_cost(F, [M["lambda_m"]])[0] == pytest.approx(4 * M["sigma_jump"] ** 2, rel=1e-9)


def test_two_scale_decomposition_inside_each_tooth():
    d = np.linspace(0, 1.0, 5001)
    envelope = 2 * math.pi * 35e6 * d / core.C
    inside = np.abs(core.wrap(2 * math.pi * d / M["lambda_m"])) < math.pi - envelope
    diff = np.abs(core.phase_cost(F, d) - core.two_scale_cost(F, d))
    assert np.max(diff[inside]) < 1e-9


def test_threshold_snr_closed_form():
    snr_total = M["n"] * core.snr_from_sigma(M["sigma_jump"])
    assert snr_total == pytest.approx((M["R"] - 1) / (2 * math.pi**2), rel=1e-12)


def test_foot_depends_only_on_R():
    x = brentq(lambda x: 8 * x * x * core.q_func(x) - 1 / M["R"], 1.0, 10.0)
    assert x == pytest.approx(4.9083, abs=1e-3)


def test_zzb_returns_to_carrier_line():
    surf = bounds.CostSurface.build(F, 0.5 * core.C / 10e6)
    s = 0.1 * M["sigma_jump"]
    assert core.db(bounds.zzb(surf, s)) == pytest.approx(core.db(core.hmse(F, s)), abs=0.01)
    assert core.db(bounds.zzb(surf, 0.083613)) == pytest.approx(-27.50, abs=0.02)


def test_lattice_limits():
    R = M["R"]
    at1 = lattice.universal_curve([1.0], R)[0]
    assert 10 * math.log10(at1) == pytest.approx(1.14, abs=0.02)
    u = 40.0
    big = lattice.universal_curve([u], R)[0]
    assert big == pytest.approx(1 + (4 / u**2) * (1 - 1 / R) ** 2 / 12, rel=1e-3)
    p = lattice.params(F, 0.5 * M["sigma_jump"])
    assert p["s"] == pytest.approx(0.5 / (2 * math.sqrt(1 - 1 / R)), rel=1e-9)
    assert lattice.ml_mse(p)["slip"] == pytest.approx(2 * core.q_func(2 * math.sqrt(1 - 1 / R)), rel=1e-6)


def test_offsets_split_into_shift_and_disagreement():
    s = 0.5 * M["sigma_jump"]
    tau = 3e-9
    p = lattice.params(F, s, 2 * math.pi * F * tau, w_rx=np.ones(8))
    assert p["beta"] == pytest.approx(0.0, abs=1e-9)
    assert p["fine_bias_m"] == pytest.approx(core.C * tau, rel=1e-9)
    eps = 0.3
    p = lattice.params(F, s, np.full(8, eps), w_rx=np.ones(8))
    assert p["beta"] == pytest.approx(-eps / (2 * math.pi), rel=2e-3)


def test_coherent_ml_recovers_range_at_high_snr():
    rng = np.random.default_rng(0)
    d0 = 12.345 + rng.random(50) * 0.1
    phi = sim.phases(rng, F, d0, 1e-3, 50)
    grid = estimators.make_grid(0, 25, M["lambda_m"])
    err = estimators.coherent_ml(phi, F, grid) - d0
    assert np.max(np.abs(err)) < 1e-4


def test_valley_mixture_matches_union_when_slips_are_rare_and_stays_bounded():
    surf = bounds.CostSurface.build(F, 15.0, step=1e-3)
    small = 0.3 * M["sigma_jump"]
    union = bounds.valley_sum_mse(surf, small)[0]
    assert bounds.valley_mixture_mse(surf, small)[0] == pytest.approx(union, rel=0.01)
    huge = 50.0
    assert bounds.valley_mixture_mse(surf, huge)[0] < 15.0**2
