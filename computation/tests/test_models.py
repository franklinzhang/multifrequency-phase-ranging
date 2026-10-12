import math
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from mfpr import bounds, core, design, lattice, multivalley, sim  # noqa: E402

F = core.equal_hz()
M = core.moments(F)
T = 0.5 * core.C / 10e6
SURF = bounds.CostSurface.build(F, T)
SJ = M["sigma_jump"]


def test_wwb_at_threshold_is_the_grid_closed_form():
    w = bounds.wwb(SURF, SJ)
    p = lattice.params(F, SJ)
    assert core.db(w / bounds.wwb_lattice(p["s"], p["lambda_eff_m"], T)) == pytest.approx(0.0, abs=0.01)
    assert core.db(w / core.mmse(F, SJ)) == pytest.approx(-0.743, abs=0.005)
    assert w > bounds.zzb(SURF, SJ)
    assert w < lattice.pm_mse(p)["mse"]


def test_wwb_with_exact_wrapped_likelihood_equals_gaussian_form_at_threshold():
    d = bounds.wwb(SURF, SJ, wrapped=True, detail=True)
    assert abs(core.db(d["wwb"] / d["wwb_gauss"])) < 1e-6


def test_wrapped_chernoff_coefficient():
    sigma = 0.08
    r = np.array([0.0, 0.05, 0.3])
    assert bounds.wrapped_log_chernoff(r, sigma, 0.5) == pytest.approx(-(r**2) / (8 * sigma**2), abs=1e-9)
    # At r = π the short and long arcs tie, so the coefficient doubles relative to the Gaussian one.
    sigma = 0.6
    got = bounds.wrapped_log_chernoff(np.array([math.pi]), sigma, 0.5)[0]
    assert got > -(math.pi**2) / (8 * sigma**2) + 0.5


def test_wwb_recovers_the_bandwidth_line_at_large_u():
    p = lattice.params(F, 40 * SJ)
    assert bounds.wwb_lattice(p["s"], p["lambda_eff_m"]) == pytest.approx(core.mmse(F, 40 * SJ) - p["hmse"], rel=1e-3)


def test_multivalley_is_the_lattice_model_on_the_equal_grid():
    for u in (0.5, 1.0):
        mv = multivalley.mse(F, u * SJ, SURF, n=100_000, seed=1)["mse"]
        assert core.db(mv / lattice.ml_mse(lattice.params(F, u * SJ))["mse"]) == pytest.approx(0.0, abs=0.06)


def test_multivalley_reduces_to_the_valley_sum_when_slips_are_rare():
    p = multivalley.mse(F, 0.25 * SJ, SURF, n=100_000, seed=2)
    assert p["mse"] == pytest.approx(p["union"], rel=0.02)


def test_multivalley_tracks_the_geometric_grid_monte_carlo():
    f = core.geometric_hz()
    surf = bounds.CostSurface.build(f, T)
    rows = [r for r in sim.load("e07_layers")["rows"] if r["design"] == "geometric" and 0.15 <= r["sigma"] <= 0.45]
    assert len(rows) >= 4
    for r in rows:
        mv = multivalley.mse(f, r["sigma"], surf, n=100_000, seed=3)["mse"]
        assert abs(core.db(mv / r["mse_coh"])) < 0.25


def test_grid_and_nongrid_valleys():
    eq = design.valley_ratios(F, T=T - 1e-3, step=1e-4)
    assert eq["rho_far"] == pytest.approx(324.0, abs=0.1)
    assert eq["rho_nongrid"] == pytest.approx(672.26, abs=0.1)
    assert eq["delta_nongrid_m"] == pytest.approx(5.725, abs=2e-3)
    assert eq["L_m"] == pytest.approx(T, rel=1e-9)
    lo = design.valley_ratios(np.array([2400, 2401, 2402, 2409, 2450, 2468, 2469, 2470]) * 1e6, step=1e-4)
    assert lo["rho_nongrid"] == pytest.approx(lo["rho_far"])
    assert lo["rho_nongrid"] == pytest.approx(62.6, abs=0.1)


def test_wwb_multi_reduces_to_single_point_at_k1():
    lam_eff = M["lambda_m"] * (1.0 - 1.0 / M["R"])
    w1 = bounds.wwb(SURF, SJ)
    wk = bounds.wwb_multi(SURF, SJ, h=np.array([lam_eff]))
    assert core.db(wk / w1) == pytest.approx(0.0, abs=1e-6)
    p = lattice.params(F, SJ)
    closed = bounds.wwb_multi_lattice(p["s"], p["lambda_eff_m"], T, k_max=1)
    assert core.db(closed["wwb"] / w1) == pytest.approx(0.0, abs=1e-4)


def test_wwb_multi_is_legal_and_tightens_at_threshold():
    p = lattice.params(F, SJ)
    pm = lattice.pm_mse(p)["mse"]
    w1 = bounds.wwb(SURF, SJ)
    wk = bounds.wwb_multi(SURF, SJ, k_grid=8)
    assert wk <= pm
    assert wk >= w1 * 0.999
    assert core.db(pm / wk) < 0.4
    assert core.db(pm / w1) > 0.6
    t = bounds.tightest_bound(SURF, SJ)
    assert t["which"] == "wwb_k"
    assert t["bound"] == pytest.approx(wk, rel=1e-12)


def test_wwb_multi_does_not_beat_zzb_on_the_downslope():
    """u ≈ 0.3–0.5: extra test points stay exponentially weaker than the ZZB integral."""
    for u in (0.3, 0.5):
        sg = u * SJ
        z = bounds.zzb(SURF, sg)
        wk = bounds.wwb_multi(SURF, sg, k_grid=8)
        pm = lattice.pm_mse(lattice.params(F, sg))["mse"]
        assert wk <= pm
        assert z > wk
        t = bounds.tightest_bound(SURF, sg)
        assert t["which"] == "zzb"


def test_zzb_record_valleys_matches_lattice_on_the_equal_grid():
    def zzb_lat(sigma):
        p = lattice.params(F, sigma)
        k = np.arange(1, 400)
        return p["hmse"] + 0.5 * p["lambda_eff_m"] ** 2 * float(np.sum((2 * k - 1) * core.q_func(k / (2 * p["s"]))))

    for u in (0.3, 1.0, 2.0):
        sg = u * SJ
        rec = bounds.zzb_record_valleys(F, sg, SURF.valleys, T)
        assert abs(core.db(rec / zzb_lat(sg))) < 0.15


def test_leave_cap_bounds_every_archived_design():
    d08 = sim.load("e08_design")
    assert design.leave_cap(F) == pytest.approx(d08["refs"]["equal"]["sigma_leave"], rel=2e-3)
    edges = np.array([2400, 2401, 2402, 2403, 2467, 2468, 2469, 2470]) * 1e6
    assert core.moments(edges)["S"] / 1e12 == pytest.approx(design.max_spread_mhz2(2400, 2470, 8))
    gain = 20 * math.log10(design.leave_cap(edges) / d08["refs"]["equal"]["sigma_leave"])
    assert gain == pytest.approx(3.60, abs=0.01)
    for r in d08["random"][:300] + [d08["best"]]:
        assert r["sigma_leave"] <= design.leave_cap(np.asarray(r["mhz"]) * 1e6) * (1 + 2e-3)


def test_lattice_matched_weights_equal_sigma_matches_scalar():
    p0 = lattice.params(F, SJ)
    p1 = lattice.params(F, np.full(F.size, SJ))
    assert p1["s"] == pytest.approx(p0["s"], rel=1e-12)
    assert p1["hmse"] == pytest.approx(p0["hmse"], rel=1e-12)
    sig = np.full(F.size, SJ)
    sig[0] = sig[-1] = 2.0 * SJ
    pw = lattice.params(F, sig)
    assert pw["s"] != pytest.approx(p0["s"], rel=1e-3)
    assert pw["hmse"] > p0["hmse"]
