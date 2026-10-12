"""Design metrics: grid and non-grid valleys, the near zone, and the lattice cap on σ_leave.

A valley of J is a grid valley when the two-scale cost explains it, J(δ_v) = J_2(δ_v): no
tone has wrapped away from the carrier comb plus envelope parabola, so the valley sits at
kλ_eff with depth k²J(λ)(1−1/R) whatever the placement. Always J ≤ J_2 (wrap picks the
shortest representative of every tone), so a non-grid valley has J_2 − J > 0.
"""

from __future__ import annotations

import math

import numpy as np
from scipy.optimize import brentq

from .core import C, moments, phase_cost, q_func, two_scale_cost

GRID_TOL = 1e-8  # rad²; J_2 − J is 0 to rounding on a grid valley


def foot_x(R: float, c: float = 8.0) -> float:
    """Root of c x² Q(x) = 1/R: c = 8 for coherent ML (math §3.3), c = 2 for the lattice ZZB."""
    return brentq(lambda x: c * x * x * q_func(x) - 1.0 / R, 0.3, 30.0)


def near_zone_m(f_hz) -> float:
    """L = c/(2Δ_max): inside it every valley with J < π²/4 is a grid valley (math §6.1)."""
    f = np.sort(np.asarray(f_hz, dtype=float))
    return C / (2.0 * float(np.max(np.diff(f))))


def sampled_valleys(f_hz, T: float = 15.0, step: float = 2e-4):
    """Valleys of J in (λ/4, T): parabola-refined (δ, J) and J_2 − J at the grid minimum."""
    f = np.asarray(f_hz, dtype=float)
    d = np.arange(step, T, step)
    J = phase_cost(f, d)
    i = np.flatnonzero((J[1:-1] <= J[:-2]) & (J[1:-1] < J[2:])) + 1
    i = i[d[i] > 0.25 * C / f.mean()]
    a, b, c = J[i - 1], J[i], J[i + 1]
    den = a - 2.0 * b + c
    ok = den > 0
    off = np.where(ok, 0.5 * (a - c) / np.where(ok, den, 1.0), 0.0)
    dv = d[i] + off * step
    jv = np.maximum(b - 0.25 * (a - c) * off, 0.0)
    return dv, jv, two_scale_cost(f, d[i]) - b


def valley_ratios(f_hz, T: float = 15.0, step: float = 2e-4) -> dict:
    """Shallowest valley over J(λ): beyond c/2B (ρ_far, old O3) and among non-grid valleys (ρ_nongrid)."""
    f = np.asarray(f_hz, dtype=float)
    jl = float(phase_cost(f, [moments(f)["lambda_m"]])[0])
    dv, jv, gap = sampled_valleys(f, T, step)
    masks = {"far": dv > C / (2.0 * (f.max() - f.min())), "nongrid": gap > GRID_TOL}
    out = {"J_lambda": jl, "L_m": near_zone_m(f), "n_valleys": int(dv.size), "n_nongrid": int(masks["nongrid"].sum())}
    for key, mask in masks.items():
        if mask.any():
            j = int(np.flatnonzero(mask)[np.argmin(jv[mask])])
            out[f"rho_{key}"], out[f"delta_{key}_m"] = float(jv[j] / jl), float(dv[j])
        else:
            out[f"rho_{key}"], out[f"delta_{key}_m"] = math.inf, None
    return out


def leave_cap(f_hz, factor: float = 2.0) -> float:
    """σ where the ±1 grid valleys alone lift the carrier line to factor × HMSE.

    Every design has the k = ±1 grid valleys (λ < c/2B ≤ L), and any other valley only adds to
    the valley sum, so σ_leave ≤ leave_cap for every design. ≈ σ_jump / x(R).
    """
    f = np.asarray(f_hz, dtype=float)
    m = moments(f)
    lam_eff = m["lambda_m"] * (1.0 - 1.0 / m["R"])
    root_j1 = math.sqrt(float(phase_cost(f, [lam_eff])[0]))
    h = m["hmse_factor"]
    return brentq(lambda s: 2.0 * q_func(root_j1 / (2.0 * s)) * lam_eff**2 - (factor - 1.0) * h * s * s, 1e-4, m["sigma_jump"])


def max_spread_mhz2(lo: int, hi: int, n: int) -> float:
    """Upper bound on S for n distinct integers in [lo, hi]: Σ(f − f̄)² ≤ Σ(f − mid)², take the n farthest."""
    mid = 0.5 * (lo + hi)
    dist = np.sort(np.abs(np.arange(lo, hi + 1) - mid))[::-1]
    return float(np.sum(dist[:n] ** 2))
