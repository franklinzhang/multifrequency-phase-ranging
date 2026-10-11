"""Two-scale lattice model of the coherent estimator.

Below the second threshold the per-tone envelope phases stay small, so the coherent
cost at δ = kλ + ε is ||e − k a − ε g||² (weighted by 1/σ_i²). Every carrier slip lies on
one line in phase space. The integer decision is then a scalar lattice decoder

    k̂ = round(z + β),   z ~ N(0, s²),

and the range error is k̂·λ_eff + (fine error), with the fine error independent of k̂.
"""

from __future__ import annotations

import math

import numpy as np

from .core import C, phi_cdf


def params(f_hz, sigma, offset_rad=None, w_rx=None) -> dict:
    """Lattice parameters for tones f with per-tone phase std sigma (scalar or vector).

    offset_rad is a deterministic per-tone phase error (calibration, multipath).
    w_rx are the weights the receiver uses; by default the matched weights 1/σ_i².
    """
    f = np.asarray(f_hz, dtype=float)
    sig = np.broadcast_to(np.asarray(sigma, dtype=float), f.shape)
    var = sig**2
    w = 1.0 / var if w_rx is None else np.broadcast_to(np.asarray(w_rx, dtype=float), f.shape)
    fbar = float(f.mean())
    lam = C / fbar
    g = 2.0 * np.pi * f / C
    a = 2.0 * np.pi * (f - fbar) / fbar
    gWg = float(np.sum(w * g * g))
    aWg = float(np.sum(w * a * g))
    u = a - (aWg / gWg) * g
    uWu = float(np.sum(w * u * u))
    out = {
        "lambda_m": lam,
        "lambda_eff_m": lam - aWg / gWg,
        "s": math.sqrt(float(np.sum(w * w * var * u * u))) / uWu,
        "hmse": float(np.sum(w * w * var * g * g)) / gWg**2,
        "beta": 0.0,
        "fine_bias_m": 0.0,
    }
    if offset_rad is not None:
        m = np.broadcast_to(np.asarray(offset_rad, dtype=float), f.shape)
        out["beta"] = float(np.sum(w * u * m)) / uWu
        out["fine_bias_m"] = float(np.sum(w * g * m)) / gWg
    return out


def integer_pmf(s: float, beta: float = 0.0, kmax: int | None = None):
    if kmax is None:
        kmax = int(math.ceil(abs(beta) + 9.0 * s)) + 3
    k = np.arange(-kmax, kmax + 1)
    p = phi_cdf((k + 0.5 - beta) / s) - phi_cdf((k - 0.5 - beta) / s)
    return k, p


def ml_mse(p: dict) -> dict:
    """MSE of coherent ML under the lattice model, with its slip probability."""
    k, pk = integer_pmf(p["s"], p["beta"])
    lam = p["lambda_eff_m"]
    b = p["fine_bias_m"]
    err = k * lam + b
    return {
        "mse": float(p["hmse"] + np.sum(pk * err**2)),
        "slip": float(1.0 - pk[k == 0][0]),
        "e_k2": float(np.sum(pk * k**2)),
    }


def pm_mse(p: dict, nz: int = 4001) -> dict:
    """Bayes MSE of the posterior mean of the integer under a flat prior on k (no bias)."""
    s = p["s"]
    zmax = 0.5 + 8.0 * s
    z = np.linspace(-zmax, zmax, nz)
    kmax = int(math.ceil(zmax + 8.0 * s)) + 2
    k = np.arange(-kmax, kmax + 1)
    # z is the decision statistic; the true integer is 0. Posterior over k given z.
    logw = -0.5 * ((z[:, None] - k[None, :]) / s) ** 2
    logw -= logw.max(axis=1, keepdims=True)
    wk = np.exp(logw)
    wk /= wk.sum(axis=1, keepdims=True)
    kbar = wk @ k
    pz = np.exp(-0.5 * (z / s) ** 2) / (s * math.sqrt(2 * math.pi))
    e2 = float(np.trapezoid(pz * kbar**2, z))
    return {"mse": float(p["hmse"] + p["lambda_eff_m"] ** 2 * e2), "e_k2": e2}


def universal_curve(u, R: float, kind: str = "ml"):
    """MSE / MMSE as a function of u = σ/σ_jump for an equal-weight design with ratio R."""
    u = np.atleast_1d(np.asarray(u, dtype=float))
    out = np.empty(u.size)
    for i, ui in enumerate(u):
        s = ui / (2.0 * math.sqrt(1.0 - 1.0 / R))
        lam_eff = 1.0 - 1.0 / R  # in units of λ
        # λ²/MMSE = 4/u², HMSE/MMSE = 1/R.
        if kind == "ml":
            k, pk = integer_pmf(s)
            ek2 = float(np.sum(pk * k**2))
        else:
            ek2 = pm_mse({"s": s, "hmse": 0.0, "lambda_eff_m": 1.0})["e_k2"]
        out[i] = 1.0 / R + (4.0 / ui**2) * lam_eff**2 * ek2
    return out
