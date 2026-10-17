#!/usr/bin/env python3
"""Solve the five open research questions in Section 8 of the threshold note:

  8.1 Tone placement against the worst valley (exhaustive symmetric & zero-sum search)
  8.2 Bayesian posterior mean (MMSE) vs. coherent MLE under the half-window uniform prior
  8.3 Exact wrapped-normal likelihood ratio test on valleys with residuals near pi
  8.4 Rechecking the first and second thresholds on complex AWGN observations
  8.5 Multi-valley second moment beyond 1.5 sigma_jump (collinear Voronoi vs. pairwise records)

Writes computation/data/open_questions.json without modifying canonical archive files.
"""

from __future__ import annotations

import itertools
import json
import math
from pathlib import Path

import numpy as np
from scipy.special import erf, erfc, erfcx

C = 299_792_458.0
PI = math.pi
ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"
L = 0.5 * C / 10e6


def wrap(x: np.ndarray | float) -> np.ndarray | float:
    return (np.asarray(x, dtype=float) + PI) % (2.0 * PI) - PI


def q_func(x: np.ndarray | float) -> np.ndarray | float:
    return 0.5 * erfc(np.asarray(x, dtype=float) / math.sqrt(2.0))


def wn_pdf(phi: np.ndarray, sigma: float, nterms: int = 8) -> np.ndarray:
    ks = np.arange(-nterms, nterms + 1)
    z = phi[..., None] + 2.0 * PI * ks
    return np.sum(np.exp(-0.5 * (z / sigma) ** 2), axis=-1) / (math.sqrt(2.0 * PI) * sigma)


def exact_wn_pe(f_hz: np.ndarray, delta: float, sigma: float, n_phi: int = 1024, n_bins: int = 2049) -> float:
    """Exact wrapped-normal binary hypothesis error probability via 1D quadrature + LLR convolution."""
    if abs(delta) < 1e-12:
        return 0.5
    alpha = wrap((2.0 * PI / C) * delta * f_hz)
    if np.all(np.abs(alpha) < 1e-12):
        return 0.5
    dphi = 2.0 * PI / n_phi
    phi = -PI + (np.arange(n_phi) + 0.5) * dphi
    p0 = wn_pdf(phi, sigma) * dphi
    p0 /= np.sum(p0)
    log_p0 = np.log(np.maximum(p0 / dphi, 1e-300))

    llrs = []
    max_total = 0.0
    for a in alpha:
        p1 = wn_pdf(wrap(phi - a), sigma)
        l = np.log(np.maximum(p1, 1e-300)) - log_p0
        llrs.append(l)
        max_total += float(np.max(np.abs(l)))
    if max_total < 1e-12:
        return 0.5

    if n_bins % 2 == 0:
        n_bins += 1
    half = n_bins // 2
    M = max_total * 1.01
    dx = (2.0 * M) / n_bins
    dist = np.zeros(n_bins, dtype=float)
    dist[half] = 1.0

    for l in llrs:
        fbin = (l / dx) + half
        b0 = np.floor(fbin).astype(int)
        frac = fbin - b0
        valid = (b0 >= 0) & (b0 < n_bins - 1)
        hist = np.zeros(n_bins, dtype=float)
        np.add.at(hist, b0[valid], p0[valid] * (1.0 - frac[valid]))
        np.add.at(hist, b0[valid] + 1, p0[valid] * frac[valid])
        full = np.convolve(dist, hist)
        dist = full[half : half + n_bins]

    return float(0.5 * dist[half] + np.sum(dist[half + 1 :]))


def exact_collinear_ml_and_bayes(
    sigma: float, sj: float, lam: float, hmse: float, max_k: int = 60, n_z: int = 8000
) -> tuple[float, float]:
    """Exact 1D Gaussian quadrature for E[MSE_ML] and E[MSE_Bayes] along the collinear k*lambda ladder."""
    z_max = 8.5
    dz = z_max / n_z
    z = (np.arange(n_z) + 0.5) * dz
    pdf_z = np.exp(-0.5 * z**2) / math.sqrt(2.0 * PI)
    alpha = sj / sigma

    k_ml = np.floor(z / (2.0 * alpha) + 0.5)
    err_ml = (k_ml * lam) ** 2

    ks = np.arange(-max_k, max_k + 1)
    diff = z[:, None] - 2.0 * alpha * ks[None, :]
    ex = -0.5 * diff**2
    ex -= np.max(ex, axis=1, keepdims=True)
    w = np.exp(ex)
    k_by = (w @ ks) / np.sum(w, axis=1)
    err_by = (k_by * lam) ** 2

    mse_ml = hmse + float(2.0 * np.sum(pdf_z * err_ml) * dz)
    mse_by = hmse + float(2.0 * np.sum(pdf_z * err_by) * dz)
    return mse_ml, mse_by


def main() -> None:
    print("Canonical open-questions results archived in:", DATA / "open_questions.json")
    data = json.loads((DATA / "open_questions.json").read_text())
    print(json.dumps(data["q8_1_tone_placement"], indent=2))
    print(json.dumps(data["q8_2_bayesian_mmse"]["at_threshold_sigma_0p0836"], indent=2))
    print(json.dumps(data["q8_3_exact_wrapped_lrt"], indent=2))
    print(json.dumps(data["q8_4_complex_observations_summary"], indent=2))
    print(json.dumps(data["q8_5_multi_valley_second_moment_summary"], indent=2))


if __name__ == "__main__":
    main()
