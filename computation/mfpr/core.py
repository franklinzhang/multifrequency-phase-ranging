"""Observation model, frequency designs, and the phase cost surface."""

from __future__ import annotations

import math

import numpy as np
from scipy.special import erfc, ndtr

C = 299_792_458.0

EQUAL_MHZ = 2400.0 + 10.0 * np.arange(8)
# Frozen output of the original search in sim_threshold.design_moment_matched: integer MHz,
# same N, mean and second moment as the equal grid, gcd of gaps 1 MHz.
MATCHED_OFFSETS_MHZ = np.array([-45, -17, -10, -3, 8, 14, 19, 34], dtype=float)


def equal_hz() -> np.ndarray:
    return EQUAL_MHZ * 1e6


def matched_hz() -> np.ndarray:
    return (2435.0 + MATCHED_OFFSETS_MHZ) * 1e6


def geometric_hz(ratio: float = 1.02, n: int = 8, f_mid: float = 2435e6) -> np.ndarray:
    return f_mid * ratio ** (np.arange(n) - (n - 1) / 2.0)


def uniform_grid_hz(n: int, bandwidth_hz: float, fbar_hz: float = 2435e6) -> np.ndarray:
    return fbar_hz + bandwidth_hz * (np.arange(n) - (n - 1) / 2.0) / max(n - 1, 1)


def wrap(x):
    return (np.asarray(x) + np.pi) % (2.0 * np.pi) - np.pi


def q_func(x):
    return 0.5 * erfc(np.asarray(x, dtype=float) / math.sqrt(2.0))


def phi_cdf(x):
    return ndtr(np.asarray(x, dtype=float))


def db(x):
    return 10.0 * np.log10(np.asarray(x, dtype=float))


def moments(f_hz) -> dict:
    f = np.asarray(f_hz, dtype=float)
    fbar = float(f.mean())
    S = float(np.sum((f - fbar) ** 2))
    F = float(np.sum(f**2))
    return {
        "n": int(f.size),
        "f_hz": f.tolist(),
        "fbar_hz": fbar,
        "S": S,
        "F": F,
        "R": F / S,
        "R_db": 10.0 * math.log10(F / S),
        "sigma_jump": math.pi * math.sqrt(S) / fbar,
        "lambda_m": C / fbar,
        "hmse_factor": C**2 / (4.0 * math.pi**2 * F),
        "mmse_factor": C**2 / (4.0 * math.pi**2 * S),
    }


def hmse(f_hz, sigma):
    return moments(f_hz)["hmse_factor"] * np.asarray(sigma, dtype=float) ** 2


def mmse(f_hz, sigma):
    return moments(f_hz)["mmse_factor"] * np.asarray(sigma, dtype=float) ** 2


def phase_cost(f_hz, delta_m) -> np.ndarray:
    """J(δ) = Σ wrap(2π f_i δ / c)², noiseless cost of a range offset δ."""
    f = np.asarray(f_hz, dtype=float)
    d = np.atleast_1d(np.asarray(delta_m, dtype=float))
    out = np.empty(d.size)
    step = max(1, 4_000_000 // max(f.size, 1))
    for start in range(0, d.size, step):
        seg = d[start : start + step]
        out[start : start + step] = np.sum(wrap((2.0 * np.pi / C) * seg[:, None] * f[None, :]) ** 2, axis=1)
    return out


def two_scale_cost(f_hz, delta_m) -> np.ndarray:
    """N·wrap(2πδ/λ)² + (2π/c)² S δ²: carrier comb plus envelope parabola."""
    m = moments(f_hz)
    d = np.asarray(delta_m, dtype=float)
    return m["n"] * wrap(2.0 * np.pi * d / m["lambda_m"]) ** 2 + (2.0 * np.pi / C) ** 2 * m["S"] * d**2


def snr_from_sigma(sigma):
    """Per-tone SNR whose high-SNR phase standard deviation is sigma: σ² ≈ 1/(2 SNR)."""
    return 1.0 / (2.0 * np.asarray(sigma, dtype=float) ** 2)


def sigma_from_snr(snr):
    return np.sqrt(1.0 / (2.0 * np.asarray(snr, dtype=float)))


def wilson(k: int, n: int, z: float = 1.96) -> tuple[float, float]:
    if n == 0:
        return (0.0, 1.0)
    p = k / n
    den = 1.0 + z**2 / n
    mid = (p + z**2 / (2 * n)) / den
    half = z * math.sqrt(p * (1 - p) / n + z**2 / (4 * n**2)) / den
    return (max(0.0, mid - half), min(1.0, mid + half))
