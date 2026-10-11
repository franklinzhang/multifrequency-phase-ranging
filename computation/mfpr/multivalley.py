"""Linearised multi-valley model: the lattice model with every valley of J as a competitor.

Below the second threshold the noiseless residual at a valley bottom δ_v is the fixed vector
A_v = wrap(g δ_v). Linearising the wrapped cost at each bottom and minimising over the fine
offset, valley v beats the true tooth when

    ΔJ_v = ‖ã_v‖² − 2 ã_vᵀ e < 0,     ã_v = A_v − (gᵀA_v / gᵀg) g,

so every competitor is a linear statistic of the same Gaussian phase noise e. The winner is
argmin_v ΔJ_v (the true tooth has ΔJ = 0), the fine error gᵀe/gᵀg is independent of it, and

    MSE = HMSE + Σ_v P(v wins) D_v²,     D_v = δ_v − gᵀA_v / gᵀg.

If all ã_v are multiples of one vector (equal spacing inside its near zone) this is round(z) of
the lattice model; while slips are rare it is the valley sum Σ 2Q(√J_v/2σ) δ_v². The win
probabilities are Gaussian polytope probabilities, estimated by Monte Carlo in R^N with a
mixture importance sampler centred on each valley's boundary point ã_v/2.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from scipy.special import logsumexp

from .bounds import CostSurface
from .core import C, moments, q_func, wrap


@dataclass
class Competitors:
    a: np.ndarray  # (V, N) projected residuals ã_v, both signs of δ
    d: np.ndarray  # (V,) displacement D_v of the winning valley
    sq: np.ndarray  # (V,) ‖ã_v‖²
    delta: np.ndarray  # (V,) valley position δ_v (signed)
    hmse_factor: float


def competitors(f_hz, valleys) -> Competitors:
    """Linearised competitors from valley bottoms (δ, J) with δ > 0; the main lobe (δ < λ/4) is dropped."""
    f = np.asarray(f_hz, dtype=float)
    m = moments(f)
    g = 2.0 * np.pi * f / C
    v = np.asarray(valleys, dtype=float).reshape(-1, 2)
    dv = v[v[:, 0] > 0.25 * m["lambda_m"], 0]
    A = wrap(dv[:, None] * g[None, :])
    shift = (A @ g) / float(g @ g)
    at = A - shift[:, None] * g[None, :]
    D = dv - shift
    at, D, dv = np.concatenate([at, -at]), np.concatenate([D, -D]), np.concatenate([dv, -dv])
    return Competitors(a=at, d=D, sq=np.sum(at * at, axis=1), delta=dv, hmse_factor=m["hmse_factor"])


def mse(f_hz, sigma: float, surface: CostSurface | None = None, T: float = 15.0, n: int = 200_000,
        seed: int = 0, chunk: int = 20_000, defensive: float = 0.2, rel_cut: float = 1e-6) -> dict:
    """MSE of coherent ML under the linearised multi-valley model at per-tone phase std sigma.

    Competitors whose union-bound share Q(‖ã‖/2σ)·max D² is below rel_cut·HMSE are dropped.
    Returns mse, slip probability, the standard error of the slip term, and win probabilities.
    """
    f = np.asarray(f_hz, dtype=float)
    surf = surface if surface is not None else CostSurface.build(f, T)
    comp = competitors(f, surf.valleys)
    hm = comp.hmse_factor * sigma**2
    ub = q_func(np.sqrt(comp.sq) / (2.0 * sigma))
    keep = ub * float(np.max(comp.d**2, initial=0.0)) >= rel_cut * hm
    out = {"sigma": sigma, "hmse": hm, "mse": hm, "slip": 0.0, "slip_se": 0.0, "n_competitors": int(keep.sum()),
           "union": hm + float(np.sum(ub * comp.d**2)), "wins": []}
    if not keep.any():
        return out
    a, D, sq, dlt = comp.a[keep], comp.d[keep], comp.sq[keep], comp.delta[keep]
    V, N = a.shape
    pi = ub[keep] * D**2
    pi = pi / pi.sum()
    log_pi = np.log(np.maximum(pi, 1e-300))
    mu = 0.5 * a
    rng = np.random.default_rng(seed)
    s_slip = s_slip2 = s_p = 0.0
    wins = np.zeros(V)
    done = 0
    while done < n:
        m_ = min(chunk, n - done)
        comp_idx = np.where(rng.random(m_) < defensive, -1, rng.choice(V, size=m_, p=pi)) if defensive < 1.0 else np.full(m_, -1)
        e = sigma * rng.standard_normal((m_, N))
        e[comp_idx >= 0] += mu[comp_idx[comp_idx >= 0]]
        proj = e @ a.T
        if defensive < 1.0:
            lr = (0.5 * proj - sq[None, :] / 8.0) / sigma**2
            log_q = np.logaddexp(math.log(defensive), math.log1p(-defensive) + logsumexp(lr + log_pi[None, :], axis=1))
            w = np.exp(-log_q)
        else:
            w = np.ones(m_)
        dJ = sq[None, :] - 2.0 * proj
        j = np.argmin(dJ, axis=1)
        won = dJ[np.arange(m_), j] < 0.0
        x = np.where(won, w * D[j] ** 2, 0.0)
        s_slip += float(x.sum())
        s_slip2 += float((x * x).sum())
        s_p += float(np.sum(w * won))
        wins += np.bincount(j[won], weights=w[won], minlength=V)
        done += m_
    slip_term = s_slip / n
    out.update({"mse": hm + slip_term, "slip": s_p / n,
                "slip_se": math.sqrt(max(s_slip2 / n - slip_term**2, 0.0) / n)})
    pw = wins / n
    order = np.argsort(-pw)
    out["wins"] = [(float(dlt[i]), float(pw[i])) for i in order[:12] if pw[i] > 0]
    out["e_d4"] = float(np.sum(pw * D**4))
    return out


def mc_band_db(pred: dict, n_mc: int, z: float = 1.96) -> tuple[float, float]:
    """Spread (dB) of an n_mc-trial sample MSE around this prediction: ±z standard errors.

    Uses the model's own error law: one trial's squared error is (D + ε)², ε ~ N(0, HMSE), so
    E[X²] = Σ P_v D_v⁴ + 6 HMSE·Σ P_v D_v² + 3 HMSE².
    """
    st, hm = pred["mse"] - pred["hmse"], pred["hmse"]
    var = max(pred.get("e_d4", 0.0) + 6.0 * hm * st + 3.0 * hm * hm - pred["mse"] ** 2, 0.0) / n_mc
    sd = math.sqrt(var) / pred["mse"]
    lo = 10.0 * math.log10(max(1.0 - z * sd, 1e-3))
    return lo, 10.0 * math.log10(1.0 + z * sd)
