"""Range estimators on wrapped phases and on complex tones.

All searches run on a grid of step λ/8 over a finite window and refine every retained
candidate with Gauss–Newton. Costs are weighted by w_i = 1/σ_i², so exp(−J/2) is the
wrapped-Gaussian likelihood in its high-SNR form.
"""

from __future__ import annotations

import numpy as np

from .core import C, wrap


def make_grid(lo: float, hi: float, lam: float, oversample: int = 8) -> np.ndarray:
    step = lam / oversample
    return np.arange(lo, hi + 0.5 * step, step)


def _weights(f, w):
    return np.ones_like(f) if w is None else np.broadcast_to(np.asarray(w, dtype=float), f.shape)


def _refine_coherent(phi, f, d, w, lo, hi, iters=12):
    """Gauss–Newton on Σ w wrap(φ − 2πfd/c)². phi: (n, N), d: (n, K)."""
    g = 2.0 * np.pi * f / C
    gw = g * w
    den = float(np.sum(gw * g))
    for _ in range(iters):
        r = wrap(phi[:, None, :] - d[..., None] * g)
        d = np.clip(d + (r @ gw) / den, lo, hi)
    r = wrap(phi[:, None, :] - d[..., None] * g)
    return d, np.sum(w * r**2, axis=-1)


def coherent_candidates(phi, f, grid, w=None, k: int = 1, block: int = 400):
    """Return the k lowest refined local minima of the coherent cost for each trial.

    phi: (n, N) wrapped phases. Returns (d, J) of shape (n, k), sorted by J.
    """
    f = np.asarray(f, dtype=float)
    w = _weights(f, w)
    g = 2.0 * np.pi * f / C
    pred = grid[:, None] * g[None, :]
    n = phi.shape[0]
    d_out = np.empty((n, k))
    j_out = np.empty((n, k))
    lo, hi = float(grid[0]), float(grid[-1])
    for s0 in range(0, n, block):
        p = phi[s0 : s0 + block]
        cost = np.sum(w * wrap(p[:, None, :] - pred[None, :, :]) ** 2, axis=-1)
        left = np.concatenate([np.full((cost.shape[0], 1), np.inf), cost[:, :-1]], axis=1)
        right = np.concatenate([cost[:, 1:], np.full((cost.shape[0], 1), np.inf)], axis=1)
        is_min = (cost <= left) & (cost < right)
        masked = np.where(is_min, cost, np.inf)
        kk = min(k, grid.size)
        idx = np.argpartition(masked, kk - 1, axis=1)[:, :kk]
        d0 = grid[idx]
        # Rows with fewer than k minima pick inf entries; their refined copies are harmless.
        d, J = _refine_coherent(p, f, d0, w, lo, hi)
        J = np.where(np.isfinite(np.take_along_axis(masked, idx, axis=1)), J, np.inf)
        order = np.argsort(J, axis=1)
        d_out[s0 : s0 + block] = np.take_along_axis(d, order, axis=1)
        j_out[s0 : s0 + block] = np.take_along_axis(J, order, axis=1)
    return d_out, j_out


def coherent_ml(phi, f, grid, w=None, k: int = 4):
    """Global minimizer of the coherent cost.

    On a λ/8 grid the sampled cost of a tooth depends on where the grid falls inside it,
    and that offset alone is about a tenth of the cost difference between neighbouring
    teeth. Several grid minima are therefore refined before the comparison.
    """
    d, _ = coherent_candidates(phi, f, grid, w, k=k)
    return d[:, 0]


def posterior_mean(phi, f, grid, w=None, k: int = 24):
    """Posterior mean under a flat prior on the grid window, Laplace on each tooth.

    Teeth share the same curvature Σ w g², so each peak is weighted by exp(−J/2).
    Also returns the ML estimate and the posterior mass of the ML tooth.
    """
    d, J = coherent_candidates(phi, f, grid, w, k=k)
    Jm = J[:, :1]
    lw = -0.5 * (J - Jm)
    wt = np.where(np.isfinite(lw), np.exp(lw), 0.0)
    wt /= wt.sum(axis=1, keepdims=True)
    return {"pm": np.sum(wt * d, axis=1), "ml": d[:, 0], "p_top": wt[:, 0], "J": J, "d": d}


def incoherent_ml(phi, f, grid, block: int = 400):
    """Maximize |Σ exp(j(φ − 2πfd/c))|: the common phase is a nuisance parameter."""
    f = np.asarray(f, dtype=float)
    g = 2.0 * np.pi * f / C
    fc = f - f.mean()
    S = float(np.sum(fc**2))
    ph = np.exp(-1j * grid[:, None] * g[None, :])
    n = phi.shape[0]
    out = np.empty(n)
    lo, hi = float(grid[0]), float(grid[-1])
    for s0 in range(0, n, block):
        z = np.exp(1j * phi[s0 : s0 + block])
        val = np.abs(z @ ph.T)
        d = grid[np.argmax(val, axis=1)]
        for _ in range(12):
            pred = d[:, None] * g[None, :]
            theta = np.angle(np.sum(np.exp(1j * (phi[s0 : s0 + block] - pred)), axis=1))
            r = wrap(phi[s0 : s0 + block] - pred - theta[:, None])
            d = np.clip(d + (C / (2.0 * np.pi)) * (r @ fc) / S, lo, hi)
        out[s0 : s0 + block] = d
    return out


def complex_coherent_ml(y, f, grid, block: int = 400, k: int = 4):
    """Maximize Re Σ y_i exp(−j2πf_i d/c): known amplitude and phase reference.

    As for the phase estimator, the k best grid maxima are refined before comparing.
    """
    f = np.asarray(f, dtype=float)
    g = 2.0 * np.pi * f / C
    ph = np.exp(-1j * grid[:, None] * g[None, :])
    n = y.shape[0]
    out = np.empty(n)
    lo, hi = float(grid[0]), float(grid[-1])
    for s0 in range(0, n, block):
        yy = y[s0 : s0 + block]
        val = np.real(yy @ ph.T)
        left = np.concatenate([np.full((val.shape[0], 1), -np.inf), val[:, :-1]], axis=1)
        right = np.concatenate([val[:, 1:], np.full((val.shape[0], 1), -np.inf)], axis=1)
        masked = np.where((val >= left) & (val > right), val, -np.inf)
        idx = np.argpartition(-masked, k - 1, axis=1)[:, :k]
        d = grid[idx]
        mag = np.abs(yy)[:, None, :]
        for _ in range(12):
            z = yy[:, None, :] * np.exp(-1j * d[..., None] * g)
            grad = np.sum(g * np.imag(z), axis=-1)
            curv = np.sum(g * g * mag, axis=-1)
            d = np.clip(d + grad / curv, lo, hi)
        obj = np.sum(np.real(yy[:, None, :] * np.exp(-1j * d[..., None] * g)), axis=-1)
        obj = np.where(np.isfinite(np.take_along_axis(masked, idx, axis=1)), obj, -np.inf)
        out[s0 : s0 + block] = np.take_along_axis(d, np.argmax(obj, axis=1)[:, None], axis=1)[:, 0]
    return out
