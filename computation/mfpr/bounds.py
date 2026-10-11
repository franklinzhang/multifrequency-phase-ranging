"""Valley-filled Ziv–Zakai bound, its attribution to valleys, the single-test-point
Weiss–Weinstein bound, and valley-sum (MIE) predictions.

Prior: true range uniform on an interval of length T. Pairwise error probability is the
Gaussian phase-vector test P_e(δ) = Q(√J(δ) / 2σ). The bound is

    ZZB = (1/T) ∫_0^T δ · max_{u≥δ} [(T−u) P_e(u)] dδ.

The integrand is resolved three ways: a geometric grid on the main lobe (width σ/β, far
below any uniform grid at high SNR), a uniform grid, and the refined bottom of every valley
of J, where (T−δ)P_e peaks. Missing a valley bottom by a fraction of a millimetre changes
the bound by decibels near the threshold.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from scipy.optimize import minimize_scalar

from .core import C, moments, phase_cost, q_func, wrap


@dataclass
class CostSurface:
    f: np.ndarray
    T: float
    delta: np.ndarray
    J: np.ndarray
    valleys: np.ndarray  # refined valley bottoms (δ, J)

    @classmethod
    def build(cls, f_hz, T: float, step: float = 2e-4) -> "CostSurface":
        f = np.asarray(f_hz, dtype=float)
        u = np.arange(step, T, step)
        Ju = phase_cost(f, u)
        is_min = np.zeros(u.size, dtype=bool)
        is_min[1:-1] = (Ju[1:-1] <= Ju[:-2]) & (Ju[1:-1] < Ju[2:])
        vals = []
        for i in np.flatnonzero(is_min):
            res = minimize_scalar(
                lambda x: float(phase_cost(f, np.array([x]))[0]),
                bounds=(u[i - 1], u[i + 1]),
                method="bounded",
                options={"xatol": 1e-9},
            )
            vals.append((res.x, res.fun))
        valleys = np.array(vals) if vals else np.zeros((0, 2))
        return cls(f=f, T=T, delta=u, J=Ju, valleys=valleys)

    def grid_for(self, sigma: float) -> tuple[np.ndarray, np.ndarray]:
        beta = 2.0 * math.pi * math.sqrt(moments(self.f)["F"]) / C
        w = sigma / beta
        lobe = np.geomspace(1e-4 * w, min(60.0 * w, self.T), 600)
        x = np.concatenate([[0.0], lobe, self.delta, self.valleys[:, 0] if self.valleys.size else []])
        x = np.unique(x[(x >= 0) & (x < self.T)])
        J = phase_cost(self.f, x)
        return x, J


def zzb(surface: CostSurface, sigma: float, detail: bool = False):
    x, J = surface.grid_for(sigma)
    T = surface.T
    pe = q_func(np.sqrt(np.maximum(J, 0.0)) / (2.0 * sigma))
    w = (T - x) * pe
    G = np.maximum.accumulate(w[::-1])[::-1]
    integrand = x * G
    total = float(np.trapezoid(integrand, x) / T)
    if not detail:
        return total
    cume = np.concatenate([[0.0], np.cumsum(0.5 * (integrand[1:] + integrand[:-1]) * np.diff(x))]) / T
    return {"zzb": total, "x": x, "J": J, "w": w, "G": G, "integrand": integrand, "cume": cume}


def zzb_pieces(surface: CostSurface, sigma: float, min_share: float = 1e-4) -> dict:
    """Split the filled integral by which valley's (T−δ)P_e sets the fill height.

    Each right-to-left record owns the segment back to the previous (nearer) record.
    Owners inside a quarter wavelength are reported together as the main lobe.
    """
    d = zzb(surface, sigma, detail=True)
    x, w, G, cume = d["x"], d["w"], d["G"], d["cume"]
    total = d["zzb"]
    lam = moments(surface.f)["lambda_m"]
    # owner[i] = index of the sample right of i where max_{u≥x_i} w(u) is attained.
    owner = np.empty(x.size, dtype=int)
    best = -1.0
    arg = x.size - 1
    for i in range(x.size - 1, -1, -1):
        if w[i] > best:
            best = w[i]
            arg = i
        owner[i] = arg
    vx = surface.valleys[:, 0] if surface.valleys.size else np.zeros(0)
    vj = surface.valleys[:, 1] if surface.valleys.size else np.zeros(0)

    def label(j: int) -> int:
        """-1 for the main lobe, else the nearest valley index."""
        if x[j] < 0.25 * lam or vx.size == 0:
            return -1
        return int(np.argmin(np.abs(vx - x[j])))

    seg_share = np.diff(cume)  # share of interval (x_i, x_{i+1}], owned by owner[i+1]
    acc: dict[int, dict] = {}
    for i in range(seg_share.size):
        if seg_share[i] <= 0:
            continue
        j = owner[i + 1]
        lab = label(j)
        item = acc.setdefault(lab, {"share": 0.0, "left_m": float(x[i]), "right_m": float(x[i + 1])})
        item["share"] += float(seg_share[i])
        item["right_m"] = float(x[i + 1])
    merged = []
    for lab, item in sorted(acc.items(), key=lambda kv: -kv[1]["share"]):
        if item["share"] < min_share * total:
            continue
        item["valley_m"] = None if lab < 0 else float(vx[lab])
        item["valley_J"] = None if lab < 0 else float(vj[lab])
        item["frac"] = item["share"] / total
        item["share_db"] = 10.0 * math.log10(item["share"])
        merged.append(item)
    k50 = int(np.searchsorted(cume, 0.5 * total))
    k90 = int(np.searchsorted(cume, 0.9 * total))
    return {
        "sigma": sigma,
        "zzb": total,
        "zzb_db": 10.0 * math.log10(total),
        "d50_m": float(x[min(k50, x.size - 1)]),
        "d90_m": float(x[min(k90, x.size - 1)]),
        "pieces": merged,
    }


def valley_sum_mse(surface: CostSurface, sigma, lam_floor: float | None = None):
    """MIE-style prediction: HMSE + Σ_valleys 2 Q(√J_k/2σ) δ_k² (union over ± sides).

    Only valleys deeper than their neighbours count. This is an approximation, not a bound.
    """
    m = moments(surface.f)
    sig = np.atleast_1d(np.asarray(sigma, dtype=float))
    v = surface.valleys
    out = m["hmse_factor"] * sig**2
    if v.size:
        p = 2.0 * q_func(np.sqrt(v[None, :, 1]) / (2.0 * sig[:, None]))
        out = out + np.sum(p * v[None, :, 0] ** 2, axis=1)
    return out


def valley_mixture_mse(surface: CostSurface, sigma):
    """Valley sum where every valley competes with the true lobe: weights p_k / (1 + Σ p_j).

    Equal to the valley sum while slips are rare; once the pairwise probabilities add up
    past one it stays a proper mixture, the true lobe keeping weight 1 / (1 + Σ p_j).
    """
    m = moments(surface.f)
    sig = np.atleast_1d(np.asarray(sigma, dtype=float))
    out = m["hmse_factor"] * sig**2
    v = surface.valleys
    if v.size:
        p = 2.0 * q_func(np.sqrt(v[None, :, 1]) / (2.0 * sig[:, None]))
        tot = 1.0 + p.sum(axis=1)
        out = out + np.sum(p * v[None, :, 0] ** 2, axis=1) / tot
    return out


def leave_sigma(surface: CostSurface, factor: float = 2.0, lo: float = 1e-3, hi: float = 3.0) -> float:
    """Smallest σ where the valley-sum prediction reaches factor × HMSE."""
    sig = np.geomspace(lo, hi, 2000)
    ratio = valley_sum_mse(surface, sig) / (moments(surface.f)["hmse_factor"] * sig**2)
    idx = np.flatnonzero(ratio >= factor)
    if idx.size == 0:
        return float("nan")
    i = int(idx[0])
    if i == 0:
        return float(sig[0])
    # Log-linear interpolation between the bracketing samples.
    r0, r1 = math.log(ratio[i - 1]), math.log(ratio[i])
    t = (math.log(factor) - r0) / (r1 - r0)
    return float(math.exp(math.log(sig[i - 1]) + t * (math.log(sig[i]) - math.log(sig[i - 1]))))


def _log_wrapped_normal(x, sigma: float, images: int):
    k = np.arange(-images, images + 1)
    z = np.asarray(x)[..., None] + 2.0 * math.pi * k
    lse = np.logaddexp.reduce(-0.5 * (z / sigma) ** 2, axis=-1)
    return lse - math.log(sigma * math.sqrt(2.0 * math.pi))


def wrapped_log_chernoff(r, sigma: float, s: float, ngrid: int | None = None) -> np.ndarray:
    """ln ∫ p(φ−r)^s p(φ)^{1−s} dφ over one turn for the wrapped normal p, elementwise in r.

    The integrand is smooth and periodic, so the midpoint rule converges spectrally.
    """
    r = np.asarray(r, dtype=float)
    n = ngrid or max(4096, int(math.ceil(40.0 * math.pi / sigma)))
    images = int(math.ceil(6.0 * sigma / (2.0 * math.pi))) + 2
    phi = -math.pi + 2.0 * math.pi * (np.arange(n) + 0.5) / n
    l0 = _log_wrapped_normal(phi, sigma, images)
    out = np.empty(r.size)
    flat = r.ravel()
    for i0 in range(0, flat.size, 64):
        rr = flat[i0 : i0 + 64]
        l1 = _log_wrapped_normal(wrap(phi[None, :] - rr[:, None]), sigma, images)
        out[i0 : i0 + 64] = np.logaddexp.reduce(s * l1 + (1.0 - s) * l0[None, :], axis=1)
    return (out + math.log(2.0 * math.pi / n)).reshape(r.shape)


def _wwb_from_eta(h, eta_s_h, eta_2s_h, eta_22s_h, eta_s_2h):
    """h² e^{2η(s,h)} / (e^{η(2s,h)} + e^{η(2−2s,−h)} − 2 e^{η(s,2h)}), in logs."""
    m = np.maximum(eta_2s_h, eta_22s_h)
    den = np.exp(eta_2s_h - m) + np.exp(eta_22s_h - m) - 2.0 * np.exp(eta_s_2h - m)
    with np.errstate(divide="ignore", invalid="ignore"):
        val = np.exp(2.0 * np.log(h) + 2.0 * eta_s_h - m) / den
    return np.where(den > 0, val, 0.0)


def wwb(surface: CostSurface, sigma: float, s=0.5, wrapped: bool = False, n_wrapped: int = 48, detail: bool = False):
    """Single-test-point Weiss–Weinstein bound under the same uniform prior of length T.

    η(s,h) = ln((T−h)/T) + ln C_s(h), with C_s(h) = ∫ p(φ|d+h)^s p(φ|d)^{1−s} dφ depending only on
    wrap(g h). Gaussian phase vector: ln C_s = −s(1−s) J(h)/2σ². wrapped=True uses the exact
    wrapped-normal coefficient tone by tone on the n_wrapped best Gaussian test points. The bound is
    the maximum over test points h (main lobe, uniform grid, valley bottoms) and over s (scalar or array).
    """
    x, J = surface.grid_for(sigma)
    T = surface.T
    h = x[x > 0]
    J1 = J[x > 0]
    J2 = phase_cost(surface.f, 2.0 * h)
    ss = np.atleast_1d(np.asarray(s, dtype=float))
    log_a = np.log((T - h) / T)
    with np.errstate(divide="ignore"):
        log_b = np.log(np.maximum(T - 2.0 * h, 0.0) / T)

    def gauss(si, Jh):
        return -si * (1.0 - si) * Jh / (2.0 * sigma**2)

    vals = np.stack([_wwb_from_eta(h, log_a + gauss(si, J1), log_a + gauss(2 * si, J1),
                                   log_a + gauss(2 - 2 * si, J1), log_b + gauss(si, J2)) for si in ss])
    best = np.unravel_index(int(np.argmax(vals)), vals.shape)
    out = {"wwb": float(vals[best]), "h_m": float(h[best[1]]), "s": float(ss[best[0]]), "wrapped": False}
    if wrapped:
        top = np.argsort(-vals.max(axis=0))[:n_wrapped]
        g = 2.0 * math.pi * surface.f / C
        r1 = wrap(h[top, None] * g[None, :])
        r2 = wrap(2.0 * h[top, None] * g[None, :])

        def lc(si, r):
            return wrapped_log_chernoff(r, sigma, si).sum(axis=1)

        wv = np.stack([_wwb_from_eta(h[top], log_a[top] + lc(si, r1), log_a[top] + lc(2 * si, r1),
                                     log_a[top] + lc(2 - 2 * si, r1), log_b[top] + lc(si, r2)) for si in ss])
        bw = np.unravel_index(int(np.argmax(wv)), wv.shape)
        out.update({"wwb": float(wv[bw]), "h_m": float(h[top][bw[1]]), "s": float(ss[bw[0]]), "wrapped": True,
                    "wwb_gauss": float(vals[best])})
    return out if detail else out["wwb"]


def wwb_lattice(s_lat: float, lam_eff: float, T: float | None = None, k: int = 1) -> float:
    """WWB at the grid test point h = kλ_eff, s = ½: J(kλ_eff) = k²J₁ and J₁/σ² = 1/s_lat².

    T = None gives the infinite-window form λ_eff² s² · t / sinh t with t = k²/(4 s²).
    """
    t = k * k / (4.0 * s_lat**2)
    if T is None:
        return lam_eff**2 * s_lat**2 * t / math.sinh(t)
    hk = k * lam_eff
    a, b = (T - hk) / T, max(T - 2 * hk, 0.0) / T
    return hk**2 * a * a * math.exp(-t) / (2.0 * a - 2.0 * b * math.exp(-2.0 * t))


def wrapped_phase_fisher(sigma: float, nterms: int = 50, ngrid: int = 8001) -> float:
    phi = np.linspace(-math.pi, math.pi, ngrid)
    dphi = phi[1] - phi[0]
    ks = np.arange(-nterms, nterms + 1)
    z = phi[None, :] + 2.0 * math.pi * ks[:, None]
    g = np.exp(-0.5 * (z / sigma) ** 2) / (sigma * math.sqrt(2.0 * math.pi))
    dp = np.sum(g * (-z / sigma**2), axis=0)
    p = np.sum(g, axis=0)
    return float(np.sum(dp**2 / np.maximum(p, 1e-300)) * dphi)
