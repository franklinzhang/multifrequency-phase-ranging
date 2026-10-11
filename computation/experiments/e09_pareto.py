#!/usr/bin/env python3
"""Leaving the carrier line late versus degrading gracefully: the Pareto front.

Objective A: σ_leave, as in e08 (larger is better).
Objective B: predicted MSE at σ_ref = 0.2 rad (about 11 dB per tone, well below every
design's first threshold), from the valley mixture of bounds.valley_mixture_mse, in which
every valley competes with the true lobe (smaller is better). At the equal grid's own
threshold (0.0836 rad) one design dominates on both objectives, so the trade-off only shows
deeper in the transition. Designs on the front are checked by Monte Carlo.
"""

from __future__ import annotations

import math
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from e08_design import BAND, N_TONE, T, fast_valleys, leave  # noqa: E402
from mfpr import core, estimators, sim  # noqa: E402

SIG_REF = 0.2
SIG_EQ = 0.083613
WEIGHTS = [0.0, 0.03, 0.06, 0.1, 0.15, 0.25, 0.4, 0.7, 1.2]


def degraded_mse(f, sigma=SIG_REF):
    m = core.moments(f)
    dv, jv = fast_valleys(f)
    p = 2.0 * core.q_func(np.sqrt(jv) / (2.0 * sigma))
    tot = float(p.sum())
    slip = float(np.sum(p * dv**2)) / (1.0 + tot)
    return m["hmse_factor"] * sigma**2 + slip


def objectives(mhz):
    f = np.asarray(mhz, dtype=float) * 1e6
    return leave(f), degraded_mse(f)


def at_equal_threshold(points):
    for p in points:
        p["mse_at_equal_threshold"] = degraded_mse(np.asarray(p["mhz"], dtype=float) * 1e6, SIG_EQ)


def anneal(args):
    w, seed = args
    rng = np.random.default_rng(seed)
    cur = np.sort(rng.choice(BAND, N_TONE, replace=False))

    def score(x):
        a, b = objectives(x)
        return math.log(a) - w * math.log(b), a, b

    cur_s, _, _ = score(cur)
    best, best_s = cur.copy(), cur_s
    iters = 2000
    for it in range(iters):
        temp = 0.15 * (1 - it / iters) + 1e-3
        cand = cur.copy()
        j = rng.integers(N_TONE)
        cand[j] = np.clip(cand[j] + rng.integers(-4, 5), BAND[0], BAND[-1]) if rng.random() < 0.6 else rng.choice(BAND)
        if len(set(cand.tolist())) < N_TONE:
            continue
        cand = np.sort(cand)
        s, _, _ = score(cand)
        if s > cur_s or rng.random() < math.exp((s - cur_s) / temp):
            cur, cur_s = cand, s
            if s > best_s:
                best, best_s = cand.copy(), s
    a, b = objectives(best)
    return {"w": w, "mhz": best.tolist(), "sigma_leave": a, "mse_ref": b}


def pareto(points):
    pts = sorted(points, key=lambda p: (-p["sigma_leave"], p["mse_ref"]))
    front, best_b = [], math.inf
    for p in pts:
        if p["mse_ref"] < best_b:
            front.append(p)
            best_b = p["mse_ref"]
    return front


def mc_curve(args):
    label, mhz, sigma, seed = args
    f = np.asarray(mhz, dtype=float) * 1e6
    m = core.moments(f)
    rng = np.random.default_rng(seed)
    n = 6000
    center = 12.345
    d0 = center + (rng.random(n) - 0.5) * m["lambda_m"]
    grid = estimators.make_grid(center - T, center + T, m["lambda_m"])
    e = estimators.coherent_ml(sim.phases(rng, f, d0, sigma, n), f, grid) - d0
    return {"label": label, "mhz": mhz, "sigma": sigma, "mse": float(np.mean(e**2)),
            "q99_abs_m": float(np.quantile(np.abs(e), 0.99)), "q90_abs_m": float(np.quantile(np.abs(e), 0.9)),
            "p_gt_half_lambda": float(np.mean(np.abs(e) > 0.5 * m["lambda_m"])),
            "hmse": float(core.hmse(f, sigma))}


def main():
    rng = np.random.default_rng(321)
    randoms = [np.sort(rng.choice(BAND, N_TONE, replace=False)).tolist() for _ in range(3000)]
    with ProcessPoolExecutor(max_workers=4) as ex:
        objs = list(ex.map(objectives, randoms, chunksize=60))
        sa = list(ex.map(anneal, [(w, 500 + i * 7 + k) for i, w in enumerate(WEIGHTS) for k in range(2)]))
    cloud = [{"mhz": r, "sigma_leave": a, "mse_ref": b} for r, (a, b) in zip(randoms, objs)]
    eq = core.EQUAL_MHZ.tolist()
    a, b = objectives(eq)
    ref = {"mhz": eq, "sigma_leave": a, "mse_ref": b}
    front = pareto(cloud + sa + [ref])
    at_equal_threshold(front + [ref])
    for p in front:
        print(f"leave {p['sigma_leave']:.4f}  mse_ref {10*math.log10(p['mse_ref']):7.2f} dB  {p['mhz']}")
    picks = {"equal": eq}
    fr = sorted(front, key=lambda p: p["sigma_leave"])
    for name, q in (("front_lo", 0.0), ("front_mid", 0.5), ("front_hi", 1.0)):
        picks[name] = fr[min(int(round(q * (len(fr) - 1))), len(fr) - 1)]["mhz"]
    sigmas = np.geomspace(0.012, 0.3, 14)
    jobs = [(lab, mhz, float(s), 99000 + 13 * i + j) for i, (lab, mhz) in enumerate(picks.items())
            for j, s in enumerate(sigmas)]
    with ProcessPoolExecutor(max_workers=4) as ex:
        mc = list(ex.map(mc_curve, jobs))
    sim.save("e09_pareto", {"sigma_ref": SIG_REF, "sigma_equal_threshold": SIG_EQ, "cloud": cloud, "anneal": sa, "equal": ref, "front": front,
                            "picks": picks, "mc": mc})


if __name__ == "__main__":
    main()
