#!/usr/bin/env python3
"""Frequency placement: can eight integer-MHz tones in 2400–2470 MHz leave the carrier line later?

Figure of merit: the leaving noise σ_leave, the smallest σ where the valley sum
HMSE + Σ_k 2Q(√J_k/2σ)δ_k² reaches twice the carrier line, over |δ| < 15 m. A larger S
raises σ_jump, but pushing tones apart also creates far valleys; σ_leave sees both.
"""

from __future__ import annotations

import math
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from mfpr import bounds, core, estimators, lattice, sim  # noqa: E402

T = 15.0
STEP = 5e-4
D = np.arange(STEP, T, STEP)
SIG = np.geomspace(3e-3, 2.0, 700)
BAND = np.arange(2400, 2471)
N_TONE = 8


def fast_valleys(f):
    J = core.phase_cost(f, D)
    i = np.flatnonzero((J[1:-1] <= J[:-2]) & (J[1:-1] < J[2:])) + 1
    a, b, c = J[i - 1], J[i], J[i + 1]
    den = a - 2 * b + c
    off = np.where(den > 0, 0.5 * (a - c) / np.where(den > 0, den, 1), 0.0)
    jmin = b - 0.25 * (a - c) * off
    dmin = D[i] + off * STEP
    lam = core.C / np.mean(f)
    keep = dmin > 0.25 * lam
    return dmin[keep], np.maximum(jmin[keep], 0.0)


def leave(f, factor=2.0):
    m = core.moments(f)
    dv, jv = fast_valleys(f)
    p = 2.0 * core.q_func(np.sqrt(jv)[None, :] / (2.0 * SIG[:, None]))
    ratio = 1.0 + np.sum(p * dv[None, :] ** 2, axis=1) / (m["hmse_factor"] * SIG**2)
    idx = np.flatnonzero(ratio >= factor)
    if idx.size == 0:
        return float(SIG[-1])
    k = int(idx[0])
    if k == 0:
        return float(SIG[0])
    r0, r1 = math.log(ratio[k - 1]), math.log(ratio[k])
    t = (math.log(factor) - r0) / (r1 - r0)
    return float(math.exp(math.log(SIG[k - 1]) + t * (math.log(SIG[k]) - math.log(SIG[k - 1]))))


def describe(mhz):
    f = np.asarray(mhz, dtype=float) * 1e6
    m = core.moments(f)
    dv, jv = fast_valleys(f)
    jl = float(core.phase_cost(f, [m["lambda_m"]])[0])
    far = dv > 2.5 * m["lambda_m"]
    worst = int(np.argmin(np.where(far, jv, np.inf))) if far.any() else None
    return {
        "mhz": [float(x) for x in mhz], "S": m["S"], "R_db": m["R_db"], "sigma_jump": m["sigma_jump"],
        "sigma_leave": leave(f), "J_lambda": jl,
        "far_valley_m": None if worst is None else float(dv[worst]),
        "far_valley_J": None if worst is None else float(jv[worst]),
    }


def anneal(seed, iters=2500):
    rng = np.random.default_rng(seed)
    cur = np.sort(rng.choice(BAND, N_TONE, replace=False))
    cur_v = leave(cur * 1e6)
    best, best_v = cur.copy(), cur_v
    for it in range(iters):
        temp = 0.15 * (1 - it / iters) + 1e-3
        cand = cur.copy()
        j = rng.integers(N_TONE)
        if rng.random() < 0.6:
            cand[j] = np.clip(cand[j] + rng.integers(-4, 5), BAND[0], BAND[-1])
        else:
            cand[j] = rng.choice(BAND)
        if len(set(cand.tolist())) < N_TONE:
            continue
        cand = np.sort(cand)
        v = leave(cand * 1e6)
        if v > cur_v or rng.random() < math.exp((math.log(v) - math.log(cur_v)) / temp):
            cur, cur_v = cand, v
            if v > best_v:
                best, best_v = cand.copy(), v
    return best.tolist(), best_v


def mc_row(args):
    label, mhz, u, seed = args
    f = np.asarray(mhz, dtype=float) * 1e6
    m = core.moments(f)
    sj_eq = core.moments(core.equal_hz())["sigma_jump"]
    sigma = u * sj_eq
    rng = np.random.default_rng(seed)
    n = 8000
    center = 12.345
    d0 = center + (rng.random(n) - 0.5) * m["lambda_m"]
    grid = estimators.make_grid(center - T, center + T, m["lambda_m"])
    e = estimators.coherent_ml(sim.phases(rng, f, d0, sigma, n), f, grid) - d0
    return {"label": label, "u_equal": u, "sigma": sigma, "mse": float(np.mean(e**2)),
            "p_gt_half_lambda": float(np.mean(np.abs(e) > 0.5 * m["lambda_m"])),
            "q99_abs_m": float(np.quantile(np.abs(e), 0.99)),
            "hmse": float(core.hmse(f, sigma)), "mmse": float(core.mmse(f, sigma))}


def main():
    refs = {
        "equal": core.EQUAL_MHZ.tolist(),
        "matched": (2435 + core.MATCHED_OFFSETS_MHZ).tolist(),
        "edges": [2400, 2401, 2402, 2403, 2467, 2468, 2469, 2470],
    }
    rng = np.random.default_rng(123)
    randoms = [np.sort(rng.choice(BAND, N_TONE, replace=False)).tolist() for _ in range(1500)]
    with ProcessPoolExecutor(max_workers=4) as ex:
        rand_desc = list(ex.map(describe, randoms, chunksize=50))
        sa = list(ex.map(anneal, range(8)))
    sa.sort(key=lambda t: -t[1])
    best = sa[0][0]
    out = {"refs": {k: describe(v) for k, v in refs.items()}, "random": rand_desc,
           "anneal": [{"mhz": s, "sigma_leave": v} for s, v in sa], "best": describe(best)}
    for k, v in out["refs"].items():
        print(f"{k:8} σ_jump {v['sigma_jump']:.4f} σ_leave {v['sigma_leave']:.4f}  R {v['R_db']:.2f} dB")
    print("best", best, out["best"]["sigma_leave"], out["best"]["sigma_jump"])

    # Monte Carlo check of equal versus the best design, and the valley-sum / ZZB curves.
    us = np.round(np.geomspace(0.15, 3.0, 14), 4)
    jobs = [(lab, mhz, float(u), 95000 + 31 * i + j) for i, (lab, mhz) in enumerate([("equal", refs["equal"]), ("best", best)])
            for j, u in enumerate(us)]
    with ProcessPoolExecutor(max_workers=4) as ex:
        mc = list(ex.map(mc_row, jobs))
    curves = {}
    for lab, mhz in (("equal", refs["equal"]), ("best", best)):
        f = np.asarray(mhz, dtype=float) * 1e6
        surf = bounds.CostSurface.build(f, T)
        sig = us * core.moments(core.equal_hz())["sigma_jump"]
        curves[lab] = {"sigma": sig, "zzb": [bounds.zzb(surf, s) for s in sig],
                       "valley_sum": bounds.valley_sum_mse(surf, sig),
                       "lattice": [lattice.ml_mse(lattice.params(f, s))["mse"] for s in sig]}
    out["mc"] = mc
    out["curves"] = curves
    sim.save("e08_design", out)


if __name__ == "__main__":
    main()
