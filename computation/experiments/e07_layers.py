#!/usr/bin/env python3
"""When the slip is not one wavelength: the geometric grid and a commensurate two-cluster grid.

Geometric (ratio 1.02): a valley at 6.211 m is as deep as the one-wavelength valley, so
the first slips are metres long and the shelf disappears.

Two clusters whose centres are integer multiples of a 55 MHz step: the synthetic
wavelength c/55 MHz is an exact multiple of both cluster wavelengths, so the cost has a
second, coarser lattice. The prediction is a staircase: carrier line, bandwidth line,
then a second jump governed by the spread inside each cluster.
"""

from __future__ import annotations

import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from mfpr import bounds, core, estimators, lattice, sim  # noqa: E402

CENTER = 12.345
HALF = 0.5 * core.C / 10e6
CLUSTER_MHZ = np.array([2414, 2418, 2422, 2426, 2469, 2473, 2477, 2481], dtype=float)
DESIGNS = {"geometric": core.geometric_hz(), "cluster": CLUSTER_MHZ * 1e6}


def run_one(args):
    name, sigma, seed = args
    f = DESIGNS[name]
    m = core.moments(f)
    rng = np.random.default_rng(seed)
    n = 8000
    d0 = CENTER + (rng.random(n) - 0.5) * m["lambda_m"]
    grid = estimators.make_grid(CENTER - HALF, CENTER + HALF, m["lambda_m"])
    phi = sim.phases(rng, f, d0, sigma, n)
    e_c = estimators.coherent_ml(phi, f, grid) - d0
    e_i = estimators.incoherent_ml(phi, f, grid) - d0
    row = {"design": name, "sigma": sigma, "u": sigma / m["sigma_jump"], "n": n,
           "mse_coh": float(np.mean(e_c**2)), "mse_inc": float(np.mean(e_i**2)),
           "slip": float(np.mean(np.abs(e_c) > 0.5 * m["lambda_m"])),
           "p_err_gt_1m": float(np.mean(np.abs(e_c) > 1.0)),
           "abs_err_q": {q: float(np.quantile(np.abs(e_c), q)) for q in (0.5, 0.9, 0.99, 0.999)},
           "hmse": float(core.hmse(f, sigma)), "mmse": float(core.mmse(f, sigma)),
           "lattice": lattice.ml_mse(lattice.params(f, sigma))["mse"]}
    big = np.abs(e_c) > 0.5 * m["lambda_m"]
    if big.any():
        row["slip_lengths_m"] = np.round(np.sort(np.abs(e_c[big]))[:: max(1, big.sum() // 200)], 3).tolist()
    print(f"{name:9} σ={sigma:.4f} u={row['u']:.3f}  coh {core.db(row['mse_coh']):7.2f}  "
          f"inc {core.db(row['mse_inc']):7.2f}  lattice {core.db(row['lattice']):7.2f}", flush=True)
    return row


def main():
    jobs = []
    seed = 90000
    for name, f in DESIGNS.items():
        for s in np.geomspace(0.01, 3.0, 22):
            seed += 1
            jobs.append((name, float(s), seed))
    with ProcessPoolExecutor(max_workers=4) as ex:
        rows = list(ex.map(run_one, jobs))
    meta = {}
    for name, f in DESIGNS.items():
        m = core.moments(f)
        surf = bounds.CostSurface.build(f, HALF)
        v = surf.valleys[np.argsort(surf.valleys[:, 1])[:8]]
        sig = np.array(sorted(r["sigma"] for r in rows if r["design"] == name))
        meta[name] = {**m, "valleys": v, "zzb": [bounds.zzb(surf, s) for s in sig],
                      "valley_sum": bounds.valley_sum_mse(surf, sig), "sigma": sig,
                      "leave_sigma": bounds.leave_sigma(surf)}
    fc = DESIGNS["cluster"]
    groups = [fc[:4], fc[4:]]
    s_within = float(sum(np.sum((g - g.mean()) ** 2) for g in groups))
    lam_syn = core.C / (groups[1].mean() - groups[0].mean())
    meta["cluster"].update({
        "S_within": s_within,
        "synthetic_wavelength_m": lam_syn,
        "mmse_within_factor": core.C**2 / (4 * np.pi**2 * s_within),
        "J_at_synthetic": float(core.phase_cost(fc, [lam_syn])[0]),
        "sigma_jump2": float(np.sqrt(core.phase_cost(fc, [lam_syn])[0]) / 2),
    })
    sim.save("e07_layers", {"designs": meta, "rows": rows})


if __name__ == "__main__":
    main()
