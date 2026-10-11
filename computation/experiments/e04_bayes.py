#!/usr/bin/env python3
"""Same prior for bound and estimators: who owns the gap at the threshold?

The true range is uniform on a window of length T = 14.99 m and every estimator searches
that window. Compared: coherent ML, posterior mean (the Bayes-optimal estimator for MSE),
the incoherent estimator, and a validated estimator that keeps the ML tooth only when its
posterior mass is at least 0.9 and otherwise falls back to the incoherent estimate.
"""

from __future__ import annotations

import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from mfpr import bounds, core, estimators, lattice, sim  # noqa: E402

T = 0.5 * core.C / 10e6
LO = 5.0
U = np.round(np.concatenate([np.geomspace(0.15, 12.0, 22), [1.0]]), 4)
CONF = 0.9


def run_one(args):
    u, seed = args
    f = core.equal_hz()
    m = core.moments(f)
    sigma = u * m["sigma_jump"]
    rng = np.random.default_rng(seed)
    n = 6000
    d0 = LO + rng.random(n) * T
    grid = estimators.make_grid(LO, LO + T, m["lambda_m"])
    phi = sim.phases(rng, f, d0, sigma, n)
    post = estimators.posterior_mean(phi, f, grid, w=1.0 / sigma**2)
    inc = estimators.incoherent_ml(phi, f, grid)
    ml, pm, ptop = post["ml"], post["pm"], post["p_top"]
    confident = ptop >= CONF
    hyb = np.where(confident, ml, inc)
    lam = m["lambda_m"]
    slip = np.abs(ml - d0) > 0.5 * lam
    row = {
        "u": float(u), "sigma": sigma, "n": n,
        "mse_ml": float(np.mean((ml - d0) ** 2)),
        "mse_pm": float(np.mean((pm - d0) ** 2)),
        "mse_inc": float(np.mean((inc - d0) ** 2)),
        "mse_hybrid": float(np.mean((hyb - d0) ** 2)),
        "slip_ml": float(np.mean(slip)),
        "confident_rate": float(np.mean(confident)),
        "undetected_slip": float(np.mean(slip & confident)),
        "hmse": float(core.hmse(f, sigma)),
        "mmse": float(core.mmse(f, sigma)),
    }
    lp = lattice.params(f, sigma)
    row["lattice_ml"] = lattice.ml_mse(lp)["mse"]
    row["lattice_pm"] = lattice.pm_mse(lp)["mse"]
    print(f"u={u:6.3f}  ML {core.db(row['mse_ml']):7.2f}  PM {core.db(row['mse_pm']):7.2f}  "
          f"inc {core.db(row['mse_inc']):7.2f}  hyb {core.db(row['mse_hybrid']):7.2f}  "
          f"undetected {row['undetected_slip']:.4f}", flush=True)
    return row


def main():
    jobs = [(float(u), 40000 + 17 * i) for i, u in enumerate(U)]
    with ProcessPoolExecutor(max_workers=4) as ex:
        rows = sorted(ex.map(run_one, jobs), key=lambda r: r["u"])
    surf = bounds.CostSurface.build(core.equal_hz(), T)
    for r in rows:
        r["zzb"] = bounds.zzb(surf, r["sigma"])
    sim.save("e04_bayes", {"T_m": T, "conf": CONF, "prior_var": T**2 / 12, "rows": rows})


if __name__ == "__main__":
    main()
