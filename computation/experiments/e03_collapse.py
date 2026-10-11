#!/usr/bin/env python3
"""Does the jump depend only on R?  Equal grids with different N and bandwidth.

Plotted as MSE/MMSE against u = σ/σ_jump, every design should follow the same curve
until its own second threshold. The lattice model says the only R dependence is the
carrier-line floor 1/R.
"""

from __future__ import annotations

import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from mfpr import core, estimators, lattice, sim  # noqa: E402

CENTER = 12.345
HALF = 15.0  # capped below the unambiguous range c/Δf of each grid
CASES = [(4, 30e6), (4, 120e6), (8, 20e6), (8, 70e6), (8, 250e6), (16, 70e6), (16, 400e6), (32, 150e6)]
U = np.round(np.geomspace(0.2, 6.0, 17), 4)


def run_one(args):
    n_tone, bw, u, seed = args
    f = core.uniform_grid_hz(n_tone, bw)
    m = core.moments(f)
    sigma = u * m["sigma_jump"]
    rng = np.random.default_rng(seed)
    n = 8000
    half = min(HALF, 0.45 * core.C / (bw / (n_tone - 1)))
    d0 = CENTER + (rng.random(n) - 0.5) * m["lambda_m"]
    grid = estimators.make_grid(CENTER - half, CENTER + half, m["lambda_m"])
    err = estimators.coherent_ml(sim.phases(rng, f, d0, sigma, n), f, grid) - d0
    row = {"N": n_tone, "B_hz": bw, "u": float(u), "sigma": sigma, "n": n, "window_half_m": half,
           "mse": float(np.mean(err**2)), "slip_rate": float(np.mean(np.abs(err) > 0.5 * m["lambda_m"]))}
    if u <= 0.5:
        nl = 1_000_000 if u < 0.3 else 400_000
        d0 = CENTER + (rng.random(nl) - 0.5) * m["lambda_m"]
        g = estimators.make_grid(CENTER - 4 * m["lambda_m"], CENTER + 4 * m["lambda_m"], m["lambda_m"])
        e2 = estimators.coherent_ml(sim.phases(rng, f, d0, sigma, nl), f, g) - d0
        row.update({"n_local": nl, "mse_local": float(np.mean(e2**2)),
                    "slip_local": float(np.mean(np.abs(e2) > 0.5 * m["lambda_m"]))})
    row["lattice_mse"] = lattice.ml_mse(lattice.params(f, sigma))["mse"]
    row["mmse"] = float(core.mmse(f, sigma))
    row["hmse"] = float(core.hmse(f, sigma))
    print(f"N={n_tone:2d} B={bw/1e6:5.0f}  u={u:5.2f}  MSE/MMSE {core.db(row['mse']/row['mmse']):6.2f}  "
          f"lattice {core.db(row['lattice_mse']/row['mmse']):6.2f}", flush=True)
    return row


def main():
    jobs = [(n, b, float(u), 20000 + 101 * i + j) for i, (n, b) in enumerate(CASES) for j, u in enumerate(U)]
    with ProcessPoolExecutor(max_workers=4) as ex:
        rows = list(ex.map(run_one, jobs))
    designs = [{"N": n, "B_hz": b, **core.moments(core.uniform_grid_hz(n, b))} for n, b in CASES]
    sim.save("e03_collapse", {"designs": designs, "rows": rows})


if __name__ == "__main__":
    main()
