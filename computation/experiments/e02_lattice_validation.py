#!/usr/bin/env python3
"""Validate the two-scale lattice model against coherent ML on the equal and matched grids.

Two Monte Carlo runs per noise level:
  full   search over a 29.98 m window, true range randomized within one wavelength;
  local  search over ±4λ, used only where slips beyond two wavelengths are negligible,
         with up to 2·10⁶ trials so the rare-slip region is resolved.
"""

from __future__ import annotations

import math
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from mfpr import core, estimators, lattice, sim  # noqa: E402

CENTER = 12.345
HALF = 0.5 * core.C / 10e6
DESIGNS = {"equal": core.equal_hz(), "matched": core.matched_hz()}
U = np.unique(np.round(np.concatenate([np.geomspace(0.15, 8.0, 25), [1.0]]), 4))


def stats(err, lam_eff):
    e2 = err**2
    k = np.rint(err / lam_eff).astype(int)
    ks, counts = np.unique(np.clip(k, -8, 8), return_counts=True)
    slips = int(np.sum(np.abs(err) > 0.5 * lam_eff))
    return {
        "n": int(err.size),
        "mse": float(e2.mean()),
        "mse_se": float(e2.std(ddof=1) / math.sqrt(err.size)),
        "slips": slips,
        "slip_rate": slips / err.size,
        "slip_ci": core.wilson(slips, err.size),
        "k_hist": {int(a): int(b) for a, b in zip(ks, counts)},
    }


def run_one(args):
    name, u, seed = args
    f = DESIGNS[name]
    m = core.moments(f)
    sigma = u * m["sigma_jump"]
    lp = lattice.params(f, sigma)
    rng = np.random.default_rng(seed)
    out = {"design": name, "u": float(u), "sigma": sigma}

    n_full = 20000 if u >= 0.35 else 8000
    d0 = CENTER + (rng.random(n_full) - 0.5) * m["lambda_m"]
    grid = estimators.make_grid(CENTER - HALF, CENTER + HALF, m["lambda_m"])
    phi = sim.phases(rng, f, d0, sigma, n_full)
    out["full"] = stats(estimators.coherent_ml(phi, f, grid) - d0, lp["lambda_eff_m"])
    out["full_incoherent"] = stats(estimators.incoherent_ml(phi, f, grid) - d0, lp["lambda_eff_m"])

    if 0.2 <= u <= 0.6:
        n_loc = 2_000_000
        errs = []
        for _ in range(n_loc // 200_000):
            d0 = CENTER + (rng.random(200_000) - 0.5) * m["lambda_m"]
            phi = sim.phases(rng, f, d0, sigma, d0.size)
            g = estimators.make_grid(CENTER - 4 * m["lambda_m"], CENTER + 4 * m["lambda_m"], m["lambda_m"])
            errs.append(estimators.coherent_ml(phi, f, g) - d0)
        out["local"] = stats(np.concatenate(errs), lp["lambda_eff_m"])

    ml = lattice.ml_mse(lp)
    k, pk = lattice.integer_pmf(lp["s"])
    out["lattice"] = {
        "mse": ml["mse"],
        "slip": ml["slip"],
        "pm_mse": lattice.pm_mse(lp)["mse"],
        "pmf": {int(a): float(b) for a, b in zip(k, pk) if abs(a) <= 8},
    }
    out["hmse"] = float(core.hmse(f, sigma))
    out["mmse"] = float(core.mmse(f, sigma))
    print(f"{name:8} u={u:6.3f}  full {core.db(out['full']['mse']):7.2f}  "
          f"lattice {core.db(ml['mse']):7.2f}  slip {out['full']['slip_rate']:.4f}/{ml['slip']:.4f}", flush=True)
    return out


def main():
    jobs = [(name, float(u), 1000 + 37 * i + (0 if name == "equal" else 5000))
            for name in DESIGNS for i, u in enumerate(U)]
    jobs.sort(key=lambda j: -int(0.2 <= j[1] <= 0.6))
    with ProcessPoolExecutor(max_workers=4) as ex:
        rows = list(ex.map(run_one, jobs))
    rows.sort(key=lambda r: (r["design"], r["u"]))
    meta = {name: core.moments(f) for name, f in DESIGNS.items()}
    sim.save("e02_lattice_validation", {"designs": meta, "rows": rows, "window_half_m": HALF})


if __name__ == "__main__":
    main()
