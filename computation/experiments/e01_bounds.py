#!/usr/bin/env python3
"""Bounds without Monte Carlo: the valley-filled ZZB, its attribution, and checks.

Replaces the archived zzb_half curve and bound_valleys.json with values computed on a
grid that resolves the main lobe and every valley bottom.
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from mfpr import bounds, core, sim  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
T = 0.5 * core.C / 10e6
DESIGNS = {"equal": core.equal_hz(), "matched": core.matched_hz(), "geometric": core.geometric_hz()}
ARCHIVE_SIGMAS = [0.008361, 0.018395, 0.083613, 0.125420, 0.403086, 0.836133, 2.341174]


def main():
    archive = json.loads((ROOT / "data" / "results.json").read_text())
    sj_eq = core.moments(DESIGNS["equal"])["sigma_jump"]
    sig_dense = sj_eq * np.geomspace(0.05, 30.0, 90)
    out = {"T_m": T, "prior_var_m2": T**2 / 12, "designs": {}}
    for name, f in DESIGNS.items():
        m = core.moments(f)
        surf = bounds.CostSurface.build(f, T)
        arch_sig = np.array(archive["sets"][name]["sigma"])
        arch_zzb = np.array(archive["sets"][name]["zzb_half"]["zzb"])
        new_at_arch = np.array([bounds.zzb(surf, s) for s in arch_sig])
        dense = np.array([bounds.zzb(surf, s) for s in sig_dense])
        pieces = [bounds.zzb_pieces(surf, s) for s in ARCHIVE_SIGMAS]
        snap = bounds.zzb(surf, 0.083613, detail=True)
        keep = np.unique(np.concatenate([np.searchsorted(snap["x"], np.geomspace(0.01, T * 0.999, 600)), [snap["x"].size - 1]]))
        v = surf.valleys
        order = np.argsort(v[:, 1])[:12]
        top = []
        for i in order:
            res = core.wrap((2 * math.pi / core.C) * v[i, 0] * f)
            top.append({"delta_m": float(v[i, 0]), "J": float(v[i, 1]), "max_abs_resid_rad": float(np.max(np.abs(res)))})
        out["designs"][name] = {
            **m,
            "archive_sigma": arch_sig, "archive_zzb": arch_zzb, "zzb_at_archive": new_at_arch,
            "sigma_dense": sig_dense, "zzb_dense": dense,
            "hmse_dense": core.hmse(f, sig_dense), "mmse_dense": core.mmse(f, sig_dense),
            "valley_sum_dense": bounds.valley_sum_mse(surf, sig_dense),
            "leave_sigma_valley_sum": bounds.leave_sigma(surf),
            "pieces": pieces,
            "cume_at_equal_jump": {"x": snap["x"][keep], "frac": snap["cume"][keep] / snap["zzb"]},
            "deepest_valleys": top,
            "J_at_lambda": float(core.phase_cost(f, [m["lambda_m"]])[0]),
        }
        print(name, "ZZB diff vs archive (dB):", np.round(core.db(new_at_arch) - core.db(arch_zzb), 2))
        for p in pieces:
            print(f"   σ={p['sigma']:.4f} ZZB {p['zzb_db']:7.2f} d50 {p['d50_m']:.3f} d90 {p['d90_m']:.3f} ",
                  [(None if q['valley_m'] is None else round(q['valley_m'], 3), round(q['frac'], 3)) for q in p["pieces"][:4]])

    f = DESIGNS["equal"]
    crb = []
    for s in sig_dense:
        i_phi = bounds.wrapped_phase_fisher(float(s))
        crb.append(1.0 / (float(np.sum((2 * math.pi * f / core.C) ** 2)) * i_phi))
    out["designs"]["equal"]["local_crb_dense"] = crb

    d = np.linspace(0.0, 3.0, 30001)
    exact = core.phase_cost(f, d)
    approx = core.two_scale_cost(f, d)
    m = core.moments(f)
    out["two_scale_check"] = {
        "max_abs_diff_within_1m": float(np.max(np.abs(exact - approx)[d <= 1.0])),
        "J_lambda_exact": float(core.phase_cost(f, [m["lambda_m"]])[0]),
        "J_lambda_formula_4sj2": 4 * m["sigma_jump"] ** 2,
        "delta": d[::30], "exact": exact[::30], "approx": approx[::30],
    }
    sim.save("e01_bounds", out)


if __name__ == "__main__":
    main()
