#!/usr/bin/env python3
"""Every number quoted in the article, recomputed from computation/data/v2 and the library.

Writes article/numbers.json and article/numbers.md. If a number in the article and this
table disagree, the table is right.
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

sys.path.insert(0, str(Path(__file__).resolve().parent))
from mfpr import core, lattice, sim  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "article"
db = core.db


def snr_db(sigma):
    return float(10 * np.log10(core.snr_from_sigma(sigma)))


def equal_grid():
    f = core.equal_hz()
    m = core.moments(f)
    R = m["R"]
    foot = brentq(lambda x: 8 * x * x * core.q_func(x) - 1 / R, 1, 10)
    sj = m["sigma_jump"]
    over = db(lattice.universal_curve(np.array([1.0]), R)[0])
    return {
        "F_MHz2": m["F"] / 1e12, "S_MHz2": m["S"] / 1e12, "R_db": m["R_db"], "sigma_jump_rad": sj,
        "lambda_m": m["lambda_m"], "threshold_snr_per_tone_db": snr_db(sj),
        "threshold_snr_total_db": snr_db(sj) + 10 * math.log10(m["n"]),
        "threshold_snr_total_closed_form_db": 10 * math.log10((R - 1) / (2 * math.pi**2)),
        "envelope_std_at_threshold_over_lambda": math.sqrt(core.mmse(f, sj)) / m["lambda_m"],
        "foot_x": foot, "foot_sigma_over_sigma_jump": 1 / foot,
        "foot_snr_per_tone_db": snr_db(sj / foot),
        "safe_rule_total_db": m["R_db"] - 10 * math.log10(2 * math.pi**2) + 20 * math.log10(foot) + 10 * math.log10(1 - 1 / R),
        "safe_rule_per_tone_db": m["R_db"] - 10 * math.log10(2 * math.pi**2) + 20 * math.log10(foot) + 10 * math.log10(1 - 1 / R) - 10 * math.log10(m["n"]),
        "safe_rule_rdb_plus_1_total_db": m["R_db"] + 1.0,
        "overshoot_at_threshold_db": over,
        "shelf_lambda2_over_12_m2": m["lambda_m"] ** 2 / 12,
    }


def lattice_vs_mc(d02):
    out = {}
    for name in ("equal", "matched"):
        rows = [r for r in d02["rows"] if r["design"] == name]
        full = [abs(db(r["full"]["mse"] / r["lattice"]["mse"])) for r in rows if r["full"]["slips"] >= 50 and r["u"] <= 6]
        loc = [abs(db(r["local"]["mse"] / r["lattice"]["mse"])) for r in rows if "local" in r and r["local"]["slips"] >= 50]
        smallest = min((r["local"]["slip_rate"] for r in rows if "local" in r and r["local"]["slips"] >= 50), default=None)
        beyond = [(r["u"], db(r["full"]["mse"] / r["lattice"]["mse"])) for r in rows if r["u"] > 6]
        out[name] = {"max_abs_db_full_u_le_6": max(full), "max_abs_db_local": max(loc) if loc else None,
                     "smallest_resolved_slip_rate": smallest, "beyond_u6_db": beyond}
    return out


def bounds_numbers(d01):
    e = d01["designs"]["equal"]
    diff = np.abs(np.array(e["zzb_at_archive"]) - np.array(e["archive_zzb"]))
    out = {"T_m": d01["T_m"], "zzb_vs_archive_max_db": float(diff.max()),
           "J_lambda_exact": d01["two_scale_check"]["J_lambda_exact"],
           "J_lambda_formula": d01["two_scale_check"]["J_lambda_formula_4sj2"]}
    for name in ("equal", "geometric"):
        des = d01["designs"][name]
        pcs = [p for p in des["pieces"] if abs(p["sigma"] - 0.083613) < 1e-4]
        if pcs:
            p = pcs[0]
            out[f"{name}_zzb_at_0.0836_db"] = p["zzb_db"]
            out[f"{name}_pieces_at_0.0836"] = [("main lobe" if q["valley_m"] is None else round(q["valley_m"], 3),
                                                round(q["frac"], 3)) for q in p["pieces"][:5]]
            out[f"{name}_d50_d90_m"] = (p["d50_m"], p["d90_m"])
        out[f"{name}_leave_sigma_valley_sum"] = des["leave_sigma_valley_sum"]
    return out


def bayes(d04):
    r = min(d04["rows"], key=lambda r: abs(r["u"] - 1))
    vals = {k: db(r[k]) for k in ("mse_ml", "mse_pm", "mse_inc", "mse_hybrid", "lattice_ml", "lattice_pm", "zzb", "mmse")}
    und = max(d04["rows"], key=lambda r: r["undetected_slip"])
    gap_hyb = [(round(x["u"], 3), db(x["mse_hybrid"] / x["mse_pm"])) for x in d04["rows"] if 0.3 <= x["u"] <= 0.9]
    return {"u": r["u"], "n": r["n"], "db_m2": vals, "ml_minus_pm_db": vals["mse_ml"] - vals["mse_pm"],
            "pm_minus_zzb_db": vals["mse_pm"] - vals["zzb"], "max_undetected_slip": und["undetected_slip"],
            "max_undetected_at_u": und["u"], "hybrid_minus_pm_db": gap_hyb, "conf": d04["conf"]}


def complex_obs(d05):
    rows = d05["rows"]
    ratio = [(r["snr_db"], r["sigma_exact"] / r["sigma_highsnr"]) for r in rows]
    gain = [(r["snr_db"], db(r["phase_ml"]["mse"] / r["complex_ml"]["mse"])) for r in rows]
    near = [g for s, g in gain if 14 <= s <= 34]
    return {"max_sigma_ratio": max(ratio, key=lambda x: x[1]), "phase_over_complex_db": gain,
            "max_abs_gap_db_14_to_34": max(abs(x) for x in near)}


def bias(d06):
    rows = d06["rows"]
    diffs = [abs(db(r["mse"] / r["lattice_mse"])) for r in rows]
    tw = [r for r in rows if r["case"] == "twopath"]
    f = core.equal_hz()
    sj = core.moments(f)["sigma_jump"]

    def at_slip(eps_deg, level):
        def g(s):
            p = lattice.params(f, s, np.full(f.size, math.radians(eps_deg)), w_rx=np.ones(f.size))
            slip = core.q_func((0.5 - p["beta"]) / p["s"]) + core.q_func((0.5 + p["beta"]) / p["s"])
            return math.log(slip) - math.log(level)
        return brentq(g, 0.1 * sj, 3 * sj)

    shift = {lvl: {e: 20 * math.log10(at_slip(0.0, lvl) / at_slip(e, lvl)) for e in (5, 10, 20, 30, 45, 90)}
             for lvl in (1e-4, 1e-2)}
    closed = {e: -20 * math.log10(1 - 2 * e / 360) for e in (5, 10, 20, 30, 45, 90)}
    tol1 = 180 * (1 - 10 ** (-1 / 20))
    return {"max_abs_db_mc_vs_lattice": max(diffs), "n_rows": len(rows),
            "twopath_max_abs_beta": max(abs(r["beta"]) for r in tw),
            "beta_common_10deg": lattice.params(f, sj, np.full(8, math.radians(10)), w_rx=np.ones(8))["beta"],
            "fine_bias_common_10deg_m": lattice.params(f, sj, np.full(8, math.radians(10)), w_rx=np.ones(8))["fine_bias_m"],
            "threshold_shift_db_at_slip_level": shift, "threshold_shift_db_closed_form": closed,
            "eps_for_1db_shift_deg_closed_form": tol1}


def collapse(d03):
    out = []
    for ds in d03["designs"]:
        rows = sorted([r for r in d03["rows"] if r["N"] == ds["N"] and r["B_hz"] == ds["B_hz"]], key=lambda r: r["u"])
        dev = [r for r in rows if r["u"] > 0.5 and db(r["mse"] / r["lattice_mse"]) > 1]
        first = dev[0] if dev else None
        valid = [r for r in rows if r["slip_rate"] * r["n"] >= 50 and (first is None or r["u"] < first["u"])]
        loc = [abs(db(r["mse_local"] / r["lattice_mse"])) for r in rows if "mse_local" in r and r["slip_local"] * r["n_local"] >= 50]
        item = {"N": ds["N"], "B_MHz": ds["B_hz"] / 1e6, "R_db": ds["R_db"],
                "max_abs_db_full": max(abs(db(r["mse"] / r["lattice_mse"])) for r in valid),
                "max_abs_db_local": max(loc) if loc else None}
        if first:
            item.update({"second_threshold_u": first["u"], "second_threshold_total_snr_db":
                         snr_db(first["sigma"]) + 10 * math.log10(ds["N"])})
        out.append(item)
    return out


def layers(d07):
    out = {}
    for name in ("geometric", "cluster"):
        des = d07["designs"][name]
        rows = sorted([r for r in d07["rows"] if r["design"] == name], key=lambda r: r["sigma"])
        left = [r for r in rows if r["mse_coh"] > 2 * r["hmse"]]
        before = [r for r in rows if r["mse_coh"] <= 2 * r["hmse"]]
        out[name] = {"R_db": des["R_db"], "sigma_jump": des["sigma_jump"], "leave_sigma_valley_sum": des["leave_sigma"],
                     "mc_leave_between": (before[-1]["sigma"], left[0]["sigma"]) if left and before else None,
                     "q99_m": [(round(r["sigma"], 4), r["abs_err_q"]["0.99"]) for r in rows]}
    c = d07["designs"]["cluster"]
    R2 = c["S"] / c["S_within"]
    rows = sorted([r for r in d07["rows"] if r["design"] == "cluster"], key=lambda r: r["sigma"])
    cmp = []
    for r in rows:
        pred = c["mmse_within_factor"] * r["sigma"] ** 2 * lattice.universal_curve(np.array([r["sigma"] / c["sigma_jump2"]]), R2)[0]
        cmp.append((round(snr_db(r["sigma"]), 1), round(db(r["mse_coh"] / pred), 2)))
    out["cluster"].update({"S_within_MHz2": c["S_within"] / 1e12, "R2_db": db(R2), "sigma_jump2": c["sigma_jump2"],
                           "synthetic_wavelength_m": c["synthetic_wavelength_m"], "J_at_synthetic": c["J_at_synthetic"],
                           "second_layer_mc_minus_model_db": cmp})
    return out


def design(d08):
    rnd = d08["random"]
    leave = np.array([r["sigma_leave"] for r in rnd])
    ratio = np.array([r["sigma_jump"] / r["sigma_leave"] for r in rnd])
    eq = d08["refs"]["equal"]
    return {"n_random": len(rnd), "leave_pct_5_50_95_100": np.percentile(leave, [5, 50, 95, 100]).tolist(),
            "sigma_jump_over_leave_pct_0_5_50_95_100": np.percentile(ratio, [0, 5, 50, 95, 100]).tolist(),
            "equal": {k: eq[k] for k in ("sigma_jump", "sigma_leave")},
            "edges": {k: d08["refs"]["edges"][k] for k in ("sigma_jump", "sigma_leave")},
            "best": d08["best"], "best_gain_db": 20 * math.log10(d08["best"]["sigma_leave"] / eq["sigma_leave"])}


def pareto(d09):
    eq = d09["equal"]
    front = sorted(d09["front"], key=lambda p: p["sigma_leave"])
    return {"sigma_ref": d09["sigma_ref"], "n_front": len(front),
            "equal": {"sigma_leave": eq["sigma_leave"], "mse_ref_db": db(eq["mse_ref"])},
            "front": [{"mhz": p["mhz"], "sigma_leave": p["sigma_leave"], "mse_ref_db": db(p["mse_ref"]),
                       "mse_at_equal_threshold_db": db(p["mse_at_equal_threshold"]) if "mse_at_equal_threshold" in p else None}
                      for p in front],
            "picks": d09["picks"]}


def main():
    load = {n: sim.load(n) for n in ("e01_bounds", "e02_lattice_validation", "e03_collapse", "e04_bayes", "e05_complex",
                                     "e06_bias", "e07_layers", "e08_design", "e09_pareto")}
    nums = {
        "equal_grid": equal_grid(),
        "lattice_vs_mc": lattice_vs_mc(load["e02_lattice_validation"]),
        "bounds": bounds_numbers(load["e01_bounds"]),
        "bayes_at_threshold": bayes(load["e04_bayes"]),
        "complex": complex_obs(load["e05_complex"]),
        "bias": bias(load["e06_bias"]),
        "collapse": collapse(load["e03_collapse"]),
        "layers": layers(load["e07_layers"]),
        "design": design(load["e08_design"]),
        "pareto": pareto(load["e09_pareto"]),
    }
    OUT.mkdir(exist_ok=True)

    def clean(x):
        if isinstance(x, dict):
            return {str(k): clean(v) for k, v in x.items()}
        if isinstance(x, (list, tuple)):
            return [clean(v) for v in x]
        if isinstance(x, (np.floating, float)):
            return float(f"{float(x):.6g}")
        if isinstance(x, np.integer):
            return int(x)
        return x

    nums = clean(nums)
    (OUT / "numbers.json").write_text(json.dumps(nums, indent=1, ensure_ascii=False) + "\n")
    lines = ["# Numbers quoted in the article", "",
             "Generated by `computation/article_numbers.py` from `computation/data/v2`. Do not edit by hand.", ""]
    for sec, val in nums.items():
        lines += [f"## {sec}", "", "```json", json.dumps(val, indent=1, ensure_ascii=False), "```", ""]
    (OUT / "numbers.md").write_text("\n".join(lines))
    print(json.dumps(nums, indent=1, ensure_ascii=False))


if __name__ == "__main__":
    main()
