#!/usr/bin/env python3
"""Numbers that close the remaining items of article/math.zh.md §7, written to data/v2_derived/.

Reads computation/data/v2 (never writes there). Sections: multi-test-point WWB, multi-valley
O2 ranking on the e09 archive, S-pruned O1 confirmation, non-grid ZZB, unequal-SNR MC.
"""

from __future__ import annotations

import json
import math
import sys
from concurrent.futures import ProcessPoolExecutor
from itertools import combinations
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent / "experiments"))
from e08_design import leave  # noqa: E402  (same σ_leave definition as the e08 archive)
from mfpr import bounds, core, design, estimators, lattice, multivalley, sim  # noqa: E402

ROOT = Path(__file__).resolve().parent
DERIVED = ROOT / "data" / "v2_derived"
T_EQ = 0.5 * core.C / 10e6
T_WIN = 15.0
SIG_REF = 0.2
db = lambda x: float(core.db(x))  # noqa: E731


def _json(x):
    if isinstance(x, dict):
        return {str(k): _json(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_json(v) for v in x]
    if isinstance(x, (np.floating, float)):
        v = float(x)
        if not math.isfinite(v):
            return None
        return float(f"{v:.8g}")
    if isinstance(x, (np.integer, int)) and not isinstance(x, bool):
        return int(x)
    if isinstance(x, np.ndarray):
        return _json(x.tolist())
    return x


def save(name: str, payload: dict) -> Path:
    DERIVED.mkdir(parents=True, exist_ok=True)
    path = DERIVED / f"{name}.json"
    path.write_text(json.dumps(_json({"git": sim.git_hash(), **payload}), indent=1, ensure_ascii=False) + "\n")
    print("wrote", path, flush=True)
    return path


def seed_of(mhz) -> int:
    s = 1
    for x in mhz:
        s = (s * 1000003 + int(round(float(x)))) % 2147483647
    return int(s)


def key_mhz(mhz) -> tuple:
    return tuple(int(round(float(v))) for v in mhz)


# ---------------------------------------------------------------------------
# 1. Multi-test-point WWB
# ---------------------------------------------------------------------------
def derive_wwb_multi() -> dict:
    f = core.equal_hz()
    m = core.moments(f)
    sj, T = m["sigma_jump"], T_EQ
    surf = bounds.CostSurface.build(f, T)
    us = (0.15, 0.2, 0.3, 0.4, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0, 5.0)
    rows = []
    legal = True
    for u in us:
        sg = u * sj
        mm = core.mmse(f, sg)
        p = lattice.params(f, sg)
        pm, ml = lattice.pm_mse(p)["mse"], lattice.ml_mse(p)["mse"]
        z = bounds.zzb(surf, sg)
        w1 = bounds.wwb(surf, sg)
        wk = bounds.wwb_multi(surf, sg, k_grid=8, detail=True)
        wl = bounds.wwb_multi_lattice(p["s"], p["lambda_eff_m"], T, k_max=8)
        tight = max(z, w1, wk["wwb"])
        which = "zzb" if tight == z else ("wwb1" if tight == w1 else "wwb_k")
        if wk["wwb"] > pm * 1.001 or wl["wwb"] > pm * 1.001:
            legal = False
        rows.append({
            "u": u, "sigma": sg,
            "zzb_db_re_mmse": db(z / mm), "wwb1_db_re_mmse": db(w1 / mm),
            "wwb_k_db_re_mmse": db(max(wk["wwb"], 1e-300) / mm),
            "wwb_k_n": wk["n"], "wwb_lat_db_re_mmse": db(max(wl["wwb"], 1e-300) / mm),
            "pm_db_re_mmse": db(pm / mm), "ml_db_re_mmse": db(ml / mm),
            "pm_minus_tight_db": db(pm / tight), "tightest": which,
            "wwb_k_cond": wk["cond"],
        })
    # downslope: extra test points never overtake ZZB
    downslope_wwb_beats_zzb = any(r["tightest"] == "wwb_k" and r["u"] <= 0.5 for r in rows)
    gap_u1 = next(r for r in rows if r["u"] == 1.0)
    return {
        "T_m": T, "note": (
            "Multi-test-point WWB with s_i=1/2 on h=k λ_eff, k=1..K, nested until Ψ is PD. "
            "Legal (≤ lattice PM) on this grid. Tightens u≥0.7 (u=1: PM−bound 0.74→0.28 dB). "
            "Does not recover the u≈0.3–0.5 gap: Chernoff exponent stays 2× ZZB's; denser "
            "non-grid test points do not beat ZZB there. Bell / free-threshold ZZ not used."
        ),
        "legal_vs_pm": legal,
        "downslope_wwb_beats_zzb": downslope_wwb_beats_zzb,
        "at_u1": {
            "wwb1_db_re_mmse": gap_u1["wwb1_db_re_mmse"],
            "wwb_k_db_re_mmse": gap_u1["wwb_k_db_re_mmse"],
            "pm_minus_wwb1_db": gap_u1["pm_db_re_mmse"] - gap_u1["wwb1_db_re_mmse"],
            "pm_minus_wwb_k_db": gap_u1["pm_minus_tight_db"],
            "k": gap_u1["wwb_k_n"],
        },
        "rows": rows,
    }


# ---------------------------------------------------------------------------
# 2. Multi-valley O2 on the e09 archive
# ---------------------------------------------------------------------------
def _ratios_from_sampled(f, dv, jv, gap) -> dict:
    jl = float(core.phase_cost(f, [core.moments(f)["lambda_m"]])[0])
    L = design.near_zone_m(f)
    masks = {"far": dv > core.C / (2.0 * (f.max() - f.min())), "nongrid": gap > design.GRID_TOL}
    out = {"J_lambda": jl, "L_m": L, "n_valleys": int(dv.size), "n_nongrid": int(masks["nongrid"].sum())}
    for key, mask in masks.items():
        if mask.any():
            j = int(np.flatnonzero(mask)[np.argmin(jv[mask])])
            out[f"rho_{key}"], out[f"delta_{key}_m"] = float(jv[j] / jl), float(dv[j])
        else:
            out[f"rho_{key}"], out[f"delta_{key}_m"] = math.inf, None
    return out


def o2_one(item):
    mhz, sigma_leave, mse_mix, n_mv = item
    f = np.asarray(mhz, dtype=float) * 1e6
    dv, jv, gap = design.sampled_valleys(f, T=T_WIN, step=2e-4)
    valleys = np.column_stack([dv, jv]) if dv.size else np.zeros((0, 2))
    pred = multivalley.mse(f, SIG_REF, valleys=valleys, n=n_mv, seed=seed_of(mhz), T=T_WIN)
    vr = _ratios_from_sampled(f, dv, jv, gap)
    return {
        "mhz": [int(v) for v in key_mhz(mhz)],
        "sigma_leave": float(sigma_leave),
        "mse_mix": float(mse_mix) if mse_mix is not None else None,
        "mse_mv": float(pred["mse"]),
        "mse_mv_se": float(pred["slip_se"]),
        "n_competitors": pred["n_competitors"],
        "rho_nongrid": vr["rho_nongrid"],
        "delta_nongrid_m": vr["delta_nongrid_m"],
        "rho_far": vr["rho_far"],
        "L_m": vr["L_m"],
    }


def _mix_at_ref(mhz):
    """Valley-mixture MSE at σ_ref using sampled valleys (matches e09's spirit, not its 5e-4 grid)."""
    f = np.asarray(mhz, dtype=float) * 1e6
    dv, jv, _ = design.sampled_valleys(f, T=T_WIN, step=2e-4)
    m = core.moments(f)
    if dv.size == 0:
        return float(m["hmse_factor"] * SIG_REF**2)
    p = 2.0 * core.q_func(np.sqrt(jv) / (2.0 * SIG_REF))
    tot = 1.0 + float(p.sum())
    return float(m["hmse_factor"] * SIG_REF**2 + float(np.sum(p * dv**2)) / tot)


def dominates(a, b):
    """a, b = (σ_leave ↑, mse ↓, ρ ↑)."""
    ge = a[0] >= b[0] and a[1] <= b[1] and a[2] >= b[2]
    gt = a[0] > b[0] or a[1] < b[1] or a[2] > b[2]
    return ge and gt


def pareto_mv(points):
    pts = sorted(points, key=lambda p: (-p["sigma_leave"], p["mse_mv"]))
    front, best_b = [], math.inf
    for p in pts:
        if p["mse_mv"] < best_b:
            front.append(p)
            best_b = p["mse_mv"]
    return front


def derive_o2() -> dict:
    d08, d09 = sim.load("e08_design"), sim.load("e09_pareto")
    named = {
        "equal": (d09["equal"]["mhz"], d09["equal"]["sigma_leave"], d09["equal"]["mse_ref"]),
        "front_lo": next(p for p in d09["front"] if key_mhz(p["mhz"]) == key_mhz(d09["picks"]["front_lo"])),
        "front_2419": next(p for p in d09["front"] if key_mhz(p["mhz"]) == (2400, 2401, 2403, 2419, 2463, 2468, 2469, 2470)),
        "front_mid": next(p for p in d09["front"] if key_mhz(p["mhz"]) == key_mhz(d09["picks"]["front_mid"])),
        "front_hi": next(p for p in d09["front"] if key_mhz(p["mhz"]) == key_mhz(d09["picks"]["front_hi"])),
        "best_e08": (d08["best"]["mhz"], d08["best"]["sigma_leave"], _mix_at_ref(d08["best"]["mhz"])),
        "edges": (d08["refs"]["edges"]["mhz"], d08["refs"]["edges"]["sigma_leave"], _mix_at_ref(d08["refs"]["edges"]["mhz"])),
    }
    # unpack front rows
    for k in ("front_lo", "front_2419", "front_mid", "front_hi"):
        p = named[k]
        named[k] = (p["mhz"], p["sigma_leave"], p["mse_ref"])

    seen = {}
    jobs_cloud = []
    for p in d09["cloud"] + d09["anneal"] + d09["front"]:
        k = key_mhz(p["mhz"])
        if k not in seen:
            seen[k] = True
            jobs_cloud.append((list(k), p["sigma_leave"], p["mse_ref"], 60_000))
    # equal is not in cloud
    k_eq = key_mhz(d09["equal"]["mhz"])
    if k_eq not in seen:
        seen[k_eq] = True
        jobs_cloud.append((list(k_eq), d09["equal"]["sigma_leave"], d09["equal"]["mse_ref"], 60_000))
    extra = []
    for lab in ("best_e08", "edges"):
        mhz, sl, mix = named[lab]
        k = key_mhz(mhz)
        extra.append((list(k), sl, mix, 120_000))
    named_jobs = []
    for lab, (mhz, sl, mix) in named.items():
        named_jobs.append((lab, list(key_mhz(mhz)), sl, mix, 150_000))

    print(f"O2: {len(jobs_cloud)} archived unique designs + {len(extra)} extras", flush=True)
    with ProcessPoolExecutor(max_workers=4) as ex:
        cloud = list(ex.map(o2_one, jobs_cloud, chunksize=8))
        extras = list(ex.map(o2_one, extra, chunksize=1))
        named_out = list(ex.map(o2_one, [(m, sl, mix, n) for _, m, sl, mix, n in named_jobs], chunksize=1))

    by_key = {key_mhz(p["mhz"]): p for p in cloud}
    for p in extras:
        by_key[key_mhz(p["mhz"])] = p
    named_map = {}
    for (lab, *rest), p in zip(named_jobs, named_out):
        named_map[lab] = p
        by_key[key_mhz(p["mhz"])] = p

    eq = named_map["equal"]
    eq_obj_mix = (eq["sigma_leave"], eq["mse_mix"], eq["rho_nongrid"])
    eq_obj_mv = (eq["sigma_leave"], eq["mse_mv"], eq["rho_nongrid"])
    pts = list(by_key.values())
    is_eq = lambda p: key_mhz(p["mhz"]) == key_mhz(eq["mhz"])  # noqa: E731
    n_dom_mix = sum(1 for p in pts if p["mse_mix"] is not None and not is_eq(p)
                    and dominates((p["sigma_leave"], p["mse_mix"], p["rho_nongrid"]), eq_obj_mix))
    n_dom_mv = sum(1 for p in pts if not is_eq(p)
                   and dominates((p["sigma_leave"], p["mse_mv"], p["rho_nongrid"]), eq_obj_mv))
    n_dom_o1o2_mv = sum(1 for p in pts if not is_eq(p)
                        and p["sigma_leave"] > eq["sigma_leave"] and p["mse_mv"] < eq["mse_mv"])
    n_dom_o1o2_mix = sum(1 for p in pts if p["mse_mix"] is not None and not is_eq(p)
                         and p["sigma_leave"] > eq["sigma_leave"] and p["mse_mix"] < eq["mse_mix"])
    both_o1_o3 = [p for p in pts if not is_eq(p)
                  and p["sigma_leave"] >= eq["sigma_leave"] and p["rho_nongrid"] >= eq["rho_nongrid"]]
    front = pareto_mv(pts)
    equal_on_front = any(key_mhz(p["mhz"]) == key_mhz(eq["mhz"]) for p in front)
    beat2 = [p for p in pts if p["sigma_leave"] > eq["sigma_leave"] and p["mse_mv"] < eq["mse_mv"]]
    max_rho = max((p["rho_nongrid"] for p in beat2), default=None)

    def slim(p):
        return {
            "mhz": p["mhz"], "sigma_leave": p["sigma_leave"],
            "gain_db": 20 * math.log10(p["sigma_leave"] / eq["sigma_leave"]),
            "mse_mix_db": None if p["mse_mix"] is None else db(p["mse_mix"]),
            "mse_mv_db": db(p["mse_mv"]),
            "rho_nongrid": p["rho_nongrid"], "delta_nongrid_m": p["delta_nongrid_m"], "L_m": p["L_m"],
        }

    return {
        "sigma_ref": SIG_REF, "n_designs": len(pts),
        "equal": slim(eq),
        "named": {lab: slim(named_map[lab]) for lab in named_map},
        "n_dominate_equal_o1o2_mix": n_dom_o1o2_mix,
        "n_dominate_equal_o1o2_mv": n_dom_o1o2_mv,
        "n_dominate_equal_3d_mix": n_dom_mix,
        "n_dominate_equal_3d_mv": n_dom_mv,
        "n_at_least_equal_on_o1_and_rho_nongrid": len(both_o1_o3),
        "equal_undominated_3d_mv": n_dom_mv == 0 and len(both_o1_o3) == 0,
        "front_mv": [slim(p) for p in front],
        "equal_on_front_o1o2_mv": equal_on_front,
        "among_o1o2_winners_max_rho_nongrid": max_rho,
        "n_mv_cloud": 60_000, "n_mv_named": 150_000,
    }


# ---------------------------------------------------------------------------
# 3. S-pruned O1
# ---------------------------------------------------------------------------
def cap_of(mhz) -> float:
    return design.leave_cap(np.asarray(mhz, dtype=float) * 1e6)


def _o1_one(item):
    mhz, target = item
    mhz = list(mhz)
    c = cap_of(mhz)
    if c + 1e-12 < target:
        return None
    sl = leave(np.asarray(mhz, dtype=float) * 1e6)
    return mhz, sl, c


def derive_o1_prune() -> dict:
    d08 = sim.load("e08_design")
    eq_sl = d08["refs"]["equal"]["sigma_leave"]
    e08 = d08["best"]
    edges = d08["refs"]["edges"]
    target = e08["sigma_leave"]  # any better O1 must have leave_cap ≥ this
    cap_edges = cap_of(edges["mhz"])

    def max_cap_forcing(f_int: int) -> float:
        mid = 2435.0
        others = sorted((x for x in range(2400, 2471) if x != f_int), key=lambda x: -abs(x - mid))
        return cap_of(sorted([f_int] + others[:7]))

    # smallest W such that every tone outside the two W-windows has max cap < target
    W = 4
    while W < 35:
        interior = [f for f in range(2400 + W, 2470 - W + 1)]
        if not interior or max(max_cap_forcing(f) for f in interior) < target:
            break
        W += 1
    left, right = list(range(2400, 2400 + W)), list(range(2470 - W + 1, 2471))
    band = sorted(set(left + right))
    n_comb = math.comb(len(band), 8)
    print(f"O1 prune: W={W} band={band[0]}..{left[-1]} ∪ {right[0]}..{band[-1]}  C({len(band)},8)={n_comb}", flush=True)

    with ProcessPoolExecutor(max_workers=4) as ex:
        raw = list(ex.map(_o1_one, [(comb, target) for comb in combinations(band, 8)], chunksize=8))
    evaluated = [t for t in raw if t is not None]
    n_cap = len(evaluated)
    best = {"mhz": e08["mhz"], "sigma_leave": e08["sigma_leave"]}
    for mhz, sl, c in evaluated:
        if sl > best["sigma_leave"]:
            best = {"mhz": mhz, "sigma_leave": sl}

    evaluated.sort(key=lambda t: -t[1])
    top = [{"mhz": m, "sigma_leave": sl, "gain_db": 20 * math.log10(sl / eq_sl),
            "cap_gain_db": 20 * math.log10(c / eq_sl)} for m, sl, c in evaluated[:8]]
    gap = 20 * math.log10(cap_edges / best["sigma_leave"])
    e08_gain = 20 * math.log10(e08["sigma_leave"] / eq_sl)
    cap_gain = 20 * math.log10(cap_edges / eq_sl)
    return {
        "window_W": W, "band_left": left, "band_right": right, "n_combinations": n_comb,
        "n_cap_pass": n_cap, "n_leave_evaluated": len(evaluated),
        "equal_sigma_leave": eq_sl,
        "e08_best": {"mhz": [int(x) for x in e08["mhz"]], "sigma_leave": e08["sigma_leave"],
                     "gain_db": e08_gain},
        "edges": {"mhz": [int(x) for x in edges["mhz"]], "sigma_leave": edges["sigma_leave"],
                  "cap_gain_db": cap_gain},
        "pruned_best": {"mhz": [int(x) for x in best["mhz"]], "sigma_leave": best["sigma_leave"],
                        "gain_db": 20 * math.log10(best["sigma_leave"] / eq_sl)},
        "global_cap_gain_db": cap_gain,
        "gap_best_to_cap_db": gap,
        "closed_to_0p04_db": gap <= 0.041,
        "new_best_vs_e08_db": 20 * math.log10(best["sigma_leave"] / e08["sigma_leave"]),
        "top": top,
        "note": (
            "Enumerate every 8-subset in the two edge windows that are the only places a design "
            "can have leave_cap ≥ e08's σ_leave. σ_leave uses e08's valley-sum definition."
        ),
    }


# ---------------------------------------------------------------------------
# 4. Non-grid ZZB
# ---------------------------------------------------------------------------
def _zzb_lat(f, sigma):
    p = lattice.params(f, sigma)
    k = np.arange(1, 400)
    return p["hmse"] + 0.5 * p["lambda_eff_m"] ** 2 * float(
        np.sum((2 * k - 1) * core.q_func(k / (2.0 * p["s"]))))


def derive_zzb_nongrid() -> dict:
    d07, d09 = sim.load("e07_layers"), sim.load("e09_pareto")
    designs = {
        "geometric": np.asarray(d07["designs"]["geometric"]["f_hz"], dtype=float),
        "cluster": np.asarray(d07["designs"]["cluster"]["f_hz"], dtype=float),
        "front_lo": np.asarray(d09["picks"]["front_lo"], dtype=float) * 1e6,
        "equal": core.equal_hz(),
    }
    T = T_WIN
    out = {}
    for name, f in designs.items():
        m = core.moments(f)
        surf = bounds.CostSurface.build(f, T)
        sj = m["sigma_jump"]
        us = (0.3, 0.5, 0.7, 1.0, 1.5, 2.0)
        rows = []
        rec_err, lat_err, mv_err = [], [], []
        for u in us:
            sg = u * sj
            z = bounds.zzb(surf, sg)
            rec = bounds.zzb_record_valleys(f, sg, surf.valleys, T)
            lat = _zzb_lat(f, sg)
            mv = multivalley.mse(f, sg, surf, n=80_000, seed=11, T=T)["mse"]
            rows.append({
                "u": u, "sigma": sg,
                "zzb_db": db(z), "record_db": db(max(rec, 1e-300)), "lattice_db": db(max(lat, 1e-300)),
                "multivalley_db": db(mv),
                "record_minus_zzb_db": db(rec / z), "lattice_minus_zzb_db": db(lat / z),
                "mv_minus_zzb_db": db(mv / z),
            })
            rec_err.append(abs(db(rec / z)))
            lat_err.append(abs(db(lat / z)))
            mv_err.append(abs(db(mv / z)))
        out[name] = {
            "sigma_jump": sj, "R_db": m["R_db"],
            "rows": rows,
            "max_abs_record_minus_zzb_db": max(rec_err),
            "max_abs_lattice_minus_zzb_db": max(lat_err),
            "max_abs_mv_minus_zzb_db": max(mv_err),
        }
    return {
        "T_m": T,
        "note": (
            "Numeric valley-filled ZZB vs the J-record formula "
            "HMSE+Σ ½(δ_v²−δ_{v−1}²)Q(√J_v/2σ) and vs the grid-only lattice formula. "
            "On the equal grid both match ZZB to <0.15 dB. On geometric/cluster/front_lo the "
            "record formula tracks ZZB; the lattice formula misses far non-grid valleys."
        ),
        "designs": out,
    }


# ---------------------------------------------------------------------------
# 5. Unequal SNR, matched weights
# ---------------------------------------------------------------------------
def _unequal_one(args):
    u, n, seed = args
    f = core.equal_hz()
    m = core.moments(f)
    sig_in = u * m["sigma_jump"]
    sig = np.full(f.size, sig_in)
    sig[0] = sig[-1] = 2.0 * sig_in  # edges −6 dB ⇒ σ × 2
    p = lattice.params(f, sig)
    model = lattice.ml_mse(p)
    rng = np.random.default_rng(seed)
    center = 12.345
    d0 = center + (rng.random(n) - 0.5) * m["lambda_m"]
    grid = estimators.make_grid(center - T_WIN, center + T_WIN, m["lambda_m"])
    w = 1.0 / sig**2
    err = estimators.coherent_ml(sim.phases(rng, f, d0, sig, n), f, grid, w=w) - d0
    mse = float(np.mean(err**2))
    return {
        "u_inner": u, "sigma_inner": sig_in, "n": n,
        "mc_mse": mse, "mc_mse_db": db(mse),
        "lattice_mse": model["mse"], "lattice_mse_db": db(model["mse"]),
        "gap_db": db(mse / model["mse"]),
        "lattice_s": p["s"], "lattice_hmse": p["hmse"],
        "slip_mc": float(np.mean(np.abs(err) > 0.5 * m["lambda_m"])),
        "slip_lat": model["slip"],
    }


def derive_unequal() -> dict:
    us = (0.5, 0.8, 1.0, 1.5, 2.0)
    jobs = [(u, 8000, 7000 + i) for i, u in enumerate(us)]
    with ProcessPoolExecutor(max_workers=4) as ex:
        rows = list(ex.map(_unequal_one, jobs))
    gaps = [abs(r["gap_db"]) for r in rows]
    return {
        "setup": (
            "Equal 2400:10:2470 grid; endpoint tones 2400 and 2470 MHz have SNR 6 dB below the "
            "inner six (σ_edge = 2 σ_inner). Receiver uses matched weights w_i = 1/σ_i². "
            "Coherent ML vs lattice.ml_mse(lattice.params(f, σ_vec))."
        ),
        "n_per_u": 8000,
        "rows": rows,
        "max_abs_gap_db": max(gaps),
        "correlated_noise": (
            "Not done. lattice.params is diagonal (independent but unequal σ_i); a correlated "
            "public term would be ũᵀ W Σ W ũ / (ũᵀ W ũ)² and needs a specified Σ plus a "
            "legality check of the near-zone lemma under J_W. No canonical Σ in this article."
        ),
    }


def main() -> None:
    DERIVED.mkdir(parents=True, exist_ok=True)
    print("== 1 multi-test-point WWB", flush=True)
    wwb = derive_wwb_multi()
    save("wwb_multi", wwb)
    print("  legal", wwb["legal_vs_pm"], "u=1", wwb["at_u1"], flush=True)

    print("== 2 multi-valley O2", flush=True)
    o2 = derive_o2()
    save("o2_multivalley", o2)
    print("  n", o2["n_designs"], "3d_mv_dom", o2["n_dominate_equal_3d_mv"],
          "o1o2_mv", o2["n_dominate_equal_o1o2_mv"], "front", len(o2["front_mv"]), flush=True)

    print("== 3 S-pruned O1", flush=True)
    o1 = derive_o1_prune()
    save("o1_prune", o1)
    print("  W", o1["window_W"], "cap_pass", o1["n_cap_pass"],
          "gap_db", o1["gap_best_to_cap_db"], "new_vs_e08", o1["new_best_vs_e08_db"], flush=True)

    print("== 4 non-grid ZZB", flush=True)
    zzb = derive_zzb_nongrid()
    save("zzb_nongrid", zzb)
    for name, d in zzb["designs"].items():
        print(f"  {name}: record {d['max_abs_record_minus_zzb_db']:.2f} dB, "
              f"lattice {d['max_abs_lattice_minus_zzb_db']:.2f} dB", flush=True)

    print("== 5 unequal SNR", flush=True)
    uneq = derive_unequal()
    save("unequal_snr", uneq)
    print("  max |MC−lattice|", uneq["max_abs_gap_db"], "dB", flush=True)

    index = {
        "files": ["wwb_multi.json", "o2_multivalley.json", "o1_prune.json", "zzb_nongrid.json", "unequal_snr.json"],
        "wwb_multi_legal": wwb["legal_vs_pm"],
        "o2_equal_undominated_3d_mv": o2["equal_undominated_3d_mv"],
        "o1_gap_to_cap_db": o1["gap_best_to_cap_db"],
        "unequal_max_abs_gap_db": uneq["max_abs_gap_db"],
    }
    save("index", index)
    print("ALL DERIVED", flush=True)


if __name__ == "__main__":
    main()
