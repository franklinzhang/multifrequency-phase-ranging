#!/usr/bin/env python3
"""Check the closed forms used in article/math.zh.md against the library and the archive.

Prints PASS/FAIL per identity and REPORT lines for the derived numbers the note quotes.
Reads computation/data/v2 and article/numbers.json; writes nothing.
"""

from __future__ import annotations

import json
import math
import sys
from functools import reduce
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

sys.path.insert(0, str(Path(__file__).resolve().parent))
from mfpr import bounds, core, lattice, sim  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
NUM = json.loads((ROOT / "article" / "numbers.json").read_text())
Q = core.q_func
db = lambda x: float(core.db(x))  # noqa: E731
FAILS: list[str] = []


def check(name: str, got: float, want: float, tol: float, rel: bool = False) -> None:
    err = abs(got - want) / (abs(want) if rel else 1.0)
    ok = err <= tol
    if not ok:
        FAILS.append(name)
    print(f"{'PASS' if ok else 'FAIL'}  {name}: got {got:.6g}, want {want:.6g} (tol {tol:g}{' rel' if rel else ''})")


def report(name: str, value) -> None:
    print(f"REPORT {name}: {value}")


def ek2_ml(s: float, beta: float = 0.0) -> float:
    k, p = lattice.integer_pmf(s, beta)
    return float(np.sum(p * k**2))


def zzb_lattice(f, sigma: float) -> float:
    """HMSE + (λ_eff²/2) Σ_{k≥1} (2k−1) Q(k/2s): valley-filled ZZB on the integer lattice."""
    p = lattice.params(f, sigma)
    k = np.arange(1, 400)
    return p["hmse"] + 0.5 * p["lambda_eff_m"] ** 2 * float(np.sum((2 * k - 1) * Q(k / (2 * p["s"]))))


def foot_x(R: float, c: float = 8.0) -> float:
    return brentq(lambda x: c * x * x * Q(x) - 1.0 / R, 0.3, 30.0)


def two_scale(f, d):
    return core.two_scale_cost(f, d)


def valleys(f, T=15.0, step=1e-4):
    d = np.arange(step, T, step)
    J = core.phase_cost(f, d)
    i = np.flatnonzero((J[1:-1] <= J[:-2]) & (J[1:-1] < J[2:])) + 1
    return d[i], J[i]


def weighted_cost(f, w, d):
    out = np.zeros(d.size)
    for fi, wi in zip(f, w):
        out += wi * core.wrap(2 * np.pi / core.C * d * fi) ** 2
    return out


def far_ratio(mhz, w=None, T=15.0, step=1e-4, lo=0.0):
    """Shallowest valley in (max(c/2B, lo), T), over J(λ). Returns (ρ, δ)."""
    f = np.asarray(mhz, dtype=float) * 1e6
    w = np.ones(f.size) if w is None else np.asarray(w, dtype=float)
    fbar = float(np.sum(w * f) / np.sum(w))
    d = np.arange(max(step, lo), T, step)
    J = weighted_cost(f, w, d)
    i = np.flatnonzero((J[1:-1] <= J[:-2]) & (J[1:-1] < J[2:])) + 1
    i = i[d[i] > core.C / (2 * (f.max() - f.min()))]
    j = i[np.argmin(J[i])]
    jl = float(np.sum(w * core.wrap(2 * np.pi * f / fbar) ** 2))
    return float(J[j] / jl), float(d[j])


def main() -> None:
    f = core.equal_hz()
    m = core.moments(f)
    N, R, sj, lam = m["n"], m["R"], m["sigma_jump"], m["lambda_m"]
    fbar, S, B = m["fbar_hz"], m["S"], float(f.max() - f.min())
    eg = NUM["equal_grid"]

    print("== §1 two local lines and R")
    check("S [MHz^2]", S / 1e12, eg["S_MHz2"], 1e-9, rel=True)
    check("R_dB", m["R_db"], eg["R_db"], 1e-4)
    check("R equal-grid formula 1+12(N-1)/(N+1)(fbar/B)^2", 1 + 12 * (N - 1) / (N + 1) * (fbar / B) ** 2, R, 1e-12, rel=True)
    check("J(λ) = 4 σ_jump^2", float(core.phase_cost(f, [lam])[0]), 4 * sj**2, 1e-12, rel=True)
    check("sqrt(MMSE(σ_jump)) / λ = 1/2", math.sqrt(core.mmse(f, sj)) / lam, 0.5, 1e-12)
    for sig in (0.01, 0.1, 0.5):
        p = lattice.params(f, sig)
        check(f"HMSE + λ_eff² s² = MMSE at σ={sig}", p["hmse"] + p["lambda_eff_m"] ** 2 * p["s"] ** 2, core.mmse(f, sig), 1e-10, rel=True)
    check("λ_eff/λ = 1 − 1/R", lattice.params(f, 0.05)["lambda_eff_m"] / lam, 1 - 1 / R, 1e-10)
    check("s(σ_jump) = 1/(2 sqrt(1−1/R))", lattice.params(f, sj)["s"], 1 / (2 * math.sqrt(1 - 1 / R)), 1e-10, rel=True)

    print("== §2 global bound (valley-filled ZZB, uniform prior of length T)")
    check("∫_0^∞ x Q(x) dx = 1/4  (ZZB → HMSE)", quad(lambda x: x * Q(x), 0, 40)[0], 0.25, 1e-10)
    T = NUM["bounds"]["T_m"]
    check("T = c/(2Δf) for the 10 MHz grid", core.C / (2 * 10e6), T, 1e-3)
    surf = bounds.CostSurface.build(f, T)
    z_thr = bounds.zzb(surf, sj)
    check("ZZB(σ_jump) [dB m^2]", db(z_thr), NUM["bounds"]["equal_zzb_at_0.0836_db"], 0.01)
    pcs = bounds.zzb_pieces(surf, 0.0836)
    shares = sorted([(round(p["valley_m"], 3), p["frac"]) for p in pcs["pieces"] if p["valley_m"]], key=lambda t: t[0])
    for (vm, fr), (vm0, fr0) in zip(shares, NUM["bounds"]["equal_pieces_at_0.0836"]):
        check(f"ZZB share of valley {vm0} m", fr, fr0, 2e-3)
    worst = 0.0
    for u in (0.15, 0.2, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0):
        worst = max(worst, abs(db(zzb_lattice(f, u * sj) / bounds.zzb(surf, u * sj))))
    check("lattice ZZB closed form vs numeric ZZB, max |dB| over u∈[0.15,3]", worst, 0.0, 0.15)
    for s in (2.0, 5.0, 20.0):
        k = np.arange(1, 4000)
        S_s = float(np.sum((2 * k - 1) * Q(k / (2 * s))))
        check(f"Σ(2k−1)Q(k/2s) ≈ 2s² − 2s/√(2π) + 1/6 at s={s}", S_s, 2 * s * s - 2 * s / math.sqrt(2 * math.pi) + 1 / 6, 2e-3, rel=True)
    report("ZZB_lat/MMSE at u=(2,3,5,10) vs 1 − 2/(√(2π)u)",
           [(u, round(zzb_lattice(f, u * sj) / core.mmse(f, u * sj), 3), round(1 - 2 / (math.sqrt(2 * math.pi) * u), 3)) for u in (2, 3, 5, 10)])

    print("== §3 knee: foot / center / overshoot")
    snr_total_c = N / (2 * sj**2)
    check("total SNR at u=1 = (R−1)/(2π²)", snr_total_c, (R - 1) / (2 * math.pi**2), 1e-10, rel=True)
    check("per-tone SNR at u=1 [dB]", db(1 / (2 * sj**2)), eg["threshold_snr_per_tone_db"], 1e-3)
    report("10log10(2π²) [dB]", round(db(2 * math.pi**2), 4))
    s1 = 1 / (2 * math.sqrt(1 - 1 / R))
    check("slip prob at u=1 ≈ 2Q(1)", 2 * Q(1 / (2 * s1)), 2 * Q(1.0), 2e-4)
    over = db(1 / R + (1 - 1 / R) * ek2_ml(s1) / s1**2)
    check("overshoot at u=1 [dB]", over, eg["overshoot_at_threshold_db"], 1e-3)
    check("overshoot series 10log10(8Q(1)+24Q(3)+40Q(5))", db(8 * Q(1) + 24 * Q(3) + 40 * Q(5)), over, 0.01)
    x = foot_x(R)
    check("foot x: 8x²Q(x)=1/R", x, eg["foot_x"], 1e-4)
    check("foot per-tone SNR [dB]", db(x * x / (2 * sj**2)), eg["foot_snr_per_tone_db"], 1e-3)
    off = 20 * math.log10(x) - db(2 * math.pi**2) + db(1 - 1 / R)
    report("foot total SNR − R_dB [dB] (equal grid)", round(off, 3))
    p = lattice.params(f, sj / x)
    report("lattice MSE/HMSE at σ_jump/x (definition gives 2)", round(lattice.ml_mse(p)["mse"] / p["hmse"], 4))
    xz = foot_x(R, 2.0)
    report("ZZB foot x_Z: 2x²Q(x)=1/R", round(xz, 4))
    report("ML foot − ZZB foot [dB SNR]", round(20 * math.log10(x / xz), 3))
    rows = []
    for rdb in (20, 26, 30, 35, 40.53, 45, 50, 60):
        Rr = 10 ** (rdb / 10)
        xx = foot_x(Rr)
        rows.append((rdb, round(xx, 3), round(20 * math.log10(xx) - db(2 * math.pi**2) + db(1 - 1 / Rr), 2), round(20 * math.log10(xx / foot_x(Rr, 2.0)), 2)))
    report("(R_dB, x, foot_total−R_dB, ML−ZZB foot dB)", rows)
    for s in (3.0, 10.0):
        check(f"E[k̂²] → s² + 1/12 at s={s}", ek2_ml(s), s * s + 1 / 12, 1e-6, rel=True)

    print("== §5 accounting at u = 1 (lattice, then archive)")
    p1 = lattice.params(f, sj)
    ml, pm, zl, mm = lattice.ml_mse(p1)["mse"], lattice.pm_mse(p1)["mse"], zzb_lattice(f, sj), core.mmse(f, sj)
    report("lattice ML−MMSE, PM−MMSE, ZZB_lat−MMSE [dB]", [round(db(ml / mm), 3), round(db(pm / mm), 4), round(db(zl / mm), 3)])
    check("PM ≈ MMSE at u=1 (Poisson summation)", db(pm / mm), 0.0, 0.02)
    report("Poisson-summation size 4π s² exp(−2π² s²) at u=1", round(4 * math.pi * s1**2 * math.exp(-2 * math.pi**2 * s1**2), 4))
    bt = NUM["bayes_at_threshold"]
    check("archive ML − PM [dB]", bt["db_m2"]["mse_ml"] - bt["db_m2"]["mse_pm"], bt["ml_minus_pm_db"], 1e-3)
    check("archive PM − ZZB [dB]", bt["db_m2"]["mse_pm"] - bt["db_m2"]["zzb"], bt["pm_minus_zzb_db"], 1e-3)
    report("archive ML − MMSE [dB]", round(bt["db_m2"]["mse_ml"] - bt["db_m2"]["mmse"], 3))
    report("lattice ML − PM, PM − ZZB_lat [dB]", [round(db(ml / pm), 3), round(db(pm / zl), 3)])
    ratios = []
    for u in (0.2, 0.3, 0.5, 1.0, 2.0):
        pp = lattice.params(f, u * sj)
        h = pp["hmse"]
        a, b, c = lattice.ml_mse(pp)["mse"] - h, lattice.pm_mse(pp)["mse"] - h, zzb_lattice(f, u * sj) - h
        ratios.append((u, round(a / c, 3), round(b / c, 3), round(a / b, 3)))
    report("slip-term ratios (u, ML/ZZB, PM/ZZB, ML/PM)", ratios)
    us = np.linspace(0.3, 2.0, 69)
    gap = [db(lattice.ml_mse(lattice.params(f, u * sj))["mse"] / lattice.pm_mse(lattice.params(f, u * sj))["mse"]) for u in us]
    report("lattice ML−PM peak (u, dB)", (round(float(us[int(np.argmax(gap))]), 3), round(max(gap), 3)))

    print("== §6.1 near-zone theorem: J<π²/4 and δ<c/(2Δ_max) ⇒ J = two-scale cost")
    rng = np.random.default_rng(7)
    designs = [np.sort(rng.choice(np.arange(2400, 2471), 8, replace=False)) * 1e6 for _ in range(150)]
    designs += [f, core.matched_hz(), np.array([2400, 2401, 2402, 2403, 2467, 2468, 2469, 2470]) * 1e6]
    worst, n_pts = 0.0, 0
    for fd in designs:
        dmax = float(np.max(np.diff(np.sort(fd))))
        d = np.arange(1e-5, core.C / (2 * dmax), 2e-5)
        J = core.phase_cost(fd, d)
        mask = J < math.pi**2 / 4
        n_pts += int(mask.sum())
        worst = max(worst, float(np.max(np.abs(J[mask] - two_scale(fd, d[mask])))) if mask.any() else 0.0)
    check(f"max |J − two-scale| over {len(designs)} designs, {n_pts} points", worst, 0.0, 1e-9)
    check("equal grid: c/(2Δ_max) equals the prior window T", core.C / (2 * 10e6), T, 1e-3)
    dv, jv = valleys(f, T=T - 1e-3)
    k = np.round(dv / (lam * (1 - 1 / R)))
    lat = jv < math.pi**2 / 4
    check("equal grid: every valley with J<π²/4 sits at kλ_eff [max |δ−kλ_eff| m]",
          float(np.max(np.abs(dv[lat] - k[lat] * lam * (1 - 1 / R)))), 0.0, 2e-4)
    check("equal grid: those valleys have depth k²J(λ)(1−1/R) [max rel err]",
          float(np.max(np.abs(jv[lat] / (k[lat] ** 2 * 4 * sj**2 * (1 - 1 / R)) - 1))), 0.0, 1e-3)
    report("equal grid: guaranteed ρ_far ≥ (R−1)/(16N)", round((R - 1) / (16 * N), 2))
    report("equal grid: smallest J among valleys with J ≥ π²/4 (over J(λ))", round(float(jv[~lat].min() / (4 * sj**2)), 2))

    print("== §6.2 design table (T = 15 m)")
    d08, d09 = sim.load("e08_design"), sim.load("e09_pareto")
    table = {
        "equal": core.EQUAL_MHZ.tolist(), "matched": (2435 + core.MATCHED_OFFSETS_MHZ).tolist(),
        "edges": [2400, 2401, 2402, 2403, 2467, 2468, 2469, 2470], "best_e08": d08["best"]["mhz"],
        "front_lo": d09["picks"]["front_lo"], "front_2419": [2400, 2401, 2403, 2419, 2463, 2468, 2469, 2470],
        "front_mid": d09["picks"]["front_mid"], "front_hi": d09["picks"]["front_hi"],
        "cluster": [2414, 2418, 2422, 2426, 2469, 2473, 2477, 2481],
    }
    sig_leave = {"equal": d08["refs"]["equal"]["sigma_leave"], "matched": d08["refs"]["matched"]["sigma_leave"],
                 "edges": d08["refs"]["edges"]["sigma_leave"], "best_e08": d08["best"]["sigma_leave"]}
    for p_ in d09["front"]:
        for name in ("front_lo", "front_2419", "front_mid", "front_hi"):
            if list(map(int, p_["mhz"])) == list(map(int, table[name])):
                sig_leave[name] = p_["sigma_leave"]
    mse_ref = {"equal": d09["equal"]["mse_ref"]}
    for p_ in d09["front"]:
        for name in ("front_lo", "front_2419", "front_mid", "front_hi"):
            if list(map(int, p_["mhz"])) == list(map(int, table[name])):
                mse_ref[name] = p_["mse_ref"]
    for name, mhz in table.items():
        fd = np.asarray(mhz, dtype=float) * 1e6
        md = core.moments(fd)
        rho, dfar = far_ratio(mhz)
        ints = [int(round(v)) for v in mhz]
        g = reduce(math.gcd, ints)
        u_far = math.sqrt(max(rho - 1.0, 0.0) / (4 * math.log(dfar / md["lambda_m"])))
        surf_d = bounds.CostSurface.build(fd, 15.0)
        vx, vj = surf_d.valleys[:, 0], surf_d.valleys[:, 1]

        def slip_owner(sig):
            c = 2 * Q(np.sqrt(vj) / (2 * sig)) * vx**2
            return round(float(vx[int(np.argmax(c))]), 3)

        sl = sig_leave.get(name, bounds.leave_sigma(surf_d))
        report(f"design {name}", {
            "S_MHz2": round(md["S"] / 1e12, 1), "R_db": round(md["R_db"], 2), "sigma_jump": round(md["sigma_jump"], 4),
            "sigma_leave": round(sl, 4), "gain_db": round(20 * math.log10(sl / sig_leave["equal"]), 2),
            "x_eff": round(md["sigma_jump"] / sl, 3), "x_lattice": round(foot_x(md["R"]), 3),
            "rho_far": round(rho, 2), "delta_far_m": round(dfar, 3), "u_far": round(u_far, 2),
            "u_ref": round(0.2 / md["sigma_jump"], 2), "slip_owner_at_leave_m": slip_owner(sl), "slip_owner_at_0.2_m": slip_owner(0.2),
            "mse_ref_db": round(db(mse_ref[name]), 2) if name in mse_ref else None,
            "gcd_MHz": g, "period_m": round(core.C / (g * 1e6), 2), "dmax_MHz": int(np.max(np.diff(np.sort(ints)))),
        })
    check("σ_leave equal (e08) vs bounds.leave_sigma", bounds.leave_sigma(surf), sig_leave["equal"], 2e-3, rel=True)
    check("foot ratio σ_jump/σ_leave of equal = x(R)", sj / sig_leave["equal"], foot_x(R), 2e-3, rel=True)
    check("S_max/S_equal for N/2 tones at each edge = 3(N−1)/(N+1)", (N * B**2 / 4) / S, 3 * (N - 1) / (N + 1), 1e-12, rel=True)
    report("ideal edge stacking σ_jump gain [dB]", round(10 * math.log10(3 * (N - 1) / (N + 1)), 3))

    print("== §6.3 anchors")
    edges = [2400, 2401, 2402, 2403, 2467, 2468, 2469, 2470]
    scan = [(fa, *far_ratio(edges + [fa])) for fa in range(2405, 2467)]
    best = max(scan, key=lambda t: t[1])
    mid = next(t for t in scan if t[0] == 2435)
    report("4+4 edges + 1 anchor: midpoint 2435 (ρ, δ)", (round(mid[1], 2), round(mid[2], 3)))
    report("4+4 edges + 1 anchor: best (f, ρ, δ)", (best[0], round(best[1], 2), round(best[2], 3)))
    report("4+4 edges + 1 anchor: 2418/2419/2451/2452 ρ", [round(t[1], 1) for t in scan if t[0] in (2418, 2419, 2451, 2452)])
    report("4+4 edges + 2435: shallowest valley in 3.5–5.5 m (ρ, δ)",
           tuple(round(v, 3) for v in far_ratio(edges + [2435], T=5.5, lo=3.5)))
    jl1 = lambda n: 2 * n * (2 * math.pi * 35e6 / 2435e6) ** 2  # noqa: E731
    lam_syn = core.C / 70e6
    orders = (1, 2, 3)
    rho0 = {mo: far_ratio([2400, 2470], [1, 1], T=mo * lam_syn + 0.5, lo=mo * lam_syn - 0.5)[0] for mo in orders}
    report("edges-only ρ0 per synthetic order m=1,2,3 (c/ΔF = %.3f m)" % lam_syn, {mo: round(v, 4) for mo, v in rho0.items()})
    for n in (64, 512):
        for alpha in (0.25, 0.2, 0.5):
            r_num, _ = far_ratio([2400, 2470, 2400 + 70 * alpha], [n, n, 1])
            pred = min(rho0[mo] + (2 * math.pi * abs(mo * alpha - round(mo * alpha))) ** 2 / jl1(n) for mo in orders)
            check(f"stiff-edge ρ ≈ min_m[ρ0_m + (2π‖mα‖)²/J(λ)], n={n}, α={alpha}", r_num, pred, 0.15, rel=True)
    rows = []
    for n in (4, 64, 256, 512):
        row = []
        for K in (0, 1, 2, 3, 4, 8):
            anc = 2400 + 70 * np.arange(1, K + 1) / (K + 1)
            row.append(round(far_ratio(np.concatenate([[2400, 2470], anc]), np.concatenate([[n, n], np.ones(K)]))[0], 3))
        rows.append((n, row))
    report("evenly spaced anchors, ρ for K=(0,1,2,3,4,8) per edge size n", rows)
    rows = []
    for n in (64, 256, 512):
        for K in (4, 16):
            anc = 2400 + 70 * ((np.arange(1, K + 1) * 0.6180339887) % 1)
            r_num = far_ratio(np.concatenate([[2400, 2470], anc]), np.concatenate([[n, n], np.ones(K)]))[0]
            rows.append((n, K, round(r_num, 2), round(K * (2435 / 70) ** 2 / (6 * n), 2)))
    report("generic anchors (n, K, ρ numeric, K(fbar/B)²/(6n))", rows)

    print("== §6.4 Pareto (archive)")
    eq = d09["equal"]
    cloud = d09["cloud"]
    dom = sum(1 for p_ in cloud if p_["sigma_leave"] > eq["sigma_leave"] and p_["mse_ref"] < eq["mse_ref"])
    report("random designs dominating equal (share)", round(dom / len(cloud), 4))
    fr = sorted(d09["front"], key=lambda p_: p_["sigma_leave"])
    check("front_lo dominates equal on both objectives", float(fr[0]["sigma_leave"] > eq["sigma_leave"] and fr[0]["mse_ref"] < eq["mse_ref"]), 1.0, 0)
    report("front σ_leave gains over equal [dB]", [round(20 * math.log10(p_["sigma_leave"] / eq["sigma_leave"]), 2) for p_ in fr])
    mc = {(r["label"], round(r["sigma"], 4)): r for r in d09["mc"]}
    for s in (0.0252, 0.087, 0.3):
        a, b = mc[("equal", s)], mc[("front_lo", s)]
        report(f"MC at σ={s} ({db(1 / (2 * s * s)):.1f} dB/tone): equal, front_lo [dB m²]", (round(db(a["mse"]), 1), round(db(b["mse"]), 1)))

    print()
    if FAILS:
        print(f"FAILED: {len(FAILS)} check(s): {FAILS}")
        sys.exit(1)
    print("ALL PASS")


if __name__ == "__main__":
    main()
