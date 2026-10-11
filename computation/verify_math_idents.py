#!/usr/bin/env python3
"""Check the closed forms used in article/math.zh.md against the library and the archive.

Prints PASS/FAIL per identity and REPORT lines for the derived numbers the note quotes.
Reads computation/data/v2 and article/numbers.json; writes nothing.
"""

from __future__ import annotations

import json
import math
import sys
from concurrent.futures import ProcessPoolExecutor
from functools import reduce
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

sys.path.insert(0, str(Path(__file__).resolve().parent))
from mfpr import bounds, core, design, estimators, lattice, multivalley, sim  # noqa: E402

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


def wrapped_pe(f, h: float, sigma: float, n: int, seed: int, d0: float = 7.0) -> tuple[float, float]:
    """Error probability of the exact wrapped-normal LRT between d0 and d0+h (equal priors), with its s.e."""
    rng = np.random.default_rng(seed)
    g = 2 * np.pi * f / core.C
    mu = (g * d0, g * (d0 + h))

    def loglik(phi, m):
        r = core.wrap(phi - m)[..., None] + 2 * np.pi * np.arange(-3, 4)
        return np.sum(np.log(np.sum(np.exp(-0.5 * (r / sigma) ** 2), axis=-1)), axis=-1)

    err = 0.0
    for t in (0, 1):
        phi = core.wrap(mu[t] + rng.normal(scale=sigma, size=(n, f.size)))
        err += float(np.mean(loglik(phi, mu[1 - t]) > loglik(phi, mu[t])))
    pe = 0.5 * err
    return pe, math.sqrt(pe * (1 - pe) / (2 * n))


def dominates(a: tuple, b: tuple) -> bool:
    """a = (σ_leave ↑, mse_ref ↓, ρ ↑) weakly better everywhere and strictly somewhere."""
    ge = a[0] >= b[0] and a[1] <= b[1] and a[2] >= b[2]
    gt = a[0] > b[0] or a[1] < b[1] or a[2] > b[2]
    return ge and gt


def _ratios(mhz):
    return design.valley_ratios(np.asarray(mhz, dtype=float) * 1e6, T=15.0, step=2e-4)


def _ml_errors(args, n: int = 6000, T: float = 15.0, center: float = 12.345):
    """Coherent-ML range errors with the e09 Monte Carlo set-up (window centre ± T)."""
    f, sigma, seed = args
    rng = np.random.default_rng(seed)
    lam = core.C / float(np.mean(f))
    d0 = center + (rng.random(n) - 0.5) * lam
    grid = estimators.make_grid(center - T, center + T, lam)
    return estimators.coherent_ml(sim.phases(rng, f, d0, sigma, n), f, grid) - d0


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

    print("== §2.1 lemma: wrapped-phase P_e(h) ≥ Q(√J(h)/2σ), so the ZZB above is valid for wrapped data")
    for sig, h, n_pe, seed in ((sj, lam * (1 - 1 / R), 200_000, 11), (0.836, 5.602, 40_000, 12)):
        jh = float(core.phase_cost(f, [h])[0])
        pe, se = wrapped_pe(f, h, sig, n_pe, seed)
        pg = float(Q(math.sqrt(jh) / (2 * sig)))
        check(f"wrapped LRT P_e ≥ Q(√J/2σ) − 3 s.e. at σ={sig:.4f}, h={h:.3f} m", float(pe >= pg - 3 * se), 1.0, 0)
        report(f"σ={sig:.4f}, h={h:.3f} m: max|wrap(gh)| [rad], P_e wrapped ± s.e., Q(√J/2σ), ratio",
               (round(float(np.max(np.abs(core.wrap(2 * np.pi * f * h / core.C)))), 2), (round(pe, 4), round(se, 4)), round(pg, 4), round(pe / pg, 3)))

    print("== §2.5 single-test-point Weiss–Weinstein bound (same prior)")
    s_grid = np.round(np.arange(0.05, 0.96, 0.05), 3)
    w1 = bounds.wwb(surf, sj, wrapped=True, detail=True)
    check("WWB(σ_jump): exact wrapped-normal η vs Gaussian η [dB]", db(w1["wwb"] / w1["wwb_gauss"]), 0.0, 1e-6)
    check("WWB(σ_jump) test point is the 1λ grid valley [m]", w1["h_m"], lam * (1 - 1 / R), 1e-4)
    p_thr = lattice.params(f, sj)
    check("WWB(σ_jump) vs grid closed form h²a²e^{−t}/(2a−2be^{−2t}) [dB]",
          db(w1["wwb"] / bounds.wwb_lattice(p_thr["s"], p_thr["lambda_eff_m"], T)), 0.0, 0.01)
    report("WWB/MMSE at u=1: numeric, closed form T→∞ (1−1/R)·t/sinh t with t=1/(4s²) [dB]",
           (round(db(w1["wwb"] / core.mmse(f, sj)), 3), round(db(bounds.wwb_lattice(p_thr["s"], p_thr["lambda_eff_m"]) / core.mmse(f, sj)), 3)))
    for u in (0.3, 1.0, 3.0):
        wo = bounds.wwb(surf, u * sj, s=s_grid, detail=True)
        check(f"s = 1/2 is optimal over s∈[0.05,0.95] at u={u}", wo["s"], 0.5, 1e-9)
    for u in (3.0, 5.0):
        ww = bounds.wwb(surf, u * sj, wrapped=True, detail=True)
        report(f"u={u}: wrapped − Gaussian WWB [dB], test point [m]", (float(f"{db(ww['wwb'] / ww['wwb_gauss']):.2e}"), round(ww["h_m"], 4)))
    rows, worst_valid = [], -np.inf
    for u in (0.15, 0.2, 0.3, 0.4, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0, 5.0):
        sg = u * sj
        mm, pp = core.mmse(f, sg), lattice.params(f, sg)
        z_, w_ = bounds.zzb(surf, sg), bounds.wwb(surf, sg)
        pm_, ml_ = lattice.pm_mse(pp)["mse"], lattice.ml_mse(pp)["mse"]
        worst_valid = max(worst_valid, db(w_ / pm_))
        rows.append((u, round(db(z_ / mm), 2), round(db(w_ / mm), 2), round(db(pm_ / mm), 2), round(db(ml_ / mm), 2),
                     round(db(pm_ / max(z_, w_)), 2)))
    report("(u, ZZB, WWB, PM_lat, ML_lat [dB re MMSE], PM − max(ZZB,WWB) [dB])", rows)
    check("WWB ≤ PM (lattice) at every u above [max dB]", float(worst_valid <= 0.0), 1.0, 0)
    u_x = brentq(lambda u: db(bounds.wwb(surf, u * sj) / bounds.zzb(surf, u * sj)), 0.5, 0.7, xtol=1e-3)
    report("u where WWB overtakes ZZB", round(u_x, 3))

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

    print("== §4.6 multi-valley linearised model (mfpr/multivalley.py)")
    worst = 0.0
    for u in (0.3, 0.5, 1.0, 2.0):
        worst = max(worst, abs(db(multivalley.mse(f, u * sj, surf, seed=3)["mse"] / lattice.ml_mse(lattice.params(f, u * sj))["mse"])))
    check("equal grid: multi-valley = lattice model, max |dB| over u∈{0.3,0.5,1,2}", worst, 0.0, 0.06)
    p_r = multivalley.mse(f, 0.25 * sj, surf, seed=4)
    check("rare slips: multi-valley = valley sum (u=0.25) [rel]", p_r["mse"], p_r["union"], 0.02, rel=True)
    d07 = sim.load("e07_layers")
    lay = {"geometric": core.geometric_hz(), "cluster": np.asarray(d07["designs"]["cluster"]["f_hz"], dtype=float)}
    for name, fd in lay.items():
        surf_d = bounds.CostSurface.build(fd, 0.5 * core.C / 10e6)
        rows, err_mv, err_mix, inside, n_pts = [], [], [], 0, 0
        for r in sorted((r for r in d07["rows"] if r["design"] == name), key=lambda r: r["sigma"]):
            sg = r["sigma"]
            if sg < 0.05 or sg > 0.45:
                continue
            pred = multivalley.mse(fd, sg, surf_d, n=100_000, seed=5)
            mix = float(bounds.valley_mixture_mse(surf_d, sg)[0])
            lo, hi = multivalley.mc_band_db(pred, r["n"])
            e_mv, e_mix = db(pred["mse"] / r["mse_coh"]), db(mix / r["mse_coh"])
            n_pts += 1
            inside += int(-hi <= e_mv <= -lo)
            if max(-lo, hi) <= 1.5:
                err_mv.append(abs(e_mv))
                err_mix.append(abs(e_mix))
            rows.append((round(sg, 4), round(r["u"], 2), round(db(r["mse_coh"]), 2), round(e_mix, 2), round(e_mv, 2),
                         (round(lo, 2), round(hi, 2))))
        report(f"e07 {name}: (σ, u, MC dB m², mixture−MC, multivalley−MC, MC 95% band of the prediction)", rows)
        report(f"e07 {name}: max |error| where the MC band is within ±1.5 dB: multivalley, mixture [dB]",
               (round(max(err_mv), 2), round(max(err_mix), 2)))
        check(f"e07 {name}: multivalley inside the MC 95% band at every σ∈[0.05,0.45] ({n_pts} points)", inside, n_pts, 0)
        if name == "geometric":
            check("e07 geometric: max |multivalley − MC| on resolved points [dB]", max(err_mv), 0.0, 0.25)
    d09_ = sim.load("e09_pareto")
    rows = []
    for lab, mhz in d09_["picks"].items():
        fd = np.asarray(mhz, dtype=float) * 1e6
        surf_d = bounds.CostSurface.build(fd, 15.0)
        mc_rows = {round(r["sigma"], 4): r for r in d09_["mc"] if r["label"] == lab}
        out = [lab]
        for sg in (0.1828, 0.2342):
            r = mc_rows[sg]
            pred = multivalley.mse(fd, r["sigma"], surf_d, n=100_000, seed=6)
            mix = float(bounds.valley_mixture_mse(surf_d, r["sigma"])[0])
            out.append((sg, round(db(r["mse"]), 2), round(db(mix) - db(r["mse"]), 2), round(db(pred["mse"]) - db(r["mse"]), 2)))
        p_ref = multivalley.mse(fd, 0.2, surf_d, n=100_000, seed=7)
        out.append(("σ_ref=0.2: mixture, multivalley [dB m²]", round(db(float(bounds.valley_mixture_mse(surf_d, 0.2)[0])), 2), round(db(p_ref["mse"]), 2)))
        rows.append(out)
    report("e09 picks: (σ, MC dB m², mixture−MC, multivalley−MC) at the two MC σ around 0.2, then the σ_ref values", rows)
    f_mid = np.asarray(d09_["picks"]["front_mid"], dtype=float) * 1e6
    with ProcessPoolExecutor(max_workers=4) as ex:
        errs = list(ex.map(_ml_errors, [(f_mid, 0.1828, 1000 + i) for i in range(8)]))
    e_all = np.concatenate(errs)
    pred = multivalley.mse(f_mid, 0.1828, bounds.CostSurface.build(f_mid, 15.0), n=200_000, seed=8)
    report("front_mid σ=0.1828 rerun: 8×6000 coherent-ML trials, MSE [dB m²], batch range, archive (6000)",
           (round(db(np.mean(e_all**2)), 2), (round(min(db(np.mean(x**2)) for x in errs), 2), round(max(db(np.mean(x**2)) for x in errs), 2)),
            round(db(next(r["mse"] for r in d09_["mc"] if r["label"] == "front_mid" and abs(r["sigma"] - 0.1828) < 1e-4)), 2)))
    check("front_mid σ=0.1828: multivalley vs 48000-trial rerun [dB]", db(pred["mse"] / np.mean(e_all**2)), 0.0, 0.3)

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
    w_db = db(w1["wwb"])
    report("WWB(σ_jump) [dB m²]", round(w_db, 3))
    report("archive PM − WWB, ML − WWB, MMSE − WWB, WWB − ZZB [dB]",
           [round(bt["db_m2"][k] - w_db, 3) for k in ("mse_pm", "mse_ml", "mmse")] + [round(w_db - bt["db_m2"]["zzb"], 3)])
    report("lattice PM − WWB, ML − WWB [dB]", [round(db(pm / w1["wwb"]), 3), round(db(ml / w1["wwb"]), 3)])
    check("archive PM − WWB at u=1 is below 1 dB", float(bt["db_m2"]["mse_pm"] - w_db < 1.0), 1.0, 0)

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
    i_big = int(np.flatnonzero(~lat)[np.argmin(jv[~lat])])
    report("  …that valley: (δ m, k, J_2 − J) — a grid valley beyond the lemma's J < π²/4",
           (round(float(dv[i_big]), 4), int(k[i_big]), float(f"{float(two_scale(f, dv[i_big:i_big + 1])[0] - jv[i_big]):.1e}")))
    vr_eq = design.valley_ratios(f, T=T - 1e-3, step=1e-4)
    report("equal grid: shallowest NON-grid valley (J < J_2) (ρ_nongrid, δ m, #non-grid valleys)",
           (round(vr_eq["rho_nongrid"], 2), round(vr_eq["delta_nongrid_m"], 3), vr_eq["n_nongrid"]))
    check("equal grid: ρ_nongrid ≥ (R−1)/(16N) (lemma over the whole window, L = T)", float(vr_eq["rho_nongrid"] >= (R - 1) / (16 * N)), 1.0, 0)

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
        vr_d = design.valley_ratios(fd, T=15.0, step=1e-4)
        i_far = int(np.argmin(np.abs(vx - dfar)))
        i_one = int(np.argmin(np.abs(vx - md["lambda_m"])))
        sg = np.geomspace(0.005, 1.0, 1600)
        far_c = Q(np.sqrt(vj[i_far]) / (2 * sg)) * vx[i_far] ** 2
        one_c = Q(np.sqrt(vj[i_one]) / (2 * sg)) * vx[i_one] ** 2
        u_cross = float(sg[np.argmax(far_c > one_c)] / md["sigma_jump"]) if (far_c > one_c).any() else float("nan")
        report(f"design {name}", {
            "hmse_vs_equal_db": round(db(md["hmse_factor"] / m["hmse_factor"]), 4), "u_cross": round(u_cross, 2),
            "S_MHz2": round(md["S"] / 1e12, 1), "R_db": round(md["R_db"], 2), "sigma_jump": round(md["sigma_jump"], 4),
            "sigma_leave": round(sl, 4), "gain_db": round(20 * math.log10(sl / sig_leave["equal"]), 2),
            "x_eff": round(md["sigma_jump"] / sl, 3), "x_lattice": round(foot_x(md["R"]), 3),
            "rho_far": round(rho, 2), "delta_far_m": round(dfar, 3), "u_far": round(u_far, 2),
            "rho_nongrid": round(vr_d["rho_nongrid"], 2), "delta_nongrid_m": round(vr_d["delta_nongrid_m"], 3),
            "L_m": round(vr_d["L_m"], 2), "cap_gain_db": round(20 * math.log10(design.leave_cap(fd) / sig_leave["equal"]), 2),
            "u_ref": round(0.2 / md["sigma_jump"], 2), "slip_owner_at_leave_m": slip_owner(sl), "slip_owner_at_0.2_m": slip_owner(0.2),
            "mse_ref_db": round(db(mse_ref[name]), 2) if name in mse_ref else None,
            "mse_ref_multivalley_db": round(db(multivalley.mse(fd, 0.2, surf_d, n=100_000, seed=7)["mse"]), 2),
            "gcd_MHz": g, "period_m": round(core.C / (g * 1e6), 2), "dmax_MHz": int(np.max(np.diff(np.sort(ints)))),
        })
    check("σ_leave equal (e08) vs bounds.leave_sigma", bounds.leave_sigma(surf), sig_leave["equal"], 2e-3, rel=True)
    check("foot ratio σ_jump/σ_leave of equal = x(R)", sj / sig_leave["equal"], foot_x(R), 2e-3, rel=True)
    check("S_max/S_equal for N/2 tones at each edge = 3(N−1)/(N+1)", (N * B**2 / 4) / S, 3 * (N - 1) / (N + 1), 1e-12, rel=True)
    report("ideal edge stacking σ_jump gain [dB]", round(10 * math.log10(3 * (N - 1) / (N + 1)), 3))

    print("== §6.2 O1 cap: σ_leave ≤ leave_cap (±1 grid valleys alone) ≈ σ_jump/x(R)")
    f_edges = np.asarray(table["edges"], dtype=float) * 1e6
    check("S_max over 8 distinct integer MHz in 2400–2470 = S(4+4 edges)", design.max_spread_mhz2(2400, 2470, 8),
          core.moments(f_edges)["S"] / 1e12, 1e-9, rel=True)
    check("leave_cap(equal) = σ_leave(equal, e08)", design.leave_cap(f), sig_leave["equal"], 2e-3, rel=True)
    cap_e = design.leave_cap(f_edges)
    me = core.moments(f_edges)
    report("4+4 edges: cap σ, cap gain over equal [dB], same from σ_jump/x(R)",
           (round(cap_e, 5), round(20 * math.log10(cap_e / sig_leave["equal"]), 3),
            round(20 * math.log10(me["sigma_jump"] / foot_x(me["R"]) / (sj / foot_x(R))), 3)))
    report("e08 best: σ_leave gain, own cap gain, gap to the global cap [dB]",
           (round(20 * math.log10(sig_leave["best_e08"] / sig_leave["equal"]), 3),
            round(20 * math.log10(design.leave_cap(np.asarray(table["best_e08"], dtype=float) * 1e6) / sig_leave["equal"]), 3),
            round(20 * math.log10(cap_e / sig_leave["best_e08"]), 3)))

    def cap_form(fb):
        s_ub = (8988.0 - 8.0 * (fb - 2435.0) ** 2) * 1e12
        return math.pi * math.sqrt(s_ub) / (fb * 1e6) / foot_x(1 + 8 * (fb * 1e6) ** 2 / s_ub)

    sup_f = max(cap_form(fb) for fb in np.linspace(2425, 2445, 401))
    check("sup over f̄ of σ_jump/x(R) under S ≤ 8988 − 8(f̄−2435)² exceeds the 4+4 value by < 0.01 dB",
          20 * math.log10(sup_f / cap_form(2435.0)), 0.0, 0.01)
    arch = [(r_["mhz"], r_["sigma_leave"]) for r_ in d08["random"]] + [(p_["mhz"], p_["sigma_leave"]) for p_ in d09["cloud"] + d09["anneal"]]
    arch += [(d08["best"]["mhz"], d08["best"]["sigma_leave"])] + [(a_["mhz"], a_["sigma_leave"]) for a_ in d08["anneal"]]
    own = max(sl_ / design.leave_cap(np.asarray(mz, dtype=float) * 1e6) for mz, sl_ in arch)
    check(f"all {len(arch)} archived designs: σ_leave ≤ own cap (max ratio − 1)", own - 1.0, 0.0, 2e-3)
    check("all archived designs: σ_leave ≤ cap of the 4+4 edges", float(max(sl_ for _, sl_ in arch) <= cap_e), 1.0, 0)

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

    print("== §6.6 O3 re-judged on the e09 archive: ρ_far (old) versus ρ_nongrid (new)")
    eq_key = tuple(int(v) for v in eq["mhz"])
    seen, pts = {eq_key}, []
    for p_ in d09["cloud"] + d09["anneal"] + d09["front"]:
        key = tuple(int(v) for v in p_["mhz"])
        if key not in seen:
            seen.add(key)
            pts.append(p_)
    with ProcessPoolExecutor(max_workers=4) as ex:
        vr = list(ex.map(_ratios, [p_["mhz"] for p_ in pts], chunksize=40))
    vr_e = _ratios(eq["mhz"])
    report(f"{len(pts)} distinct archived designs; equal: ρ_far, ρ_nongrid, L [m]",
           (round(vr_e["rho_far"], 1), round(vr_e["rho_nongrid"], 1), round(vr_e["L_m"], 2)))
    obj = lambda p_, r_, key: (p_["sigma_leave"], p_["mse_ref"], r_[key])  # noqa: E731
    old = [(p_, r_) for p_, r_ in zip(pts, vr) if dominates(obj(p_, r_, "rho_far"), obj(eq, vr_e, "rho_far"))]
    new = [(p_, r_) for p_, r_ in zip(pts, vr) if dominates(obj(p_, r_, "rho_nongrid"), obj(eq, vr_e, "rho_nongrid"))]
    report("designs dominating equal in (σ_leave, MSE_ref, ρ_far)", len(old))
    report("  e.g. (MHz, σ_leave, MSE_ref dB, ρ_far @ δ, ρ_nongrid)",
           [(list(map(int, p_["mhz"])), round(p_["sigma_leave"], 5), round(db(p_["mse_ref"]), 2),
             f"{r_['rho_far']:.0f}@{r_['delta_far_m']:.2f}", round(r_["rho_nongrid"], 1)) for p_, r_ in old[:3]])
    check("designs dominating equal in (σ_leave, MSE_ref, ρ_nongrid)", len(new), 0, 0)
    both = [(p_, r_) for p_, r_ in zip(pts, vr) if p_["sigma_leave"] >= eq["sigma_leave"] and r_["rho_nongrid"] >= vr_e["rho_nongrid"]]
    check("designs at least as good as equal on σ_leave and ρ_nongrid alone (so any O2 proxy gives the same verdict)", len(both), 0, 0)
    hi_ng = [(p_, r_) for p_, r_ in zip(pts, vr) if r_["rho_nongrid"] >= vr_e["rho_nongrid"]]
    report("designs with ρ_nongrid ≥ equal's: count, share, best σ_leave gain among them [dB], median B [MHz]",
           (len(hi_ng), round(len(hi_ng) / len(pts), 3),
            round(20 * math.log10(max(p_["sigma_leave"] for p_, _ in hi_ng) / eq["sigma_leave"]), 2),
            float(np.median([max(p_["mhz"]) - min(p_["mhz"]) for p_, _ in hi_ng]))))
    report("designs with ρ_far ≥ equal's", sum(1 for r_ in vr if r_["rho_far"] >= vr_e["rho_far"]))
    beat2 = [(p_, r_) for p_, r_ in zip(pts, vr) if p_["sigma_leave"] > eq["sigma_leave"] and p_["mse_ref"] < eq["mse_ref"]]
    report("designs beating equal on (σ_leave, MSE_ref): count, largest ρ_nongrid among them (@ δ)",
           (len(beat2), max((round(r_["rho_nongrid"], 1), round(r_["delta_nongrid_m"], 2)) for _, r_ in beat2)))
    lb = [r_["L_m"] * 2 * (max(p_["mhz"]) - min(p_["mhz"])) * 1e6 / (core.C * (N - 1)) for p_, r_ in zip(pts, vr)]
    check("Corollary 2 on the archive: L < c(N−1)/(2B) for every non-equal design", float(max(lb) < 1.0 - 1e-9), 1.0, 0)
    report("largest L among the other designs (all bandwidths) [m]", round(max(r_["L_m"] for r_ in vr), 2))
    l_good = [r_["L_m"] for p_, r_ in zip(pts, vr) if p_["sigma_leave"] >= eq["sigma_leave"]]
    check("no design with σ_leave ≥ equal's has L ≥ equal's", float(max(l_good) < vr_e["L_m"]), 1.0, 0)
    report("largest L among designs with σ_leave ≥ equal's [m]", round(max(l_good), 2))
    report("front (σ_leave ↑): (MHz, ρ_far, ρ_nongrid @ δ, L)",
           [(list(map(int, p_["mhz"])), round(r_["rho_far"], 2), f"{r_['rho_nongrid']:.2f}@{r_['delta_nongrid_m']:.2f}", round(r_["L_m"], 2))
            for p_, r_ in sorted(((p_, _ratios(p_["mhz"])) for p_ in d09["front"]), key=lambda t: t[0]["sigma_leave"])])

    print("== deep §12.1: foot per tone from x(R) (not threshold + 13.8 dB)")
    rows = []
    for f0 in (900, 2400, 5800, 24000):
        fd = (f0 + 10.0 * np.arange(8)) * 1e6
        md = core.moments(fd)
        thr = db(1 / (2 * md["sigma_jump"] ** 2))
        xx = foot_x(md["R"])
        rows.append((f0, round(md["R_db"], 1), round(md["sigma_jump"], 4), round(thr, 1), round(xx, 3),
                     round(thr + 20 * math.log10(xx), 1), round(thr + 13.8, 1)))
    report("(f0 MHz, R_dB, σ_jump, centre dB/tone, x(R), foot dB/tone, old centre+13.8)", rows)

    print()
    if FAILS:
        print(f"FAILED: {len(FAILS)} check(s): {FAILS}")
        sys.exit(1)
    print("ALL PASS")


if __name__ == "__main__":
    main()
