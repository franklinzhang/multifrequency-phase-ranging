#!/usr/bin/env python3
"""Figures for the article, in Chinese and English. Reads computation/data/v2/*.json."""

from __future__ import annotations

import math
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import font_manager

sys.path.insert(0, str(Path(__file__).resolve().parent))
from mfpr import bounds, core, lattice, sim  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "article" / "figures"
HALF_E07 = 0.5 * core.C / 10e6
CJK_CANDIDATES = ["Source Han Sans SC", "Noto Sans CJK SC", "PingFang SC", "Microsoft YaHei",
                  "WenQuanYi Micro Hei", "Droid Sans Fallback"]

COL = {"coh": "#1f4e79", "inc": "#c45911", "pm": "#2e7d32", "hyb": "#7b1fa2", "zzb": "#2e7d32",
       "theory": "#d62728", "line": "#333333", "grey": "#888888"}

T = {
    "zh": {
        "snr": "每个频点的信噪比 / dB", "mse": "均方误差 / m²", "u": r"$\sigma\,/\,\sigma_{\mathrm{jump}}$",
        "coh": "相干估计（蒙特卡洛）", "inc": "非相干估计（蒙特卡洛）", "hmse": "载频直线", "mmse": "带宽直线",
        "theory": "格点模型（闭式）", "floor": "搜索窗地板", "jump": "门限", "delta": "距离差 δ / m",
        "cost": "相位代价 J(δ)", "exact": "精确代价", "comb": "载波梳子", "parab": "包络抛物线",
        "lik": "归一化似然", "dist": "距离 / m", "slip": "滑周概率", "ratio": "均方误差 / 带宽直线",
        "ml": "极大似然", "pm": "后验均值", "hyb": "带检验的估计", "zzb": "Ziv–Zakai 下界",
        "lat_ml": "格点模型：极大似然", "lat_pm": "格点模型：后验均值", "prior": "先验方差",
        "undetected": "未被发现的滑周", "confident": "判为可信的比例",
        "complex": "复观测相干估计", "phase": "只用相位的相干估计",
        "eps": "公共相位偏差", "twopath": "反射多径的超额路径 Δ / m", "pred": "格点模型预测",
        "frac": "下界累计占比", "equal": "等间距", "matched": "矩匹配", "geometric": "几何间隔",
        "valley": "谷求和预测", "cluster": "两簇（公度）", "sj": r"$\sigma_{\mathrm{jump}}$ / rad",
        "leave": "离开载频直线的相位噪声 / rad", "random": "随机排布", "best": "搜索得到的最优",
        "edges": "两端各四个", "k": "滑周整数 k", "freq": "频率", "count": "概率",
        "mc_full": "蒙特卡洛（全窗）", "mc_local": "蒙特卡洛（局部窗，$2\\times10^6$ 次）",
        "within": "簇内带宽直线", "q99": "99% 分位误差 / m",
        "beta": "包络与载波的分歧 β / 波长", "rms_cm": "均方根误差 / cm", "layer2": "第二层格点模型",
        "foot": "σ_jump / 4.91（等间距的门限脚）",
        "anneal": "退火搜索", "front": "帕累托前沿", "mse_ref": "σ = {s:g} rad 处的预测均方误差 / dB(m²)",
        "front_lo": "前沿：平台最低", "front_mid": "前沿：折中", "front_hi": "前沿：最晚离开",
    },
    "en": {
        "snr": "per-tone SNR / dB", "mse": "mean-square error / m²", "u": r"$\sigma\,/\,\sigma_{\mathrm{jump}}$",
        "coh": "coherent (Monte Carlo)", "inc": "incoherent (Monte Carlo)", "hmse": "carrier line", "mmse": "bandwidth line",
        "theory": "lattice model (closed form)", "floor": "search-window floor", "jump": "threshold", "delta": "range offset δ / m",
        "cost": "phase cost J(δ)", "exact": "exact cost", "comb": "carrier comb", "parab": "envelope parabola",
        "lik": "normalized likelihood", "dist": "range / m", "slip": "slip probability", "ratio": "MSE / bandwidth line",
        "ml": "maximum likelihood", "pm": "posterior mean", "hyb": "validated estimator", "zzb": "Ziv–Zakai bound",
        "lat_ml": "lattice model: ML", "lat_pm": "lattice model: posterior mean", "prior": "prior variance",
        "undetected": "undetected slips", "confident": "fraction declared reliable",
        "complex": "coherent ML on complex tones", "phase": "coherent ML on phases only",
        "eps": "common phase offset", "twopath": "excess path of the reflection Δ / m", "pred": "lattice prediction",
        "frac": "cumulative share of the bound", "equal": "equal spacing", "matched": "moment matched", "geometric": "geometric",
        "valley": "valley-sum prediction", "cluster": "two clusters (commensurate)", "sj": r"$\sigma_{\mathrm{jump}}$ / rad",
        "leave": "phase noise at which the carrier line is left / rad", "random": "random placements", "best": "best found",
        "edges": "four at each edge", "k": "slip integer k", "freq": "frequency", "count": "probability",
        "mc_full": "Monte Carlo (full window)", "mc_local": "Monte Carlo (local window, $2\\times10^6$ trials)",
        "within": "within-cluster bandwidth line", "q99": "99th-percentile error / m",
        "beta": "envelope–carrier disagreement β / wavelengths", "rms_cm": "RMS error / cm",
        "layer2": "second-layer lattice model", "foot": "σ_jump / 4.91 (foot of the equal grid)",
        "anneal": "annealing", "front": "Pareto front", "mse_ref": "predicted MSE at σ = {s:g} rad / dB(m²)",
        "front_lo": "front: lowest shelf", "front_mid": "front: compromise", "front_hi": "front: latest departure",
    },
}


def setup_fonts():
    names = {f.name for f in font_manager.fontManager.ttflist}
    for cand in CJK_CANDIDATES:
        if cand in names:
            plt.rcParams["font.sans-serif"] = [cand, "DejaVu Sans"]
            break
    plt.rcParams.update({"font.family": "sans-serif", "axes.unicode_minus": False, "font.size": 12,
                         "axes.titlesize": 13, "axes.labelsize": 13, "legend.fontsize": 10.5,
                         "xtick.labelsize": 11, "ytick.labelsize": 11, "axes.grid": True,
                         "grid.color": "#dddddd", "grid.linewidth": 0.7, "axes.axisbelow": True,
                         "savefig.dpi": 200, "figure.dpi": 100})


def fig(w=8.0, h=4.8, ncols=1):
    f, ax = plt.subplots(1, ncols, figsize=(w, h), constrained_layout=True)
    return f, ax


def save(f, lang, name):
    path = OUT / lang / name
    path.parent.mkdir(parents=True, exist_ok=True)
    f.savefig(path)
    plt.close(f)


def u_axis(ax, ticks=(0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10)):
    from matplotlib.ticker import FixedLocator, FuncFormatter, NullFormatter

    lo, hi = ax.get_xlim()
    tk = [x for x in ticks if lo <= x <= hi]
    ax.xaxis.set_major_locator(FixedLocator(tk))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
    ax.xaxis.set_minor_formatter(NullFormatter())


def snr_db_of(sigma):
    return 10 * np.log10(core.snr_from_sigma(np.asarray(sigma)))


_SURF = {}


def mixture(f_hz, half, sigma):
    key = (tuple(np.round(np.asarray(f_hz) / 1e3).astype(int)), round(half, 6))
    if key not in _SURF:
        _SURF[key] = bounds.CostSurface.build(np.asarray(f_hz, dtype=float), half)
    return bounds.valley_mixture_mse(_SURF[key], sigma)


# --------------------------------------------------------------------------- figures


def f01_hero(lang, d02, d01):
    t = T[lang]
    rows = [r for r in d02["rows"] if r["design"] == "equal"]
    m = d02["designs"]["equal"]
    sig = np.array([r["sigma"] for r in rows])
    x = snr_db_of(sig)
    sd = np.geomspace(0.004, 1.2, 400)
    xd = snr_db_of(sd)
    f, ax = fig(8.4, 5.0)
    ax.semilogy(xd, m["hmse_factor"] * sd**2, "--", color=COL["line"], lw=1.6, label=t["hmse"])
    ax.semilogy(xd, m["mmse_factor"] * sd**2, "-.", color=COL["line"], lw=1.6, label=t["mmse"])
    lat = [lattice.ml_mse(lattice.params(core.equal_hz(), s))["mse"] for s in sd]
    keep = sd < 0.45
    ax.semilogy(xd[keep], np.array(lat)[keep], color=COL["theory"], lw=2.0, label=t["theory"])
    ax.semilogy(x, [r.get("local", r["full"])["mse"] for r in rows], "o", color=COL["coh"], ms=6, label=t["coh"])
    ax.semilogy(x, [r["full_incoherent"]["mse"] for r in rows], "s", color=COL["inc"], ms=5, mfc="none",
                label=t["inc"])
    ax.axhline((2 * d02["window_half_m"]) ** 2 / 12, color=COL["grey"], ls=":", lw=1.4, label=t["floor"])
    xj = float(snr_db_of(m["sigma_jump"]))
    ax.axvline(xj, color=COL["inc"], ls="--", lw=1.0)
    ax.annotate(f"{t['jump']}  {xj:.1f} dB", (xj, 3e-8), xytext=(xj + 1.0, 3e-8), fontsize=10.5, color=COL["inc"])
    ax.annotate("", xy=(37.0, m["mmse_factor"] * core.sigma_from_snr(10 ** 3.7) ** 2),
                xytext=(37.0, m["hmse_factor"] * core.sigma_from_snr(10 ** 3.7) ** 2),
                arrowprops=dict(arrowstyle="<->", color=COL["line"], lw=1.2))
    ax.text(36.3, 2e-6, f"R = {m['R_db']:.1f} dB", fontsize=10.5)
    ax.set_xlim(-2, 48)
    ax.set_ylim(1e-9, 300)
    ax.invert_xaxis()
    ax.set_xlabel(t["snr"])
    ax.set_ylabel(t["mse"])
    ax.legend(loc="upper left", framealpha=0.95)
    save(f, lang, "f01_hero.png")


def f02_two_rulers(lang):
    t = T[lang]
    F = core.equal_hz()
    m = core.moments(F)
    d = np.linspace(-0.62, 0.62, 6001)
    exact = core.phase_cost(F, np.abs(d))
    comb = m["n"] * core.wrap(2 * np.pi * d / m["lambda_m"]) ** 2
    par = (2 * np.pi / core.C) ** 2 * m["S"] * d**2
    f, ax = fig(8.4, 4.4)
    ax.plot(d, comb, color="#9ecae1", lw=1.2, label=t["comb"])
    ax.plot(d, exact, color=COL["coh"], lw=2.0, label=t["exact"])
    ax.plot(d, par, color=COL["inc"], lw=2.0, ls="--", label=t["parab"])
    for k in range(-4, 5):
        ax.plot([k * m["lambda_m"]], [par[np.argmin(np.abs(d - k * m["lambda_m"]))]], "o", color=COL["inc"], ms=5)
    ax.set_ylim(0, 1.0)
    ax.set_xlim(d[0], d[-1])
    ax.set_xlabel(t["delta"])
    ax.set_ylabel(t["cost"])
    ax.legend(loc="upper center", ncol=3, framealpha=0.95)
    save(f, lang, "f02_two_rulers.png")


def likelihood_panel(ax, phi, sigma, F, d0, grid_fine, t, title):
    J = np.sum(core.wrap(phi[None, :] - (2 * np.pi / core.C) * grid_fine[:, None] * F[None, :]) ** 2, axis=1)
    L = np.exp(-(J - J.min()) / (2 * sigma**2))
    ax.plot(grid_fine - d0, L, color=COL["coh"], lw=1.0)
    ax.axvline(0, color=COL["inc"], lw=1.0, ls="--")
    ax.set_title(title, fontsize=11.5)
    ax.set_xlim(grid_fine[0] - d0, grid_fine[-1] - d0)
    ax.set_ylim(0, 1.08)


def f03_likelihoods(lang):
    t = T[lang]
    F = core.equal_hz()
    m = core.moments(F)
    rng = np.random.default_rng(11)
    d0 = 12.345
    grid = np.linspace(d0 - 0.75, d0 + 0.75, 12001)
    e = rng.normal(size=F.size)
    f, axes = plt.subplots(3, 1, figsize=(8.4, 6.6), constrained_layout=True, sharex=True)
    for ax, u in zip(axes, (0.15, 1.0, 4.0)):
        s = u * m["sigma_jump"]
        phi = core.wrap((2 * np.pi / core.C) * d0 * F + s * e)
        snr = float(snr_db_of(s))
        likelihood_panel(ax, phi, s, F, d0, grid, t, f"u = {u:g}   ({snr:.1f} dB)")
        ax.set_ylabel(t["lik"])
    axes[-1].set_xlabel(t["delta"])
    save(f, lang, "f03_likelihoods.png")


def f03_gif(lang):
    from matplotlib.animation import FuncAnimation, PillowWriter

    t = T[lang]
    F = core.equal_hz()
    m = core.moments(F)
    rng = np.random.default_rng(11)
    d0 = 12.345
    grid = np.linspace(d0 - 0.75, d0 + 0.75, 6001)
    e = rng.normal(size=F.size)
    us = np.concatenate([np.geomspace(0.1, 6.0, 60), np.geomspace(6.0, 0.1, 60)])
    f, ax = plt.subplots(figsize=(7.2, 3.6), constrained_layout=True)
    line, = ax.plot([], [], color=COL["coh"], lw=1.2)
    ax.axvline(0, color=COL["inc"], lw=1.0, ls="--")
    ax.set_xlim(-0.75, 0.75)
    ax.set_ylim(0, 1.08)
    ax.set_xlabel(t["delta"])
    ax.set_ylabel(t["lik"])
    title = ax.set_title("")

    def update(i):
        s = us[i] * m["sigma_jump"]
        phi = core.wrap((2 * np.pi / core.C) * d0 * F + s * e)
        J = np.sum(core.wrap(phi[None, :] - (2 * np.pi / core.C) * grid[:, None] * F[None, :]) ** 2, axis=1)
        line.set_data(grid - d0, np.exp(-(J - J.min()) / (2 * s**2)))
        title.set_text(f"u = {us[i]:.2f}   ({float(snr_db_of(s)):.1f} dB)")
        return line, title

    anim = FuncAnimation(f, update, frames=len(us), blit=False)
    path = OUT / lang / "f03_likelihood.gif"
    path.parent.mkdir(parents=True, exist_ok=True)
    anim.save(path, writer=PillowWriter(fps=12), dpi=90)
    plt.close(f)


def f04_validation(lang, d02):
    t = T[lang]
    rows = [r for r in d02["rows"] if r["design"] == "equal"]
    R = d02["designs"]["equal"]["R"]
    u = np.array([r["u"] for r in rows])
    ud = np.geomspace(0.15, 8, 300)
    f, (a1, a2) = plt.subplots(1, 2, figsize=(11.5, 4.6), constrained_layout=True)
    slip_th = [lattice.ml_mse(lattice.params(core.equal_hz(), x * d02["designs"]["equal"]["sigma_jump"]))["slip"] for x in ud]
    a1.loglog(ud, slip_th, color=COL["theory"], lw=2, label=t["theory"])
    full = np.array([r["full"]["slip_rate"] for r in rows])
    a1.loglog(u[full > 0], full[full > 0], "o", color=COL["coh"], label=t["mc_full"])
    loc = [r for r in rows if "local" in r]
    ul = np.array([r["u"] for r in loc])
    pl = np.array([r["local"]["slip_rate"] for r in loc])
    ci = np.array([r["local"]["slip_ci"] for r in loc])
    a1.errorbar(ul, pl, yerr=[pl - ci[:, 0], ci[:, 1] - pl], fmt="s", color=COL["inc"], mfc="none", capsize=3,
                label=t["mc_local"])
    a1.set_xlabel(t["u"])
    u_axis(a1)
    a1.set_ylabel(t["slip"])
    a1.set_ylim(1e-7, 1.5)
    a1.legend(loc="lower right")
    a2.semilogx(ud, 10 * np.log10(lattice.universal_curve(ud, R)), color=COL["theory"], lw=2, label=t["theory"])
    mm = np.array([r["mmse"] for r in rows])
    a2.semilogx(u, 10 * np.log10([r["full"]["mse"] for r in rows] / mm), "o", color=COL["coh"], label=t["mc_full"])
    a2.semilogx(ul, 10 * np.log10([r["local"]["mse"] / r["mmse"] for r in loc]), "s", color=COL["inc"], mfc="none",
                label=t["mc_local"])
    a2.axhline(0, color=COL["line"], ls="-.", lw=1.2)
    a2.axhline(-10 * np.log10(R), color=COL["line"], ls="--", lw=1.2)
    a2.set_ylim(-44, 18)
    a2.set_xlabel(t["u"])
    u_axis(a2)
    a2.set_ylabel(t["ratio"] + " / dB")
    a2.legend(loc="upper left")
    save(f, lang, "f04_validation.png")


def f05_kpmf(lang, d02):
    t = T[lang]
    rows = {round(r["u"], 3): r for r in d02["rows"] if r["design"] == "equal"}
    picks = [min(rows, key=lambda x: abs(x - target)) for target in (1.0, 2.5)]
    f, axes = plt.subplots(1, 2, figsize=(10.5, 3.8), constrained_layout=True, sharey=True)
    for ax, key in zip(axes, picks):
        r = rows[key]
        h = {int(k): v for k, v in r["full"]["k_hist"].items()}
        n = r["full"]["n"]
        ks = np.arange(-6, 7)
        ax.bar(ks - 0.18, [h.get(int(k), 0) / n for k in ks], width=0.36, color=COL["coh"], label=t["mc_full"])
        pmf = {int(k): v for k, v in r["lattice"]["pmf"].items()}
        ax.bar(ks + 0.18, [pmf.get(int(k), 0) for k in ks], width=0.36, color=COL["theory"], label=t["theory"])
        ax.set_title(f"u = {r['u']:.2f}")
        ax.set_xlabel(t["k"])
        ax.set_yscale("log")
        ax.set_ylim(1e-4, 1.2)
    axes[0].set_ylabel(t["count"])
    axes[0].legend(loc="upper left")
    save(f, lang, "f05_k_pmf.png")


def f06_collapse(lang, d03):
    t = T[lang]
    designs = d03["designs"]
    f, (a1, a2) = plt.subplots(1, 2, figsize=(11.5, 4.6), constrained_layout=True)
    cmap = plt.get_cmap("viridis")
    for i, dsg in enumerate(designs):
        rows = [r for r in d03["rows"] if r["N"] == dsg["N"] and r["B_hz"] == dsg["B_hz"]]
        u = np.array([r["u"] for r in rows])
        mse = np.array([r.get("mse_local", r["mse"]) for r in rows])
        mm = np.array([r["mmse"] for r in rows])
        c = cmap(i / max(len(designs) - 1, 1))
        lab = f"N={dsg['N']}, B={dsg['B_hz']/1e6:.0f} MHz, R={dsg['R_db']:.1f} dB"
        a1.semilogx(u, 10 * np.log10(mse / mm), "o", color=c, ms=4.5, label=lab)
        ud = np.geomspace(0.2, 6, 200)
        a1.semilogx(ud, 10 * np.log10(lattice.universal_curve(ud, dsg["R"])), color=c, lw=0.9, alpha=0.8)
        sig = np.array([r["sigma"] for r in rows])
        a2.semilogy(snr_db_of(sig), mse, "o-", color=c, ms=4, lw=1.0)
    a1.set_ylim(-55, 15)
    a1.set_xlabel(t["u"])
    u_axis(a1)
    a1.set_ylabel(t["ratio"] + " / dB")
    a1.legend(fontsize=8.5, loc="lower right")
    a2.invert_xaxis()
    a2.set_xlabel(t["snr"])
    a2.set_ylabel(t["mse"])
    save(f, lang, "f06_collapse.png")


def f07_bayes(lang, d04):
    t = T[lang]
    rows = d04["rows"]
    u = np.array([r["u"] for r in rows])
    f, (a1, a2) = plt.subplots(1, 2, figsize=(11.5, 4.6), constrained_layout=True, gridspec_kw={"width_ratios": [1.6, 1]})
    mm = np.array([r["mmse"] for r in rows])
    expected = np.array([r["n"] * lattice.ml_mse(lattice.params(core.equal_hz(), r["sigma"]))["slip"] for r in rows])
    ok = expected >= 3
    for key, lab, sty in (("mse_ml", t["ml"], dict(marker="o", color=COL["coh"])),
                          ("mse_pm", t["pm"], dict(marker="^", color=COL["pm"])),
                          ("mse_hybrid", t["hyb"], dict(marker="D", color=COL["hyb"], mfc="none")),
                          ("mse_inc", t["inc"], dict(marker="s", color=COL["inc"], mfc="none"))):
        a1.semilogx(u[ok], 10 * np.log10(np.array([r[key] for r in rows])[ok] / mm[ok]), ls="none", ms=5, label=lab,
                    **sty)
    a1.semilogx(u, 10 * np.log10(np.array([r["lattice_ml"] for r in rows]) / mm), color=COL["coh"], lw=1.2, label=t["lat_ml"])
    a1.semilogx(u, 10 * np.log10(np.array([r["lattice_pm"] for r in rows]) / mm), color=COL["pm"], lw=1.2, ls="--",
                label=t["lat_pm"])
    a1.semilogx(u, 10 * np.log10(np.array([r["zzb"] for r in rows]) / mm), color=COL["zzb"], lw=2.4, alpha=0.5, label=t["zzb"])
    a1.semilogx(u, 10 * np.log10(np.array([r["hmse"] for r in rows]) / mm), color=COL["line"], ls="--", lw=1)
    a1.axhline(0, color=COL["line"], ls="-.", lw=1)
    a1.set_ylim(-45, 10)
    a1.set_xlim(0.14, 13)
    a1.set_xlabel(t["u"])
    u_axis(a1)
    a1.set_ylabel(t["ratio"] + " / dB")
    a1.legend(fontsize=9, loc="lower right", ncol=2)
    a2.loglog(u, [max(r["slip_ml"], 1e-5) for r in rows], "o-", color=COL["coh"], ms=4, label=t["slip"] + " (ML)")
    und = np.array([r["undetected_slip"] for r in rows])
    a2.loglog(u[und > 0], und[und > 0], "D-", color=COL["hyb"], ms=4, label=t["undetected"])
    a2.loglog(u, [1 - r["confident_rate"] + 1e-5 for r in rows], "s--", color=COL["grey"], ms=4,
              label=("1 − " + t["confident"]))
    a2.set_ylim(1e-4, 1.5)
    a2.set_xlim(0.14, 13)
    a2.set_xlabel(t["u"])
    u_axis(a2)
    a2.legend(fontsize=9, loc="lower right")
    save(f, lang, "f07_bayes.png")


def f08_complex(lang, d05):
    t = T[lang]
    rows = d05["rows"]
    x = np.array([r["snr_db"] for r in rows])
    f, ax = fig(8.4, 4.8)
    ax.semilogy(x, [r["complex_ml"]["mse"] for r in rows], "o", color=COL["coh"], label=t["complex"])
    ax.semilogy(x, [r["phase_ml"]["mse"] for r in rows], "^", color=COL["pm"], mfc="none", label=t["phase"])
    ax.semilogy(x, [r["phase_inc"]["mse"] for r in rows], "s", color=COL["inc"], mfc="none", label=t["inc"])
    xd = np.linspace(x.min(), x.max(), 400)
    sd = core.sigma_from_snr(10 ** (xd / 10))
    ax.semilogy(xd, [lattice.ml_mse(lattice.params(core.equal_hz(), s))["mse"] for s in sd], color=COL["theory"], lw=1.8,
                label=t["theory"])
    ax.semilogy(x, [r["hmse"] for r in rows], "--", color=COL["line"], lw=1.2, label=t["hmse"])
    ax.semilogy(x, [r["mmse"] for r in rows], "-.", color=COL["line"], lw=1.2, label=t["mmse"])
    ax.axhline(d05["floor_m2"], color=COL["grey"], ls=":", lw=1.3, label=t["floor"])
    ax.axvline(d05["snr_threshold_per_tone_db"], color=COL["inc"], ls="--", lw=1)
    ax.invert_xaxis()
    ax.set_ylim(1e-9, 300)
    ax.set_xlabel(t["snr"])
    ax.set_ylabel(t["mse"])
    ax.legend(loc="upper left", fontsize=9.5)
    save(f, lang, "f08_complex.png")


def f09_bias(lang, d06):
    t = T[lang]
    rows = d06["rows"]
    f, (a1, a2) = plt.subplots(1, 2, figsize=(11.5, 4.5), constrained_layout=True)
    cmap = plt.get_cmap("plasma")
    eps_list = sorted({r["eps_deg"] for r in rows if r["case"] == "common"})
    sj = core.moments(core.equal_hz())["sigma_jump"]
    ud = np.geomspace(0.28, 1.6, 100)
    for i, eps in enumerate(eps_list):
        c = cmap(i / (len(eps_list) + 0.5))
        rr = sorted([r for r in rows if r["case"] == "common" and r["eps_deg"] == eps], key=lambda r: r["u"])
        a1.loglog([r["u"] for r in rr], [max(r["slip"], 1e-5) for r in rr], "o", color=c, ms=5)
        th = [lattice.ml_mse(lattice.params(core.equal_hz(), x * sj, np.full(8, math.radians(eps)), w_rx=np.ones(8)))["slip"]
              for x in ud]
        a1.loglog(ud, th, color=c, lw=1.4, label=f"{t['eps']} {eps:g}°")
    a1.set_xlabel(t["u"])
    u_axis(a1)
    a1.set_ylabel(t["slip"])
    a1.set_ylim(5e-4, 1)
    a1.legend(fontsize=9, loc="lower right")
    rr = sorted([r for r in rows if r["case"] == "twopath" and abs(r["snr_db"] - d06["snr_u05_db"]) < 1e-9],
                key=lambda r: r["delta_m"])
    dm = np.array([r["delta_m"] for r in rr])
    a3 = a2.twinx()
    a3.bar(dm, [r["beta"] for r in rr], width=0.06, color="#cccccc", zorder=0)
    a3.axhline(0, color="#999999", lw=0.8)
    a3.set_ylabel(t["beta"], color="#777777")
    a3.set_ylim(-4.5, 4.5)
    a3.grid(False)
    a2.set_zorder(a3.get_zorder() + 1)
    a2.patch.set_visible(False)
    a2.semilogy(dm, 100 * np.sqrt([r["mse"] for r in rr]), "o", color=COL["coh"], mfc="none", ms=6, label="MC")
    a2.semilogy(dm, 100 * np.sqrt([r["lattice_mse"] for r in rr]), color=COL["theory"], lw=1.4, label=t["pred"])
    a2.axhline(100 * math.sqrt(rr[0]["lattice_mse"]), color=COL["grey"], ls=":", lw=1.2)
    a2.set_xlabel(t["twopath"])
    a2.set_ylabel(t["rms_cm"])
    a2.set_title(f"a = {rr[0]['a']:g},  u = {rr[0]['u_nominal']:.2f}", fontsize=11.5)
    a2.legend(fontsize=9, loc="upper left")
    save(f, lang, "f09_bias.png")


def f10_geometric(lang, d07, d01):
    t = T[lang]
    g = d07["designs"]["geometric"]
    rows = sorted([r for r in d07["rows"] if r["design"] == "geometric"], key=lambda r: r["sigma"])
    sig = np.array([r["sigma"] for r in rows])
    u = sig / g["sigma_jump"]
    f, (a1, a2) = plt.subplots(1, 2, figsize=(11.5, 4.6), constrained_layout=True)
    a1.loglog(u, [r["mse_coh"] for r in rows], "o", color=COL["coh"], label=t["coh"])
    a1.loglog(u, [r["mse_inc"] for r in rows], "s", color=COL["inc"], mfc="none", label=t["inc"])
    a1.loglog(u, [r["lattice"] for r in rows], color=COL["theory"], lw=1.4, label=t["theory"])
    sd = np.geomspace(sig[0], sig[-1], 300)
    a1.loglog(sd / g["sigma_jump"], mixture(g["f_hz"], HALF_E07, sd), color=COL["hyb"], lw=1.8, ls="--",
              label=t["valley"])
    a1.loglog(np.array(g["sigma"]) / g["sigma_jump"], g["zzb"], color=COL["zzb"], lw=2.2, alpha=0.6, label=t["zzb"])
    a1.axvline(g["leave_sigma"] / g["sigma_jump"], color=COL["grey"], ls=":", lw=1.2)
    a1.loglog(u, [r["hmse"] for r in rows], "--", color=COL["line"], lw=1)
    a1.loglog(u, [r["mmse"] for r in rows], "-.", color=COL["line"], lw=1)
    a1.set_ylim(1e-9, 300)
    a1.set_xlabel(t["u"])
    u_axis(a1)
    a1.set_ylabel(t["mse"])
    a1.legend(fontsize=9, loc="upper left")
    for key, ls in (("equal", "-"), ("matched", "--"), ("geometric", "-")):
        c = d01["designs"][key]["cume_at_equal_jump"]
        x = np.array(c["x"])
        y = np.array(c["frac"])
        a2.semilogx(x[x > 0.02], y[x > 0.02], ls=ls, lw=2.2 if key != "matched" else 1.6,
                    color={"equal": COL["coh"], "matched": COL["inc"], "geometric": COL["pm"]}[key], label=t[key])
    a2.set_xlabel(t["delta"])
    a2.set_ylabel(t["frac"])
    a2.set_ylim(0, 1.05)
    a2.legend(loc="upper left")
    save(f, lang, "f10_geometric.png")


def f11_cluster(lang, d07):
    t = T[lang]
    c = d07["designs"]["cluster"]
    rows = sorted([r for r in d07["rows"] if r["design"] == "cluster"], key=lambda r: r["sigma"])
    sig = np.array([r["sigma"] for r in rows])
    x = snr_db_of(sig)
    f, ax = fig(8.4, 4.8)
    ax.semilogy(x, [r["mse_coh"] for r in rows], "o", color=COL["coh"], label=t["coh"])
    ax.semilogy(x, [r["mse_inc"] for r in rows], "s", color=COL["inc"], mfc="none", label=t["inc"])
    sd = np.geomspace(sig[0], sig[-1], 300)
    ax.semilogy(snr_db_of(sd), mixture(c["f_hz"], HALF_E07, sd), color=COL["hyb"], ls="--", lw=1.6, label=t["valley"])
    ax.semilogy(x, [r["lattice"] for r in rows], color=COL["theory"], lw=1.2, alpha=0.7, label=t["theory"])
    ax.semilogy(x, [r["hmse"] for r in rows], "--", color=COL["line"], lw=1.2, label=t["hmse"])
    ax.semilogy(x, [r["mmse"] for r in rows], "-.", color=COL["line"], lw=1.2, label=t["mmse"])
    ax.semilogy(x, c["mmse_within_factor"] * sig**2, ":", color=COL["line"], lw=1.6, label=t["within"])
    # second layer: carrier line → full-band line, bandwidth line → within-cluster line
    R2 = c["S"] / c["S_within"]
    l2 = c["mmse_within_factor"] * sd**2 * lattice.universal_curve(sd / c["sigma_jump2"], R2)
    ax.semilogy(snr_db_of(sd), l2, color=COL["pm"], lw=1.8, label=t["layer2"])
    for s in (c["sigma_jump"], c["sigma_jump2"]):
        ax.axvline(float(snr_db_of(s)), color=COL["inc"], ls="--", lw=0.9)
    ax.invert_xaxis()
    ax.set_ylim(1e-9, 300)
    ax.set_xlabel(t["snr"])
    ax.set_ylabel(t["mse"])
    ax.legend(loc="upper left", fontsize=9.5)
    save(f, lang, "f11_cluster.png")


def f12_design(lang, d08):
    t = T[lang]
    f, (a1, a2) = plt.subplots(1, 2, figsize=(11.5, 4.6), constrained_layout=True)
    rnd = d08["random"]
    a1.scatter([r["sigma_jump"] for r in rnd], [r["sigma_leave"] for r in rnd], s=7, color="#bbbbbb", label=t["random"])
    sj = np.linspace(0.025, 0.135, 50)
    a1.plot(sj, sj / 4.9083, color=COL["theory"], lw=1.2, ls="--", label=t["foot"])
    a1.set_xlim(0.025, 0.135)
    a1.set_ylim(0.004, 0.029)
    for key, mk, col in (("equal", "o", COL["coh"]), ("edges", "s", COL["inc"])):
        r = d08["refs"][key]
        a1.plot(r["sigma_jump"], r["sigma_leave"], mk, color=col, ms=9, label=t[key])
    b = d08["best"]
    a1.plot(b["sigma_jump"], b["sigma_leave"], "*", color=COL["hyb"], ms=15, label=t["best"])
    a1.set_xlabel(t["sj"])
    a1.set_ylabel(t["leave"])
    a1.legend(loc="upper left", fontsize=9.5)
    sj_eq = core.moments(core.equal_hz())["sigma_jump"]
    for lab, col, mk in (("equal", COL["coh"], "o"), ("best", COL["hyb"], "*")):
        rr = sorted([r for r in d08["mc"] if r["label"] == lab], key=lambda r: r["sigma"])
        x = snr_db_of([r["sigma"] for r in rr])
        a2.semilogy(x, [r["mse"] for r in rr], mk, color=col, ms=7 if mk == "*" else 5, label=t[lab] + " MC")
        mhz = core.EQUAL_MHZ if lab == "equal" else np.array(d08["best"]["mhz"], dtype=float)
        sd = np.geomspace(min(r["sigma"] for r in rr), max(r["sigma"] for r in rr), 200)
        a2.semilogy(snr_db_of(sd), mixture(mhz * 1e6, 15.0, sd), color=col, lw=1.2, ls="--")
        a2.semilogy(x, [r["hmse"] for r in rr], ":", color=col, lw=1)
    a2.invert_xaxis()
    a2.set_xlabel(t["snr"])
    a2.set_ylabel(t["mse"])
    a2.legend(loc="upper left", fontsize=9.5)
    save(f, lang, "f12_design.png")


def f13_pareto(lang, d09):
    t = T[lang]
    f, (a1, a2) = plt.subplots(1, 2, figsize=(11.5, 4.6), constrained_layout=True)
    cl = d09["cloud"]
    a1.scatter([p["sigma_leave"] for p in cl], [10 * np.log10(p["mse_ref"]) for p in cl], s=6, color="#bbbbbb",
               label=t["random"])
    sa = d09["anneal"]
    a1.scatter([p["sigma_leave"] for p in sa], [10 * np.log10(p["mse_ref"]) for p in sa], s=18, color=COL["hyb"],
               marker="x", label=t["anneal"])
    fr = sorted(d09["front"], key=lambda p: p["sigma_leave"])
    a1.plot([p["sigma_leave"] for p in fr], [10 * np.log10(p["mse_ref"]) for p in fr], "-", color=COL["theory"], lw=1.6,
            label=t["front"])
    styles = {"equal": ("o", COL["coh"]), "front_lo": ("s", "#1b9e77"), "front_mid": ("D", "#d95f02"),
              "front_hi": ("*", COL["hyb"])}
    for lab, mhz in d09["picks"].items():
        mk, col = styles[lab]
        p = d09["equal"] if lab == "equal" else next(q for q in fr if q["mhz"] == mhz)
        a1.plot(p["sigma_leave"], 10 * np.log10(p["mse_ref"]), mk, color=col, ms=10 if mk == "*" else 8, mec="k",
                label=t.get(lab, lab))
    a1.set_xlabel(t["leave"])
    a1.set_ylabel(t["mse_ref"].format(s=d09["sigma_ref"]))
    a1.legend(fontsize=8.5, loc="upper left")
    for lab, mhz in d09["picks"].items():
        mk, col = styles[lab]
        rr = sorted([r for r in d09["mc"] if r["label"] == lab], key=lambda r: r["sigma"])
        x = snr_db_of([r["sigma"] for r in rr])
        a2.semilogy(x, [r["mse"] for r in rr], mk, color=col, ms=8 if mk == "*" else 5, label=t.get(lab, lab))
        sd = np.geomspace(rr[0]["sigma"], rr[-1]["sigma"], 200)
        a2.semilogy(snr_db_of(sd), mixture(np.array(mhz, dtype=float) * 1e6, 15.0, sd), color=col, lw=1.0, ls="--")
    a2.invert_xaxis()
    a2.set_xlabel(t["snr"])
    a2.set_ylabel(t["mse"])
    a2.legend(fontsize=9, loc="upper left")
    save(f, lang, "f13_pareto.png")


def main(which=None):
    setup_fonts()
    data = {}
    for name in ("e01_bounds", "e02_lattice_validation", "e03_collapse", "e04_bayes", "e05_complex", "e06_bias",
                 "e07_layers", "e08_design", "e09_pareto"):
        try:
            data[name] = sim.load(name)
        except FileNotFoundError:
            print("missing", name)
    jobs = {
        "f01": (f01_hero, ("e02_lattice_validation", "e01_bounds")),
        "f02": (f02_two_rulers, ()),
        "f03": (f03_likelihoods, ()),
        "f03gif": (f03_gif, ()),
        "f04": (f04_validation, ("e02_lattice_validation",)),
        "f05": (f05_kpmf, ("e02_lattice_validation",)),
        "f06": (f06_collapse, ("e03_collapse",)),
        "f07": (f07_bayes, ("e04_bayes",)),
        "f08": (f08_complex, ("e05_complex",)),
        "f09": (f09_bias, ("e06_bias",)),
        "f10": (f10_geometric, ("e07_layers", "e01_bounds")),
        "f11": (f11_cluster, ("e07_layers",)),
        "f12": (f12_design, ("e08_design",)),
        "f13": (f13_pareto, ("e09_pareto",)),
    }
    for key, (fn, deps) in jobs.items():
        if which and key not in which:
            continue
        if any(d not in data for d in deps):
            print("skip", key)
            continue
        for lang in ("zh", "en"):
            fn(lang, *[data[d] for d in deps])
        print("done", key, flush=True)


if __name__ == "__main__":
    main(sys.argv[1:] or None)
