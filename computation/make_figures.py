#!/usr/bin/env python3
"""Figures for the multi-frequency phase-ranging note.

All axis labels and legends are English. Canvas width 2000 px.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import font_manager

ROOT = Path(__file__).resolve().parent
DATA = json.loads((ROOT / "data" / "results.json").read_text())
FIG = ROOT.parent / "notes" / "figures"
FIG.mkdir(parents=True, exist_ok=True)
C = 299_792_458.0

# Portable default face (English labels only).
SANS = font_manager.FontProperties(size=46)

DPI = 100
MAIN = 46
TICK = 42
LEGEND = 40

X_SIGMA = r"$\sigma\,/\,\sigma_{\mathrm{jump}}$"


def style_ax(ax, xlabel: str, ylabel: str, math_x: bool = False) -> None:
    if math_x:
        ax.set_xlabel(xlabel, fontsize=MAIN, labelpad=16)
    else:
        ax.set_xlabel(xlabel, fontproperties=SANS, fontsize=MAIN, labelpad=16)
    ax.set_ylabel(ylabel, fontproperties=SANS, fontsize=MAIN, labelpad=16)
    ax.tick_params(axis="both", which="major", labelsize=TICK, width=2.2, length=10, pad=8)
    ax.tick_params(axis="both", which="minor", width=1.4, length=6)
    for lab in ax.get_xticklabels() + ax.get_yticklabels():
        lab.set_fontproperties(SANS)
        lab.set_fontsize(TICK)
    ax.grid(True, which="major", color="#d0d0d0", linewidth=1.2)
    ax.grid(True, which="minor", color="#e6e6e6", linewidth=0.8)
    for spine in ax.spines.values():
        spine.set_linewidth(2.0)
    ax.set_axisbelow(True)


def legend(ax, loc="best", ncol: int = 1, above: bool = False) -> None:
    kwargs = dict(prop=SANS, fontsize=LEGEND, frameon=True, fancybox=False, edgecolor="#222222", ncol=ncol)
    if above:
        leg = ax.legend(loc="lower left", bbox_to_anchor=(0.0, 1.02), borderaxespad=0.0, **kwargs)
    else:
        leg = ax.legend(loc=loc, **kwargs)
    leg.get_frame().set_linewidth(1.6)


def new_fig(height: float = 11.6, top: float = 0.955):
    fig, ax = plt.subplots(figsize=(20, height))
    fig.subplots_adjust(left=0.155, right=0.965, bottom=0.175, top=top)
    return fig, ax


def save(fig, name: str) -> None:
    path = FIG / name
    fig.savefig(path, dpi=DPI)
    plt.close(fig)
    print(path, "px", fig.get_size_inches() * DPI)


def plot_cost(ax, delta, cost, y_cap: float, **kwargs):
    shown = np.where(cost <= y_cap, cost, np.nan)
    ax.plot(delta, shown, **kwargs)
    ax.set_ylim(0, y_cap * 1.16)


def wrap(x):
    return (x + np.pi) % (2 * np.pi) - np.pi


def cost_curve(f_hz: np.ndarray, d_max: float, n: int = 2500):
    delta = np.linspace(0.0, d_max, n)
    phase = (2 * np.pi / C) * delta[:, None] * f_hz[None, :]
    return delta, np.sum(wrap(phase) ** 2, axis=1)


def mse_ax(ax, ratios, mse, label, color, marker="o"):
    ax.loglog(ratios, mse, marker=marker, color=color, ms=14, lw=3.2, mew=1.5, mfc=color, label=label)


def main() -> None:
    eq = DATA["sets"]["equal"]
    md = DATA["sets"]["matched"]
    geo = DATA["sets"]["geometric"]
    f_eq = np.array(eq["f_hz"])
    f_md = np.array(md["f_hz"])
    f_geo = np.array(geo["f_hz"])
    ratio = np.array(eq["sigma"]) / eq["sigma_jump"]
    ratio_g = np.array(geo["sigma"]) / geo["sigma_jump"]
    floor = eq["floor_L2_12"]

    d, j = cost_curve(f_eq, 1.05, 1800)
    fig, ax = new_fig(11.2, top=0.86)
    plot_cost(ax, d, j, 2.0, color="#1f4e79", lw=3.6)
    ax.axvline(eq["lambda_m"], color="#c45911", lw=2.6, ls="--", label=r"$\lambda$")
    style_ax(ax, "Range offset / m", "Phase cost")
    legend(ax, above=True, ncol=1)
    ax.set_xlim(0, 1.05)
    save(fig, "fig01_equal_cost.png")

    fig, ax = new_fig()
    mse_ax(ax, ratio, eq["mse_coh"], "Coherent", "#1f4e79")
    ax.loglog(ratio, eq["hmse"], color="#222222", lw=3.0, ls="--", label="HMSE")
    ax.loglog(ratio, eq["mmse"], color="#555555", lw=3.0, ls="-.", label="MMSE")
    ax.loglog(ratio, eq["zzb_half"]["zzb"], color="#2e7d32", lw=3.2, label="Half-window bound")
    ax.axhline(floor, color="#8a8a8a", lw=2.6, ls=":", label="Search-window floor")
    ax.axvline(1.0, color="#c45911", lw=2.4, ls="--", label="Threshold")
    style_ax(ax, X_SIGMA, r"MSE / m$^2$", math_x=True)
    legend(ax, loc="upper left")
    ax.set_xlim(0.08, 40)
    ax.set_ylim(2e-9, 200)
    save(fig, "fig02_equal_mse.png")

    d, je = cost_curve(f_eq, 2.0, 2400)
    _, jm = cost_curve(f_md, 2.0, 2400)
    fig, ax = new_fig(11.2, top=0.86)
    plot_cost(ax, d, je, 6.5, color="#1f4e79", lw=3.8, label="Equal spacing")
    plot_cost(ax, d, jm, 6.5, color="#c45911", lw=3.0, ls="--", label="Moment-matched")
    style_ax(ax, "Range offset / m", "Phase cost")
    legend(ax, above=True, ncol=2)
    ax.set_xlim(0, 2.0)
    save(fig, "fig03_matched_cost.png")

    fig, ax = new_fig()
    mse_ax(ax, ratio, eq["mse_coh"], "Equal spacing", "#1f4e79")
    mse_ax(ax, ratio, md["mse_coh"], "Moment-matched", "#c45911", marker="s")
    ax.loglog(ratio, eq["hmse"], color="#222222", lw=3.0, ls="--", label="HMSE")
    ax.loglog(ratio, eq["mmse"], color="#555555", lw=3.0, ls="-.", label="MMSE")
    ax.axhline(floor, color="#8a8a8a", lw=2.6, ls=":", label="Search-window floor")
    ax.axvline(1.0, color="#c45911", lw=2.4, ls="--", label="Threshold")
    style_ax(ax, X_SIGMA, r"MSE / m$^2$", math_x=True)
    legend(ax, loc="upper left")
    ax.set_xlim(0.08, 40)
    ax.set_ylim(2e-9, 200)
    save(fig, "fig04_matched_mse.png")

    d, jg = cost_curve(f_geo, 14.0, 4200)
    fig, ax = new_fig(11.2)
    plot_cost(ax, d, jg, 10.0, color="#1f4e79", lw=3.2)
    ax.axvline(geo["lambda_m"], color="#c45911", lw=2.6, ls="--", label=r"$\lambda$")
    ax.axvline(6.211, color="#2e7d32", lw=2.6, ls="--", label="6.21 m")
    style_ax(ax, "Range offset / m", "Phase cost")
    leg = ax.legend(
        loc="center",
        bbox_to_anchor=(0.20, 0.62),
        prop=SANS,
        fontsize=LEGEND,
        frameon=True,
        fancybox=False,
        edgecolor="#222222",
    )
    leg.get_frame().set_linewidth(1.6)
    ax.set_xlim(0, 14)
    save(fig, "fig05_geometric_cost.png")

    fig, ax = new_fig()
    mse_ax(ax, ratio_g, geo["mse_coh"], "Coherent", "#1f4e79")
    mse_ax(ax, ratio_g, geo["mse_inc"], "Incoherent", "#c45911", marker="s")
    ax.loglog(ratio_g, geo["hmse"], color="#222222", lw=3.0, ls="--", label="HMSE")
    ax.loglog(ratio_g, geo["mmse"], color="#555555", lw=3.0, ls="-.", label="MMSE")
    ax.loglog(ratio_g, geo["zzb_half"]["zzb"], color="#2e7d32", lw=3.2, label="Half-window bound")
    ax.axvline(1.0, color="#c45911", lw=2.4, ls="--", label="Threshold")
    style_ax(ax, X_SIGMA, r"MSE / m$^2$", math_x=True)
    legend(ax, loc="upper left")
    ax.set_xlim(0.015, 8)
    ax.set_ylim(1e-9, 200)
    save(fig, "fig06_geometric_mse.png")

    fig, ax = new_fig()
    mse_ax(ax, ratio, eq["mse_coh"], "Coherent", "#1f4e79")
    mse_ax(ax, ratio, eq["mse_inc"], "Incoherent", "#c45911", marker="s")
    ax.loglog(ratio, eq["hmse"], color="#222222", lw=3.0, ls="--", label="HMSE")
    ax.loglog(ratio, eq["mmse"], color="#555555", lw=3.0, ls="-.", label="MMSE")
    ax.axhline(floor, color="#8a8a8a", lw=2.6, ls=":", label="Search-window floor")
    ax.axvline(1.0, color="#c45911", lw=2.4, ls="--", label="Threshold")
    style_ax(ax, X_SIGMA, r"MSE / m$^2$", math_x=True)
    legend(ax, loc="upper left")
    ax.set_xlim(0.08, 40)
    ax.set_ylim(2e-9, 200)
    save(fig, "fig07_coherent_incoherent.png")

    bound = json.loads((ROOT / "data" / "bound_valleys.json").read_text())
    fig, ax = new_fig(11.6, top=0.86)
    styles = (
        ("equal", "Equal spacing", "#1f4e79", "-", 4.0),
        ("matched", "Moment-matched", "#c45911", "--", 3.2),
        ("geometric", "Geometric spacing", "#2e7d32", "-", 3.4),
    )
    for key, label, color, ls, lw in styles:
        block = bound["sets"][key]
        x = np.asarray(block["snap_delta_m"], dtype=float)
        y = np.asarray(block["snap_cume_frac"], dtype=float)
        keep = x >= 0.03
        ax.semilogx(x[keep], y[keep], color=color, ls=ls, lw=lw, label=label)
    style_ax(ax, "Range offset / m", "Cumulative fraction")
    legend(ax, ncol=3, above=True)
    ax.set_xlim(0.03, 15)
    ax.set_ylim(0.0, 1.18)
    save(fig, "fig08_bound_cumulative.png")


if __name__ == "__main__":
    main()
