#!/usr/bin/env python3
"""Attribute the half-window Ziv–Zakai integral to valleys. No Monte Carlo."""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
from scipy.special import erfc

from sim_threshold import C, design_moment_matched, geometric_hz, j_of_delta, moments

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "data" / "bound_valleys.json"
L = 0.5 * C / 10e6


def q_func(x: np.ndarray) -> np.ndarray:
    return 0.5 * erfc(np.asarray(x, dtype=float) / math.sqrt(2.0))


def valley_list(f: np.ndarray, delta: np.ndarray, cost: np.ndarray) -> list[dict]:
    """Minima of J on (0.05 m, L), refined already by the fine grid."""
    found = []
    for i in range(1, cost.size - 1):
        if delta[i] < 0.05 or delta[i] > L - 1e-6:
            continue
        if cost[i] <= cost[i - 1] and cost[i] < cost[i + 1] and cost[i] < cost[i - 1] - 1e-6:
            found.append((float(delta[i]), float(cost[i])))
    merged: list[dict] = []
    for x, y in found:
        if merged and x - merged[-1]["delta_m"] < 0.02:
            if y < merged[-1]["J"]:
                merged[-1] = {"delta_m": x, "J": y}
        else:
            merged.append({"delta_m": x, "J": y})
    return merged


def zzb_parts(delta: np.ndarray, cost: np.ndarray, sigma: float, beta: float) -> dict:
    j = cost.copy()
    j[0] = 0.0
    norm = np.sqrt(np.maximum(j, 0.0))
    pe = q_func(norm / (2.0 * sigma))
    pe[0] = 0.5
    w = (L - delta) * pe
    # G(δ) = max_{u ≥ δ} w(u)
    g = np.maximum.accumulate(w[::-1])[::-1]
    integrand = delta * g
    total = float(np.trapezoid(integrand, delta) / L)
    # Right-to-left records: a sample is a record if it sets the accumulate.
    rev_max = np.maximum.accumulate(w[::-1])
    is_record = np.zeros(delta.size, dtype=bool)
    is_record[-1] = True
    # A point is the location of a new record when the reversed max increases.
    inc = np.diff(rev_max, prepend=rev_max[0])
    record_from_right = inc > 0
    is_record = record_from_right[::-1]
    record_idx = np.flatnonzero(is_record)
    # Keep only records that are local peaks of w, plus the first interior sample.
    valleys = valley_list(np.array([]), delta, cost)
    # Map each record to the nearest valley within 1.5 cm; otherwise main lobe if < 2 cm.
    pieces = []
    order = record_idx  # increasing δ
    # The fill interval of record k (from the left) is (prev_record, this_record].
    edges = [0.0]
    for idx in order:
        edges.append(float(delta[idx]))
    # Drop the trivial δ=0 record.
    useful = [(idx, edges[k]) for k, idx in enumerate(order) if delta[idx] > 1e-4]
    for idx, left in useful:
        right = float(delta[idx])
        mask = (delta > left) & (delta <= right + 1e-15)
        if not np.any(mask):
            continue
        piece = float(np.trapezoid(integrand[mask], delta[mask]) / L)
        pieces.append(
            {
                "delta_m": right,
                "J": float(j[idx]),
                "pe": float(pe[idx]),
                "w": float(w[idx]),
                "left_m": left,
                "share": piece,
            }
        )
    # Cumulative fraction.
    cume = np.cumsum((integrand[:-1] + integrand[1:]) * 0.5 * np.diff(delta)) / L
    cume = np.concatenate([[0.0], cume])
    def at_frac(frac: float) -> float:
        k = int(np.searchsorted(cume, frac * total))
        k = min(max(k, 0), delta.size - 1)
        return float(delta[k])

    hmse = (sigma / beta) ** 2
    return {
        "sigma": sigma,
        "zzb": total,
        "zzb_db": 10.0 * math.log10(total),
        "hmse_db": 10.0 * math.log10(hmse),
        "d50_m": at_frac(0.5),
        "d90_m": at_frac(0.9),
        "pieces": pieces,
        "delta": delta,
        "cume_frac": cume / total,
        "integrand": integrand,
    }


def single_valley(delta_a: float, j_a: float, sigma: float, beta: float) -> float:
    """Fill from δ_× = √J / β to δ_a at height (L-δ_a) P_e. Returns the piece in m²."""
    pe = float(q_func(math.sqrt(j_a) / (2.0 * sigma)))
    delta_cross = math.sqrt(j_a) / beta
    if delta_cross >= delta_a:
        return 0.0
    return (L - delta_a) / L * pe * (delta_a**2 - delta_cross**2) / 2.0


def main() -> None:
    equal = (2400.0 + 10.0 * np.arange(8)) * 1e6
    matched = (2435.0 + design_moment_matched(4200)) * 1e6
    geometric = geometric_hz()
    sets = {"equal": equal, "matched": matched, "geometric": geometric}
    delta = np.linspace(0.0, L, 30001)
    sigmas = [
        0.008361,
        0.018395,
        0.083613,
        0.125420,
        0.403086,
        0.836133,
        2.341174,
    ]
    payload: dict = {"L_m": L, "sigmas": sigmas, "sets": {}}
    # Store a decimated cumulative curve at the equal-spacing threshold.
    snap_sigma = 0.083613
    for name, f in sets.items():
        meta = moments(f)
        beta = 2.0 * math.pi * math.sqrt(meta["F"]) / C
        cost = j_of_delta(f, delta)
        rows = []
        snap = None
        for sigma in sigmas:
            part = zzb_parts(delta, cost, sigma, beta)
            top = sorted(part["pieces"], key=lambda item: -item["share"])[:6]
            row = {
                "sigma": sigma,
                "zzb_db": part["zzb_db"],
                "hmse_db": part["hmse_db"],
                "d50_m": part["d50_m"],
                "d90_m": part["d90_m"],
                "top": [
                    {
                        "delta_m": round(item["delta_m"], 4),
                        "J": item["J"],
                        "pe": item["pe"],
                        "share_db": 10.0 * math.log10(item["share"]) if item["share"] > 0 else None,
                        "share_frac": item["share"] / (10 ** (part["zzb_db"] / 10.0)),
                    }
                    for item in top
                ],
            }
            rows.append(row)
            print(
                f"{name:10} σ={sigma:.4f}  ZZB {part['zzb_db']:7.2f}  HMSE {part['hmse_db']:7.2f}  "
                f"d50 {part['d50_m']:.3f}  d90 {part['d90_m']:.3f}  "
                f"top {row['top'][0]['delta_m'] if row['top'] else None}"
            )
            if abs(sigma - snap_sigma) < 1e-9:
                snap = part
        # Analytic pieces at the equal-spacing threshold.
        lam = meta["lambda_m"]
        j_lam = float(j_of_delta(f, np.array([lam]))[0])
        j_far = float(j_of_delta(f, np.array([6.21105]))[0])
        analytic = {
            "lambda_m": lam,
            "J_lambda": j_lam,
            "beta": beta,
            "piece_lambda_m2": single_valley(lam, j_lam, snap_sigma, beta),
            "piece_6p211_m2": single_valley(6.21105, j_far, snap_sigma, beta),
        }
        step = 30  # 1.5 cm
        payload["sets"][name] = {
            "fbar_hz": meta["fbar_hz"],
            "sigma_jump": meta["sigma_jump"],
            "R_db": meta["R_db"],
            "rows": rows,
            "analytic_at_equal_jump": analytic,
            "snap_sigma": snap_sigma,
            "snap_delta_m": snap["delta"][::step].tolist(),
            "snap_cume_frac": snap["cume_frac"][::step].tolist(),
        }
        a = analytic
        print(
            f"  analytic λ piece {10*math.log10(a['piece_lambda_m2']):.2f} dB, "
            f"6.211 piece {10*math.log10(max(a['piece_6p211_m2'], 1e-30)):.2f} dB"
        )
    OUT.write_text(json.dumps(payload), encoding="utf-8")
    print("wrote", OUT)


if __name__ == "__main__":
    main()
