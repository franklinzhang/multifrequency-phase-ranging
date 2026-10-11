#!/usr/bin/env python3
"""Multi-frequency phase ranging: ambiguity, MSE, and a Ziv–Zakai integral.

Writes computation/data/results.json and computation/data/curves.npz.

Archived files under computation/data/ are authoritative for the notes.
Do not silently overwrite them; if you re-run, state any numerical differences.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
from scipy.special import erfc

C = 299_792_458.0
OUT = Path(__file__).resolve().parent / "data"
RNG_SEED = 7


def wrap(x: np.ndarray) -> np.ndarray:
    return (x + np.pi) % (2.0 * np.pi) - np.pi


def moments(f_hz: np.ndarray) -> dict:
    f = np.asarray(f_hz, dtype=float)
    n = f.size
    fbar = float(f.mean())
    S = float(np.sum((f - fbar) ** 2))
    F = float(np.sum(f**2))
    sig_jump = math.pi * math.sqrt(S) / fbar
    R = F / S
    return {
        "n": int(n),
        "f_hz": f.tolist(),
        "fbar_hz": fbar,
        "S": S,
        "F": F,
        "sigma_jump": sig_jump,
        "R": R,
        "R_db": 10.0 * math.log10(R),
        "lambda_m": C / fbar,
    }


def design_moment_matched(s_target: int = 4200) -> np.ndarray:
    """Integer-MHz offsets, sum 0, sum of squares = s_target, gcd of gaps = 1."""
    rng = np.random.default_rng(1)
    best = None
    best_score = 1e9
    equal = [-35, -25, -15, -5, 5, 15, 25, 35]
    for _ in range(30):
        raw = rng.normal(0.0, 16.0, size=(200_000, 8))
        raw -= raw.mean(axis=1, keepdims=True)
        ss = np.sum(raw**2, axis=1, keepdims=True)
        raw *= math.sqrt(s_target) / np.sqrt(ss)
        cand = np.rint(raw).astype(np.int64)
        ordered = np.sort(cand, axis=1)
        mask = (
            (cand.sum(axis=1) == 0)
            & (np.sum(cand**2, axis=1) == s_target)
            & np.all(np.diff(ordered, axis=1) > 0, axis=1)
        )
        span = ordered[:, -1] - ordered[:, 0]
        mask &= (span >= 60) & (span <= 110)
        for row in ordered[mask]:
            values = row.tolist()
            if values == equal:
                continue
            gaps = [int(gap) for gap in np.diff(row)]
            if min(gaps) < 5:
                continue
            g = 0
            for gap in gaps:
                g = math.gcd(g, gap)
            if g != 1:
                continue
            score = abs(int(span[0] if False else row[-1] - row[0]) - 84)
            if score < best_score:
                best_score = score
                best = row.copy()
            if best_score <= 6:
                return best
    if best is None:
        raise RuntimeError("moment-matched frequency set not found")
    return best


def geometric_hz() -> np.ndarray:
    """Eight tones in geometric progression, carrier near 2.435 GHz."""
    ratio = 1.02
    exponents = np.arange(8) - 3.5
    f_mid = 2435e6
    return f_mid * ratio**exponents


def j_of_delta(f: np.ndarray, delta: np.ndarray) -> np.ndarray:
    phase = (2.0 * np.pi / C) * delta[:, None] * f[None, :]
    return np.sum(wrap(phase) ** 2, axis=1)


def local_minima(delta: np.ndarray, cost: np.ndarray, f: np.ndarray, min_delta: float) -> list[dict]:
    found = []
    for i in range(1, cost.size - 1):
        if delta[i] < min_delta:
            continue
        if cost[i] <= cost[i - 1] and cost[i] < cost[i + 1]:
            # Parabolic refine.
            x1, x2, x3 = delta[i - 1], delta[i], delta[i + 1]
            y1, y2, y3 = cost[i - 1], cost[i], cost[i + 1]
            denom = (x1 - x2) * (x1 - x3) * (x2 - x3)
            if denom == 0:
                x = x2
                y = y2
            else:
                a = (x3 * (y2 - y1) + x2 * (y1 - y3) + x1 * (y3 - y2)) / denom
                b = (x3**2 * (y1 - y2) + x2**2 * (y3 - y1) + x1**2 * (y2 - y3)) / denom
                x = x2 if a == 0 else -b / (2 * a)
                y = float(np.interp(x, delta, cost))
            found.append({"delta_m": float(delta[i]), "J": float(cost[i]), "norm": float(math.sqrt(max(float(cost[i]), 0.0)))})
    # Deduplicate neighbours closer than 2 cm, keep the lower one.
    found.sort(key=lambda item: item["delta_m"])
    merged: list[dict] = []
    for item in found:
        if merged and abs(item["delta_m"] - merged[-1]["delta_m"]) < 0.02:
            if item["J"] < merged[-1]["J"]:
                merged[-1] = item
        else:
            merged.append(item)
    # Re-evaluate on a fine grid. The coarse sample can sit off the bottom of a narrow dip.
    refined = []
    for item in merged:
        xs = np.linspace(item["delta_m"] - 0.008, item["delta_m"] + 0.008, 321)
        js = j_of_delta(f, xs)
        k = int(np.argmin(js))
        y = float(js[k])
        x = float(xs[k])
        refined.append({"delta_m": x, "J": y, "norm": float(math.sqrt(max(y, 0.0)))})
    refined.sort(key=lambda item: item["J"])
    return refined


def q_func(x: np.ndarray | float) -> np.ndarray | float:
    arr = np.asarray(x, dtype=float)
    out = 0.5 * erfc(arr / math.sqrt(2.0))
    if np.isscalar(x):
        return float(out)
    return out


def predict_leave(j_side: float, delta_side: float, hmse_factor: float) -> dict:
    """Solve 2 Q(sqrt(J)/(2σ)) * δ^2 = HMSE(σ).

    HMSE = hmse_factor * σ^2, hmse_factor = c^2 / (4 π^2 F).
    """
    norm = math.sqrt(j_side)
    # Scan σ because the equation is transcendental.
    sigmas = np.geomspace(1e-3, 2.0, 4000)
    p = 2.0 * q_func(norm / (2.0 * sigmas))
    pred = p * delta_side**2
    hmse = hmse_factor * sigmas**2
    # First crossing where the slip term exceeds HMSE, from the quiet side.
    idx = np.where(pred >= hmse)[0]
    if idx.size == 0:
        return {"sigma_leave": None}
    i = int(idx[0])
    return {
        "sigma_leave": float(sigmas[i]),
        "x": float(norm / (2.0 * sigmas[i])),
        "J": j_side,
        "delta_m": delta_side,
        "norm": norm,
    }


def wrapped_phase_fisher(sigma: float, nterms: int = 50, ngrid: int = 8001) -> float:
    phi = np.linspace(-math.pi, math.pi, ngrid)
    dphi = phi[1] - phi[0]
    ks = np.arange(-nterms, nterms + 1)
    z = phi[None, :] + 2.0 * math.pi * ks[:, None]
    inv = 1.0 / sigma
    g = np.exp(-0.5 * (z * inv) ** 2) * inv / math.sqrt(2.0 * math.pi)
    dp = np.sum(g * (-z / sigma**2), axis=0)
    p = np.sum(g, axis=0)
    return float(np.sum((dp**2) / np.maximum(p, 1e-300)) * dphi)


def zzb_curve(f: np.ndarray, delta: np.ndarray, cost: np.ndarray, sigmas: np.ndarray, length: float) -> dict:
    """Valley-filled ZZB on a uniform cell of length `length`."""
    # Use δ in (0, L]. cost was computed on a symmetric grid; fold by interpolation.
    grid = np.linspace(0.0, length, 4000)
    j = np.interp(grid, delta, cost, left=0.0, right=cost[-1])
    j[0] = 0.0
    out = []
    for sigma in sigmas:
        norm = np.sqrt(np.maximum(j, 0.0))
        pe = q_func(norm / (2.0 * sigma))
        pe = np.asarray(pe, dtype=float)
        pe[0] = 0.5  # indistinguishable at δ=0; limit of Q(0)
        # At δ=0 the two hypotheses coincide, Pe=1/2. Q(0)=1/2. Good.
        h = (length - grid) * pe
        # Valley fill from the right.
        g = np.maximum.accumulate(h[::-1])[::-1]
        integrand = grid * g
        val = float(np.trapezoid(integrand, grid) / length)
        out.append(val)
    return {"L": length, "sigma": sigmas.tolist(), "zzb": out}


def simulate_set(
    name: str,
    f: np.ndarray,
    sigmas: np.ndarray,
    d0: float,
    half_window: float,
    n_each: np.ndarray,
    rng: np.random.Generator,
) -> dict:
    meta = moments(f)
    step = meta["lambda_m"] / 8.0
    dgrid = np.arange(d0 - half_window, d0 + half_window, step)
    # Drop a perfect alias if it lands on the grid.
    pred0 = (2.0 * np.pi / C) * (dgrid - d0)[:, None] * f[None, :]
    alias = np.sum(wrap(pred0) ** 2, axis=1) < 1e-8
    alias[np.argmin(np.abs(dgrid - d0))] = False
    dgrid = dgrid[~alias]

    f64 = f.astype(float)
    F = float(np.sum(f64**2))
    f_center = f64 - f64.mean()
    S = float(np.sum(f_center**2))
    chunk = 1500
    mse_coh = []
    mse_inc = []
    slip = []
    for k, sigma in enumerate(sigmas):
        n_use = int(n_each[k])
        phi = wrap((2.0 * np.pi / C) * f64 * d0 + rng.normal(0.0, sigma, size=(n_use, f64.size)))
        best_c = np.full(n_use, np.inf)
        arg_c = np.zeros(n_use, dtype=int)
        best_i = np.full(n_use, np.inf)
        arg_i = np.zeros(n_use, dtype=int)
        for start in range(0, dgrid.size, chunk):
            grid = dgrid[start : start + chunk]
            pred = (2.0 * np.pi / C) * grid[:, None] * f64[None, :]
            resid = wrap(phi[:, None, :] - pred[None, :, :])
            cost_c = np.sum(resid**2, axis=-1)
            z = np.exp(1j * (phi[:, None, :] - pred[None, :, :]))
            cost_i = -np.abs(np.sum(z, axis=-1)) ** 2
            i_c = np.argmin(cost_c, axis=1)
            v_c = cost_c[np.arange(n_use), i_c]
            take_c = v_c < best_c
            best_c[take_c] = v_c[take_c]
            arg_c[take_c] = i_c[take_c] + start
            i_i = np.argmin(cost_i, axis=1)
            v_i = cost_i[np.arange(n_use), i_i]
            take_i = v_i < best_i
            best_i[take_i] = v_i[take_i]
            arg_i[take_i] = i_i[take_i] + start
            offset = start
        d_c = dgrid[arg_c].copy()
        d_i = dgrid[arg_i].copy()
        for _ in range(16):
            resid = wrap(phi - (2.0 * np.pi / C) * d_c[:, None] * f64[None, :])
            d_c = d_c + (C / (2.0 * np.pi)) * (resid @ f64) / F
            pred = (2.0 * np.pi / C) * d_i[:, None] * f64[None, :]
            theta = np.angle(np.sum(np.exp(1j * (phi - pred)), axis=1))
            resid_i = wrap(phi - pred - theta[:, None])
            d_i = d_i + (C / (2.0 * np.pi)) * (resid_i @ f_center) / S
        d_c = np.clip(d_c, dgrid[0], dgrid[-1])
        d_i = np.clip(d_i, dgrid[0], dgrid[-1])
        err_c = d_c - d0
        err_i = d_i - d0
        mse_coh.append(float(np.mean(err_c**2)))
        mse_inc.append(float(np.mean(err_i**2)))
        slip.append(float(np.mean(np.abs(err_c) > 0.5 * meta["lambda_m"])))
        print(
            f"{name} σ={sigma:.5f} n={n_use}  coh {10*math.log10(mse_coh[-1]):7.2f} dBm²  "
            f"inc {10*math.log10(mse_inc[-1]):7.2f} dBm²  slip {slip[-1]:.4f}",
            flush=True,
        )
    meta["mse_coh"] = mse_coh
    meta["mse_inc"] = mse_inc
    meta["slip_rate"] = slip
    meta["sigma"] = sigmas.tolist()
    meta["half_window_m"] = half_window
    meta["grid_step_m"] = step
    meta["ntrial"] = [int(n) for n in n_each]
    meta["d0_m"] = d0
    return meta


def ambiguity_bundle(f: np.ndarray, d_max: float, min_delta: float) -> dict:
    delta = np.linspace(0.0, d_max, int(d_max / 0.002) + 1)
    cost = j_of_delta(f, delta)
    minima = local_minima(delta, cost, f, min_delta=min_delta)
    # Carrier-wavelength candidate.
    fbar = float(f.mean())
    lam = C / fbar
    j_lam = float(j_of_delta(f, np.array([lam]))[0])
    return {
        "delta_m": delta,
        "J": cost,
        "minima": minima[:80],
        "J_at_lambda": j_lam,
        "lambda_m": lam,
    }


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    equal_mhz = 2400.0 + 10.0 * np.arange(8)
    matched_mhz = design_moment_matched(4200)
    print("matched offsets MHz", matched_mhz.tolist(), flush=True)
    sets = {
        "equal": equal_mhz * 1e6,
        "matched": (2435.0 + matched_mhz) * 1e6,
        "geometric": geometric_hz(),
    }
    # Ambiguity out to one equal-spacing cell, plus a longer look for geometry.
    d_cell = C / 10e6
    bundles = {}
    for name, f in sets.items():
        d_max = d_cell if name != "geometric" else max(d_cell, 80.0)
        meta_m = moments(f)
        min_delta = 0.25 * meta_m["lambda_m"]
        bundles[name] = ambiguity_bundle(f, d_max=d_max, min_delta=min_delta)
        top = bundles[name]["minima"][:5]
        print(name, "R_dB", round(meta_m["R_db"], 3), "sj", round(meta_m["sigma_jump"], 5))
        print("  J(λ)", bundles[name]["J_at_lambda"], "top minima", top)

    sig_ref = moments(sets["equal"])["sigma_jump"]
    ratios = np.array(
        [0.10, 0.16, 0.22, 0.25, 0.28, 0.32, 0.40, 0.50, 0.65, 0.85, 1.0, 1.5, 2.2, 3.5, 6.0, 10.0, 16.0, 28.0]
    )
    sigmas = sig_ref * ratios
    n_each = np.array([8000 if 0.22 <= r <= 0.50 else 1600 for r in ratios])
    rng = np.random.default_rng(RNG_SEED)
    d0 = 12.345
    half = 0.5 * d_cell
    results = {}
    curves = {}
    for name, f in sets.items():
        results[name] = simulate_set(name, f, sigmas, d0, half, n_each, rng)
        amb = bundles[name]
        # Store a decimated ambiguity for the figure (every ~1 cm already 2 mm; decimate to 4 mm).
        step = max(1, int(round(0.004 / (amb["delta_m"][1] - amb["delta_m"][0]))))
        curves[f"{name}_delta"] = amb["delta_m"][::step]
        curves[f"{name}_J"] = amb["J"][::step]
        results[name]["J_at_lambda"] = amb["J_at_lambda"]
        results[name]["minima"] = amb["minima"]
        hmse_factor = C**2 / (4.0 * math.pi**2 * results[name]["F"])
        mmse_factor = C**2 / (4.0 * math.pi**2 * results[name]["S"])
        results[name]["hmse_factor"] = hmse_factor
        results[name]["mmse_factor"] = mmse_factor
        inside = [item for item in amb["minima"] if item["delta_m"] < half]
        if inside:
            near = min(inside, key=lambda item: item["delta_m"])
            low = min(inside, key=lambda item: item["J"])
            results[name]["leave_lowest"] = predict_leave(low["J"], low["delta_m"], hmse_factor)
            results[name]["leave_nearest"] = predict_leave(near["J"], near["delta_m"], hmse_factor)
            results[name]["nearest_minimum"] = near
            results[name]["lowest_minimum"] = low
            results[name]["minima_in_window"] = inside
        # λ-shift prediction even if it is not a minimum.
        results[name]["leave_lambda"] = predict_leave(
            amb["J_at_lambda"], results[name]["lambda_m"], hmse_factor
        )

    # Local wrapped-normal CRB for the equal set, carrier information only.
    f_eq = sets["equal"]
    crb = []
    for sigma in sigmas:
        i_phi = wrapped_phase_fisher(float(sigma))
        i_d = float(np.sum((2.0 * math.pi * f_eq / C) ** 2) * i_phi)
        crb.append(1.0 / i_d)
    results["equal"]["crb"] = crb

    zzb_sigmas = sigmas
    for name in ("equal", "matched"):
        delta = bundles[name]["delta_m"]
        cost = bundles[name]["J"]
        # Ambiguity array starts at 0 and goes to d_cell. Mirror is even, J(δ)=J(-δ).
        results[name]["zzb"] = zzb_curve(sets[name], delta, cost, zzb_sigmas, d_cell)

    # Analytic reference lines at the simulated sigmas.
    for name in results:
        if name not in sets:
            continue
        sig = np.array(results[name]["sigma"])
        results[name]["hmse"] = (results[name]["hmse_factor"] * sig**2).tolist()
        results[name]["mmse"] = (results[name]["mmse_factor"] * sig**2).tolist()
        sj = results[name]["sigma_jump"]
        # One-cycle mixture against λ.
        norm = math.sqrt(results[name]["J_at_lambda"])
        p = 2.0 * np.asarray(q_func(norm / (2.0 * sig)))
        lam = results[name]["lambda_m"]
        mix = (1.0 - p) * results[name]["hmse_factor"] * sig**2 + p * lam**2
        results[name]["mixture_lambda"] = mix.tolist()
        floor = (d_cell**2) / 12.0
        results[name]["floor_L2_12"] = floor

    # JSON cannot hold numpy. results are plain.
    payload = {
        "c": C,
        "d_cell_m": d_cell,
        "seed": RNG_SEED,
        "sigma_ref_equal": sig_ref,
        "sets": results,
    }
    (OUT / "results.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    np.savez_compressed(OUT / "curves.npz", **{k: np.asarray(v) for k, v in curves.items()})
    print("wrote", OUT / "results.json", flush=True)


if __name__ == "__main__":
    main()
