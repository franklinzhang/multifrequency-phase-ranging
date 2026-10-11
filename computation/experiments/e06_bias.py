#!/usr/bin/env python3
"""Deterministic phase errors: which ones move the jump?

A per-tone phase error m splits into three parts: along g (a consistent range shift),
along the slip direction u (an envelope–carrier disagreement of β wavelengths, which
moves the jump), and the rest (second order). Three cases:

  common   m_i = ε, an uncalibrated common phase;
  delay    m_i = 2π f_i τ, an uncalibrated group delay;
  twopath  h_i = 1 + a·exp(j(2π f_i Δ/c + ψ)), a single reflection with excess path Δ.

The receiver is the equal-weight coherent wrapped least squares on arg y_i.
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
F = core.equal_hz()
M = core.moments(F)
LAM = M["lambda_m"]


def mc_phase(rng, sigma_vec, offset, n=8000):
    d0 = CENTER + (rng.random(n) - 0.5) * LAM
    grid = estimators.make_grid(CENTER - HALF, CENTER + HALF, LAM)
    phi = sim.phases(rng, F, d0, sigma_vec, n, offset=offset)
    err = estimators.coherent_ml(phi, F, grid) - d0
    return err


def summarize(err, lp):
    bias = lp["fine_bias_m"]
    ml = lattice.ml_mse(lp)
    return {
        "mse": float(np.mean(err**2)),
        "slip": float(np.mean(np.abs(err - bias) > 0.5 * LAM)),
        "lattice_mse": ml["mse"],
        "lattice_slip": ml["slip"],
        "beta": lp["beta"],
        "fine_bias_m": bias,
    }


def job_common(args):
    eps_deg, u, seed = args
    rng = np.random.default_rng(seed)
    sigma = u * M["sigma_jump"]
    off = np.full(F.size, math.radians(eps_deg))
    lp = lattice.params(F, sigma, off, w_rx=np.ones(F.size))
    row = {"case": "common", "eps_deg": eps_deg, "u": u, **summarize(mc_phase(rng, sigma, off), lp)}
    print(f"common ε={eps_deg:4.0f}° u={u:4.2f} slip {row['slip']:.4f}/{row['lattice_slip']:.4f}", flush=True)
    return row


def job_delay(args):
    tau_ns, u, seed = args
    rng = np.random.default_rng(seed)
    sigma = u * M["sigma_jump"]
    off = 2 * math.pi * F * tau_ns * 1e-9
    lp = lattice.params(F, sigma, off, w_rx=np.ones(F.size))
    row = {"case": "delay", "tau_ns": tau_ns, "u": u, **summarize(mc_phase(rng, sigma, off), lp)}
    print(f"delay τ={tau_ns} ns u={u:4.2f} slip {row['slip']:.4f}/{row['lattice_slip']:.4f}", flush=True)
    return row


def twopath(a, delta_m, psi=math.pi):
    return 1.0 + a * np.exp(1j * (2 * math.pi * F * delta_m / core.C + psi))


def job_twopath(args):
    a, delta_m, snr_db, seed = args
    rng = np.random.default_rng(seed)
    h = twopath(a, delta_m)
    snr = 10 ** (snr_db / 10)
    n = 8000
    d0 = CENTER + (rng.random(n) - 0.5) * LAM
    grid = estimators.make_grid(CENTER - HALF, CENTER + HALF, LAM)
    y = sim.tones(rng, F, d0, snr, n, channel=h)
    err = estimators.coherent_ml(np.angle(y), F, grid) - d0
    sig_vec = core.sigma_from_snr(snr * np.abs(h) ** 2)
    lp = lattice.params(F, sig_vec, np.angle(h), w_rx=np.ones(F.size))
    row = {"case": "twopath", "a": a, "delta_m": delta_m, "snr_db": snr_db,
           "u_nominal": float(core.sigma_from_snr(snr) / M["sigma_jump"]), **summarize(err, lp)}
    print(f"twopath a={a} Δ={delta_m:5.2f} SNR={snr_db:5.1f}  slip {row['slip']:.4f}/{row['lattice_slip']:.4f} "
          f"β={lp['beta']:+.3f}", flush=True)
    return row


def main():
    seed = 70000
    jobs = []
    u_list = [0.3, 0.4, 0.5, 0.7, 1.0, 1.5]
    for eps in (0, 10, 20, 45, 90):
        for u in u_list:
            seed += 1
            jobs.append(("common", (float(eps), u, seed)))
    for tau in (0.0, 1.0, 5.0):
        for u in (0.4, 0.7, 1.0):
            seed += 1
            jobs.append(("delay", (tau, u, seed)))
    snr_u05 = float(10 * math.log10(core.snr_from_sigma(0.5 * M["sigma_jump"])))
    snr_u035 = float(10 * math.log10(core.snr_from_sigma(0.35 * M["sigma_jump"])))
    for snr_db in (snr_u05, snr_u035):
        for d in np.round(np.linspace(0.0, 3.0, 31), 3):
            seed += 1
            jobs.append(("twopath", (0.3, float(d), snr_db, seed)))
    for d in (0.25, 1.0):
        for u in (0.25, 0.3, 0.35, 0.4, 0.5, 0.7, 1.0, 1.5, 2.2):
            seed += 1
            snr_db = float(10 * math.log10(core.snr_from_sigma(u * M["sigma_jump"])))
            jobs.append(("twopath", (0.3, d, snr_db, seed)))
    fn = {"common": job_common, "delay": job_delay, "twopath": job_twopath}
    with ProcessPoolExecutor(max_workers=4) as ex:
        futs = [ex.submit(fn[k], a) for k, a in jobs]
        rows = [fu.result() for fu in futs]
    sim.save("e06_bias", {"rows": rows, "snr_u05_db": snr_u05, "snr_u035_db": snr_u035})


if __name__ == "__main__":
    main()
