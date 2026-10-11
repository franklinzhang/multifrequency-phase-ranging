#!/usr/bin/env python3
"""Complex tones instead of phases: does the jump stay where σ² = 1/(2 SNR) puts it?

y_i = exp(j2πf_i d/c) + n_i with E|n_i|² = 1/SNR. Estimators: coherent ML on the complex
tones (known amplitude and phase reference), coherent wrapped least squares on arg y_i,
and the incoherent estimator on arg y_i.
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
SNR_DB = np.arange(-6.0, 46.0, 2.0)


def phase_std_exact(snr: float, n: int = 400_000, seed: int = 3) -> float:
    rng = np.random.default_rng(seed)
    z = 1.0 + (rng.normal(size=n) + 1j * rng.normal(size=n)) * math.sqrt(0.5 / snr)
    return float(np.sqrt(np.mean(np.angle(z) ** 2)))


def run_one(args):
    snr_db, seed = args
    f = core.equal_hz()
    m = core.moments(f)
    snr = 10 ** (snr_db / 10)
    rng = np.random.default_rng(seed)
    n = 8000
    d0 = CENTER + (rng.random(n) - 0.5) * m["lambda_m"]
    grid = estimators.make_grid(CENTER - HALF, CENTER + HALF, m["lambda_m"])
    y = sim.tones(rng, f, d0, snr, n)
    ph = np.angle(y)
    lam = m["lambda_m"]
    res = {}
    for key, est in (
        ("complex_ml", estimators.complex_coherent_ml(y, f, grid)),
        ("phase_ml", estimators.coherent_ml(ph, f, grid)),
        ("phase_inc", estimators.incoherent_ml(ph, f, grid)),
    ):
        e = est - d0
        res[key] = {"mse": float(np.mean(e**2)), "slip": float(np.mean(np.abs(e) > 0.5 * lam))}
    sig_hi = float(core.sigma_from_snr(snr))
    sig_ex = phase_std_exact(snr)
    row = {
        "snr_db": float(snr_db), "n": n, "sigma_highsnr": sig_hi, "sigma_exact": sig_ex,
        "u": sig_hi / m["sigma_jump"], **res,
        "hmse": float(core.hmse(f, sig_hi)), "mmse": float(core.mmse(f, sig_hi)),
        "lattice": lattice.ml_mse(lattice.params(f, sig_hi))["mse"],
        "lattice_exact_sigma": lattice.ml_mse(lattice.params(f, sig_ex))["mse"],
    }
    print(f"SNR {snr_db:5.1f} dB  complex {core.db(res['complex_ml']['mse']):7.2f}  "
          f"phase {core.db(res['phase_ml']['mse']):7.2f}  lattice {core.db(row['lattice']):7.2f}", flush=True)
    return row


def main():
    jobs = [(float(s), 60000 + i) for i, s in enumerate(SNR_DB)]
    with ProcessPoolExecutor(max_workers=4) as ex:
        rows = sorted(ex.map(run_one, jobs), key=lambda r: r["snr_db"])
    m = core.moments(core.equal_hz())
    sim.save("e05_complex", {
        "rows": rows,
        "snr_threshold_per_tone_db": float(10 * math.log10(1 / (2 * m["sigma_jump"] ** 2))),
        "floor_m2": (2 * HALF) ** 2 / 12,
    })


if __name__ == "__main__":
    main()
