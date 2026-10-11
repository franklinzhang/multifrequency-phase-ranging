"""Monte Carlo helpers shared by the experiments."""

from __future__ import annotations

import json
import math
import subprocess
from pathlib import Path

import numpy as np

from .core import C, wrap

DATA = Path(__file__).resolve().parents[1] / "data" / "v2"


def phases(rng, f, d0, sigma, n, offset=None):
    """Wrapped phases for n trials. d0 scalar or (n,), sigma scalar or (N,)."""
    f = np.asarray(f, dtype=float)
    d0 = np.broadcast_to(np.asarray(d0, dtype=float), (n,))
    e = rng.normal(0.0, 1.0, size=(n, f.size)) * np.broadcast_to(np.asarray(sigma, dtype=float), f.shape)
    ph = (2.0 * math.pi / C) * d0[:, None] * f[None, :] + e
    if offset is not None:
        ph = ph + np.asarray(offset, dtype=float)
    return wrap(ph)


def tones(rng, f, d0, snr, n, channel=None):
    """Complex tones y_i = h_i exp(j2πf_i d/c) + n_i with E|n|² = 1/snr (per tone, unit |h|)."""
    f = np.asarray(f, dtype=float)
    d0 = np.broadcast_to(np.asarray(d0, dtype=float), (n,))
    h = np.ones(f.size, dtype=complex) if channel is None else np.asarray(channel, dtype=complex)
    s = h[None, :] * np.exp(1j * (2.0 * math.pi / C) * d0[:, None] * f[None, :])
    noise = (rng.normal(size=(n, f.size)) + 1j * rng.normal(size=(n, f.size))) * math.sqrt(0.5 / snr)
    return s + noise


def git_hash() -> str:
    try:
        return subprocess.check_output(["git", "rev-parse", "--short", "HEAD"], text=True).strip()
    except Exception:
        return "unknown"


def save(name: str, payload: dict) -> Path:
    DATA.mkdir(parents=True, exist_ok=True)
    payload = {"git": git_hash(), **payload}
    path = DATA / f"{name}.json"
    path.write_text(json.dumps(payload, indent=1, default=_np), encoding="utf-8")
    print("wrote", path, flush=True)
    return path


def load(name: str) -> dict:
    return json.loads((DATA / f"{name}.json").read_text(encoding="utf-8"))


def _np(o):
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    raise TypeError(type(o))
