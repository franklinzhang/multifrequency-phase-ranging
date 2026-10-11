# Archived trials

This directory holds the **authoritative** numbers cited by the Chinese note and the interactive notebook.

## Canonical data (do not silently overwrite)

| File | Role |
| --- | --- |
| [`data/results.json`](data/results.json) | Per-design moments, trial MSE, slip rates, HMSE/MMSE, full-cell `zzb`, and half-window `zzb_half` |
| [`data/curves.npz`](data/curves.npz) | Cost-surface samples (`*_delta`, `*_J`) for equal / matched / geometric sets |
| [`data/bound_valleys.json`](data/bound_valleys.json) | Half-window valley fill segments and cumulative fractions |

- Random seed used for the archive: **7** (also recorded in `results.json`).
- The notes and [`notes/多频相位测距的门限.ipynb`](../notes/多频相位测距的门限.ipynb) **read** these files. They must not rewrite them.
- A later run of `sim_threshold.py` or `bound_valleys.py` is **not** a substitute for the archived numbers unless every difference is stated in the notes.

## Half-window bound (`zzb_half`)

Figures and the notebook use `sets[*].zzb_half` inside `results.json` and the cumulative fractions in `bound_valleys.json`. Treat those archived fields as authoritative. The current `sim_threshold.py` may not regenerate every half-window field in exactly the same post-processed form; do not assume a clean re-run reproduces the published curves without checking.

## Scripts

| Script | Purpose |
| --- | --- |
| `sim_threshold.py` | Monte Carlo MSE / slip rates / bound integrals (writes under `data/`) |
| `bound_valleys.py` | Valley-fill segment breakdown for the half-window bound |
| `make_figures.py` | Static PNGs for the Markdown twin → `notes/figures/` |

Install dependencies from the repository root: `pip install -r requirements.txt`.
