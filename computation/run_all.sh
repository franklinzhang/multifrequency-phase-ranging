#!/usr/bin/env bash
# Regenerate derived numbers, the article number table, and run tests + identities.
# Does not rewrite computation/data/v2 (the frozen archive). Set RUN_ARCHIVE=1 to
# rebuild e01–e09 and the figures; that is about an hour on four cores.
set -euo pipefail
cd "$(dirname "$0")"
python3 -m pytest -q tests
if [ "${RUN_ARCHIVE:-}" = 1 ]; then
  for e in e01_bounds e02_lattice_validation e03_collapse e04_bayes e05_complex e06_bias e07_layers e08_design e09_pareto; do
    echo "== $e"
    python3 "experiments/$e.py"
  done
  python3 figures_article.py
fi
python3 derive_remaining.py
python3 article_numbers.py
python3 verify_math_idents.py
