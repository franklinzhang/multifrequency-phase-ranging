#!/usr/bin/env bash
# Regenerate every number and figure of the article. About an hour on four cores.
set -euo pipefail
cd "$(dirname "$0")"
python3 -m pytest -q tests
for e in e01_bounds e02_lattice_validation e03_collapse e04_bayes e05_complex e06_bias e07_layers e08_design e09_pareto; do
  echo "== $e"
  python3 "experiments/$e.py"
done
python3 figures_article.py
python3 article_numbers.py
