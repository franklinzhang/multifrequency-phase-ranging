# Changelog

## [v0.2.0] — 2026-10-16

Master-edition overhaul: complete closure of all five open research questions, self-contained inline mathematical derivations, annotated log-log overview diagram, and four interactive figures with one-click bilingual switching.

- **Long Abstract & Annotated Master Overview Diagram:** Added bilingual extended abstract (`## 摘要` / `## Abstract`) paired with `figures/fig00_overview.png` (`2000×1240` px) illustrating the four regimes, vertical–horizontal coupling law \(R = 1 + N\pi^2/\sigma_{\mathrm{jump}}^2\), leaving criterion \(2P_e(\delta_m)\delta_m^2 = \mathrm{HMSE}\), three-layer threshold risk split, and collinear Voronoi plateau.
- **Closed All Five Open Questions ([`computation/open_questions.py`](computation/open_questions.py), [`computation/data/open_questions.json`](computation/data/open_questions.json)):**
  1. **Worst-valley tone placement (8.1):** Exposed the \(4.494\,\mathrm{m}\) bimodal beat trap (`2400..2403, 2467..2470 MHz`) and proved by exhaustive search over all \(\binom{34}{3}=5984\) symmetric sets in \([2400,2470]\,\mathrm{MHz}\) that shifting the innermost pair inward by \(2\,\mathrm{MHz}\) (`2400, 2401, 2402, 2405, 2465, 2468, 2469, 2470 MHz`) delays \(\sigma_{\mathrm{leave}}\) by **\(+49.1\%\)** (from \(0.0170\) to \(0.0254\,\mathrm{rad}\)) while keeping every far valley \(\ge 4.00\,J(\lambda)\).
  2. **Three-way attribution of the \(4.82\,\mathrm{dB}\) threshold gap (8.2):** Using exact piecewise-parabolic quadrature of \(C(d)\), decomposed the gap at \(\sigma_{\mathrm{jump}}\) into prior boundary truncation (\(0.24\,\mathrm{dB}\)), hard-rounding overshoot of MLE vs. smooth Bayesian posterior mean (\(1.20\,\mathrm{dB}\), with \(\hat d_{\mathrm{Bayes}} = -24.22\,\mathrm{dB}\) landing on \(\mathrm{MMSE}\)), and pairwise-test slack of ZZB (\(3.38\,\mathrm{dB}\)).
  3. **Exact wrapped-normal likelihood ratio on \(\pi\)-residual valleys (8.3):** 8-fold LLR convolution shows wrapping across \(\pm\pi\) raises \(P_e^{\mathrm{WN}}\) by \(21\%\text{–}24\%\) on the \(10.6\text{–}11.3\,\mathrm{m}\) cluster at \(\sigma=0.836\,\mathrm{rad}\), widening the moment-matched vs. equal-grid bound gap from \(1.38\,\mathrm{dB}\) to \(1.70\,\mathrm{dB}\).
  4. **Complex baseband verification & envelope weighting (8.4):** Verified that the first threshold is unchanged up to \(3.5\,\sigma_{\mathrm{jump}}\), while at \(6\,\sigma_{\mathrm{jump}}\) (\(\mathrm{SNR}\approx 3\,\mathrm{dB}\)) envelope weighting \(|y_i|\) in complex MLE saves \(5.37\,\mathrm{dB}\) over unweighted phase least squares.
  5. **Collinear Voronoi multi-wavelength series (8.5):** Proved that \(2\lambda^2\sum_{k\ge 1}(2k-1)Q\bigl((2k-1)\sigma_{\mathrm{jump}}/\sigma\bigr)\) matches the \(1.5\text{–}6\,\sigma_{\mathrm{jump}}\) plateau within \(0.07\,\mathrm{dB}\) and converges analytically to \(\mathrm{MMSE}\).
- **Four Interactive Figures & One-Click Bilingual UX (`site/`):** Added live interactive canvases for (1) Four-Regime Log-Log Anatomy, (2) 8-Tone Synthesizer & Cost Surface Lab, (3) Bayesian Sigmoid vs. Hard MLE Staircase, and (4) Bellini–Tartara Valley Fill & Exact Wrapped LRT, with section-preserving `中文 / English` switching.

## [v0.1.0] — 2026-10-11

First research-style GitHub Release of the threshold account (notes + archived trials).

- **Read:** typographic HTML in [`site/`](site/) (Release asset `threshold-reading-site.zip`; CI → `gh-pages`)
- **Explore:** no-install [`site/explore.html`](site/explore.html) on frozen equal-grid archive
- **Interact:** Binder + [`environment.yml`](environment.yml); Chinese and English notebooks
- Chinese storyteller note: [`notes/门限问题的进展.md`](notes/门限问题的进展.md)
- English digest: [`README.md`](README.md)
- Frozen trials: [`computation/data/`](computation/data/) (seed 7)
- License CC BY 4.0; [`CITATION.cff`](CITATION.cff); gate checklist [`QUALITY.md`](QUALITY.md)
