# Pre-public quality gate

Confirmed for repository visibility **public** (2026-10-11).

## 1. Storyteller writing (JIT, not obscure)

- Chinese canonical narrative: `notes/门限问题的进展.md` / Chinese notebook / `site/zh.html`
- Long abstract states what was done and answered (numbers, not adjectives)
- No Background / Preliminaries / Related Work dumping; wrap, CRB, ZZB, Weiss–Weinstein factor 12, Rife–Boorstyn appear only when first needed
- No metaphors (脚/膝/旋钮/创纪录/没有落点 removed)
- Problem-before-formula section titles; open questions have completion tests
- English digest + English notebook follow the same architecture

## 2. Refined typography

- `site/` HTML: Source Serif 4 + Noto Serif SC, IBM Plex Sans UI, measure ~42rem, sticky TOC, CN/EN switch, KaTeX, English figure labels
- Rebuild: `python3 site/build_site.py`
- Offline via Release asset `threshold-reading-site.zip`
- GitHub Pages workflow: `.github/workflows/pages.yml`

## 3. Interaction

- **No install:** `site/explore.html` — slider over frozen equal-grid archive
- **Full notebooks:** Binder links on `site/index.html` and README (`environment.yml`)
- Local: `pip install -r requirements.txt` then Jupyter
- GitHub’s built-in notebook preview is intentionally *not* treated as the interactive edition

## 4. Easy citation

- `LICENSE` = CC BY 4.0 (attribution required)
- `CITATION.cff` with `preferred-citation`, version `0.1.0`, release URL
- README BibTeX under License and citation
- Archive reuse rule: cite the repo and state values come from `computation/data/`

## Owner one-time switches (agent token cannot)

The cloud agent can push, tag, and attach Release assets, but **cannot** change
visibility or enable Pages (GitHub returns 403 for those admin APIs). After review:

1. **Public:** Settings → General → Danger Zone → Change repository visibility → Public  
   (or `gh repo edit franklinzhang/multifrequency-phase-ranging --visibility public`)
2. **Pages:** Settings → Pages → Build and deployment → Source: **GitHub Actions**  
   Then re-run workflow `Deploy reading site` (or push any commit on `main` / this branch).
3. **Merge PR #2** into `main` so the default branch matches the release tag tip.
