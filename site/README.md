# Reading edition

Typographic HTML for the threshold notes (Chinese + English).

```bash
# from repo root
python3 site/build_site.py
python3 -m http.server -d site 8765
# open http://127.0.0.1:8765/
```

| File | Role |
| --- | --- |
| `index.html` | Landing: CN / EN / explorer / Binder |
| `zh.html` / `en.html` | Typeset reading (KaTeX, sticky TOC, English figure labels) |
| `explore.html` | No-install slider on frozen equal-grid archive |
| `data/results.json` | Slim archive for the explorer |
| `assets/reading.css` | Typography |
| `build_site.py` | Rebuild from `notes/门限问题的进展.md` and `README.md` |

This folder is zipped into the GitHub Release asset `threshold-reading-site.zip`, published to the `gh-pages` branch by CI, and (after Pages is enabled) served at GitHub Pages.
