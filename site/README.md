# Reading edition

Typographic HTML for the threshold notes (Chinese + English).

```bash
# from repo root
python3 site/build_site.py
open site/index.html   # or any static file server
```

- `index.html` — language picker  
- `zh.html` / `en.html` — full reading pages (KaTeX, TOC, English figure labels)  
- `build_site.py` — rebuild from `notes/门限问题的进展.md` and `README.md`

This folder is what gets zipped into the GitHub Release asset `threshold-reading-site.zip`.
