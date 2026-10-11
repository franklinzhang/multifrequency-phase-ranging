#!/usr/bin/env python3
"""Build typographic HTML reading pages from the Markdown notes."""

from __future__ import annotations

import re
import shutil
from pathlib import Path

import markdown
from markdown.extensions.toc import slugify_unicode

ROOT = Path(__file__).resolve().parents[1]
SITE = ROOT / "site"
NOTES = ROOT / "notes"
FIG_SRC = NOTES / "figures"
FIG_DST = SITE / "figures"


def md_to_html(text: str) -> str:
    # Normalize display math delimiters for better KaTeX coverage.
    text = re.sub(r"\\\[(.*?)\\\]", r"$$\1$$", text, flags=re.S)
    text = re.sub(r"(?<!\\)\\\((.*?)(?<!\\)\\\)", r"$\1$", text, flags=re.S)
    return markdown.markdown(
        text,
        extensions=[
            "pymdownx.superfences",
            "pymdownx.arithmatex",
            "tables",
            "toc",
            "sane_lists",
            "smarty",
        ],
        extension_configs={
            "pymdownx.arithmatex": {"generic": True},
            "toc": {
                "permalink": False,
                "slugify": slugify_unicode,
                "toc_depth": "2-3",
            },
        },
        output_format="html5",
    )


def extract_toc(html: str) -> str:
    m = re.search(r'<div class="toc">(.*?)</div>', html, flags=re.S)
    if not m:
        return "<strong>Contents</strong><ol></ol>"
    inner = m.group(1).strip()
    # markdown toc wraps a ul; keep it.
    return "<strong>Contents</strong>" + inner


def strip_toc(html: str) -> str:
    return re.sub(r'<div class="toc">.*?</div>\s*', "", html, count=1, flags=re.S)


def retarget_figures(html: str) -> str:
    html = html.replace("figures/", "figures/")
    html = html.replace("(figures/", "(figures/")
    # Markdown image paths like figures/fig02...
    html = html.replace('src="figures/', 'src="figures/')
    return html


PAGE = """<!DOCTYPE html>
<html lang="{lang}">
<head>
  <meta charset="utf-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1" />
  <title>{title}</title>
  <link rel="preconnect" href="https://fonts.googleapis.com" />
  <link rel="preconnect" href="https://fonts.gstatic.com" crossorigin />
  <link href="https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@400;500&family=IBM+Plex+Sans:wght@400;500;600&family=Noto+Sans+SC:wght@400;500;600&family=Noto+Serif+SC:wght@500;600&family=Source+Serif+4:opsz,wght@8..60,400;8..60,600&display=swap" rel="stylesheet" />
  <link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/katex@0.16.11/dist/katex.min.css" />
  <link rel="stylesheet" href="assets/reading.css" />
</head>
<body>
  <div class="shell">
    <aside class="toc">{toc}</aside>
    <main class="article">
      <div class="topbar">
        <div class="lang-switch">
          <a href="zh.html" class="{zh_active}">中文</a>
          <a href="en.html" class="{en_active}">English</a>
        </div>
        <div class="meta"><a href="index.html">Home</a> · <a href="../README.md">Repo README</a></div>
      </div>
      {body}
      <div class="footnote-bar">
        CC BY 4.0 · Archived numbers in <code>computation/data/</code> · Interactive notebooks remain available in <code>notes/</code>
      </div>
    </main>
  </div>
  <script defer src="https://cdn.jsdelivr.net/npm/katex@0.16.11/dist/katex.min.js"></script>
  <script defer src="https://cdn.jsdelivr.net/npm/katex@0.16.11/dist/contrib/auto-render.min.js"
    onload="renderMathInElement(document.body, {{
      delimiters: [
        {{left: '$$', right: '$$', display: true}},
        {{left: '$', right: '$', display: false}},
        {{left: '\\\\(', right: '\\\\)', display: false}},
        {{left: '\\\\[', right: '\\\\]', display: true}}
      ],
      throwOnError: false
    }});"></script>
</body>
</html>
"""

HOME = """<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1" />
  <title>Multi-Frequency Phase Ranging — Threshold Notes</title>
  <link rel="preconnect" href="https://fonts.googleapis.com" />
  <link rel="preconnect" href="https://fonts.gstatic.com" crossorigin />
  <link href="https://fonts.googleapis.com/css2?family=IBM+Plex+Sans:wght@400;500;600&family=Noto+Sans+SC:wght@400;500;600&family=Noto+Serif+SC:wght@600&family=Source+Serif+4:opsz,wght@8..60,600&display=swap" rel="stylesheet" />
  <link rel="stylesheet" href="assets/reading.css" />
</head>
<body>
  <div class="home">
    <div class="home-card">
      <h1>The Threshold of Multi-Frequency Phase Ranging</h1>
      <p>Read the consolidation in the browser with typeset math and English figure labels. Notebooks remain available for archive interaction.</p>
      <div class="choices">
        <a class="choice" href="zh.html">
          <strong>中文阅读版</strong>
          <span>终局感叙事 · 精排 HTML · KaTeX</span>
        </a>
        <a class="choice" href="en.html">
          <strong>English reading edition</strong>
          <span>Storyteller digest · typographic HTML · KaTeX</span>
        </a>
        <a class="choice" href="../notes/多频相位测距的门限.ipynb">
          <strong>中文 Jupyter（交互）</strong>
          <span>滑条查看归档噪声比；需本地或支持 widgets 的环境</span>
        </a>
        <a class="choice" href="../notes/threshold.ipynb">
          <strong>English Jupyter (interactive)</strong>
          <span>Archive slider; run locally with requirements.txt</span>
        </a>
      </div>
      <div class="footnote-bar">CC BY 4.0 · If GitHub links 404 while logged out, the repository is private — open the Release zip or clone locally.</div>
    </div>
  </div>
</body>
</html>
"""


def prepare_zh_source(text: str) -> str:
    # Point image paths at site/figures and keep captions.
    text = text.replace("](figures/", "](figures/")
    # Drop the repo-path logistics table clutter slightly? Keep — useful.
    return text


def prepare_en_source(text: str) -> str:
    # Prefer a reading-focused English page from README, without license bulk at top.
    # Keep Open the notes, abstract, and technical sections; trim license/bibtex to short note.
    cut = text.find("## License and citation")
    if cut >= 0:
        refs = text.find("## References")
        head = text[:cut].rstrip() + "\n\n"
        if refs >= 0:
            text = head + text[refs:]
        else:
            text = head
    text = "# The Threshold of Multi-Frequency Phase Ranging\n\n" + text.split("\n", 2)[-1]
    # Remove duplicate H1 if README already starts with it after strip — handled above.
    # Fix relative notebook links for site root.
    text = text.replace("(notes/", "(../notes/")
    text = text.replace("](LICENSE)", "](../LICENSE)")
    text = text.replace("](CITATION.cff)", "](../CITATION.cff)")
    text = text.replace("](computation/", "](../computation/")
    return text


def main() -> None:
    FIG_DST.mkdir(parents=True, exist_ok=True)
    for png in FIG_SRC.glob("*.png"):
        shutil.copy2(png, FIG_DST / png.name)

    zh_md = prepare_zh_source((NOTES / "门限问题的进展.md").read_text(encoding="utf-8"))
    en_md = prepare_en_source((ROOT / "README.md").read_text(encoding="utf-8"))

    zh_body = retarget_figures(md_to_html(zh_md))
    en_body = retarget_figures(md_to_html(en_md))

    # Force TOC generation by ensuring toc extension runs; markdown places toc where [TOC] is.
    if '[TOC]' not in zh_md:
        zh_body = md_to_html("[TOC]\n\n" + zh_md)
        zh_body = retarget_figures(zh_body)
    if '[TOC]' not in en_md:
        en_body = md_to_html("[TOC]\n\n" + en_md)
        en_body = retarget_figures(en_body)

    zh_toc = extract_toc(zh_body)
    en_toc = extract_toc(en_body)
    zh_body = strip_toc(zh_body)
    en_body = strip_toc(en_body)

    (SITE / "zh.html").write_text(
        PAGE.format(
            lang="zh-CN",
            title="多频相位测距的门限",
            toc=zh_toc,
            body=zh_body,
            zh_active="active",
            en_active="",
        ),
        encoding="utf-8",
    )
    (SITE / "en.html").write_text(
        PAGE.format(
            lang="en",
            title="The Threshold of Multi-Frequency Phase Ranging",
            toc=en_toc,
            body=en_body,
            zh_active="",
            en_active="active",
        ),
        encoding="utf-8",
    )
    (SITE / "index.html").write_text(HOME, encoding="utf-8")
    print("wrote", SITE / "index.html")
    print("wrote", SITE / "zh.html")
    print("wrote", SITE / "en.html")


if __name__ == "__main__":
    main()
