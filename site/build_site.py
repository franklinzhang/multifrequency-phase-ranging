#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Build typographic HTML reading pages and interactive figure laboratory from the Markdown notes."""

from __future__ import annotations

import html
import json
import re
import shutil
import unicodedata
from pathlib import Path

try:
    import markdown
    from markdown.extensions.toc import slugify_unicode
    HAS_MARKDOWN_PKG = True
except ImportError:
    HAS_MARKDOWN_PKG = False

ROOT = Path(__file__).resolve().parents[1]
SITE = ROOT / "site"
NOTES = ROOT / "notes"
FIG_SRC = NOTES / "figures"
FIG_DST = SITE / "figures"
DATA_DIR = ROOT / "computation" / "data"
DATA_DST = SITE / "data"
REPO = "https://github.com/franklinzhang/multifrequency-phase-ranging"
REPO_BLOB = f"{REPO}/blob/main"
BINDER_BASE = "https://mybinder.org/v2/gh/franklinzhang/multifrequency-phase-ranging/main"
BINDER_EN = f"{BINDER_BASE}?labpath=notes%2Fthreshold.ipynb"
BINDER_ZH = f"{BINDER_BASE}?labpath=notes%2F%E5%A4%9A%E9%A2%91%E7%9B%B8%E4%BD%8D%E6%B5%8B%E8%B7%9D%E7%9A%84%E9%97%A8%E9%99%90.ipynb"


def slugify_heading(text: str) -> str:
    s = re.sub(r"\\\((.*?)\\\)", r"\1", text)
    s = re.sub(r"[\\{}`*_]", "", s)
    s = unicodedata.normalize("NFKC", s).strip().lower()
    out = []
    for ch in s:
        cat = unicodedata.category(ch)
        if cat.startswith("L") or cat.startswith("N"):
            out.append(ch)
        elif ch in (" ", "-", "\t"):
            out.append("-")
    res = re.sub(r"-+", "-", "".join(out)).strip("-")
    return res or "section"


def inline_format(text: str) -> str:
    math_spans = []
    def save_math(m):
        idx = len(math_spans)
        math_spans.append(f'<span class="arithmatex">\\({m.group(1)}\\)</span>')
        return f"\x00MATH{idx}\x00"

    text = re.sub(r"\\\((.*?)\\\)", save_math, text)

    code_spans = []
    def save_code(m):
        idx = len(code_spans)
        code_spans.append(f"<code>{html.escape(m.group(1))}</code>")
        return f"\x00CODE{idx}\x00"

    text = re.sub(r"`([^`]+)`", save_code, text)
    text = re.sub(r"!\[([^\]]*)\]\(([^)]+)\)", r'<img alt="\1" src="\2">', text)
    text = re.sub(r"\[([^\]]+)\]\(([^)]+)\)", r'<a href="\2">\1</a>', text)
    text = re.sub(r"\*\*(.+?)\*\*", r"<strong>\1</strong>", text)
    text = re.sub(r"\*(.+?)\*", r"<em>\1</em>", text)

    for idx, val in enumerate(code_spans):
        text = text.replace(f"\x00CODE{idx}\x00", val)
    for idx, val in enumerate(math_spans):
        text = text.replace(f"\x00MATH{idx}\x00", val)
    return text


def stdlib_md_to_html_and_toc(md_text: str) -> tuple[str, str]:
    blocks = re.split(r"\n\s*\n", md_text.strip())
    html_parts = []
    toc_items = []

    for blk in blocks:
        b = blk.strip()
        if not b:
            continue
        if b.startswith(r"\[") and b.endswith(r"\]"):
            inner = b[2:-2].strip()
            html_parts.append(f'<div class="arithmatex">\\[\n{inner}\n\\]</div>')
            continue
        m_h = re.match(r"^(#{1,3})\s+(.*)$", b)
        if m_h:
            level = len(m_h.group(1))
            raw_title = m_h.group(2).strip()
            slug = slugify_heading(raw_title)
            title_html = inline_format(raw_title)
            toc_title = re.sub(r'<span class="arithmatex">(.*?)</span>', r"\1", title_html)
            if level in (2, 3):
                toc_items.append((level, slug, toc_title))
            html_parts.append(f'<h{level} id="{slug}">{title_html}</h{level}>')
            continue
        lines = [ln.strip() for ln in b.splitlines() if ln.strip()]
        if len(lines) >= 2 and lines[0].startswith("|") and re.match(r"^\|[\s\-:|]+\|$", lines[1]):
            headers = [c.strip() for c in lines[0].strip("|").split("|")]
            t_lines = ["<table>", "<thead>", "<tr>"]
            for h_cell in headers:
                t_lines.append(f"<th>{inline_format(h_cell)}</th>")
            t_lines.extend(["</tr>", "</thead>", "<tbody>"])
            for row_ln in lines[2:]:
                if not row_ln.startswith("|"):
                    continue
                cells = [c.strip() for c in row_ln.strip("|").split("|")]
                t_lines.append("<tr>")
                for c_cell in cells:
                    t_lines.append(f"<td>{inline_format(c_cell)}</td>")
                t_lines.append("</tr>")
            t_lines.extend(["</tbody>", "</table>"])
            html_parts.append("\n".join(t_lines))
            continue
        if all(ln.startswith("- ") for ln in lines):
            ul = ["<ul>"]
            for ln in lines:
                ul.append(f"<li>{inline_format(ln[2:].strip())}</li>")
            ul.append("</ul>")
            html_parts.append("\n".join(ul))
            continue
        if re.match(r"^1\.\s+", lines[0]):
            items = re.split(r"\n(?=\d+\.\s+)", b)
            ol = ["<ol>"]
            for it in items:
                it_clean = re.sub(r"^\d+\.\s+", "", it.strip())
                it_clean = re.sub(
                    r"\\\[(.*?)\\\]",
                    lambda m: f'<div class="arithmatex">\\[\n{m.group(1).strip()}\n\\]</div>',
                    it_clean,
                    flags=re.S,
                )
                ol.append(f"<li>{inline_format(it_clean)}</li>")
            ol.append("</ol>")
            html_parts.append("\n".join(ol))
            continue
        html_parts.append(f"<p>{inline_format(b)}</p>")

    toc_lines = ["<strong>Contents</strong><ul>"]
    in_sub = False
    for i, (lvl, slug, title) in enumerate(toc_items):
        nxt_lvl = toc_items[i + 1][0] if i + 1 < len(toc_items) else 2
        if lvl == 2:
            if in_sub:
                toc_lines.append("</ul></li>")
                in_sub = False
            if nxt_lvl == 3:
                toc_lines.append(f'<li><a href="#{slug}">{title}</a><ul>')
                in_sub = True
            else:
                toc_lines.append(f'<li><a href="#{slug}">{title}</a></li>')
        elif lvl == 3:
            toc_lines.append(f'<li><a href="#{slug}">{title}</a></li>')
    if in_sub:
        toc_lines.append("</ul></li>")
    toc_lines.append("</ul>")
    return "\n".join(html_parts), "\n".join(toc_lines)


def inject_interactive_figures(body_html: str, lang: str) -> str:
    """Inject the 4 interactive figure containers into the article HTML at their natural sections."""
    # 1. Overview interactive explorer after fig02_equal_mse.png paragraph
    m1 = re.search(r'(<p><img[^>]*fig02_equal_mse\.png[^>]*></p>\s*<p>.*?</p>)', body_html, flags=re.S)
    if m1:
        card1 = f'\n<div class="ifig-card" data-ifig="overview" data-lang="{lang}"></div>\n'
        body_html = body_html[: m1.end()] + card1 + body_html[m1.end() :]

    # 2. Bayes vs MLE staircase explorer at end of Section 5 (before <h2 id="6-...)
    m2 = re.search(r'(<h2 id="6-[^"]*">)', body_html)
    if m2:
        card2 = f'<div class="ifig-card" data-ifig="bayes" data-lang="{lang}"></div>\n'
        body_html = body_html[: m2.start()] + card2 + body_html[m2.start() :]

    # 3. Valley fill & Exact Wrapped LRT explorer at end of Section 6 (before <h2 id="7-...)
    m3 = re.search(r'(<h2 id="7-[^"]*">)', body_html)
    if m3:
        card3 = f'<div class="ifig-card" data-ifig="valleys" data-lang="{lang}"></div>\n'
        body_html = body_html[: m3.start()] + card3 + body_html[m3.start() :]

    # 4. Live 8-tone synthesizer & worst-valley optimizer at end of Section 7 (before <h2 id="8-...)
    m4 = re.search(r'(<h2 id="8-[^"]*">)', body_html)
    if m4:
        card4 = f'<div class="ifig-card" data-ifig="synthesizer" data-lang="{lang}"></div>\n'
        body_html = body_html[: m4.start()] + card4 + body_html[m4.start() :]

    return body_html


def rewrite_repo_links(html_str: str) -> str:
    replacements = [
        ('href="../README.md"', f'href="{REPO_BLOB}/README.md"'),
        ('href="../LICENSE"', f'href="{REPO_BLOB}/LICENSE"'),
        ('href="../CITATION.cff"', f'href="{REPO_BLOB}/CITATION.cff"'),
        ('href="../computation/README.md"', f'href="{REPO_BLOB}/computation/README.md"'),
        ('href="../computation/open_questions.py"', f'href="{REPO_BLOB}/computation/open_questions.py"'),
        ('href="../computation/data/open_questions.json"', f'href="{REPO_BLOB}/computation/data/open_questions.json"'),
        ('href="../computation/data/"', f'href="{REPO_BLOB}/computation/data"'),
        ('href="../computation/"', f'href="{REPO_BLOB}/computation"'),
        ('href="../notes/', f'href="{REPO_BLOB}/notes/'),
        ('href="README.md"', f'href="{REPO_BLOB}/README.md"'),
        ('href="LICENSE"', f'href="{REPO_BLOB}/LICENSE"'),
        ('href="CITATION.cff"', f'href="{REPO_BLOB}/CITATION.cff"'),
        ('href="computation/README.md"', f'href="{REPO_BLOB}/computation/README.md"'),
        ('href="computation/open_questions.py"', f'href="{REPO_BLOB}/computation/open_questions.py"'),
        ('href="computation/data/open_questions.json"', f'href="{REPO_BLOB}/computation/data/open_questions.json"'),
        ('href="computation/data/"', f'href="{REPO_BLOB}/computation/data"'),
        ('href="computation/"', f'href="{REPO_BLOB}/computation"'),
        ('href="notes/', f'href="{REPO_BLOB}/notes/'),
    ]
    for old, new in replacements:
        html_str = html_str.replace(old, new)
    html_str = html_str.replace(
        "https://franklinzhang.github.io/multifrequency-phase-ranging/figures/",
        "figures/",
    )
    return html_str


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
          <a href="zh.html" class="{zh_active}" title="切换为中文（保持当前章节位置）">中文</a>
          <a href="en.html" class="{en_active}" title="Switch to English (preserving section)">English</a>
        </div>
        <div class="meta"><a href="index.html">{nav_home}</a> · <a href="explore.html">{nav_lab}</a> · <a href="https://github.com/franklinzhang/multifrequency-phase-ranging">GitHub</a></div>
      </div>
      {body}
      <div class="footnote-bar">
        CC BY 4.0 · Cite via <a href="https://github.com/franklinzhang/multifrequency-phase-ranging/blob/main/CITATION.cff">CITATION.cff</a> · Archived numbers in <code>computation/data/</code>
      </div>
    </main>
  </div>
  <script defer src="assets/interactive.js"></script>
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
<html lang="zh-CN">
<head>
  <meta charset="utf-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1" />
  <title>多频相位测距的门限 · The Threshold of Multi-Frequency Phase Ranging</title>
  <meta name="description" content="Coherent range estimation from wrapped multi-frequency phases: the first threshold σ_jump = π√S/f̄, collinear Voronoi series, half-window Ziv–Zakai bound, piecewise-parabolic Bayesian risk, and worst-valley tone placement. In Chinese and English." />
  <link rel="canonical" href="https://franklinzhang.github.io/multifrequency-phase-ranging/" />
  <link rel="preconnect" href="https://fonts.googleapis.com" />
  <link rel="preconnect" href="https://fonts.gstatic.com" crossorigin />
  <link href="https://fonts.googleapis.com/css2?family=IBM+Plex+Sans:wght@400;500;600&family=Noto+Sans+SC:wght@400;500;600&family=Noto+Serif+SC:wght@600&family=Source+Serif+4:opsz,wght@8..60,600&display=swap" rel="stylesheet" />
  <link rel="stylesheet" href="assets/reading.css" />
</head>
<body>
  <div class="home">
    <div class="home-card">
      <div class="topbar" style="position:static;margin-bottom:1.25rem;">
        <div class="lang-switch">
          <button type="button" class="active" data-set-lang="zh">中文</button>
          <button type="button" data-set-lang="en">English</button>
        </div>
        <div class="meta"><a href="explore.html" data-i18n-zh="交互实验室 (4 组交互图)" data-i18n-en="Interactive Lab (4 Figures)">交互实验室 (4 组交互图)</a> · <a href="https://github.com/franklinzhang/multifrequency-phase-ranging">GitHub</a></div>
      </div>
      <p class="edition-label" data-i18n-zh="多频相位测距的门限 · 双语定本" data-i18n-en="The Threshold of Multi-Frequency Phase Ranging · Bilingual Monograph">多频相位测距的门限 · 双语定本</p>
      <h1 data-i18n-zh="多频相位测距的门限" data-i18n-en="The Threshold of Multi-Frequency Phase Ranging">多频相位测距的门限</h1>
      <p class="lede" data-i18n-zh="从包裹多频相位观测模型出发，用参数估计的经典数学工具将双对数图上的平行双线、第一道门限尺度 σ_jump = π√S/f̄、共线多波长 Voronoi 级数、半窗谷填充 Ziv–Zakai 下界、分段抛物线贝叶斯风险与同带宽最坏谷排频彻底算穿。" data-i18n-en="Starting from wrapped multi-frequency phase observations, this monograph applies classical parameter-estimation tools to completely resolve the parallel carrier/bandwidth lines, the first threshold scale σ_jump = π√S/f̄, the collinear Voronoi series, the half-window Ziv–Zakai bound, piecewise-parabolic Bayesian risk, and worst-valley tone placement.">从包裹多频相位观测模型出发，用参数估计的经典数学工具将双对数图上的平行双线、第一道门限尺度 σ_jump = π√S/f̄、共线多波长 Voronoi 级数、半窗谷填充 Ziv–Zakai 下界、分段抛物线贝叶斯风险与同带宽最坏谷排频彻底算穿。</p>

      <div class="contrib-grid">
        <div class="contrib-item">
          <strong data-i18n-zh="1. 纵横尺度锁定律" data-i18n-en="1. Vertical–Horizontal Coupling">1. 纵横尺度锁定律</strong>
          <span data-i18n-zh="导出门限相位尺度 σ_jump = π√S/f̄，证明其与两线落差 R = F/S 由同一无量纲矩锁定：R = 1 + Nπ²/σ_jump²。" data-i18n-en="Derives σ_jump = π√S/f̄ and proves it is locked to the vertical gap R = F/S via R = 1 + Nπ²/σ_jump².">导出门限相位尺度 σ_jump = π√S/f̄，证明其与两线落差 R = F/S 由同一无量纲矩锁定：R = 1 + Nπ²/σ_jump²。</span>
        </div>
        <div class="contrib-item">
          <strong data-i18n-zh="2. 共线 Voronoi 级数" data-i18n-en="2. Collinear Voronoi Bridge">2. 共线 Voronoi 级数</strong>
          <span data-i18n-zh="证明近处 kλ 谷沿 f_i − f̄ 共线，Voronoi 级数在带宽平台上误差 ≤ 0.07 dB，黎曼和解析收敛于 MMSE。" data-i18n-en="Proves near kλ valleys are collinear along f_i − f̄; their Voronoi series matches the plateau within 0.07 dB and converges to MMSE.">证明近处 kλ 谷沿 f_i − f̄ 共线，Voronoi 级数在带宽平台上误差 ≤ 0.07 dB，黎曼和解析收敛于 MMSE。</span>
        </div>
        <div class="contrib-item">
          <strong data-i18n-zh="3. 门限 4.82 dB 三方拆解" data-i18n-en="3. 4.82 dB Three-Way Split">3. 门限 4.82 dB 三方拆解</strong>
          <span data-i18n-zh="利用等曲率分段抛物线算出精确贝叶斯后验均值，拆清先验截断 (0.24 dB)、MLE 硬取整过冲 (1.20 dB) 与 ZZB 松弛 (3.38 dB)。" data-i18n-en="Uses exact piecewise-parabolic quadrature to split the 4.82 dB gap into prior truncation (0.24 dB), MLE overshoot (1.20 dB), and ZZB slack (3.38 dB).">利用等曲率分段抛物线算出精确贝叶斯后验均值，拆清先验截断 (0.24 dB)、MLE 硬取整过冲 (1.20 dB) 与 ZZB 松弛 (3.38 dB)。</span>
        </div>
        <div class="contrib-item">
          <strong data-i18n-zh="4. 同带宽最坏谷排频 (+49.1%)" data-i18n-en="4. Worst-Valley Tone Placement (+49.1%)">4. 同带宽最坏谷排频 (+49.1%)</strong>
          <span data-i18n-zh="揭示极端双簇的拍频远谷陷阱，给出内侧微错开 2 MHz 的穷举最优解，推迟门限 49.1% 且消除拍频远谷。" data-i18n-en="Exposes the bimodal beat trap and proves a 2 MHz inner stagger delays σ_leave by 49.1% while keeping all far valleys ≥ 4 J(λ).">揭示极端双簇的拍频远谷陷阱，给出内侧微错开 2 MHz 的穷举最优解，推迟门限 49.1% 且消除拍频远谷。</span>
        </div>
      </div>

      <div class="choices">
        <a class="choice primary" href="zh.html">
          <strong data-i18n-zh="阅读中文版全文 (含 4 组交互图)" data-i18n-en="Read Chinese Edition (with 4 Interactive Figures)">阅读中文版全文 (含 4 组交互图)</strong>
          <span data-i18n-zh="长摘要总览图 · 自包含数学推导 · 三组频点归档试验 · 完整闭环定本" data-i18n-en="Annotated Master Abstract · Self-Contained Derivations · Archived Trials">长摘要总览图 · 自包含数学推导 · 三组频点归档试验 · 完整闭环定本</span>
        </a>
        <a class="choice primary" href="en.html">
          <strong data-i18n-zh="阅读英文版全文 (English Edition)" data-i18n-en="Read English Edition (with 4 Interactive Figures)">阅读英文版全文 (English Edition)</strong>
          <span data-i18n-zh="Extended Abstract · Self-Contained Mathematical Tools · Complete Closed Account" data-i18n-en="Extended Abstract · Self-Contained Mathematical Tools · Complete Closed Account">Extended Abstract · Self-Contained Mathematical Tools · Complete Closed Account</span>
        </a>
        <a class="choice" href="explore.html">
          <strong data-i18n-zh="交互式图表实验室 (Interactive Lab)" data-i18n-en="Interactive Figure Laboratory">交互式图表实验室 (Interactive Lab)</strong>
          <span data-i18n-zh="集中体验全部 4 组交互式图表：双对数四区解剖、八频点排频合成器、贝叶斯软硬阶梯、谷填充与包裹似然比" data-i18n-en="Explore all 4 interactive figures: Log-Log Anatomy, 8-Tone Synthesizer Lab, Bayes vs MLE Staircase & Valley-Fill LRT">集中体验全部 4 组交互式图表：双对数四区解剖、八频点排频合成器、贝叶斯软硬阶梯、谷填充与包裹似然比</span>
        </a>
        <a class="choice" href="{binder_zh}">
          <strong data-i18n-zh="Jupyter 交互式笔记本 (Binder)" data-i18n-en="Jupyter Notebooks on Binder">Jupyter 交互式笔记本 (Binder)</strong>
          <span data-i18n-zh="在浏览器中直接运行归档数据与复现脚本（中英双语 Notebook）" data-i18n-en="Run the archived trials and reproduction scripts live in your browser">在浏览器中直接运行归档数据与复现脚本（中英双语 Notebook）</span>
        </a>
      </div>
      <div class="footnote-bar">CC BY 4.0 · Cite <a href="https://github.com/franklinzhang/multifrequency-phase-ranging/blob/main/CITATION.cff">CITATION.cff</a> · Source <a href="https://github.com/franklinzhang/multifrequency-phase-ranging">GitHub</a> · Release <a href="https://github.com/franklinzhang/multifrequency-phase-ranging/releases">Releases</a></div>
    </div>
  </div>
  <script>
    (function () {
      const btns = document.querySelectorAll("[data-set-lang]");
      function applyLang(lang) {
        document.documentElement.lang = lang === "zh" ? "zh-CN" : "en";
        btns.forEach(b => b.classList.toggle("active", b.getAttribute("data-set-lang") === lang));
        document.querySelectorAll("[data-i18n-zh]").forEach(el => {
          el.textContent = el.getAttribute(lang === "zh" ? "data-i18n-zh" : "data-i18n-en");
        });
      }
      btns.forEach(b => b.addEventListener("click", () => applyLang(b.getAttribute("data-set-lang"))));
    })();
  </script>
</body>
</html>
""".replace("{binder_zh}", BINDER_ZH).replace("{binder_en}", BINDER_EN)

EXPLORE = """<!DOCTYPE html>
<html lang="zh-CN">
<head>
  <meta charset="utf-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1" />
  <title>交互式图表实验室 · Interactive Figure Laboratory</title>
  <link rel="preconnect" href="https://fonts.googleapis.com" />
  <link rel="preconnect" href="https://fonts.gstatic.com" crossorigin />
  <link href="https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@400;500&family=IBM+Plex+Sans:wght@400;500;600&family=Noto+Sans+SC:wght@400;500;600&family=Source+Serif+4:opsz,wght@8..60,600&display=swap" rel="stylesheet" />
  <link rel="stylesheet" href="assets/reading.css" />
  <style>
    .explore-wrap { max-width: 64rem; margin: 0 auto; padding: 1.75rem 1rem 4rem; }
  </style>
</head>
<body>
  <div class="explore-wrap">
    <div class="topbar">
      <div class="lang-switch">
        <button type="button" class="active" data-lab-lang="zh">中文</button>
        <button type="button" data-lab-lang="en">English</button>
      </div>
      <div class="meta">
        <a href="index.html">Home</a> · <a href="zh.html">中文全文</a> · <a href="en.html">English Text</a> · <a href="https://github.com/franklinzhang/multifrequency-phase-ranging">GitHub</a>
      </div>
    </div>
    <div class="article" style="padding: 1.75rem 2rem 2.25rem;">
      <h1 id="lab-title" style="margin-bottom:0.3rem;">多频相位测距 · 交互式图表实验室</h1>
      <p class="lede" id="lab-sub" style="margin-bottom:1.5rem;">基于 <code>computation/data/</code> 归档试验（随机种子 7）与闭式理论公式的 4 组实时交互图表。点击右上角「中文 / English」可一键实时切换全页语言。</p>
      <div class="ifig-card" data-ifig="overview" data-lang="zh"></div>
      <div class="ifig-card" data-ifig="synthesizer" data-lang="zh"></div>
      <div class="ifig-card" data-ifig="bayes" data-lang="zh"></div>
      <div class="ifig-card" data-ifig="valleys" data-lang="zh"></div>
    </div>
  </div>
  <script src="assets/interactive.js"></script>
  <script>
    (function () {
      const btns = document.querySelectorAll("[data-lab-lang]");
      btns.forEach(btn => {
        btn.addEventListener("click", () => {
          const lang = btn.getAttribute("data-lab-lang");
          btns.forEach(b => b.classList.toggle("active", b === btn));
          document.documentElement.lang = lang === "zh" ? "zh-CN" : "en";
          document.getElementById("lab-title").textContent = lang === "zh"
            ? "多频相位测距 · 交互式图表实验室"
            : "Multi-Frequency Phase Ranging · Interactive Figure Lab";
          document.getElementById("lab-sub").innerHTML = lang === "zh"
            ? "基于 <code>computation/data/</code> 归档试验（随机种子 7）与闭式理论公式的 4 组实时交互图表。点击右上角「中文 / English」可一键实时切换全页语言。"
            : "Four live interactive figures driven by the archived trials in <code>computation/data/</code> (seed 7) and closed-form theory. Click 中文 / English above to switch language instantly.";
          document.querySelectorAll("[data-ifig]").forEach(el => el.setAttribute("data-lang", lang));
          // Re-dispatch DOMContentLoaded or re-boot interactive figures
          window.dispatchEvent(new Event("DOMContentLoaded"));
        });
      });
    })();
  </script>
</body>
</html>
"""


def write_interactive_data() -> None:
    DATA_DST.mkdir(parents=True, exist_ok=True)
    full = json.loads((DATA_DIR / "results.json").read_text(encoding="utf-8"))
    bval = json.loads((DATA_DIR / "bound_valleys.json").read_text(encoding="utf-8"))
    open_q = json.loads((DATA_DIR / "open_questions.json").read_text(encoding="utf-8"))

    sets_slim = {}
    for k in ("equal", "matched", "geometric"):
        s = full["sets"][k]
        sets_slim[k] = {
            "f_hz": s["f_hz"],
            "fbar_hz": s["fbar_hz"],
            "S": s["S"],
            "F": s["F"],
            "R": s["R"],
            "R_db": s["R_db"],
            "sigma_jump": s["sigma_jump"],
            "lambda_m": s["lambda_m"],
            "hmse_factor": s["hmse_factor"],
            "mmse_factor": s["mmse_factor"],
            "sigma": s["sigma"],
            "mse_coh": s["mse_coh"],
            "mse_inc": s["mse_inc"],
            "hmse": s["hmse"],
            "mmse": s["mmse"],
            "slip_rate": s["slip_rate"],
            "floor_L2_12": s["floor_L2_12"],
            "zzb_half": s["zzb_half"],
        }

    payload = {
        "seed": full["seed"],
        "sets": sets_slim,
        "bound_valleys": {
            "L_m": bval["L_m"],
            "sets": {
                k: {
                    "snap_delta_m": bval["sets"][k]["snap_delta_m"],
                    "snap_cume_frac": bval["sets"][k]["snap_cume_frac"],
                }
                for k in ("equal", "matched", "geometric")
            },
        },
        "open_questions": open_q,
    }
    (DATA_DST / "results.json").write_text(json.dumps(payload, separators=(",", ":")), encoding="utf-8")
    (SITE / "explore.html").write_text(EXPLORE, encoding="utf-8")


def main() -> None:
    FIG_DST.mkdir(parents=True, exist_ok=True)
    for ext in ("*.png", "*.svg"):
        for fig in FIG_SRC.glob(ext):
            shutil.copy2(fig, FIG_DST / fig.name)

    write_interactive_data()

    zh_md = (NOTES / "门限问题的进展.md").read_text(encoding="utf-8")
    en_md = (ROOT / "README.md").read_text(encoding="utf-8")
    cut = en_md.find("## License and citation")
    if cut >= 0:
        en_md = en_md[:cut].rstrip() + "\n"
    en_md = re.sub(r"^# .*\n+", "", en_md, count=1)
    en_md = re.sub(r"^\[\*\*Read online\*\*\].*\n+", "", en_md, count=1)
    en_md = re.sub(r"^October 2026\n+", "", en_md, count=1)
    en_md = "# The Threshold of Multi-Frequency Phase Ranging\n\nOctober 2026\n\n" + en_md.lstrip()

    zh_body, zh_toc = stdlib_md_to_html_and_toc(zh_md)
    en_body, en_toc = stdlib_md_to_html_and_toc(en_md)

    zh_body = inject_interactive_figures(rewrite_repo_links(zh_body), "zh")
    en_body = inject_interactive_figures(rewrite_repo_links(en_body), "en")

    (SITE / "zh.html").write_text(
        PAGE.format(
            lang="zh-CN",
            title="多频相位测距的门限",
            toc=zh_toc,
            body=zh_body,
            zh_active="active",
            en_active="",
            nav_home="首页",
            nav_lab="交互实验室",
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
            nav_home="Home",
            nav_lab="Interactive Lab",
        ),
        encoding="utf-8",
    )
    (SITE / "index.html").write_text(HOME, encoding="utf-8")
    print("wrote", SITE / "index.html")
    print("wrote", SITE / "zh.html")
    print("wrote", SITE / "en.html")
    print("wrote", SITE / "explore.html")


if __name__ == "__main__":
    main()
