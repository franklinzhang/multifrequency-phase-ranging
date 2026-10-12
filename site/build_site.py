#!/usr/bin/env python3
"""Build typographic HTML reading pages from the Markdown notes."""

from __future__ import annotations

import json
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
DATA_SRC = ROOT / "computation" / "data" / "results.json"
DATA_DST = SITE / "data"
REPO = "https://github.com/franklinzhang/multifrequency-phase-ranging"
REPO_BLOB = f"{REPO}/blob/main"
BINDER_BASE = (
    "https://mybinder.org/v2/gh/franklinzhang/multifrequency-phase-ranging/v0.1.0"
)
BINDER_EN = f"{BINDER_BASE}?labpath=notes%2Fthreshold.ipynb"
BINDER_ZH = (
    f"{BINDER_BASE}?labpath=notes%2F%E5%A4%9A%E9%A2%91%E7%9B%B8%E4%BD%8D"
    "%E6%B5%8B%E8%B7%9D%E7%9A%84%E9%97%A8%E9%99%90.ipynb"
)

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


def rewrite_repo_links(html: str) -> str:
    """Point repository paths at GitHub so the Pages edition does not 404."""
    replacements = [
        ('href="../README.md"', f'href="{REPO_BLOB}/README.md"'),
        ('href="../LICENSE"', f'href="{REPO_BLOB}/LICENSE"'),
        ('href="../CITATION.cff"', f'href="{REPO_BLOB}/CITATION.cff"'),
        ('href="../computation/README.md"', f'href="{REPO_BLOB}/computation/README.md"'),
        ('href="../computation/data/"', f'href="{REPO_BLOB}/computation/data"'),
        ('href="../computation/"', f'href="{REPO_BLOB}/computation"'),
        ('href="../notes/', f'href="{REPO_BLOB}/notes/'),
        ('href="README.md"', f'href="{REPO_BLOB}/README.md"'),
        ('href="LICENSE"', f'href="{REPO_BLOB}/LICENSE"'),
        ('href="CITATION.cff"', f'href="{REPO_BLOB}/CITATION.cff"'),
        ('href="computation/README.md"', f'href="{REPO_BLOB}/computation/README.md"'),
        ('href="computation/data/"', f'href="{REPO_BLOB}/computation/data"'),
        ('href="computation/"', f'href="{REPO_BLOB}/computation"'),
        ('href="notes/', f'href="{REPO_BLOB}/notes/'),
    ]
    for old, new in replacements:
        html = html.replace(old, new)
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
        <div class="meta"><a href="index.html">Home</a> · <a href="explore.html">Explore</a> · <a href="https://github.com/franklinzhang/multifrequency-phase-ranging">GitHub</a></div>
      </div>
      {body}
      <div class="footnote-bar">
        CC BY 4.0 · Cite via <a href="https://github.com/franklinzhang/multifrequency-phase-ranging/blob/main/CITATION.cff">CITATION.cff</a> · Archived numbers in <code>computation/data/</code>
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
  <title>The Threshold of Multi-Frequency Phase Ranging</title>
  <meta name="description" content="Coherent range estimation from wrapped multi-frequency phases: the first threshold σ_jump = π√S/f̄, the half-window Ziv–Zakai bound, unequal spacing, and engineering rules. In Chinese and English." />
  <meta property="og:title" content="The Threshold of Multi-Frequency Phase Ranging" />
  <meta property="og:description" content="Coherent range estimation from wrapped multi-frequency phases: σ_jump, the half-window Ziv–Zakai bound, unequal spacing, and engineering rules." />
  <meta property="og:type" content="website" />
  <meta property="og:url" content="https://franklinzhang.github.io/multifrequency-phase-ranging/" />
  <link rel="canonical" href="https://franklinzhang.github.io/multifrequency-phase-ranging/" />
  <link rel="preconnect" href="https://fonts.googleapis.com" />
  <link rel="preconnect" href="https://fonts.gstatic.com" crossorigin />
  <link href="https://fonts.googleapis.com/css2?family=IBM+Plex+Sans:wght@400;500;600&family=Noto+Sans+SC:wght@400;500;600&family=Noto+Serif+SC:wght@600&family=Source+Serif+4:opsz,wght@8..60,600&display=swap" rel="stylesheet" />
  <link rel="stylesheet" href="assets/reading.css" />
</head>
<body>
  <div class="home">
    <div class="home-card">
      <p class="edition-label">多频相位测距的门限</p>
      <h1>The Threshold of Multi-Frequency Phase Ranging</h1>
      <p class="lede">A self-contained parameter-estimation account of coherent range from wrapped multi-frequency phase: two Cramér–Rao lines, the phase scale of the jump, the one-wavelength mixture, and what the half-window Ziv–Zakai bound does and does not say.</p>
      <div class="choices">
        <a class="choice primary" href="zh.html">
          <strong>中文</strong>
          <span>观测模型 · 两条克拉美–罗直线 · 门限尺度 · 下界 · 频点与起始频率</span>
        </a>
        <a class="choice primary" href="en.html">
          <strong>English</strong>
          <span>Observation model · two Cramér–Rao lines · threshold scale · bound · tones and the start of the band</span>
        </a>
        <a class="choice" href="explore.html">
          <strong>Equal-grid curves</strong>
          <span>Step through the 18 archived noise levels: trial MSE, slip rate, carrier and bandwidth lines, half-window bound, mixture</span>
        </a>
        <a class="choice" href="{binder_zh}">
          <strong>中文 notebook</strong>
          <span>在 Binder 上运行；图由归档试验重绘，带噪声点滑块</span>
        </a>
        <a class="choice" href="{binder_en}">
          <strong>English notebook</strong>
          <span>Runs on Binder; figures redrawn from the archived trials, with a noise-level slider</span>
        </a>
      </div>
      <div class="footnote-bar">CC BY 4.0 · Cite <a href="https://github.com/franklinzhang/multifrequency-phase-ranging/blob/main/CITATION.cff">CITATION.cff</a> · Source <a href="https://github.com/franklinzhang/multifrequency-phase-ranging">GitHub</a> · Offline copy <a href="https://github.com/franklinzhang/multifrequency-phase-ranging/releases/tag/v0.1.0">v0.1.0</a></div>
    </div>
  </div>
</body>
</html>
""".replace("{binder_zh}", BINDER_ZH).replace("{binder_en}", BINDER_EN)

EXPLORE = """<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1" />
  <title>Archive explorer — threshold notes</title>
  <link rel="preconnect" href="https://fonts.googleapis.com" />
  <link rel="preconnect" href="https://fonts.gstatic.com" crossorigin />
  <link href="https://fonts.googleapis.com/css2?family=IBM+Plex+Sans:wght@400;500;600&family=Source+Serif+4:opsz,wght@8..60,600&display=swap" rel="stylesheet" />
  <link rel="stylesheet" href="assets/reading.css" />
  <script src="https://cdn.jsdelivr.net/npm/chart.js@4.4.6/dist/chart.umd.min.js"></script>
  <style>
    .explore-wrap { max-width: 52rem; margin: 0 auto; padding: 1.5rem 1rem 3rem; }
    .explore-panel {
      background: var(--panel); border: 1px solid var(--line); border-radius: 1rem;
      box-shadow: 0 18px 50px var(--shade); padding: 1.25rem 1.35rem 1.5rem;
    }
    .explore-panel h1 { margin: 0 0 0.4rem; font-size: 1.7rem; }
    .explore-panel .sub { color: var(--muted); font-family: var(--sans); margin: 0 0 1rem; }
    .controls { font-family: var(--sans); display: grid; gap: 0.55rem; margin: 0.75rem 0 1rem; }
    .controls input[type=range] { width: 100%; }
    .readout { font-family: var(--sans); font-size: 0.95rem; color: var(--ink-soft); min-height: 3.2rem; }
    canvas { width: 100% !important; max-height: 420px; }
  </style>
</head>
<body>
  <div class="explore-wrap">
    <div class="topbar" style="margin-bottom:1rem;font-family:var(--sans);font-size:0.92rem;">
      <div class="lang-switch">
        <a href="index.html">Home</a>
        <a href="zh.html">中文</a>
        <a href="en.html">English</a>
      </div>
      <div class="meta">Frozen archive · seed 7 · English labels</div>
    </div>
    <div class="explore-panel">
      <h1>Archive explorer</h1>
      <p class="sub">Equal-spacing coherent MSE, HMSE, MMSE, half-window bound, and mixture prediction from <code>computation/data/results.json</code>. Read-only — nothing is recomputed.</p>
      <div class="controls">
        <label for="idx">Noise ratio point on the archived grid</label>
        <input id="idx" type="range" min="0" max="1" step="1" value="0" />
        <div class="readout" id="readout"></div>
      </div>
      <canvas id="chart" height="220"></canvas>
    </div>
  </div>
  <script>
  (async function () {
    const payload = await fetch('data/results.json').then(r => r.json());
    const eq = payload.sets.equal;
    const ratio = eq.sigma.map((s, i) => s / eq.sigma_jump);
    const mse = eq.mse_coh;
    const hmse = eq.hmse;
    const mmse = eq.mmse;
    const zzb = eq.zzb_half.zzb;
    const slip = eq.slip_rate;
    const lam2 = eq.lambda_m * eq.lambda_m;
    const mix = slip.map((p, i) => (1 - p) * hmse[i] + p * lam2);
    const db = x => 10 * Math.log10(x);
    const idxDefault = ratio.findIndex(r => Math.abs(r - 1) < 1e-9);
    const slider = document.getElementById('idx');
    slider.max = String(ratio.length - 1);
    slider.value = String(idxDefault >= 0 ? idxDefault : 0);
    const readout = document.getElementById('readout');
    const ctx = document.getElementById('chart');
    const chart = new Chart(ctx, {
      type: 'line',
      data: {
        datasets: [
          { label: 'Coherent MSE', data: ratio.map((r,i)=>({x:r,y:mse[i]})), borderColor:'#1f4e79', pointRadius:3, tension:0 },
          { label: 'HMSE', data: ratio.map((r,i)=>({x:r,y:hmse[i]})), borderColor:'#222', borderDash:[6,4], pointRadius:0 },
          { label: 'MMSE', data: ratio.map((r,i)=>({x:r,y:mmse[i]})), borderColor:'#666', borderDash:[2,3], pointRadius:0 },
          { label: 'Half-window bound', data: ratio.map((r,i)=>({x:r,y:zzb[i]})), borderColor:'#2e7d32', pointRadius:0 },
          { label: 'Mixture', data: ratio.map((r,i)=>({x:r,y:mix[i]})), borderColor:'#6a1b9a', borderDash:[1,2], pointRadius:0 },
        ]
      },
      options: {
        responsive: true,
        scales: {
          x: { type: 'logarithmic', title: { display: true, text: 'σ / σ_jump' } },
          y: { type: 'logarithmic', title: { display: true, text: 'MSE / m²' }, min: 1e-9, max: 200 }
        },
        plugins: { legend: { position: 'bottom' } }
      }
    });
    let marker;
    function paint() {
      const i = Number(slider.value);
      const r = ratio[i];
      readout.innerHTML = `Selected <strong>${r.toFixed(3)}</strong> σ/σ_jump · sample <strong>${db(mse[i]).toFixed(1)} dB</strong> · HMSE ${db(hmse[i]).toFixed(1)} · MMSE ${db(mmse[i]).toFixed(1)} · bound ${db(zzb[i]).toFixed(1)} · slip ${slip[i].toFixed(4)} · mixture ${db(mix[i]).toFixed(1)} dB`;
      if (marker) chart.data.datasets = chart.data.datasets.filter(d => d.label !== 'Selected');
      chart.data.datasets.push({
        label: 'Selected', data: [{x:r,y:mse[i]}], showLine:false,
        pointRadius: 7, pointBackgroundColor: '#c45911', borderColor: '#c45911'
      });
      marker = true;
      chart.update('none');
    }
    slider.addEventListener('input', paint);
    paint();
  })();
  </script>
</body>
</html>
"""

def prepare_zh_source(text: str) -> str:
    return text


def prepare_en_source(text: str) -> str:
    # Science body + references only; license/reproduce stay in the repo README.
    cut = text.find("## License and citation")
    if cut >= 0:
        refs = text.find("## References")
        head = text[:cut].rstrip() + "\n\n"
        # Append References only when they sit after License (older layout).
        if refs > cut:
            text = head + text[refs:]
        else:
            text = head
    # Drop the repo title and the one-line Pages link; site chrome covers navigation.
    text = re.sub(r"^# .*\n+", "", text, count=1)
    text = re.sub(r"^\[\*\*Read online\*\*\].*\n+", "", text, count=1)
    text = re.sub(r"^October 2026\n+", "", text, count=1)
    text = "# The Threshold of Multi-Frequency Phase Ranging\n\n" + text.lstrip()
    return text


def write_explore() -> None:
    DATA_DST.mkdir(parents=True, exist_ok=True)
    # Slim copy: keep only fields the explorer needs, still from the archive.
    full = json.loads(DATA_SRC.read_text(encoding="utf-8"))
    eq = full["sets"]["equal"]
    slim = {
        "seed": full["seed"],
        "sets": {
            "equal": {
                "sigma_jump": eq["sigma_jump"],
                "lambda_m": eq["lambda_m"],
                "sigma": eq["sigma"],
                "mse_coh": eq["mse_coh"],
                "hmse": eq["hmse"],
                "mmse": eq["mmse"],
                "slip_rate": eq["slip_rate"],
                "zzb_half": eq["zzb_half"],
            }
        },
    }
    (DATA_DST / "results.json").write_text(json.dumps(slim), encoding="utf-8")
    (SITE / "explore.html").write_text(EXPLORE, encoding="utf-8")


def main() -> None:
    FIG_DST.mkdir(parents=True, exist_ok=True)
    for png in FIG_SRC.glob("*.png"):
        shutil.copy2(png, FIG_DST / png.name)
    write_explore()

    zh_md = prepare_zh_source((NOTES / "门限问题的进展.md").read_text(encoding="utf-8"))
    en_md = prepare_en_source((ROOT / "README.md").read_text(encoding="utf-8"))
    zh_body = rewrite_repo_links(retarget_figures(md_to_html(zh_md)))
    en_body = rewrite_repo_links(retarget_figures(md_to_html(en_md)))

    # Force TOC generation by ensuring toc extension runs; markdown places toc where [TOC] is.
    if "[TOC]" not in zh_md:
        zh_body = rewrite_repo_links(retarget_figures(md_to_html("[TOC]\n\n" + zh_md)))
    if "[TOC]" not in en_md:
        en_body = rewrite_repo_links(retarget_figures(md_to_html("[TOC]\n\n" + en_md)))

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
    print("wrote", SITE / "explore.html")


if __name__ == "__main__":
    main()
