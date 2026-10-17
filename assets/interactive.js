/* Interactive figures & bilingual switching for The Threshold of Multi-Frequency Phase Ranging */
(function () {
  "use strict";

  const C_LIGHT = 299792458.0;
  const PI = Math.PI;
  const L_HALF = 0.5 * C_LIGHT / 10e6; // 14.98962 m

  function erfc(x) {
    // Abramowitz and Stegun 7.1.26 approximation (max error 1.5e-7)
    const z = Math.abs(x);
    const t = 1.0 / (1.0 + 0.3275911 * z);
    const a1 = 0.254829592, a2 = -0.284496736, a3 = 1.421413741, a4 = -1.453152027, a5 = 1.061405429;
    const poly = ((((a5 * t + a4) * t + a3) * t + a2) * t + a1) * t;
    const ans = poly * Math.exp(-z * z);
    return x >= 0 ? ans : 2.0 - ans;
  }

  function qFunc(x) {
    return 0.5 * erfc(x / Math.SQRT2);
  }

  function wrapRad(x) {
    let y = (x + PI) % (2.0 * PI);
    if (y < 0) y += 2.0 * PI;
    return y - PI;
  }

  function toDb(val) {
    return 10.0 * Math.log10(Math.max(val, 1e-30));
  }

  function setupHiDPICanvas(canvas, logicalW, logicalH) {
    const dpr = Math.min(window.devicePixelRatio || 1, 2.5);
    canvas.width = Math.round(logicalW * dpr);
    canvas.height = Math.round(logicalH * dpr);
    const ctx = canvas.getContext("2d");
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    return ctx;
  }

  // Seamless section-preserving language switch between zh.html and en.html
  function initLangSwitchPreserver() {
    const links = document.querySelectorAll(".lang-switch a[href$='.html']");
    if (!links.length) return;
    links.forEach((a) => {
      a.addEventListener("click", function (e) {
        const href = a.getAttribute("href");
        if (!href || a.classList.contains("active")) return;
        const headings = Array.from(document.querySelectorAll("main.article h2, main.article h3"));
        let activeIdx = -1;
        const scrollY = window.scrollY;
        for (let i = 0; i < headings.length; i++) {
          if (headings[i].getBoundingClientRect().top <= 140) {
            activeIdx = i;
          } else {
            break;
          }
        }
        if (activeIdx >= 0) {
          e.preventDefault();
          window.location.href = href + "#sec-idx-" + activeIdx;
        }
      });
    });

    // Restore section index if hash is #sec-idx-N
    const m = window.location.hash.match(/^#sec-idx-(\d+)$/);
    if (m) {
      const idx = parseInt(m[1], 10);
      const headings = Array.from(document.querySelectorAll("main.article h2, main.article h3"));
      if (headings[idx]) {
        setTimeout(() => {
          headings[idx].scrollIntoView({ behavior: "instant", block: "start" });
        }, 40);
      }
    }
  }

  // =========================================================================
  // INTERACTIVE FIGURE 1: Log-Log Master Anatomy & Multi-Estimator Explorer
  // =========================================================================
  function mountOverviewExplorer(container, data, lang) {
    const zh = lang === "zh";
    let curSetKey = "equal";
    let curIdx = 10; // sigma / sigma_jump = 1.00
    const showCurves = {
      coh: true,
      inc: true,
      vor: true,
      bayes: true,
      zzb: true,
      lines: true,
    };

    container.innerHTML = `
      <div class="ifig-header">
        <div>
          <span class="ifig-badge">${zh ? "交互式图 1 · 双对数解剖" : "Interactive Fig 1 · Log-Log Anatomy"}</span>
          <h4 class="ifig-title" style="display:inline;">${
            zh
              ? "四区门限结构、相干/非相干估计、共线 Voronoi 级数与半窗下界"
              : "Four-Regime Threshold Structure, Coherent/Incoherent MLE, Voronoi Series & ZZB"
          }</h4>
        </div>
        <div class="ifig-btn-group" data-role="set-btns">
          <button class="ifig-btn active" data-set="equal">${zh ? "等间距 (Equal)" : "Equal Grid"}</button>
          <button class="ifig-btn" data-set="matched">${zh ? "矩匹配 (Matched)" : "Moment-Matched"}</button>
          <button class="ifig-btn" data-set="geometric">${zh ? "几何间隔 (Geometric)" : "Geometric (1.02)"}</button>
        </div>
      </div>
      <p class="ifig-sub">${
        zh
          ? "拖动下方噪声滑块或直接点击画布上的数据点，实时查看当前 σ/σ_jump 下的相干均方、非相干均方、共线多波长 Voronoi 级数、贝叶斯后验均值与半窗 Ziv–Zakai 下界。"
          : "Drag the noise slider or click directly on the plot to inspect Coherent MSE, Incoherent MSE, Collinear Voronoi Series, Bayesian MMSE, and the Half-Window Ziv–Zakai Bound at each σ/σ_jump."
      }</p>
      <div class="ifig-slider-row">
        <strong>${zh ? "归一化噪声点 σ / σ_jump：" : "Normalized noise σ / σ_jump:"}</strong>
        <input type="range" min="0" max="17" step="1" value="10" data-role="idx-slider" />
        <span data-role="slider-val" style="font-family:var(--mono);font-weight:600;min-width:9.5rem;text-align:right;">1.000 σ_jump</span>
      </div>
      <div class="ifig-readout" data-role="readout"></div>
      <div class="ifig-canvas-wrap">
        <canvas class="ifig-canvas" width="920" height="470"></canvas>
      </div>
      <div class="ifig-checks">
        <label><input type="checkbox" data-curve="coh" checked /> ${zh ? "相干极大似然 (Coherent MLE)" : "Coherent MLE"}</label>
        <label><input type="checkbox" data-curve="inc" checked /> ${zh ? "非相干估计 (Incoherent)" : "Incoherent MLE"}</label>
        <label><input type="checkbox" data-curve="vor" checked /> ${zh ? "共线 Voronoi 级数 (Voronoi Series)" : "Voronoi Series"}</label>
        <label><input type="checkbox" data-curve="bayes" checked /> ${zh ? "同先验贝叶斯 E[d|φ]" : "Bayesian MMSE E[d|φ]"}</label>
        <label><input type="checkbox" data-curve="zzb" checked /> ${zh ? "半窗下界 (Half-Window ZZB)" : "Half-Window ZZB"}</label>
      </div>
    `;

    const canvas = container.querySelector("canvas");
    const slider = container.querySelector('[data-role="idx-slider"]');
    const sliderVal = container.querySelector('[data-role="slider-val"]');
    const readout = container.querySelector('[data-role="readout"]');
    const setBtns = container.querySelectorAll("[data-set]");
    const curveChecks = container.querySelectorAll("[data-curve]");

    function computeVoronoiSeries(sigmas, sj, lam, hmse) {
      return sigmas.map((s, idx) => {
        const p1 = 2.0 * qFunc(sj / s);
        let sumV = 0.0;
        for (let k = 1; k < 40; k++) {
          sumV += (2 * k - 1) * qFunc(((2 * k - 1) * sj) / s);
        }
        return (1.0 - p1) * hmse[idx] + 2.0 * lam * lam * sumV;
      });
    }

    function render() {
      const sData = data.sets[curSetKey];
      const sj = sData.sigma_jump;
      const lam = sData.lambda_m;
      const sigmas = sData.sigma;
      const ratios = sigmas.map((s) => s / sj);
      const vor = computeVoronoiSeries(sigmas, sj, lam, sData.hmse);
      const bayesDbList = data.open_questions.q8_2_bayesian_mmse.curve.map((r) => r.expected_bayes_collinear_db);

      const rCur = ratios[curIdx];
      const sigCur = sigmas[curIdx];
      const snrDb = 10.0 * Math.log10(1.0 / (2.0 * sigCur * sigCur));

      sliderVal.textContent = `${rCur.toFixed(2)} σ_jump (${sigCur.toFixed(4)} rad)`;

      let regimeText = "";
      if (rCur <= 0.22) {
        regimeText = zh ? "I. 载频主瓣区（整周无误）" : "I. Carrier Main-Lobe Regime";
      } else if (rCur <= 1.0) {
        regimeText = curSetKey === "geometric"
          ? (zh ? "II. 远谷起跳前夜 / 跃升区" : "II. Pre-Jump / Far-Valley Onset")
          : (zh ? "II. 第一道门限起跳（1λ 滑移主导）" : "II. First Threshold Jump (1λ Slips)");
      } else if (rCur <= 6.0) {
        regimeText = curSetKey === "geometric"
          ? (zh ? "III. 6.21 m 远谷主导（无带宽平台）" : "III. 6.21 m Far-Valley Jump (No Plateau)")
          : (zh ? "III. 带宽直线平台（共线 kλ Voronoi）" : "III. Bandwidth Plateau (Collinear kλ)");
      } else {
        regimeText = zh ? "IV. 第二道门限（趋向搜索窗地板）" : "IV. Second Threshold (Toward Floor)";
      }

      const bayesStr = curSetKey === "equal" ? `${bayesDbList[curIdx].toFixed(2)} dB` : "—";

      readout.innerHTML = `
        <div class="ifig-metric"><span>${zh ? "当前工作区" : "Regime"}</span><strong style="font-family:var(--sans);font-size:0.84rem;color:var(--accent-ink);">${regimeText}</strong></div>
        <div class="ifig-metric"><span>${zh ? "单频点复 SNR" : "Per-Tone SNR"}</span><strong>${snrDb.toFixed(1)} dB</strong></div>
        <div class="ifig-metric"><span>${zh ? "相干试验 MSE" : "Coherent MSE"}</span><strong style="color:#1f4e79;">${toDb(sData.mse_coh[curIdx]).toFixed(2)} dB</strong></div>
        <div class="ifig-metric"><span>${zh ? "非相干试验 MSE" : "Incoherent MSE"}</span><strong style="color:#d96b27;">${toDb(sData.mse_inc[curIdx]).toFixed(2)} dB</strong></div>
        <div class="ifig-metric"><span>${zh ? "Voronoi / 贝叶斯" : "Voronoi / Bayes"}</span><strong style="color:#6a1b9a;">${toDb(vor[curIdx]).toFixed(2)} / ${bayesStr}</strong></div>
        <div class="ifig-metric"><span>${zh ? "半窗下界 ZZB" : "Half-Window ZZB"}</span><strong style="color:#2e7d32;">${toDb(sData.zzb_half.zzb[curIdx]).toFixed(2)} dB</strong></div>
        <div class="ifig-metric"><span>${zh ? "HMSE / MMSE" : "HMSE / MMSE"}</span><strong>${toDb(sData.hmse[curIdx]).toFixed(1)} / ${toDb(sData.mmse[curIdx]).toFixed(1)} dB</strong></div>
        <div class="ifig-metric"><span>${zh ? "样本滑移率" : "Slip Rate"}</span><strong>${(sData.slip_rate[curIdx] * 100).toFixed(2)}%</strong></div>
      `;

      const W = 920, H = 470;
      const ctx = setupHiDPICanvas(canvas, W, H);
      const L = 74, R = 28, T = 36, B = 52;
      const PW = W - L - R, PH = H - T - B;

      const rMin = curSetKey === "geometric" ? 0.016 : 0.08;
      const rMax = curSetKey === "geometric" ? 8.0 : 35.0;
      const lxMin = Math.log10(rMin), lxMax = Math.log10(rMax);
      const lyMin = -90.0, lyMax = 25.0; // in dB

      const px = (r) => L + ((Math.log10(r) - lxMin) / (lxMax - lxMin)) * PW;
      const py = (db) => T + PH - ((Math.max(Math.min(db, lyMax), lyMin) - lyMin) / (lyMax - lyMin)) * PH;

      ctx.fillStyle = "#ffffff";
      ctx.fillRect(0, 0, W, H);

      // Regime shading for equal/matched
      if (curSetKey !== "geometric") {
        const xLv = px(0.204), xJp = px(1.0), xSc = px(6.0);
        ctx.fillStyle = "rgba(11, 110, 99, 0.045)";
        ctx.fillRect(L, T, xLv - L, PH);
        ctx.fillStyle = "rgba(196, 89, 17, 0.055)";
        ctx.fillRect(xLv, T, xJp - xLv, PH);
        ctx.fillStyle = "rgba(31, 78, 121, 0.04)";
        ctx.fillRect(xJp, T, xSc - xJp, PH);
      }

      // Grid lines
      ctx.strokeStyle = "#e4e8ec";
      ctx.lineWidth = 1;
      for (let db = -80; db <= 20; db += 20) {
        const yy = py(db);
        ctx.beginPath();
        ctx.moveTo(L, yy);
        ctx.lineTo(L + PW, yy);
        ctx.stroke();
        ctx.fillStyle = "#5a6772";
        ctx.font = "12px -apple-system, BlinkMacSystemFont, sans-serif";
        ctx.textAlign = "right";
        ctx.fillText(`${db} dB`, L - 8, yy + 4);
      }

      const xTicks = curSetKey === "geometric" ? [0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5] : [0.1, 0.2, 0.5, 1, 2, 5, 10, 20];
      xTicks.forEach((xt) => {
        const xx = px(xt);
        ctx.beginPath();
        ctx.moveTo(xx, T);
        ctx.lineTo(xx, T + PH);
        ctx.stroke();
        ctx.fillStyle = "#5a6772";
        ctx.font = "12px -apple-system, BlinkMacSystemFont, sans-serif";
        ctx.textAlign = "center";
        ctx.fillText(String(xt), xx, T + PH + 18);
      });

      // Threshold vertical line at sigma = sigma_jump (ratio = 1.0)
      const xOne = px(1.0);
      ctx.save();
      ctx.strokeStyle = "#c45911";
      ctx.lineWidth = 1.8;
      ctx.setLineDash([6, 4]);
      ctx.beginPath();
      ctx.moveTo(xOne, T);
      ctx.lineTo(xOne, T + PH);
      ctx.stroke();
      ctx.restore();

      // Search-window floor
      const yFl = py(toDb(sData.floor_L2_12));
      ctx.save();
      ctx.strokeStyle = "#888";
      ctx.lineWidth = 1.5;
      ctx.setLineDash([3, 4]);
      ctx.beginPath();
      ctx.moveTo(L, yFl);
      ctx.lineTo(L + PW, yFl);
      ctx.stroke();
      ctx.restore();

      function drawLine(arrDb, color, lw, dash) {
        ctx.save();
        ctx.strokeStyle = color;
        ctx.lineWidth = lw;
        if (dash) ctx.setLineDash(dash);
        ctx.beginPath();
        ratios.forEach((r, i) => {
          const x = px(r), y = py(arrDb[i]);
          if (i === 0) ctx.moveTo(x, y);
          else ctx.lineTo(x, y);
        });
        ctx.stroke();
        ctx.restore();
      }

      if (showCurves.lines) {
        drawLine(sData.hmse.map(toDb), "#1c232b", 2.0, [7, 5]);
        drawLine(sData.mmse.map(toDb), "#606a74", 2.0, [9, 4, 2, 4]);
      }
      if (showCurves.zzb) {
        drawLine(sData.zzb_half.zzb.map(toDb), "#2e7d32", 2.6, null);
      }
      if (showCurves.vor && curSetKey !== "geometric") {
        drawLine(vor.slice(0, 15).map(toDb), "#6a1b9a", 2.2, [4, 4]);
      }
      if (showCurves.bayes && curSetKey === "equal") {
        drawLine(bayesDbList.slice(0, 15), "#0b6e63", 2.2, [2, 3]);
      }
      if (showCurves.inc) {
        const incDb = sData.mse_inc.map(toDb);
        drawLine(incDb, "#d96b27", 2.0, null);
        ratios.forEach((r, i) => {
          const x = px(r), y = py(incDb[i]);
          ctx.fillStyle = "#fff";
          ctx.strokeStyle = "#d96b27";
          ctx.lineWidth = 1.8;
          ctx.fillRect(x - 4, y - 4, 8, 8);
          ctx.strokeRect(x - 4, y - 4, 8, 8);
        });
      }
      if (showCurves.coh) {
        const cohDb = sData.mse_coh.map(toDb);
        drawLine(cohDb, "#1f4e79", 2.8, null);
        ratios.forEach((r, i) => {
          const x = px(r), y = py(cohDb[i]);
          ctx.beginPath();
          ctx.arc(x, y, i === curIdx ? 6.5 : 4.2, 0, 2 * PI);
          ctx.fillStyle = i === curIdx ? "#c45911" : "#1f4e79";
          ctx.fill();
          ctx.strokeStyle = "#fff";
          ctx.lineWidth = 1.5;
          ctx.stroke();
        });
      }

      // Highlight vertical cursor at curIdx
      const xCur = px(rCur);
      ctx.save();
      ctx.strokeStyle = "rgba(196, 89, 17, 0.55)";
      ctx.lineWidth = 1.5;
      ctx.beginPath();
      ctx.moveTo(xCur, T);
      ctx.lineTo(xCur, T + PH);
      ctx.stroke();
      ctx.restore();

      // On-canvas mathematical annotations
      ctx.fillStyle = "#182026";
      ctx.font = "600 12px -apple-system, BlinkMacSystemFont, sans-serif";
      ctx.textAlign = "left";
      const rDb = sData.R_db.toFixed(2);
      ctx.fillText(
        zh
          ? `平行落差 R = F/S = ${rDb} dB | 门限尺度 σ_jump = ${sj.toFixed(4)} rad`
          : `Parallel Gap R = F/S = ${rDb} dB | Threshold σ_jump = ${sj.toFixed(4)} rad`,
        L + 12,
        T + 20
      );

      // Frame & Axis labels
      ctx.strokeStyle = "#182026";
      ctx.lineWidth = 1.5;
      ctx.strokeRect(L, T, PW, PH);

      ctx.fillStyle = "#182026";
      ctx.font = "600 13px -apple-system, BlinkMacSystemFont, sans-serif";
      ctx.textAlign = "center";
      ctx.fillText(zh ? "归一化相位噪声标准差  σ / σ_jump" : "Normalized Phase Noise  σ / σ_jump", L + 0.5 * PW, H - 12);

      ctx.save();
      ctx.translate(18, T + 0.5 * PH);
      ctx.rotate(-PI / 2);
      ctx.fillText(zh ? "距离均方误差  MSE (dB rel. 1 m²)" : "Mean-Square Error  MSE (dB rel. 1 m²)", 0, 0);
      ctx.restore();
    }

    slider.addEventListener("input", () => {
      curIdx = parseInt(slider.value, 10);
      render();
    });

    setBtns.forEach((btn) => {
      btn.addEventListener("click", () => {
        setBtns.forEach((b) => b.classList.remove("active"));
        btn.classList.add("active");
        curSetKey = btn.getAttribute("data-set");
        render();
      });
    });

    curveChecks.forEach((chk) => {
      chk.addEventListener("change", () => {
        showCurves[chk.getAttribute("data-curve")] = chk.checked;
        render();
      });
    });

    canvas.addEventListener("click", (e) => {
      const rect = canvas.getBoundingClientRect();
      const relX = ((e.clientX - rect.left) / rect.width) * 920;
      const sData = data.sets[curSetKey];
      const sj = sData.sigma_jump;
      const ratios = sData.sigma.map((s) => s / sj);
      const rMin = curSetKey === "geometric" ? 0.016 : 0.08;
      const rMax = curSetKey === "geometric" ? 8.0 : 35.0;
      const lxMin = Math.log10(rMin), lxMax = Math.log10(rMax);
      const px = (r) => 74 + ((Math.log10(r) - lxMin) / (lxMax - lxMin)) * (920 - 74 - 28);
      let bestI = 0, bestDist = 1e9;
      ratios.forEach((r, i) => {
        const d = Math.abs(px(r) - relX);
        if (d < bestDist) {
          bestDist = d;
          bestI = i;
        }
      });
      curIdx = bestI;
      slider.value = String(bestI);
      render();
    });

    render();
  }

  // =========================================================================
  // INTERACTIVE FIGURE 2: Cost Surface J(delta) & Live 8-Tone Synthesizer Lab
  // =========================================================================
  function mountSynthesizerLab(container, data, lang) {
    const zh = lang === "zh";
    const presets = {
      equal: {
        label: zh ? "等间距基准 (Equal)" : "Equal Grid",
        freqs: [2400, 2410, 2420, 2430, 2440, 2450, 2460, 2470],
        maxDelta: 2.0,
        yCap: 6.0,
      },
      matched: {
        label: zh ? "矩匹配 (Matched)" : "Moment-Matched",
        freqs: [2390, 2418, 2425, 2432, 2443, 2449, 2454, 2469],
        maxDelta: 2.0,
        yCap: 6.0,
      },
      geometric: {
        label: zh ? "几何间隔 (6.21m 远谷)" : "Geometric (6.21m Valley)",
        freqs: [2271.95, 2317.39, 2363.73, 2411.01, 2459.23, 2508.41, 2558.58, 2609.75],
        maxDelta: 7.0,
        yCap: 6.0,
      },
      bimodal: {
        label: zh ? "极端双簇 (4.49m 拍频陷阱)" : "Extreme Bimodal (4.49m Trap)",
        freqs: [2400, 2401, 2402, 2403, 2467, 2468, 2469, 2470],
        maxDelta: 5.2,
        yCap: 2.0,
      },
      optimal: {
        label: zh ? "同带宽微错开最优 (+49.1%)" : "Edge-Staggered Optimum (+49.1%)",
        freqs: [2400, 2401, 2402, 2405, 2465, 2468, 2469, 2470],
        maxDelta: 5.2,
        yCap: 2.0,
      },
    };

    let curPreset = "optimal";
    let customInnerOffset = 30; // corresponds to 2435 +/- 30 -> 2405 / 2465 MHz

    container.innerHTML = `
      <div class="ifig-header">
        <div>
          <span class="ifig-badge">${zh ? "交互式图 2 · 排频实验室" : "Interactive Fig 2 · Tone Placement Lab"}</span>
          <h4 class="ifig-title" style="display:inline;">${
            zh
              ? "包裹相位代价曲面 J(δ)、双簇拍频陷阱与最坏谷离开点优化"
              : "Wrapped Phase Cost J(δ), Bimodal Beat Trap & Worst-Valley Optimization"
          }</h4>
        </div>
        <div class="ifig-btn-group" data-role="preset-btns">
          ${Object.keys(presets)
            .map(
              (k) =>
                `<button class="ifig-btn ${k === curPreset ? "active" : ""}" data-preset="${k}">${presets[k].label}</button>`
            )
            .join("")}
        </div>
      </div>
      <p class="ifig-sub">${
        zh
          ? "切换预设或拖动下方「最内侧频偏 u₁」滑块（固定外侧 ±35, ±34, ±33 MHz），实时观察 4.49 m 双簇拍频远谷如何随内侧频点微移而消长，以及最早离开点 σ_leave 由谁决定。"
          : "Switch presets or drag the innermost offset slider u₁ (with outer tones at ±35, ±34, ±33 MHz) to watch the 4.49 m inter-cluster beat valley appear or vanish in real time."
      }</p>
      <div class="ifig-slider-row">
        <strong>${zh ? "微调内侧对称频偏 ±u₁ (2435 ± u₁ MHz)：" : "Adjust inner symmetric offset ±u₁ (2435 ± u₁ MHz):"}</strong>
        <input type="range" min="5" max="32" step="1" value="30" data-role="u1-slider" />
        <span data-role="u1-val" style="font-family:var(--mono);font-weight:600;min-width:11rem;text-align:right;">±30 MHz (2405, 2465)</span>
      </div>
      <div class="ifig-readout" data-role="synth-readout"></div>
      <div class="ifig-canvas-wrap">
        <canvas class="ifig-canvas" width="920" height="400"></canvas>
      </div>
    `;

    const canvas = container.querySelector("canvas");
    const readout = container.querySelector('[data-role="synth-readout"]');
    const u1Slider = container.querySelector('[data-role="u1-slider"]');
    const u1Val = container.querySelector('[data-role="u1-val"]');
    const presetBtns = container.querySelectorAll("[data-preset]");

    function solveSigmaLeave(jVal, delta, hmseFactor) {
      if (jVal <= 1e-12 || delta <= 1e-6) return 0.0;
      const norm = Math.sqrt(jVal);
      let lo = 1e-4, hi = 1.5;
      for (let it = 0; it < 40; it++) {
        const mid = 0.5 * (lo + hi);
        const lhs = 2.0 * qFunc(norm / (2.0 * mid)) * delta * delta;
        const rhs = hmseFactor * mid * mid;
        if (lhs >= rhs) hi = mid;
        else lo = mid;
      }
      return 0.5 * (lo + hi);
    }

    function evaluateFreqs(freqsMhz, maxDelta, yCap) {
      const N = freqsMhz.length;
      const fbar = freqsMhz.reduce((a, b) => a + b, 0) / N;
      const S = freqsMhz.reduce((acc, f) => acc + (f - fbar) * (f - fbar), 0);
      const F = freqsMhz.reduce((acc, f) => acc + f * f, 0);
      const R = F / S;
      const sj = (PI * Math.sqrt(S)) / fbar;
      const lam = C_LIGHT / (fbar * 1e6);
      const hmseFactor = (C_LIGHT * C_LIGHT) / (4.0 * PI * PI * F * 1e12);

      // Evaluate J(lambda)
      let jLam = 0.0;
      for (let k = 0; k < N; k++) {
        const r = wrapRad(((2.0 * PI) / C_LIGHT) * lam * freqsMhz[k] * 1e6);
        jLam += r * r;
      }
      const leaveLam = solveSigmaLeave(jLam, lam, hmseFactor);

      // Scan valleys up to 7.0 m to find earliest leaving point & deepest far valley
      const nScan = 3500;
      const dMaxScan = 7.0;
      const step = dMaxScan / nScan;
      let earliestSig = leaveLam;
      let earliestDelta = lam;
      let farMinJ = 1e9;
      let farMinDelta = 0.0;

      let prev2 = 0.0, prev1 = 0.0;
      for (let i = 0; i <= nScan; i++) {
        const d = i * step;
        let s = 0.0;
        const sc = ((2.0 * PI) / C_LIGHT) * d * 1e6;
        for (let k = 0; k < N; k++) {
          const r = wrapRad(sc * freqsMhz[k]);
          s += r * r;
        }
        if (i >= 2) {
          const dMid = (i - 1) * step;
          if (dMid >= 0.6 * lam && prev1 <= prev2 && prev1 < s) {
            // Refine minimum locally
            let bestD = dMid, bestJ = prev1;
            for (let m = -8; m <= 8; m++) {
              const dm = dMid + m * (step / 8.0);
              let sm = 0.0;
              const scm = ((2.0 * PI) / C_LIGHT) * dm * 1e6;
              for (let k = 0; k < N; k++) {
                const r = wrapRad(scm * freqsMhz[k]);
                sm += r * r;
              }
              if (sm < bestJ) {
                bestJ = sm;
                bestD = dm;
              }
            }
            const sigV = solveSigmaLeave(bestJ, bestD, hmseFactor);
            if (sigV < earliestSig) {
              earliestSig = sigV;
              earliestDelta = bestD;
            }
            if (bestD >= 1.8 * lam && bestJ < farMinJ) {
              farMinJ = bestJ;
              farMinDelta = bestD;
            }
          }
        }
        prev2 = prev1;
        prev1 = s;
      }

      // Build plot curve on [0, maxDelta]
      const nPlot = 1400;
      const xs = [], ys = [];
      for (let i = 0; i <= nPlot; i++) {
        const d = (i / nPlot) * maxDelta;
        let s = 0.0;
        const sc = ((2.0 * PI) / C_LIGHT) * d * 1e6;
        for (let k = 0; k < N; k++) {
          const r = wrapRad(sc * freqsMhz[k]);
          s += r * r;
        }
        xs.push(d);
        ys.push(s);
      }

      return {
        freqsMhz, fbar, S, R, sj, lam, jLam, leaveLam,
        earliestSig, earliestDelta, farMinJ, farMinDelta,
        xs, ys, maxDelta, yCap
      };
    }

    function render(freqsMhz, maxDelta, yCap) {
      const m = evaluateFreqs(freqsMhz, maxDelta, yCap);
      const gainPct = ((m.earliestSig - 0.01703) / 0.01703) * 100.0;
      const gainStr = `${gainPct >= 0 ? "+" : ""}${gainPct.toFixed(1)}%`;

      readout.innerHTML = `
        <div class="ifig-metric"><span>${zh ? "频点集合 (MHz)" : "Tones (MHz)"}</span><strong style="font-size:0.78rem;">${m.freqsMhz.map((f) => Math.round(f)).join(", ")}</strong></div>
        <div class="ifig-metric"><span>${zh ? "二阶矩 S" : "Moment S"}</span><strong>${Math.round(m.S)} MHz²</strong></div>
        <div class="ifig-metric"><span>${zh ? "两线落差 R" : "Parallel Gap R"}</span><strong>${(10 * Math.log10(m.R)).toFixed(2)} dB</strong></div>
        <div class="ifig-metric"><span>${zh ? "门限尺度 σ_jump" : "Scale σ_jump"}</span><strong>${m.sj.toFixed(4)} rad</strong></div>
        <div class="ifig-metric"><span>${zh ? "最早离开点 σ_leave" : "Earliest σ_leave"}</span><strong style="color:var(--accent-ink);">${m.earliestSig.toFixed(4)} rad (${gainStr})</strong></div>
        <div class="ifig-metric"><span>${zh ? "决定离开点的谷" : "Deciding Valley"}</span><strong style="color:${m.earliestDelta > 0.2 ? "#c45911" : "#0b6e63"};">${m.earliestDelta.toFixed(3)} m</strong></div>
        <div class="ifig-metric"><span>${zh ? "最深远谷 J / J(λ)" : "Deepest Far Valley"}</span><strong>${(m.farMinJ / m.jLam).toFixed(2)}× J(λ) (@ ${m.farMinDelta.toFixed(2)}m)</strong></div>
      `;

      const W = 920, H = 400;
      const ctx = setupHiDPICanvas(canvas, W, H);
      const L = 68, R = 26, T = 32, B = 48;
      const PW = W - L - R, PH = H - T - B;

      const px = (d) => L + (d / m.maxDelta) * PW;
      const py = (j) => T + PH - (Math.min(j, m.yCap) / m.yCap) * PH;

      ctx.fillStyle = "#ffffff";
      ctx.fillRect(0, 0, W, H);

      // Grid
      ctx.strokeStyle = "#e6eaee";
      ctx.lineWidth = 1;
      for (let k = 0; k <= 4; k++) {
        const jv = (k / 4) * m.yCap;
        const yy = py(jv);
        ctx.beginPath();
        ctx.moveTo(L, yy);
        ctx.lineTo(L + PW, yy);
        ctx.stroke();
        ctx.fillStyle = "#5a6772";
        ctx.font = "12px -apple-system, BlinkMacSystemFont, sans-serif";
        ctx.textAlign = "right";
        ctx.fillText(jv.toFixed(1), L - 8, yy + 4);
      }

      const nXTicks = 6;
      for (let k = 0; k <= nXTicks; k++) {
        const dv = (k / nXTicks) * m.maxDelta;
        const xx = px(dv);
        ctx.beginPath();
        ctx.moveTo(xx, T);
        ctx.lineTo(xx, T + PH);
        ctx.stroke();
        ctx.fillStyle = "#5a6772";
        ctx.font = "12px -apple-system, BlinkMacSystemFont, sans-serif";
        ctx.textAlign = "center";
        ctx.fillText(`${dv.toFixed(1)} m`, xx, T + PH + 18);
      }

      // Horizontal reference line at J(lambda)
      const yJlam = py(m.jLam);
      ctx.save();
      ctx.strokeStyle = "#c45911";
      ctx.lineWidth = 1.6;
      ctx.setLineDash([5, 4]);
      ctx.beginPath();
      ctx.moveTo(L, yJlam);
      ctx.lineTo(L + PW, yJlam);
      ctx.stroke();
      ctx.restore();
      ctx.fillStyle = "#c45911";
      ctx.font = "600 12px -apple-system, BlinkMacSystemFont, sans-serif";
      ctx.textAlign = "left";
      ctx.fillText(`J(λ) = ${m.jLam.toFixed(4)} rad²`, L + 10, Math.max(T + 16, yJlam - 6));

      // Vertical marker at lambda
      const xLam = px(m.lam);
      ctx.save();
      ctx.strokeStyle = "#0b6e63";
      ctx.lineWidth = 1.8;
      ctx.setLineDash([4, 4]);
      ctx.beginPath();
      ctx.moveTo(xLam, T);
      ctx.lineTo(xLam, T + PH);
      ctx.stroke();
      ctx.restore();

      // Plot J(delta) clipped at yCap
      ctx.save();
      ctx.strokeStyle = "#1f4e79";
      ctx.lineWidth = 2.2;
      ctx.beginPath();
      let penDown = false;
      for (let i = 0; i < m.xs.length; i++) {
        if (m.ys[i] <= m.yCap * 1.02) {
          const x = px(m.xs[i]), y = py(m.ys[i]);
          if (!penDown) {
            ctx.moveTo(x, y);
            penDown = true;
          } else {
            ctx.lineTo(x, y);
          }
        } else {
          penDown = false;
        }
      }
      ctx.stroke();
      ctx.restore();

      // Highlight deciding valley
      if (m.earliestDelta <= m.maxDelta) {
        const xDec = px(m.earliestDelta);
        const yDec = py(m.earliestDelta > 0.2 ? m.farMinJ : m.jLam);
        ctx.beginPath();
        ctx.arc(xDec, yDec, 6.5, 0, 2 * PI);
        ctx.fillStyle = m.earliestDelta > 0.2 ? "#c45911" : "#0b6e63";
        ctx.fill();
        ctx.strokeStyle = "#fff";
        ctx.lineWidth = 1.8;
        ctx.stroke();
      }

      // Frame & labels
      ctx.strokeStyle = "#182026";
      ctx.lineWidth = 1.5;
      ctx.strokeRect(L, T, PW, PH);

      ctx.fillStyle = "#182026";
      ctx.font = "600 13px -apple-system, BlinkMacSystemFont, sans-serif";
      ctx.textAlign = "center";
      ctx.fillText(zh ? "距离偏移  δ (m)" : "Range Offset  δ (m)", L + 0.5 * PW, H - 10);

      ctx.save();
      ctx.translate(18, T + 0.5 * PH);
      ctx.rotate(-PI / 2);
      ctx.fillText(zh ? "包裹相位代价  J(δ) (rad²)" : "Wrapped Phase Cost  J(δ) (rad²)", 0, 0);
      ctx.restore();
    }

    presetBtns.forEach((btn) => {
      btn.addEventListener("click", () => {
        presetBtns.forEach((b) => b.classList.remove("active"));
        btn.classList.add("active");
        curPreset = btn.getAttribute("data-preset");
        const p = presets[curPreset];
        if (curPreset === "bimodal") {
          u1Slider.value = "32";
          u1Val.textContent = "±32 MHz (2403, 2467)";
        } else if (curPreset === "optimal") {
          u1Slider.value = "30";
          u1Val.textContent = "±30 MHz (2405, 2465)";
        }
        render(p.freqs, p.maxDelta, p.yCap);
      });
    });

    u1Slider.addEventListener("input", () => {
      presetBtns.forEach((b) => b.classList.remove("active"));
      const u1 = parseInt(u1Slider.value, 10);
      customInnerOffset = u1;
      const fLow = 2435 - u1, fHigh = 2435 + u1;
      u1Val.textContent = `±${u1} MHz (${fLow}, ${fHigh})`;
      const customFreqs = [2400, 2401, 2402, fLow, fHigh, 2468, 2469, 2470].sort((a, b) => a - b);
      render(customFreqs, 5.2, 2.0);
    });

    const p0 = presets[curPreset];
    render(p0.freqs, p0.maxDelta, p0.yCap);
  }

  // =========================================================================
  // INTERACTIVE FIGURE 3: Hard MLE Staircase vs Smooth Bayesian Sigmoid (8.2)
  // =========================================================================
  function mountBayesExplorer(container, data, lang) {
    const zh = lang === "zh";
    let curRatio = 0.85; // sigma / sigma_jump

    container.innerHTML = `
      <div class="ifig-header">
        <div>
          <span class="ifig-badge">${zh ? "交互式图 3 · 贝叶斯软阶梯 vs 极大似然硬取整" : "Interactive Fig 3 · Bayes Sigmoid vs Hard MLE Staircase"}</span>
          <h4 class="ifig-title" style="display:inline;">${
            zh
              ? "分段抛物线重心加权如何消除门限过冲，以及 4.82 dB 差距的三方归属"
              : "How Boltzmann Averaging Eliminates Threshold Overshoot & the 4.82 dB Attribution"
          }</h4>
        </div>
      </div>
      <p class="ifig-sub">${
        zh
          ? "沿共线 kλ 谷方向，相干极大似然对投影噪声 Z 执行硬取整阶梯 round(Z·σ/(2σ_jump))，而贝叶斯后验均值 E[d|φ] 执行平滑 Boltzmann 加权，在过渡带省下 1.15–1.30 dB 并恰好落在带宽直线 MMSE 上。"
          : "Along the collinear kλ ladder, Coherent MLE applies a hard-rounding staircase round(Z·σ/(2σ_jump)), whereas the Bayesian posterior mean E[d|φ] applies a smooth Boltzmann sigmoid, saving 1.15–1.30 dB in the transition strip."
      }</p>
      <div class="ifig-slider-row">
        <strong>${zh ? "调节噪声比 σ / σ_jump：" : "Noise ratio σ / σ_jump:"}</strong>
        <input type="range" min="0.20" max="2.50" step="0.05" value="0.85" data-role="bayes-slider" />
        <span data-role="bayes-val" style="font-family:var(--mono);font-weight:600;min-width:8rem;text-align:right;">0.85 σ_jump</span>
      </div>
      <div class="ifig-readout" data-role="bayes-readout"></div>
      <div class="ifig-canvas-wrap">
        <canvas class="ifig-canvas" width="920" height="360"></canvas>
      </div>
    `;

    const canvas = container.querySelector("canvas");
    const slider = container.querySelector('[data-role="bayes-slider"]');
    const valSpan = container.querySelector('[data-role="bayes-val"]');
    const readout = container.querySelector('[data-role="bayes-readout"]');

    function render() {
      valSpan.textContent = `${curRatio.toFixed(2)} σ_jump`;
      const alpha = 1.0 / curRatio; // sigma_jump / sigma

      // Compute expected MSE_ML and MSE_Bayes via 1D Gaussian quadrature on Z in [0, 7]
      const lam = data.sets.equal.lambda_m;
      const sj = data.sets.equal.sigma_jump;
      const sig = curRatio * sj;
      const hmse = data.sets.equal.hmse_factor ? data.sets.equal.hmse_factor * sig * sig : (data.sets.equal.hmse[10] * curRatio * curRatio);
      const mmse = data.sets.equal.mmse[10] * curRatio * curRatio;

      const nZ = 600;
      const zMax = 6.5;
      const dz = zMax / nZ;
      let sumMl = 0.0, sumBy = 0.0;
      for (let i = 0; i < nZ; i++) {
        const z = (i + 0.5) * dz;
        const pdf = Math.exp(-0.5 * z * z) / Math.sqrt(2.0 * PI);
        const kMl = Math.floor(z / (2.0 * alpha) + 0.5);
        let maxEx = -1e30;
        const exps = [];
        for (let k = -12; k <= 12; k++) {
          const diff = z - 2.0 * alpha * k;
          const ex = -0.5 * diff * diff;
          if (ex > maxEx) maxEx = ex;
          exps.push({ k, ex });
        }
        let num = 0.0, den = 0.0;
        exps.forEach((item) => {
          const w = Math.exp(item.ex - maxEx);
          den += w;
          num += item.k * w;
        });
        const kBy = num / den;
        sumMl += 2.0 * pdf * (kMl * lam) * (kMl * lam) * dz;
        sumBy += 2.0 * pdf * (kBy * lam) * (kBy * lam) * dz;
      }
      const mseMlDb = toDb(hmse + sumMl);
      const mseByDb = toDb(hmse + sumBy);
      const mmseDb = toDb(mmse);

      readout.innerHTML = `
        <div class="ifig-metric"><span>${zh ? "第一跳阈值 Z₁ = σ_jump/σ" : "1st Step Threshold Z₁"}</span><strong>Z = ±${alpha.toFixed(2)}</strong></div>
        <div class="ifig-metric"><span>${zh ? "极大似然期望 E[MSE_ML]" : "Expected E[MSE_ML]"}</span><strong style="color:#1f4e79;">${mseMlDb.toFixed(2)} dB</strong></div>
        <div class="ifig-metric"><span>${zh ? "贝叶斯期望 E[MSE_Bayes]" : "Expected E[MSE_Bayes]"}</span><strong style="color:#0b6e63;">${mseByDb.toFixed(2)} dB</strong></div>
        <div class="ifig-metric"><span>${zh ? "带宽直线 MMSE" : "Bandwidth Line MMSE"}</span><strong>${mmseDb.toFixed(2)} dB</strong></div>
        <div class="ifig-metric"><span>${zh ? "软加权节省 (ML − Bayes)" : "Bayes Saving (ML − Bayes)"}</span><strong style="color:#c45911;">+${(mseMlDb - mseByDb).toFixed(2)} dB</strong></div>
        <div class="ifig-metric"><span>${zh ? "门限 4.82 dB 归属拆解" : "4.82 dB Split at σ_jump"}</span><strong>0.24 + 1.20 + 3.38 dB</strong></div>
      `;

      const W = 920, H = 360;
      const ctx = setupHiDPICanvas(canvas, W, H);
      const L = 66, R = 28, T = 32, B = 48;
      const PW = W - L - R, PH = H - T - B;

      const zPlotMax = 4.5;
      const kPlotMax = 3.2;
      const px = (z) => L + ((z + zPlotMax) / (2.0 * zPlotMax)) * PW;
      const py = (k) => T + PH - ((k + kPlotMax) / (2.0 * kPlotMax)) * PH;

      ctx.fillStyle = "#ffffff";
      ctx.fillRect(0, 0, W, H);

      // Grid
      ctx.strokeStyle = "#e6eaee";
      ctx.lineWidth = 1;
      for (let k = -3; k <= 3; k++) {
        const yy = py(k);
        ctx.beginPath();
        ctx.moveTo(L, yy);
        ctx.lineTo(L + PW, yy);
        ctx.stroke();
        ctx.fillStyle = "#5a6772";
        ctx.font = "12px -apple-system, BlinkMacSystemFont, sans-serif";
        ctx.textAlign = "right";
        ctx.fillText(`${k}λ`, L - 8, yy + 4);
      }

      for (let z = -4; z <= 4; z++) {
        const xx = px(z);
        ctx.beginPath();
        ctx.moveTo(xx, T);
        ctx.lineTo(xx, T + PH);
        ctx.stroke();
        ctx.fillStyle = "#5a6772";
        ctx.font = "12px -apple-system, BlinkMacSystemFont, sans-serif";
        ctx.textAlign = "center";
        ctx.fillText(`Z=${z}`, xx, T + PH + 18);
      }

      // Draw standard normal density bell curve in background
      ctx.save();
      ctx.fillStyle = "rgba(11, 110, 99, 0.08)";
      ctx.beginPath();
      ctx.moveTo(px(-zPlotMax), py(0));
      for (let i = 0; i <= 300; i++) {
        const z = -zPlotMax + (i / 300) * (2.0 * zPlotMax);
        const bell = Math.exp(-0.5 * z * z) * 2.4;
        ctx.lineTo(px(z), py(bell));
      }
      ctx.lineTo(px(zPlotMax), py(0));
      ctx.closePath();
      ctx.fill();
      ctx.restore();

      // Draw MLE hard staircase
      ctx.save();
      ctx.strokeStyle = "#1f4e79";
      ctx.lineWidth = 2.6;
      ctx.beginPath();
      for (let i = 0; i <= 600; i++) {
        const z = -zPlotMax + (i / 600) * (2.0 * zPlotMax);
        const kMl = Math.max(-3, Math.min(3, Math.floor(z / (2.0 * alpha) + 0.5)));
        const x = px(z), y = py(kMl);
        if (i === 0) ctx.moveTo(x, y);
        else ctx.lineTo(x, y);
      }
      ctx.stroke();
      ctx.restore();

      // Draw Bayesian smooth Boltzmann curve
      ctx.save();
      ctx.strokeStyle = "#c45911";
      ctx.lineWidth = 3.0;
      ctx.beginPath();
      for (let i = 0; i <= 400; i++) {
        const z = -zPlotMax + (i / 400) * (2.0 * zPlotMax);
        let maxEx = -1e30;
        const exps = [];
        for (let k = -10; k <= 10; k++) {
          const diff = z - 2.0 * alpha * k;
          const ex = -0.5 * diff * diff;
          if (ex > maxEx) maxEx = ex;
          exps.push({ k, ex });
        }
        let num = 0.0, den = 0.0;
        exps.forEach((item) => {
          const w = Math.exp(item.ex - maxEx);
          den += w;
          num += item.k * w;
        });
        const kBy = Math.max(-3.1, Math.min(3.1, num / den));
        const x = px(z), y = py(kBy);
        if (i === 0) ctx.moveTo(x, y);
        else ctx.lineTo(x, y);
      }
      ctx.stroke();
      ctx.restore();

      // Legend on canvas
      ctx.fillStyle = "#1f4e79";
      ctx.font = "600 12px -apple-system, BlinkMacSystemFont, sans-serif";
      ctx.textAlign = "left";
      ctx.fillText(zh ? "── 相干极大似然硬阶梯 δ_ML(Z)" : "── Coherent MLE Hard Staircase δ_ML(Z)", L + 14, T + 22);
      ctx.fillStyle = "#c45911";
      ctx.fillText(zh ? "── 贝叶斯后验均值软阶梯 δ_Bayes(Z)" : "── Bayesian Posterior Mean Sigmoid δ_Bayes(Z)", L + 14, T + 42);
      ctx.fillStyle = "#0b6e63";
      ctx.fillText(zh ? "▒ 背景阴影：标准正态投影噪声密度 φ(Z)" : "▒ Shaded Bell: Standard Normal Noise Density φ(Z)", L + 14, T + 62);

      ctx.strokeStyle = "#182026";
      ctx.lineWidth = 1.5;
      ctx.strokeRect(L, T, PW, PH);

      ctx.fillStyle = "#182026";
      ctx.font = "600 13px -apple-system, BlinkMacSystemFont, sans-serif";
      ctx.textAlign = "center";
      ctx.fillText(zh ? "共线方向标准正态投影噪声  Z = uᵀe / σ" : "Projected Standard Normal Noise  Z = uᵀe / σ", L + 0.5 * PW, H - 10);

      ctx.save();
      ctx.translate(18, T + 0.5 * PH);
      ctx.rotate(-PI / 2);
      ctx.fillText(zh ? "整周滑移偏移量  δ / λ" : "Estimator Shift  δ / λ", 0, 0);
      ctx.restore();
    }

    slider.addEventListener("input", () => {
      curRatio = parseFloat(slider.value);
      render();
    });

    render();
  }

  // =========================================================================
  // INTERACTIVE FIGURE 4: Cumulative Ziv-Zakai Bound & Exact Wrapped LRT (8.3)
  // =========================================================================
  function mountValleyCumulativeExplorer(container, data, lang) {
    const zh = lang === "zh";
    let mode = "cume"; // "cume" or "lrt"

    container.innerHTML = `
      <div class="ifig-header">
        <div>
          <span class="ifig-badge">${zh ? "交互式图 4 · 谷填充积分与精确包裹似然比" : "Interactive Fig 4 · Valley Fill & Exact Wrapped LRT"}</span>
          <h4 class="ifig-title" style="display:inline;">${
            zh
              ? "半窗下界沿距离的累计占比 vs 大噪声远谷簇的精确包裹似然比检验"
              : "Cumulative Half-Window Bound Fraction & Exact Wrapped LRT on Far Valleys"
          }</h4>
        </div>
        <div class="ifig-btn-group">
          <button class="ifig-btn active" data-mode="cume">${zh ? "累计积分占比 (σ = 0.0836 rad)" : "Cumulative Fraction (σ = 0.0836 rad)"}</button>
          <button class="ifig-btn" data-mode="lrt">${zh ? "精确包裹似然比核验 (σ = 0.836 rad)" : "Exact Wrapped LRT Check (σ = 0.836 rad)"}</button>
        </div>
      </div>
      <div class="ifig-readout" data-role="v-readout"></div>
      <div class="ifig-canvas-wrap">
        <canvas class="ifig-canvas" width="920" height="370"></canvas>
      </div>
    `;

    const canvas = container.querySelector("canvas");
    const readout = container.querySelector('[data-role="v-readout"]');
    const modeBtns = container.querySelectorAll("[data-mode]");

    function render() {
      const W = 920, H = 370;
      const ctx = setupHiDPICanvas(canvas, W, H);
      const L = 70, R = 28, T = 32, B = 48;
      const PW = W - L - R, PH = H - T - B;

      ctx.fillStyle = "#ffffff";
      ctx.fillRect(0, 0, W, H);

      if (mode === "cume") {
        readout.innerHTML = `
          <div class="ifig-metric"><span>${zh ? "等间距 1λ / 2λ 占比" : "Equal 1λ / 2λ Share"}</span><strong style="color:#1f4e79;">67% / 29% (0.25m 达 96%)</strong></div>
          <div class="ifig-metric"><span>${zh ? "矩匹配 0.4m 内占比" : "Matched ≤ 0.4m Share"}</span><strong style="color:#c45911;">100% (与等间距重合)</strong></div>
          <div class="ifig-metric"><span>${zh ? "几何间隔 6.21m 远谷段" : "Geometric 6.21m Segment"}</span><strong style="color:#2e7d32;">96% (延伸至 6.211 m)</strong></div>
          <div class="ifig-metric"><span>${zh ? "几何间隔积分中位数 d₅₀" : "Geometric Median d₅₀"}</span><strong>4.34 m (d₉₀ = 5.91 m)</strong></div>
        `;

        const lxMin = Math.log10(0.03), lxMax = Math.log10(15.0);
        const px = (d) => L + ((Math.log10(Math.max(d, 0.03)) - lxMin) / (lxMax - lxMin)) * PW;
        const py = (f) => T + PH - f * PH;

        ctx.strokeStyle = "#e6eaee";
        ctx.lineWidth = 1;
        for (let k = 0; k <= 5; k++) {
          const fv = k * 0.2;
          const yy = py(fv);
          ctx.beginPath();
          ctx.moveTo(L, yy);
          ctx.lineTo(L + PW, yy);
          ctx.stroke();
          ctx.fillStyle = "#5a6772";
          ctx.font = "12px -apple-system, BlinkMacSystemFont, sans-serif";
          ctx.textAlign = "right";
          ctx.fillText(`${Math.round(fv * 100)}%`, L - 8, yy + 4);
        }

        [0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10].forEach((dt) => {
          const xx = px(dt);
          ctx.beginPath();
          ctx.moveTo(xx, T);
          ctx.lineTo(xx, T + PH);
          ctx.stroke();
          ctx.fillStyle = "#5a6772";
          ctx.font = "12px -apple-system, BlinkMacSystemFont, sans-serif";
          ctx.textAlign = "center";
          ctx.fillText(`${dt} m`, xx, T + PH + 18);
        });

        const styles = [
          { key: "equal", col: "#1f4e79", lw: 3.0, dash: null },
          { key: "matched", col: "#c45911", lw: 2.4, dash: [6, 4] },
          { key: "geometric", col: "#2e7d32", lw: 3.0, dash: null },
        ];
        styles.forEach((st) => {
          const blk = data.bound_valleys.sets[st.key];
          ctx.save();
          ctx.strokeStyle = st.col;
          ctx.lineWidth = st.lw;
          if (st.dash) ctx.setLineDash(st.dash);
          ctx.beginPath();
          let started = false;
          blk.snap_delta_m.forEach((d, idx) => {
            if (d < 0.03) return;
            const x = px(d), y = py(blk.snap_cume_frac[idx]);
            if (!started) {
              ctx.moveTo(x, y);
              started = true;
            } else {
              ctx.lineTo(x, y);
            }
          });
          ctx.stroke();
          ctx.restore();
        });

        ctx.strokeStyle = "#182026";
        ctx.lineWidth = 1.5;
        ctx.strokeRect(L, T, PW, PH);

        ctx.fillStyle = "#182026";
        ctx.font = "600 13px -apple-system, BlinkMacSystemFont, sans-serif";
        ctx.textAlign = "center";
        ctx.fillText(zh ? "距离偏移  δ (m，对数坐标)" : "Range Offset  δ (m, log scale)", L + 0.5 * PW, H - 10);
      } else {
        const lrt = data.open_questions.q8_3_exact_wrapped_lrt;
        readout.innerHTML = `
          <div class="ifig-metric"><span>${zh ? "高斯检验下界 (Equal / Matched)" : "Gaussian ZZB (Eq / Mat)"}</span><strong>−3.80 dB / −2.42 dB (Δ = 1.38 dB)</strong></div>
          <div class="ifig-metric"><span>${zh ? "10.6–11.3m 换精确包裹似然比" : "Exact WN LRT on 10.6–11.3m"}</span><strong style="color:#c45911;">−3.76 dB / −2.05 dB (Δ = 1.70 dB)</strong></div>
          <div class="ifig-metric"><span>${zh ? "全半窗 J≤30 换精确包裹似然比" : "Exact WN LRT on All J≤30"}</span><strong style="color:#0b6e63;">−3.17 dB / −1.83 dB (Δ = 1.34 dB)</strong></div>
          <div class="ifig-metric"><span>${zh ? "物理机制结论" : "Physical Mechanism"}</span><strong>${zh ? "±π 折回使远谷 Pe 升高 21–24%" : "±π wrapping raises valley Pe by 21–24%"}</strong></div>
        `;

        // Draw bar comparison of Gaussian vs Exact Cluster vs Exact Full
        const groups = [
          {
            label: zh ? "高斯向量检验 (Gaussian Test)" : "Gaussian Vector Test",
            eq: lrt.equal.zzb_gaussian_db,
            ma: lrt.matched.zzb_gaussian_db,
            gap: lrt.gap_matched_minus_equal_db.gaussian_db,
          },
          {
            label: zh ? "10.6–11.3m 谷簇换精确包裹似然比" : "Exact Wrapped LRT on 10.6–11.3m Cluster",
            eq: lrt.equal.zzb_exact_cluster_10p6_11p3_db,
            ma: lrt.matched.zzb_exact_cluster_10p6_11p3_db,
            gap: lrt.gap_matched_minus_equal_db.exact_cluster_only_db,
          },
          {
            label: zh ? "半窗全部 J≤30 谷换精确包裹似然比" : "Exact Wrapped LRT Across All Valleys",
            eq: lrt.equal.zzb_exact_all_valleys_db,
            ma: lrt.matched.zzb_exact_all_valleys_db,
            gap: lrt.gap_matched_minus_equal_db.exact_all_valleys_db,
          },
        ];

        const yMin = -4.5, yMax = 0.0;
        const py = (db) => T + PH - ((db - yMin) / (yMax - yMin)) * PH;

        ctx.strokeStyle = "#e6eaee";
        ctx.lineWidth = 1;
        [-4, -3, -2, -1, 0].forEach((db) => {
          const yy = py(db);
          ctx.beginPath();
          ctx.moveTo(L, yy);
          ctx.lineTo(L + PW, yy);
          ctx.stroke();
          ctx.fillStyle = "#5a6772";
          ctx.font = "12px -apple-system, BlinkMacSystemFont, sans-serif";
          ctx.textAlign = "right";
          ctx.fillText(`${db.toFixed(1)} dB`, L - 8, yy + 4);
        });

        const gw = PW / 3;
        groups.forEach((g, idx) => {
          const cx = L + (idx + 0.5) * gw;
          const bw = 56;
          const yEq = py(g.eq), yMa = py(g.ma), yBot = py(yMin);

          ctx.fillStyle = "#1f4e79";
          ctx.fillRect(cx - bw - 6, yEq, bw, yBot - yEq);
          ctx.fillStyle = "#c45911";
          ctx.fillRect(cx + 6, yMa, bw, yBot - yMa);

          ctx.fillStyle = "#1f4e79";
          ctx.font = "600 12px -apple-system, BlinkMacSystemFont, sans-serif";
          ctx.textAlign = "center";
          ctx.fillText(`${g.eq.toFixed(2)} dB`, cx - 0.5 * bw - 6, yEq - 8);

          ctx.fillStyle = "#c45911";
          ctx.fillText(`${g.ma.toFixed(2)} dB`, cx + 0.5 * bw + 6, yMa - 8);

          ctx.fillStyle = "#085248";
          ctx.font = "600 13px -apple-system, BlinkMacSystemFont, sans-serif";
          ctx.fillText(`Δ = +${g.gap.toFixed(2)} dB`, cx, Math.min(yEq, yMa) - 28);

          ctx.fillStyle = "#182026";
          ctx.font = "500 12px -apple-system, BlinkMacSystemFont, sans-serif";
          ctx.fillText(g.label, cx, T + PH + 22);
        });

        ctx.strokeStyle = "#182026";
        ctx.lineWidth = 1.5;
        ctx.strokeRect(L, T, PW, PH);
      }
    }

    modeBtns.forEach((btn) => {
      btn.addEventListener("click", () => {
        modeBtns.forEach((b) => b.classList.remove("active"));
        btn.classList.add("active");
        mode = btn.getAttribute("data-mode");
        render();
      });
    });

    render();
  }

  async function boot() {
    initLangSwitchPreserver();

    const containers = document.querySelectorAll("[data-ifig]");
    if (!containers.length) return;

    try {
      const resp = await fetch("data/results.json");
      const data = await resp.json();
      const pageLang = (document.documentElement.lang || "en").startsWith("zh") ? "zh" : "en";

      containers.forEach((el) => {
        const kind = el.getAttribute("data-ifig");
        const lang = el.getAttribute("data-lang") || pageLang;
        if (kind === "overview") mountOverviewExplorer(el, data, lang);
        else if (kind === "synthesizer") mountSynthesizerLab(el, data, lang);
        else if (kind === "bayes") mountBayesExplorer(el, data, lang);
        else if (kind === "valleys") mountValleyCumulativeExplorer(el, data, lang);
      });
    } catch (err) {
      console.error("Failed to initialize interactive figures:", err);
    }
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", boot);
  } else {
    boot();
  }
})();
