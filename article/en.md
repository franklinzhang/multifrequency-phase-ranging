# The jump between two rulers

*A story for engineers that estimation theorists can cite. Full derivations, every experiment, and three design Q&As: [deep dive](deep.en.md). 中文：[README.md](README.md). Numbers: [numbers.md](numbers.md).*

---

Open an error curve from Bluetooth 6.0 Channel Sounding, a stepped-frequency radar, or a multi-wavelength interferometer, and you almost always see the same picture:

![Figure 1. Error rides a fine line, then jumps onto a coarse one](figures/en/f01_hero.png)

**Figure 1.** Eight equally spaced tones (2400–2470 MHz, 10 MHz apart). Blue dots are simulation; the red curve is the closed-form model below — no fitted parameters.

At high SNR the error follows a very low straight line. At some point it jumps almost vertically by about 40 dB, then follows a second, parallel line. Engineers call this a threshold. The Cramér–Rao bound in the textbooks does not jump here at all.

This article tells one story: **how two rulers work together, how they fall out at a precise moment, why that moment has a formula, and which cell on Monday’s link budget you should change.**

You do not need a book on parameter estimation or signal processing. Every formula is a block of ice that freezes the story into a number — skip the ice, and the story still runs.

---

## Four acts

| Act | The sentence you take home |
| --- | --- |
| [I. Two rulers](#act-i-two-rulers) | Carrier is the fine ruler, bandwidth the coarse one; when the fine ruler picks the wrong tick, error jumps onto the coarse line. |
| [II. The falling-out](#act-ii-the-falling-out) | Falling out = the estimator miscounts whole cycles; miscounting is rounding a Gaussian to an integer. |
| [III. One budget card](#act-iii-one-budget-card) | Threshold centre ≈ total SNR $R_{\mathrm{dB}}-13$; safe zone ≈ $R_{\mathrm{dB}}+1$ (for $R\approx40$ dB). |
| [IV. How you place tones](#act-iv-how-you-place-tones) | Equal spacing minimizes ambiguity scale; to postpone the threshold, use edges plus enough anchors. |

---

## Act I. Two rulers

Suppose you have phases and no clock. Over a distance $d$, a carrier at frequency $f$ turns by $2\pi f d/c$ — but the receiver only sees the leftover after whole turns. One turn at 2.4 GHz is 12.3 cm. So **one tone only tells you which 12.3 cm tick you are on**; which tick number, it does not know.

Measure a few tones and two rulers appear (Figure 2):

![Figure 2. Fine ruler as a comb, coarse ruler as a parabola](figures/en/f02_two_rulers.png)

**Figure 2.** Horizontal axis: how far the guess is wrong. Blue is the true cost; the light-blue comb returns to zero every wavelength — the fine ruler; the orange parabola rises slowly — the coarse ruler.

- **Fine ruler (carrier).** All tones sit near 2.4 GHz and turn together with period 12.3 cm. Very fine marks; repeats every 12.3 cm.
- **Coarse ruler (bandwidth).** Tones differ by tens of megahertz, so phase differences drift slowly with range. Coarse marks; no repeat over tens of metres.

They cooperate intuitively: the coarse ruler first says “roughly which tick”, the fine ruler reads the leftover. As long as the coarse error is well below half a wavelength, cooperation holds and you get fine-ruler precision. Once the coarse ruler wobbles by half a wavelength, it points at the wrong tick and the error becomes “off by a whole stick” — millimetres to decimetres in one step.

The gap between the two error lines is a purely geometric number of the frequency set:

\[
R = \frac{\sum f_i^2}{\sum (f_i-\bar f)^2}.
\]

On our eight tones, $R = 40.5\,\mathrm{dB}$: the fine ruler is about 106× more precise than the coarse one. At 35 dB per tone, about 0.09 mm versus 9.3 mm.

> **Takeaway:** The threshold is not noise shredding the phase. It is the **coarse ruler losing the right to pick a tick for the fine one**.

---

## Act II. The falling-out

Spread out every possible range and the likelihood looks like a row of teeth (Figure 3). The truth is the middle tooth.

![Figure 3. As noise grows, the teeth become equally tall](figures/en/f03_likelihoods.png)

**Figure 3.** One noise sample at three noise levels. [Animation](figures/en/f03_likelihood.gif).

Low noise: the middle tooth towers over its neighbours — the estimator bites it. Medium noise: neighbours are almost as tall — a small push and it bites the next one, off by exactly 12.3 cm. Still more noise: the fine teeth smear out and only the coarse envelope remains.

That is the threshold of multimodal maximum likelihood. Usually there is no closed form. **Here the teeth line up: the next tooth, and the one after, all lie on one line.** The whole story collapses to one sentence in plain language:

> The estimator is really doing one thing — **rounding a Gaussian random number to an integer $k$**. That $k$ is how many wavelengths it miscounted.

The probability of a miscount, and the mean-square error after it, are entirely set by that rounding. We call this the **lattice model**. No fitted parameters; against simulation the gap is under 0.2 dB, checked down to slip probabilities of $10^{-5}$ (Figure 4).

![Figure 4. Model and simulation coincide](figures/en/f04_validation.png)

**Figure 4.** Left: slip probability. Right: error relative to the coarse line. Red is the model.

More important: change the tone count, change the bandwidth, normalize by “noise / falling-out scale”, and the curves stack (Figure 6). The shape of the falling-out depends on only one number, $R$.

![Figure 6. Eight designs collapse onto one curve](figures/en/f06_collapse.png)

<details>
<summary>Open for the formulas (skippable)</summary>

The phase scale of the falling-out is $\sigma_{\mathrm{jump}} = \pi\sqrt{S}/\bar f$ (0.084 rad in the example). Normalized noise $u=\sigma/\sigma_{\mathrm{jump}}$. The integer decision has std $s = u\big/\bigl(2\sqrt{1-1/R}\bigr)$; MSE is the fine-ruler error plus $\lambda^2\,\mathbb{E}[\hat k^2]$. Derivation in the [deep dive, Appendix A](deep.en.md#appendix-a-derivation-of-the-lattice-model).

</details>

> **Takeaway:** You do not need “ambiguity fixing” or “Ziv–Zakai bounds” to use this model. It is one sentence: wrong tick = rounding.

---

## Act III. One budget card

Finish the rounding and engineering is a few rows. For any tone set, compute $R$, then read:

| What you want | Answer (total SNR) | Eight tones, 70 MHz (per tone) |
| --- | --- | --- |
| Safe zone: budget the fine ruler | $\gtrsim R_{\mathrm{dB}} + 1$ (only for $R\approx40$–45 dB; in general $R_{\mathrm{dB}}-12.95+20\log_{10}x(R)$, see the [math note](math.zh.md) §3.3, in Chinese) | $\gtrsim 32.5\,\mathrm{dB}$ |
| Centre of the falling-out | $\approx R_{\mathrm{dB}} - 13$ | $18.5\,\mathrm{dB}$ |
| What the centre is physically | coarse RMS error $= \lambda/2$ | 6.2 cm |
| Chance of the wrong tick there | about 32% | — |
| Shelf after the falling-out | coarse line, about +1.1 dB | — |
| Even the coarse line is lost | total SNR about 8–12 dB | a second threshold; use complex samples |

**Move from 2.4 GHz to 5.8 GHz with the same spacing?** Coarse precision unchanged, but the centre rises from 18.5 to 26 dB — higher carrier, narrower relative bandwidth, harder to keep the fine ruler. Details in [deep §12.1](deep.en.md#121-does-changing-the-starting-frequency-matter).

**Group delay, common phase, multipath — who moves the threshold?**

- Group delay: a global range shift; **threshold unchanged**.
- A phase offset common to every tone: the fine ruler thinks you moved, the coarse one saw nothing — the rulers quarrel — **threshold moves**. Within about 20°, the cost stays under 1 dB.
- One reflection: as excess path sweeps a few metres, the disagreement can swing by several wavelengths; past half a wavelength, **the wrong tick is chosen even with no noise** (Figure 9).

![Figure 9. Who makes the two rulers quarrel](figures/en/f09_bias.png)

**Figure 9.** Left: common phase raises the miscount rate. Right: multipath swings the disagreement $\beta$; error peaks follow.

> **Takeaway:** On the budget sheet, write $R$ first, then two SNR lines (safe zone / centre). Calibrate group delay for range accuracy; calibrate inter-tone phase disagreement so the tick is not pointed wrong.

---

## Act IV. How you place tones

The last act is a choice.

**Equal spacing is harder than people think.** All its deep valleys sit at $1\lambda$, $2\lambda$, $3\lambda$, …. When it falls out, you are off by 12.3 cm, not four metres. That is **minimal ambiguity scale** — if your only fear is “how far wrong when it fails”, equal spacing is the structural answer.

**Drawing a few tones uniformly at random from many candidates is not optimal.** The median usually falls out earlier and has a higher shelf than equal spacing.

**To postpone the falling-out, push tones to the edges** — $S$ grows, $\sigma_{\mathrm{jump}}$ grows. Pure edges leave a valley near 4.5 m as deep as $1\lambda$: later threshold, ruined ambiguity scale. Keep **anchors** in the middle to fill that far valley.

An anchor is not another ruler. At the far valley the edge residuals are near zero; a middle frequency is off an integer cycle, the cost rises, and the valley becomes unreachable. On eight tones, the set

```
2400, 2401, 2402, 2409, 2450, 2468, 2469, 2470 MHz
```

keeps the $1\lambda$ ambiguity scale, yet falls out later and shelves lower than equal spacing (Figure 13).

![Figure 13. Later threshold, shelf kept](figures/en/f13_pareto.png)

**Figure 13.** Left: Pareto front of “late threshold × low shelf”. Right: Monte Carlo check.

**Many tones at the edges and only one in the middle — drowned? Your intuition is right.** Not because that tone’s SNR is low, but because thick edges make $J(1\lambda)$ grow linearly while one anchor’s lift of the far valley barely grows (the valley also relocates to dodge it). With about 1000 edge tones, one middle tone is wasted; the number of anchors must scale with the edges. The metric is far-valley cost / $J(1\lambda)$; about 10 is hard enough. See [deep §12.2](deep.en.md#122-drawing-n-out-of-m-equal-spaced-candidates-is-uniform-random-best).

> **Takeaway:** Ask whether you fear “falling out earlier” or “falling out by four metres”. Fear the first → edges. Fear the second → equal spacing, or edges plus enough anchors. Do not let a random draw make that choice for you.

---

## Curtain: why this picture arrived so late

Multi-frequency phase ranging is not a small problem, yet “where the jump is, and how high” is rarely written in closed form. It sat in the gaps between communities: threshold theory mostly faces continuous waveforms; GNSS cares about fixing success rates, not the whole error curve; Chinese-remainder work optimizes unambiguous range — a different knob; engineering uses synthetic-wavelength rules of thumb and paints the threshold region “does not work”.

What we filled in is the middle piece: **write the falling-out as rounding, then fold SNR, calibration, and tone placement into one card.**

Limits, stated plainly: designs with shallow far valleys need another approximation; several multipaths and correlated phase noise are not yet systematic; a tighter global bound for the downslope of the threshold is still missing. Details, the complex-observation check, and the Bayes split (1.2 dB for hard decisions at the threshold; 3.2 dB of Ziv–Zakai looseness, 2.5 dB of which a single-test-point Weiss–Weinstein bound recovers) live in the [deep dive](deep.en.md) and the [math note](math.zh.md).

---

## If you will spend two more minutes

1. Compute your $R$.  
2. Safe budget: total SNR $\gtrsim R_{\mathrm{dB}}+1$ (valid for $R\approx40$–45 dB; for other $R$ use the $x(R)$ table in §3.3 of the [math note](math.zh.md)); centre $\approx R_{\mathrm{dB}}-13$.  
3. Leave about 20° for common phase; treat multipath as “will the two rulers point at different ticks”, not only as a noise margin.  
4. Placement: equal spacing = smallest ambiguity scale; later threshold → edges + anchors in proportion.  

Code: [`computation/run_all.sh`](../computation/run_all.sh). References and errata: [end of the deep dive](deep.en.md#references).
