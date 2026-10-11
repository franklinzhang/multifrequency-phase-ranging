# Range Estimation from Wrapped Multi-Frequency Phase Observations

Wrapped phases at a finite set of frequencies determine a range only up to a lattice of aliases. The coherent maximum-likelihood estimator uses those phases against one common phase reference. On a log-log plot its mean-square error follows the carrier line, then jumps, at the same slope, onto the bandwidth line.

**Start here:** the short story [*The jump between two rulers*](article/en.md) ([中文](article/README.md)) — four acts, one budget card, no textbook required. Full derivations and experiments: [deep.en.md](article/deep.en.md) / [deep.zh.md](article/deep.zh.md). Code: `computation/mfpr`, `computation/experiments/e01`–`e09`, `computation/run_all.sh`.

The Chinese notes in this repository are the original account. This page states the model, the two lines, the first threshold, what unequal spacing actually changes, and what the bound does not claim.

## Model

The observation at frequency \(f_i\) is

\[
\varphi_i = \mathrm{wrap}(2\pi f_i d/c + e_i), \qquad e_i \sim \mathcal{N}(0,\sigma^2),
\]

with independent phase errors and \(c = 299\,792\,458\,\mathrm{m/s}\). The coherent estimator minimizes the sum of squared wrapped residuals. If every tone shares an unknown phase, the estimator maximizes the modulus of the phase sum, and the carrier frequency drops out of the range information.

This phase model is the high-SNR edge of complex observations \(y_i = \exp(\mathrm{j}\,2\pi f_i d/c) + n_i\). The first threshold derived below sits near \(0.08\,\mathrm{rad}\) for the equal grid used here, where that edge is still valid. Phase noise of order \(1\,\mathrm{rad}\) is a second threshold: the phase edge itself fails, and the budget has to return to the complex observations.

## Two parallel lines

With \(F = \sum_i f_i^2\) and \(S = \sum_i (f_i-\bar f)^2\),

\[
\mathrm{HMSE} = \frac{c^2\sigma^2}{4\pi^2 F}, \qquad \mathrm{MMSE} = \frac{c^2\sigma^2}{4\pi^2 S}.
\]

The ratio \(R = F/S\) does not depend on \(\sigma\), so the lines are parallel. For \(N\) equally spaced tones of bandwidth \(B\),

\[
R = 1 + \frac{12(N-1)}{N+1}\left(\frac{\bar f}{B}\right)^2.
\]

The coefficient tends to 12 for large \(N\). That is the factor Weiss and Weinstein obtained for a continuous flat spectrum. The observation model is different, and a discrete grid is not that spectrum. At \(N = 8\) the coefficient is 9.33, about 1.1 dB lower than the continuous-spectrum formula.

On eight tones from 2400 to 2470 MHz in steps of 10 MHz, \(\bar f = 2435\,\mathrm{MHz}\), \(S = 4200\,\mathrm{MHz}^2\), and \(R = 40.53\,\mathrm{dB}\).

Li Wei and Qi Wangdong give these two heights (arXiv:1202.1424, equations 30 and 42) by perturbation. They do not write the noise level of the jump as a function of the frequency set.

## Where the jump sits

The phase scale of a one-wavelength error is

\[
\sigma_{\mathrm{jump}} = \pi\sqrt{S}/\bar f.
\]

For the equal grid this is \(0.0836\,\mathrm{rad}\), about 4.8 degrees. A single phase sample is nowhere near \(\pm\pi\) there. The wrap happens because the range hypothesis has crossed one carrier cycle.

A Gaussian test between two range hypotheses separated by \(\delta\) has error probability \(P_e(\delta) = Q(\sqrt{J(\delta)}/(2\sigma))\). At one mean wavelength the argument is \(\sigma_{\mathrm{jump}}/\sigma\). This \(P_e\) is the Gaussian phase-vector test, not the fully wrapped likelihood ratio.

The mixture \((1-p)\mathrm{HMSE} + p\lambda^2\) with \(p = 2Q(\sigma_{\mathrm{jump}}/\sigma)\) leaves the carrier line where \(p\lambda^2 = \mathrm{HMSE}\). On the equal grid that foot is \(0.0170\,\mathrm{rad}\), or \(0.204\,\sigma_{\mathrm{jump}}\). It is the solution of the mixture equation, not a knee read off the measured curve. In 8000 trials the first slip appears at \(0.25\,\sigma_{\mathrm{jump}}\), and the sample mean square there equals one wavelength error in 8000 trials. From \(0.25\) to \(1.0\,\sigma_{\mathrm{jump}}\) the sample mean square equals the slip rate times \(\lambda^2\), within 0.1 dB. At the threshold, \(\lambda^2 = 4\,\mathrm{MMSE}\), so \(p = 2Q(1) = 0.317\) predicts +1.0 dB above the bandwidth line. The measured slip rate is 0.344 and the sample mean square is −22.8 dB against a bandwidth line of −24.2 dB, which is +1.4 dB.

Local Fisher information of a wrapped-normal phase is smooth in \(\sigma\). The jump is the global estimator leaving the main lobe.

An estimator that does not use absolute carrier phase never sees this jump. The incoherent estimator, which maximizes the modulus of the phase sum, lies on the bandwidth line from the lowest noise in the trial (−44.4 dB against a line of −44.2 dB). \(\sigma_{\mathrm{jump}}\) is not a threshold for every estimator.

A round trip replaces \(c\) by \(c/2\). Both mean squares divide by 4. \(\sigma_{\mathrm{jump}}\) in radians is unchanged.

## Three knobs, not substitutes

\(S/\bar f^2\) sets the phase scale of the first threshold and the gap between the lines. Spreading the tones shortens the jump and moves it to a larger phase noise.

Each deep valley’s cost \(J(\delta)\) and its squared distance \(\delta^2\) decide whether a bandwidth shelf exists, and whether a blunder costs one wavelength squared or several metres squared. A valley with \(J(\delta)\approx J(\lambda)\) is as easy to slip into as a one-wavelength error.

The greatest common divisor of the tones sets the exact period on which every phase returns to an integer cycle, which is the unambiguous range. A coprime or Chinese-remainder layout moves only that period. It does not postpone the first threshold.

Keeping \(N\), \(\bar f\), and \(S\) fixed, and only breaking the equal spacing onto integer megahertz, leaves the first threshold where it was. At \(\sigma = 0.0836\,\mathrm{rad}\) the half-window bound splits into the same pieces as the equal grid. The two designs diverge only at larger noise.

A geometric progression with ratio 1.02, centered near 2435 MHz, has its own second moment: \(R = 26.87\,\mathrm{dB}\) and \(\sigma_{\mathrm{jump}} = 0.403\,\mathrm{rad}\). The valley at 6.211 m has the same cost as one mean wavelength, \(J = 0.650\), and phase residuals no larger than 0.44 rad. The squared error is about 34 dB worse than one wavelength. That valley does not fill all the way back to the main lobe, because the prior weight \((L-\delta)\) is smaller there. It owns the segment from one wavelength out to 6.211 m, and that segment is 96% of the half-window integral. The mass is spread along the segment: half of it lies beyond 4.3 m, and nine tenths lie inside 5.9 m. It is not a spike at 6.211 m. At the equal-grid threshold, which is only 0.21 of the geometric threshold, the half-window bound is already −51 dB while the carrier line is still −64.8 dB. The bandwidth shelf is gone.

At \(\sigma = 0.836\,\mathrm{rad}\), ten times the equal-grid threshold, the moment-matched bound is −2.4 dB and the equal-grid bound is −3.8 dB. A cluster of valleys from 10.6 to 11.3 m, with cost about 14 against about 25 on the equal grid, accounts for the gap. The record at 10.65 m carries 29% of the matched integral. Removing that segment drops the bound from −2.4 dB to about −3.9 dB. Those valleys are not exact integer-cycle aliases.

## What the bound is

The half-window Ziv–Zakai bound, with the Gaussian phase-vector \(P_e\) and valley fill, returns to the carrier line at high SNR. The identity \(\int_0^\infty u\,Q(u)\,\mathrm{d}u = 1/4\) makes the numerical difference 0.00 dB. The bound is the Bayes risk of a uniform prior on a window of length 14.990 m. The trials are conditional at the window center. At the equal-grid threshold the bound is −27.6 dB, the bandwidth line is −24.2 dB, and the coherent trial is −22.8 dB. The bound rises before the trials do. It is not the knee.

The half-window prior variance is 12.7 dB. The trial floor is the variance of the full search window, 18.74 dB. At \(\sigma = 2.34\,\mathrm{rad}\) the equal-grid bound is still 8.1 dB, below its own prior variance. The full-period bound is inflated by a replica valley and is not the curve reported here.

A right-to-left record of \((L-\delta)P_e(\delta)\) owns the segment from the previous higher record to itself. If the left edge meets the main lobe, that edge is \(\sqrt{J}/\beta\). On the equal grid at its own threshold the three contributing pieces are one, two, and three wavelengths: −29.2, −32.9, and −43.0 dB, summing to the stored integral −27.6 dB.

## Engineering rules

With a phase reference, budget the carrier line only below about \(0.2\,\sigma_{\mathrm{jump}}\). On this equal grid that is a phase standard deviation of about one degree. Between about \(0.5\,\sigma_{\mathrm{jump}}\) and the onset of phase wrap, budget the bandwidth line, which is 40.53 dB worse here. Without a common phase, budget the bandwidth line directly. A far valley with \(J\approx J(\lambda)\) still lifts the error after the carrier line has been removed.

When placing tones, solve \(P_e(\delta)\,\delta^2 = \mathrm{HMSE}\) at every deep valley and keep the earliest \(\sigma\). Quoting only \(\sigma_{\mathrm{jump}}\) misses the 6.211 m valley. Do not use a longer unambiguous range to push the first threshold back.

## Limits

There is no closed form for the global wrapped least-squares mean square at finite \(N\). The mixture foot is earlier than the rise this trial size can see. The bound is a bound. The incoherent estimator is the counterexample to “every estimator jumps at \(\sigma_{\mathrm{jump}}\)”.

## Open questions

Five questions are ready to be closed. Each has a completion test.

1. Place tones under a band limit, a fixed \(N\), and an integer-megahertz synthesizer step so that the earliest leaving noise is as large as possible, and no valley with \(J(\delta)\approx J(\lambda)\) sits at \(\delta\gg\lambda\). Done when a design beats the equal grid’s first threshold and has no such valley. Shuffling gaps at fixed \(S\) does not.
2. Compute the posterior-mean mean square under the same uniform prior as the bound, and compare it with coherent maximum likelihood. Done when the two curves, on that prior, attribute the roughly 5 dB gap at the threshold (−27.6 against −22.8 dB) either to a loose bound or to a suboptimal estimator. The same comparison says whether integer least squares, or marginalizing the ambiguity, is worth doing.
3. Replace the Gaussian-vector \(P_e\) by the wrapped likelihood ratio only on the 10.6–11.3 m cluster, where residuals reach about 2.5 rad. The first threshold and the 6.211 m valley are not the suspicious place. Done when it is clear whether the 1.4 dB gap at \(\sigma = 0.836\,\mathrm{rad}\) moves.
4. Repeat the coherent estimator on the complex observations, and only check the first threshold. The prediction is that the first jump stays and the floor moves. Amplitude fading, unequal tone SNR, and shared-oscillator phase noise are budget corrections after that check.
5. Past about \(1.5\,\sigma_{\mathrm{jump}}\), weight the record valleys by \(\delta^2\) instead of charging every slip as \(\lambda^2\). Done when that weighting tracks the sample mean square more closely on that side of the threshold, with the remaining gap stated. A closed form for the global mean square is not the goal.

Two tasks are not worth the effort. One is a proof that every estimator jumps at \(\sigma_{\mathrm{jump}}\). The incoherent estimator already does not. The other is using the Chinese remainder theorem to postpone the first threshold. That theorem moves the unambiguous range.

## Documents

The derivation, with the cost surfaces and both estimators drawn separately, is the note *多频相位测距的平行门限*. The consolidation of the conclusions, the engineering rules, and the open questions is *多频相位测距的门限*, also stored as Markdown under the same Chinese title. The three curves in the consolidation are the equal-spacing mean square, the cumulative fraction of the half-window bound at \(0.0836\,\mathrm{rad}\), and the geometric-spacing mean square.

The computation directory holds the frozen trials behind those curves and the programs that produced them. The numbers in the notes are the archived trials. A later run is not a substitute for those numbers unless the difference is stated.

## References

Wei Li and Wangdong Qi, Optimization in Multi-Frequency Interferometry Ranging: Theory and Experiment, arXiv:1202.1424, 2012.

A. Weiss and E. Weinstein, Fundamental Limitations in Passive Time-Delay Estimation — Part I, IEEE Transactions on Acoustics, Speech, and Signal Processing, 1983; Part II, 1984.

J. Ziv and M. Zakai, Some Lower Bounds on Signal Parameter Estimation, IEEE Transactions on Information Theory, 1969.

D. C. Rife and R. R. Boorstyn, Single-Tone Parameter Estimation from Discrete-Time Observations, IEEE Transactions on Information Theory, vol. 20, no. 5, pp. 591–598, 1974.
