# Range Estimation from Wrapped Multi-Frequency Phase Observations

The Chinese note *多频相位测距的门限* (`notes/门限问题的进展.md`) is the full account. This page is the English digest: what was established, what the numbers say, what to budget, and what remains open.

## Long abstract: what was done, what was answered

Wrapped phases at a finite set of frequencies determine a range only up to a lattice of aliases. For the coherent maximum-likelihood estimator, mean-square error on a log-log plot follows the carrier line, then jumps, at the same slope, onto the bandwidth line. The phase scale of that first threshold is a function of the frequency set,

\[
\sigma_{\mathrm{jump}} = \pi\sqrt{S}/\bar f,
\]

with \(S\) the second moment of the tones about their mean frequency \(\bar f\). On eight equal-spaced tones from 2400 to 2470 MHz in 10 MHz steps, \(\sigma_{\mathrm{jump}} = 0.0836\,\mathrm{rad}\) (about 4.8 degrees). The mixture \((1-p)\,\mathrm{HMSE} + p\lambda^2\) leaves the carrier line at \(0.204\,\sigma_{\mathrm{jump}}\). In 8000 trials the first slip appears at \(0.25\,\sigma_{\mathrm{jump}}\); from there to \(1.0\,\sigma_{\mathrm{jump}}\) the sample mean square equals the slip rate times \(\lambda^2\), within 0.1 dB. At the threshold the slip rate is 0.344, the sample mean square is −22.8 dB, and the bandwidth line is −24.2 dB.

The half-window Ziv–Zakai bound at the same threshold is −27.6 dB. It is the Bayes risk under a uniform prior on the half-window; the trials are conditional at the window center. The bound leaves the carrier line before the trials do. Local information does not jump. An incoherent estimator that never uses absolute carrier phase lies on the bandwidth line from the lowest noise in the trial, so \(\sigma_{\mathrm{jump}}\) is not a threshold for every estimator.

Unequal spacing does not automatically cancel the jump. Holding \(N\), \(\bar f\), and \(S\) fixed and only breaking equal spacing leaves the first threshold where it was. A geometric progression with ratio 1.02 raises \(S\) by about a factor of 23 and \(\sigma_{\mathrm{jump}}\) to 0.403 rad; a valley at 6.211 m with the same cost as one mean wavelength removes the bandwidth shelf and makes a blunder cost metres squared, not one wavelength squared. At ten times the equal-grid threshold, a cluster of valleys from 10.6 to 11.3 m lifts the moment-matched half-window bound about 1.4 dB above the equal-grid bound.

With a phase reference, budget the carrier line only below about \(0.2\,\sigma_{\mathrm{jump}}\); between about \(0.5\,\sigma_{\mathrm{jump}}\) and the failure of the phase model, budget the bandwidth line. Without a common phase, budget the bandwidth line directly. When placing tones, solve \(P_e(\delta)\,\delta^2 = \mathrm{HMSE}\) at every deep valley and keep the earliest \(\sigma\). A longer unambiguous range does not postpone the first threshold.

Five questions remain open and closable: place tones under band, \(N\), and integer-megahertz constraints against the worst valley; compare posterior-mean and coherent ML mean squares under the same prior as the bound; replace the Gaussian-vector \(P_e\) by the wrapped likelihood ratio only on the 10.6–11.3 m cluster; repeat the first-threshold check on complex observations; past about \(1.5\,\sigma_{\mathrm{jump}}\), weight record segments by \(\delta^2\). Two tasks are not worth doing: proving that every estimator jumps at \(\sigma_{\mathrm{jump}}\), and using the Chinese remainder theorem to postpone the first threshold.

## Observation and two estimators

The observation at frequency \(f_i\) is

\[
\varphi_i = \mathrm{wrap}(2\pi f_i d/c + e_i), \qquad e_i \sim \mathcal{N}(0,\sigma^2),
\]

with independent phase errors and \(c = 299\,792\,458\,\mathrm{m/s}\). The wrap map folds the angle into \((-\pi,\pi]\): the true phase may differ by any multiple of \(2\pi\), and the reading keeps only the principal value. A finite set of phases therefore determines range only up to aliases.

The coherent estimator minimizes the sum of squared wrapped residuals. It uses absolute phase relative to one local oscillator, so the carrier enters the range information. If every tone shares an unknown common phase, the estimator maximizes the modulus of the phase sum; the carrier drops out and only the spacings remain. That second estimator is the incoherent estimator below.

This phase model is the high-SNR edge of complex observations \(y_i = \exp(\mathrm{j}\,2\pi f_i d/c) + n_i\). The first threshold derived below sits near \(0.08\,\mathrm{rad}\) on the equal grid, where that edge is still valid. Phase noise of order \(1\,\mathrm{rad}\) is a second threshold: the phase edge itself fails, and the budget returns to the complex observations.

A round trip replaces \(c\) by \(c/2\). Both mean squares divide by 4. \(\sigma_{\mathrm{jump}}\) in radians is unchanged.

## Two parallel lines

With \(F = \sum_i f_i^2\) and \(S = \sum_i (f_i-\bar f)^2\),

\[
\mathrm{HMSE} = \frac{c^2\sigma^2}{4\pi^2 F}, \qquad \mathrm{MMSE} = \frac{c^2\sigma^2}{4\pi^2 S}.
\]

The ratio \(R = F/S\) does not depend on \(\sigma\), so the lines are parallel. Li Wei and Qi Wangdong give these two heights (arXiv:1202.1424, equations 30 and 42) by perturbation. They do not write the noise level of the jump as a function of the frequency set.

For \(N\) equally spaced tones of bandwidth \(B\),

\[
R = 1 + \frac{12(N-1)}{N+1}\left(\frac{\bar f}{B}\right)^2.
\]

The coefficient tends to 12 for large \(N\). That is the factor Weiss and Weinstein obtained for continuous-time delay estimation with a flat spectrum. The observation model is different; a discrete grid of wrapped phases is not that spectrum. At \(N = 8\) the coefficient is 9.33, about 1.1 dB lower than the continuous-spectrum formula at the same \(\bar f/B\).

On the equal grid used here, \(\bar f = 2435\,\mathrm{MHz}\), \(S = 4200\,\mathrm{MHz}^2\), \(R = 40.53\,\mathrm{dB}\), and \(\lambda = 0.12312\,\mathrm{m}\).

## Where the jump sits

The phase scale of a one-wavelength error is

\[
\sigma_{\mathrm{jump}} = \pi\sqrt{S}/\bar f.
\]

For the equal grid this is \(0.0836\,\mathrm{rad}\). A single phase sample is nowhere near \(\pm\pi\) there. The wrap happens because the range hypothesis has crossed one carrier cycle.

A Gaussian test between two range hypotheses separated by \(\delta\) has error probability \(P_e(\delta) = Q(\sqrt{J(\delta)}/(2\sigma))\). At one mean wavelength the argument is \(\sigma_{\mathrm{jump}}/\sigma\). This \(P_e\) is the Gaussian phase-vector test, not the fully wrapped likelihood ratio. Near the main lobe, where residuals are far below \(\pi\), the two are close.

The mixture \((1-p)\,\mathrm{HMSE} + p\lambda^2\) with \(p = 2Q(\sigma_{\mathrm{jump}}/\sigma)\) leaves the carrier line where \(p\lambda^2 = \mathrm{HMSE}\). On the equal grid that solution is \(0.0170\,\mathrm{rad}\), or \(0.204\,\sigma_{\mathrm{jump}}\). It is the solution of the mixture equation, not a turning point read off the measured curve. In 8000 trials one expects a first slip near \(0.26\,\sigma_{\mathrm{jump}}\); the observed first slip is at \(0.25\,\sigma_{\mathrm{jump}}\), and the sample mean square there equals one wavelength error in 8000 trials. From \(0.25\) to \(1.0\,\sigma_{\mathrm{jump}}\) the sample mean square equals the slip rate times \(\lambda^2\), within 0.1 dB. At the threshold, \(\lambda^2 = 4\,\mathrm{MMSE}\), so \(p = 2Q(1) = 0.317\) predicts +1.0 dB above the bandwidth line. The measured slip rate is 0.344 and the sample mean square is −22.8 dB against a bandwidth line of −24.2 dB (+1.4 dB).

Local Fisher information of a wrapped-normal phase is smooth in \(\sigma\). The local Cramér–Rao bound for range follows the carrier line and does not jump; on this equal grid it stays within 0.03 dB of the carrier line out to \(\sigma = 0.84\,\mathrm{rad}\). The jump is the global estimator leaving the main lobe. The same skeleton appears in single-tone frequency estimation (Rife and Boorstyn): high-SNR variance tracks the Cramér–Rao bound, then leaves it when the estimator locks onto a wrong spectral peak. There the parameter is frequency and the observations are not wrapped phases.

An estimator that does not use absolute carrier phase never sees this jump. The incoherent estimator lies on the bandwidth line from the lowest noise in the trial (−44.4 dB against a line of −44.2 dB).

## What the bound is

The half-window Ziv–Zakai bound, with the Gaussian phase-vector \(P_e\) and valley fill, returns to the carrier line at high SNR. The identity \(\int_0^\infty u\,Q(u)\,\mathrm{d}u = 1/4\) makes the numerical difference 0.00 dB. The bound is the Bayes risk of a uniform prior on a window of half-length 14.990 m. The trials are conditional at the window center. At the equal-grid threshold the bound is −27.6 dB, the bandwidth line is −24.2 dB, and the coherent trial is −22.8 dB. The bound rises before the trials do. It constrains mean square; it is not the turning point of the trial curve.

The half-window prior variance is 12.7 dB. The trial floor is the variance of the full search window, 18.74 dB. At \(\sigma = 2.34\,\mathrm{rad}\) the equal-grid bound is still 8.1 dB, below its own prior variance. The full-period bound is inflated by a replica valley and is not the curve reported here.

A right-to-left running high of \((L-\delta)P_e(\delta)\) accounts for the segment from the previous higher point to itself. If the left edge meets the main lobe, that edge is \(\sqrt{J}/\beta\). On the equal grid at its own threshold the three contributing pieces are one, two, and three wavelengths: −29.2, −32.9, and −43.0 dB, summing to the stored integral −27.6 dB.

## Three design quantities, not substitutes

\(S/\bar f^2\) sets the phase scale of the first threshold and the gap between the lines. Spreading the tones shortens the jump and moves it to a larger phase noise.

Each deep valley’s cost \(J(\delta)\) and its squared distance \(\delta^2\) decide whether a bandwidth shelf exists, and whether a blunder costs one wavelength squared or several metres squared. A valley with \(J(\delta)\approx J(\lambda)\) is as easy to slip into as a one-wavelength error.

The greatest common divisor of the tones sets the exact period on which every phase returns to an integer cycle, which is the unambiguous range. A coprime or Chinese-remainder layout moves only that period. It does not postpone the first threshold.

Keeping \(N\), \(\bar f\), and \(S\) fixed, and only breaking the equal spacing onto integer megahertz, leaves the first threshold where it was. At \(\sigma = 0.0836\,\mathrm{rad}\) the half-window bound splits into the same pieces as the equal grid. The two designs diverge only at larger noise.

A geometric progression with ratio 1.02, centered near 2435 MHz, has \(R = 26.87\,\mathrm{dB}\) and \(\sigma_{\mathrm{jump}} = 0.403\,\mathrm{rad}\). The valley at 6.211 m has the same cost as one mean wavelength, \(J = 0.650\), and phase residuals no larger than 0.44 rad. The squared error is about 34 dB worse than one wavelength. That valley does not fill all the way back to the main lobe, because the prior weight \((L-\delta)\) is smaller there. It accounts for the segment from one wavelength out to 6.211 m, and that segment is 96% of the half-window integral. The probability mass is spread along the segment: half of it lies beyond 4.3 m, and nine tenths lie inside 5.9 m. It is not a spike at 6.211 m. At the equal-grid threshold, which is only 0.21 of the geometric threshold, the half-window bound is already −51 dB while the carrier line is still −64.8 dB. The bandwidth shelf is gone.

At \(\sigma = 0.836\,\mathrm{rad}\), ten times the equal-grid threshold, the moment-matched bound is −2.4 dB and the equal-grid bound is −3.8 dB. A cluster of valleys from 10.6 to 11.3 m, with cost about 14 against about 25 on the equal grid, accounts for the gap. The right-to-left high at 10.65 m carries 29% of the matched integral. Removing that segment drops the bound from −2.4 dB to about −3.9 dB. Those valleys are not exact integer-cycle aliases.

## Engineering rules

With a phase reference, budget the carrier line only below about \(0.2\,\sigma_{\mathrm{jump}}\). On this equal grid that is a phase standard deviation of about one degree (0.017 rad). Between about \(0.5\,\sigma_{\mathrm{jump}}\) and the onset of phase wrap, budget the bandwidth line, which is 40.53 dB worse here. Without a common phase, budget the bandwidth line directly. A far valley with \(J\approx J(\lambda)\) still lifts the error after the carrier line has been removed.

When placing tones, solve \(P_e(\delta)\,\delta^2 = \mathrm{HMSE}\) at every deep valley and keep the earliest \(\sigma\). Quoting only \(\sigma_{\mathrm{jump}}\) misses the 6.211 m valley. Do not use a longer unambiguous range to push the first threshold back.

On this equal grid \(\sigma_{\mathrm{jump}} = 0.0836\,\mathrm{rad}\). Phase noise of order 1 rad is the second threshold: return the budget to complex observations. Do not mix the two thresholds.

## Closed and still open

What is closed: the two line heights; \(\sigma_{\mathrm{jump}}\) as the phase scale of a one-wavelength Gaussian test; the mixture leaving point at \(0.204\,\sigma_{\mathrm{jump}}\) and the trial first slip at \(0.25\,\sigma_{\mathrm{jump}}\) on the equal grid; local information that does not jump; the incoherent counterexample; the half-window bound as prior risk, not a trial turning point; the fact that moment-matched unequal spacing does not cancel the first jump, while a far valley with \(J\approx J(\lambda)\) removes the shelf; the separation of unambiguous range from the first threshold.

There is no closed form for the global wrapped least-squares mean square at finite \(N\). The numbers in the notes are the archived trials in `computation/`.

Five questions are ready to be closed. Each has a completion test.

1. Place tones under a band limit, a fixed \(N\), and an integer-megahertz synthesizer step so that the earliest leaving noise is as large as possible, and no valley with \(J(\delta)\approx J(\lambda)\) sits at \(\delta\gg\lambda\). Done when a design beats the equal grid’s first threshold and has no such valley. Shuffling gaps at fixed \(S\) does not.
2. Compute the posterior-mean mean square under the same uniform prior as the bound, and compare it with coherent maximum likelihood. Done when the two curves, on that prior, attribute the roughly 5 dB gap at the threshold (−27.6 against −22.8 dB) either to a loose bound or to a suboptimal estimator.
3. Replace the Gaussian-vector \(P_e\) by the wrapped likelihood ratio only on the 10.6–11.3 m cluster, where residuals reach about 2.5 rad. Done when it is clear whether the 1.4 dB gap at \(\sigma = 0.836\,\mathrm{rad}\) moves.
4. Repeat the coherent estimator on the complex observations, and only check the first threshold. The prediction is that the first jump stays and the floor moves.
5. Past about \(1.5\,\sigma_{\mathrm{jump}}\), weight the right-to-left high segments by \(\delta^2\) instead of charging every slip as \(\lambda^2\). Done when that weighting tracks the sample mean square more closely on that side of the threshold, with the remaining gap stated.

Two tasks are not worth the effort. One is a proof that every estimator jumps at \(\sigma_{\mathrm{jump}}\). The incoherent estimator already does not. The other is using the Chinese remainder theorem to postpone the first threshold. That theorem moves the unambiguous range.

## Documents

The derivation, with the cost surfaces and both estimators drawn separately, is the note *多频相位测距的平行门限*. The consolidation of the conclusions, the engineering rules, and the open questions is *多频相位测距的门限*, stored as Markdown at `notes/门限问题的进展.md`. The three curves in the consolidation are the equal-spacing mean square, the cumulative fraction of the half-window bound at \(0.0836\,\mathrm{rad}\), and the geometric-spacing mean square.

The computation directory holds the frozen trials behind those curves and the programs that produced them. A later run is not a substitute for those numbers unless the difference is stated.

## References

Wei Li and Wangdong Qi, Optimization in Multi-Frequency Interferometry Ranging: Theory and Experiment, arXiv:1202.1424, 2012. Used in the two-line heights; no \(\sigma_{\mathrm{jump}}\) as a function of the frequency set.

A. Weiss and E. Weinstein, Fundamental Limitations in Passive Time-Delay Estimation — Part I, IEEE Transactions on Acoustics, Speech, and Signal Processing, 1983; Part II, 1984. Used at the large-\(N\) limit of equal-grid \(R\); different observation model.

J. Ziv and M. Zakai, Some Lower Bounds on Signal Parameter Estimation, IEEE Transactions on Information Theory, 1969. Used for the half-window bound with valley fill.

D. C. Rife and R. R. Boorstyn, Single-Tone Parameter Estimation from Discrete-Time Observations, IEEE Transactions on Information Theory, vol. 20, no. 5, pp. 591–598, 1974. Used to mark the same local-bound / global-jump skeleton for a different parameter.
