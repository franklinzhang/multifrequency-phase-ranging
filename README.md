# Range Estimation from Wrapped Multi-Frequency Phase Observations

[**Read online**](https://franklinzhang.github.io/multifrequency-phase-ranging/)

October 2026

## 1 Observation model

A phase is measured at each of \(N\) frequencies \(f_1,\dots,f_N\), and the range \(d\) is estimated from these phases. The reading at the \(i\)-th frequency is

\[
\varphi_i = \mathrm{wrap}(2\pi f_i d/c + e_i), \qquad e_i \sim \mathcal{N}(0,\sigma^2),
\]

with independent \(e_i\) and \(c = 299\,792\,458\,\mathrm{m/s}\). The map \(\mathrm{wrap}\) reduces an angle modulo \(2\pi\) to \((-\pi,\pi]\).

The true phase \(2\pi f_i d/c\) is a number of whole cycles plus a fraction of a cycle, and the reading keeps only the fraction. When the range changes by \(\delta\), the \(i\)-th reading turns by \(2\pi f_i\delta/c\); if every frequency turns by nearly a whole number of cycles, all the readings stay nearly the same, and \(d\) and \(d+\delta\) are hard to tell apart. Write \(\bar f\) for the mean frequency and \(\lambda = c/\bar f\) for the mean wavelength. At \(\delta = \lambda\) the \(i\)-th reading turns by \(2\pi + 2\pi(f_i-\bar f)/\bar f\), only a small angle more than a whole cycle, and the closer the tone lies to the mean frequency, the smaller that angle; when every \(f_i\delta/c\) is exactly an integer, the readings do not change at all. A finite set of wrapped phases therefore determines the range only up to a set of aliases: approximate aliases near integer multiples of the mean wavelength, and an exact period set by the greatest common divisor of the frequencies.

The coherent estimator minimizes the sum of squared wrapped residuals:

\[
\hat d = \arg\min_{d} \sum_i \mathrm{wrap}\bigl(\varphi_i - 2\pi f_i d/c\bigr)^2 .
\]

When \(\sigma\) is far below \(\pi\), the wrapped normal density on the principal interval differs from the normal density only by negligible terms, and this minimization is the maximum-likelihood estimator. It uses the absolute phase of each tone relative to a single local oscillator, so the carrier enters the range information.

If the tones also share an unknown phase (for example, when transmitter and receiver have no common phase reference), that phase must be eliminated. The incoherent estimator maximizes the modulus of the phasor sum:

\[
\hat d_{\mathrm{inc}} = \arg\max_{d} \Bigl|\sum_i e^{\mathrm{j}(\varphi_i - 2\pi f_i d/c)}\Bigr| .
\]

Multiplying every phasor by the same \(e^{\mathrm{j}\theta}\) does not change the modulus, so the common phase has no effect on it. The carrier term \(2\pi\bar f d/c\) is also the same constant at every tone, cannot be separated from the common phase, and is removed with it; the range information that remains comes only from \(f_i-\bar f\), that is, from the tone spacings.

The sections below work out how the mean-square error of the coherent estimator under this model changes as \(\sigma\) grows. All the tools are standard ones from parameter estimation: wrapped maximum likelihood, Fisher information and the Cramér–Rao bound, the two-hypothesis test and the Ziv–Zakai bound, the mixture approximation for slips, and archived trials on three tone sets.

## 2 The phenomenon: two parallel lines and a jump

Take eight tones from 2400 to 2470 MHz in 10 MHz steps (the equal grid below), repeat trials at a series of \(\sigma\), and record the sample mean-square error of the coherent estimator. The result is Figure 1. The horizontal axis is \(\sigma\) divided by the scale \(\sigma_{\mathrm{jump}} = 0.0836\,\mathrm{rad}\) derived in Section 3. The lower black dashed line is the carrier line, HMSE (high-SNR mean-square error); the upper grey dash-dot line is the bandwidth line, MMSE (moderate-SNR mean-square error, not the minimum mean-square error).

![Figure 1: mean-square error on the equal grid](https://franklinzhang.github.io/multifrequency-phase-ranging/figures/fig02_equal_mse.png)

Figure 1. Mean-square error of the coherent estimator on eight equally spaced tones. Blue circles: trials, 1600 or 8000 per point, random seed 7. Black dashed: carrier line (HMSE); grey dash-dot: bandwidth line (MMSE); green: half-window Ziv–Zakai bound (Section 5); grey dotted: search-window floor; orange vertical line: \(\sigma = \sigma_{\mathrm{jump}}\).

Read the figure from left to right. Up to \(0.22\,\sigma_{\mathrm{jump}}\) the points lie on the carrier line: the whole cycles at every frequency are resolved correctly, and the carrier itself measures the range. As \(\sigma\) grows by less than another factor of three, to about \(0.65\,\sigma_{\mathrm{jump}}\), the points rise almost vertically onto the bandwidth line; from there to \(6\,\sigma_{\mathrm{jump}}\) they stay at most 1.4 dB above it. In this stretch the whole cycles are no longer reliable, and only the tone spacings remain usable. The two lines are parallel and 40.53 dB apart. As \(\sigma\) grows further, the mean square leaves the bandwidth line and rises to the search-window floor, the mean square when the estimate carries no information and the error is uniform over the whole search window. Counted in the direction of growing noise, the jump is the first threshold and this rise is the second, at phase noise of order \(1\,\mathrm{rad}\). The phase model above is the high-SNR limit of complex observations \(y_i = \exp(\mathrm{j}\,2\pi f_i d/c) + n_i\); at phase noise of order \(1\,\mathrm{rad}\) that limit itself no longer holds, and a budget for the second threshold has to return to the complex observations. Only its position is recorded here.

This is the phenomenon to be worked out: the coherent mean square first follows the carrier line, then, within a narrow range of noise, jumps onto a bandwidth line parallel to it. The heights of the two lines are known results. At what \(\sigma\) the jump occurs, which property of the frequency set decides it, and whether tone placement can postpone it are the subject of the sections that follow.

The heights of the two lines come from linear least squares inside the main lobe, the stretch around the true range that has not yet crossed into a neighbouring alias. While the error stays in the main lobe the wrap can be dropped: with phase plotted against frequency, the coherent estimator fits a line through the origin whose slope gives the range, and the incoherent estimator also fits an intercept, the common phase. With \(F = \sum_i f_i^2\) and \(S = \sum_i (f_i-\bar f)^2\), the two mean squares are

\[
\mathrm{HMSE} = \frac{c^2\sigma^2}{4\pi^2 F}, \qquad \mathrm{MMSE} = \frac{c^2\sigma^2}{4\pi^2 S}.
\]

\(F\) comes mostly from the carrier, while \(S\) depends only on how the tones are spread about their mean; hence the names carrier line and bandwidth line. Both are proportional to \(\sigma^2\), and their ratio \(R = F/S\) does not depend on \(\sigma\), so the two lines are parallel on a log-log plot. Fitting an intercept amounts to subtracting the mean frequency first: \(F = S + N\bar f^2\), so \(R = 1 + N\bar f^2/S\). On the equal grid \(\bar f = 2435\,\mathrm{MHz}\), \(S = 4200\,\mathrm{MHz}^2\), and \(R = 11295\), or 40.53 dB; the mean wavelength is \(\lambda = 0.12312\,\mathrm{m}\).

These heights are equations (30) and (42) of Wei Li and Wangdong Qi, obtained by perturbation analysis; the names HMSE and MMSE are also theirs. They call the phenomenon the double-threshold effect and number thresholds in the direction of increasing SNR, so the transition from MMSE down to HMSE is their second threshold, which is the first threshold here. They do not write the noise level of this transition as a function of the frequency set.

Continuous-time passive time-delay estimation has the same three regions. With an improved Ziv–Zakai bound, Weinstein and Weiss found that at high SNR the bound coincides with the Cramér–Rao bound, at moderate SNR it is higher by a factor \(12(\omega_0/W)^2\) (\(\omega_0\) the centre frequency, \(W\) the bandwidth), and at low SNR it tends to a constant set only by the prior interval. For \(N\) equally spaced tones spanning \(B = (N-1)\Delta f\),

\[
R = 1 + \frac{12(N-1)}{N+1}\left(\frac{\bar f}{B}\right)^2 ,
\]

and the coefficient tends to 12 as \(N\) grows, the same as their factor. That is the same limit reached from another observation model, and discrete wrapped phases are not a flat continuous spectrum: at \(N = 8\) the coefficient is 9.33, about 1.1 dB below the continuous-spectrum value at the same \(\bar f/B\).

## 3 The phase scale of the first threshold

### 3.1 The phase scale \(\sigma_{\mathrm{jump}}\)

Two range hypotheses \(\delta\) apart differ in phase by \(2\pi f_i\delta/c\) at the \(i\)-th frequency. Write the squared length of the wrapped phase-difference vector as

\[
J(\delta) = \sum_i \mathrm{wrap}(2\pi f_i\delta/c)^2 .
\]

\(J(0) = 0\). A local minimum of \(J\) at \(\delta \neq 0\) is called a valley; the lower the bottom of a valley, the harder it is to tell apart two ranges \(\delta\) apart. Figure 2 shows \(J(\delta)\) for the equal grid.

![Figure 2: J(δ) for the equal grid](https://franklinzhang.github.io/multifrequency-phase-ranging/figures/fig01_equal_cost.png)

Figure 2. \(J(\delta)\) for the equal grid (in rad²), drawn only where \(J \le 2\); the orange vertical line marks the mean wavelength \(\lambda\). The valleys lie near integer multiples of \(\lambda\), and \(J\) at the bottom of the \(k\)-th valley is about \(k^2 J(\lambda)\). At the mean wavelength the bottom is \(J = 0.028\), almost on zero at this scale.

Take \(\delta = \lambda\). The \(i\)-th component is \(2\pi f_i/\bar f = 2\pi + 2\pi(f_i-\bar f)/\bar f\); once the whole cycle is removed, \(2\pi(f_i-\bar f)/\bar f\) remains, so \(\sqrt{J(\lambda)} = 2\pi\sqrt{S}/\bar f\). Define

\[
\sigma_{\mathrm{jump}} = \pi\sqrt{S}/\bar f ;
\]

then \(\sqrt{J(\lambda)} = 2\sigma_{\mathrm{jump}}\).

Treating the phase-difference vector as a Gaussian observation, the error probability of deciding between the two hypotheses is

\[
P_e(\delta) = Q\bigl(\sqrt{J(\delta)}/(2\sigma)\bigr),
\]

where \(Q\) is the right-tail probability of the standard normal distribution. At \(\delta = \lambda\) the argument is exactly \(\sigma_{\mathrm{jump}}/\sigma\): \(\sigma_{\mathrm{jump}}\) is the noise scale for telling apart two ranges one mean wavelength apart. At \(\sigma = \sigma_{\mathrm{jump}}\) the probability of a one-wavelength error to either side is \(2Q(1) = 0.317\) in total. The valley near the mean wavelength has its bottom at \(\delta = \lambda(1-1/R)\), where \(J = (1-1/R)\,J(\lambda)\), differing from \(J(\lambda)\) only by the fraction \(1/R\).

Since \(F = S + N\bar f^2\), also \(R = 1 + N\pi^2/\sigma_{\mathrm{jump}}^2\): the phase scale of the threshold and the gap between the two lines are set by one quantity, \(S/\bar f^2\).

On the equal grid \(\sigma_{\mathrm{jump}} = 0.0836\,\mathrm{rad}\), about 4.8°. At that noise no single reading comes near \(\pm\pi\); the wrap matters because the range hypothesis has crossed a carrier cycle, not because a reading has been folded back.

This \(P_e\) comes from the Gaussian vector test, not from the likelihood ratio after wrapping. The two are close when the residuals are far below \(\pi\); valleys with residuals near \(\pi\) appear in Section 6.3.

For a round trip the phase is \(4\pi f d/c\), which amounts to replacing \(c\) by \(c/2\): the wavelength and the exact period halve, both mean-square lines divide by 4, and \(\sigma_{\mathrm{jump}}\), which depends only on frequency ratios, keeps its value in radians.

### 3.2 Where the mixture leaves the carrier line

Suppose that each trial slips by one wavelength with probability \(p = 2Q(\sigma_{\mathrm{jump}}/\sigma)\) and otherwise stays in the main lobe. The mixture approximation of the mean square is

\[
\mathrm{MSE} \approx (1-p)\,\mathrm{HMSE} + p\lambda^2 .
\]

Setting \(p\lambda^2 = \mathrm{HMSE}\) gives \(\sigma = 0.0170\,\mathrm{rad}\) on the equal grid, or \(0.204\,\sigma_{\mathrm{jump}}\); this is where the mixture leaves the carrier line. At \(\sigma = \sigma_{\mathrm{jump}}\), \(\mathrm{MMSE} = c^2/(4\bar f^2) = \lambda^2/4\); with \(p = 0.317\), \(p\lambda^2\) is 1.0 dB above the bandwidth line. At the threshold, one-wavelength slips alone are enough to lift the mean square to the height of the bandwidth line.

### 3.3 The equal-grid trials

The coherent estimator takes the global minimum over a search window of half-length 14.990 m centred on the true range, 12.345 m: first the minimum over a grid with step \(\lambda/8\), then iterative refinement around that point. The window length, 29.979 m, equals the equal grid's exact period \(c/\Delta f\). The six noise levels from 0.22 to 0.50 \(\sigma_{\mathrm{jump}}\) have 8000 trials each and the others 1600; the random seed is 7, and the results are in [`computation/data/`](computation/data/). A slip is an error larger than half a mean wavelength in magnitude. Decibel values from here on are relative to 1 m².

At the three noise levels from 0.10 to 0.22 \(\sigma_{\mathrm{jump}}\) there are no slips, and the sample mean square is within 0.35 dB of the carrier line. At 0.25 \(\sigma_{\mathrm{jump}}\) one slip occurs in 8000 trials and the sample mean square rises to −57.2 dB; a single one-wavelength error averaged over 8000 trials contributes exactly −57.2 dB. Where a trial curve starts to rise also depends on the number of trials: one expected slip in 8000 trials needs \(p = 1/8000\), which occurs at \(0.26\,\sigma_{\mathrm{jump}}\). The rise at this point comes from a single large error, not from another asymptote between the two lines.

From 0.25 to 1.0 \(\sigma_{\mathrm{jump}}\) the sample mean square equals the slip rate times \(\lambda^2\) within 0.1 dB: in this stretch nearly all of the mean square comes from one-wavelength slips. At the threshold the slip rate is 0.344 and the sample mean square −22.8 dB, 1.4 dB above the bandwidth line (−24.2 dB); the mixture with \(p = 0.317\) predicts 1.0 dB above, in the same direction. From 1.5 \(\sigma_{\mathrm{jump}}\) on, some slips are longer than one wavelength, and the sample mean square starts to exceed the slip rate times \(\lambda^2\).

## 4 Local information does not jump; the incoherent counterexample

The Fisher information of a wrapped-normal phase about its mean is a smooth function of \(\sigma\). The Cramér–Rao bound for range built from it follows the carrier line and does not jump; on the equal grid it is still within 0.03 dB of the carrier line at \(\sigma = 0.84\,\mathrm{rad}\). Li and Qi likewise find, on complex observations, that HMSE coincides with the Cramér–Rao bound. The jump comes from the global estimator leaving the main lobe. Fisher information describes only the curvature of the likelihood near the true range and does not involve the valleys outside the main lobe; on the jump, that is all it can say.

Maximum-likelihood estimation of a single tone's frequency has the same structure: Rife and Boorstyn showed that at high SNR the variance coincides with the Cramér–Rao bound, and that below a threshold SNR the estimate falls on a wrong spectral peak. There the parameter is frequency, and the observations are not a set of wrapped phases.

The incoherent estimator shows whose jump this is. It uses the same observations at the same tones as the coherent estimator, without the absolute carrier phase. In Figure 3 it lies on the bandwidth line from the lowest noise on: −44.4 dB at 0.10 \(\sigma_{\mathrm{jump}}\), against −44.2 dB for the bandwidth line. The common phase removes the carrier from the start, so there is no lower line for it to leave, and no first threshold. \(\sigma_{\mathrm{jump}}\) is therefore a threshold of estimators that use the absolute carrier phase, not a property of the observations; there is no theorem that every estimator jumps at \(\sigma_{\mathrm{jump}}\).

![Figure 3: coherent and incoherent estimators on the equal grid](https://franklinzhang.github.io/multifrequency-phase-ranging/figures/fig07_coherent_incoherent.png)

Figure 3. Mean-square error of the coherent estimator (blue circles) and the incoherent estimator (orange squares) on the equal grid, against \(\sigma/\sigma_{\mathrm{jump}}\). The incoherent estimator lies on the bandwidth line from the lowest noise on; from about \(2\,\sigma_{\mathrm{jump}}\) the two curves nearly coincide, and at \(10\,\sigma_{\mathrm{jump}}\) they leave the bandwidth line together.

## 5 The half-window Ziv–Zakai bound

The green curve in Figure 1 leaves the carrier line before the trials do. What does it constrain, and what does it not?

Let the true range be uniform on an interval of length \(L\). Ziv and Zakai bounded the estimation error through the error probability of deciding between two hypotheses; Bellini and Tartara added valley filling, and the mean-square bound becomes

\[
\mathrm{ZZB} = \frac{1}{L}\int_0^L \delta\, G\bigl[(L-\delta)\,P_e(\delta)\bigr]\,\mathrm{d}\delta,
\qquad G[h](\delta) = \max_{\eta\ge\delta} h(\eta),
\]

that is, \(G\) replaces a function by its running maximum from the right. The pairwise error probability \(P_e\) is still the Gaussian vector test of Section 3.1. Averaged over this prior, the mean square of any estimator is no lower than ZZB. Here \(L = 14.990\,\mathrm{m}\), half the length of the search window, hence the name half-window bound: no trial error can exceed half the window length, and if \(L\) were the whole period, the replica valley at the far end of the period would raise the bound.

At high SNR only the main lobe contributes: \(\sqrt{J(\delta)} \approx \beta\delta\) with \(\beta = 2\pi\sqrt{F}/c\), and the integral becomes \((2\sigma/\beta)^2\int_0^\infty u\,Q(u)\,\mathrm{d}u\). Since \(\int_0^\infty u\,Q(u)\,\mathrm{d}u = 1/4\), the bound equals \(\sigma^2/\beta^2 = \mathrm{HMSE}\). On the equal grid at \(\sigma = 0.008\,\mathrm{rad}\) the two differ by 0.00 dB, and only range differences up to about 0.5 mm contribute to the integral.

### 5.1 Valley-by-valley decomposition

Valley filling does not add up the error probabilities of individual valleys. By the definition of \(G\), read the weighted error probability \(h(\delta) = (L-\delta)P_e(\delta)\) from right to left: each new running maximum at \(\delta_a\) sets the integrand over the whole stretch to its left, down to the first point \(\delta_b\) where \(h\) is higher. That stretch contributes

\[
\frac{L-\delta_a}{L}\, P_e(\delta_a)\, \frac{\delta_a^2-\delta_b^2}{2}
\]

to the bound. When the left end falls in the main lobe, \(\delta_b \approx \sqrt{J(\delta_a)}/\beta\).

On the equal grid at the threshold \(\sigma = 0.0836\,\mathrm{rad}\), the new running maxima from right to left fall on the \(k\lambda\) valleys. The three largest segments come from the valleys at one, two, and three mean wavelengths and contribute −29.2, −32.9, and −43.0 dB; together they give −27.6 dB, the value of the numerical half-window integral. The left end of the one-wavelength segment falls in the main lobe, at only 1.2 mm. Split by the numerical integral, the one-wavelength valley carries 67% and the two-wavelength valley 29%. Valleys further out at \(k\lambda\) barely enter the integral: \(J\) at the bottom of the \(k\)-th valley grows as \(k^2\), and \(P_e\) falls with it much faster than \(\delta^2\) grows.

### 5.2 What it constrains and what it does not

The bound is for a prior that makes the true range uniform on the half-window; the trials fix the true range at the centre of the search window and give a sample estimate of the conditional mean square. At the threshold the bound is −27.6 dB, the bandwidth line −24.2 dB, and the coherent trials −22.8 dB.

At 0.22 \(\sigma_{\mathrm{jump}}\) the bound is already about 3 dB above the carrier line (−74.8 dB), yet the 8000 trials there have no slip, and the sample mean square, −77.8 dB, lies below the bound. The two do not contradict each other: the bound constrains the expected mean square. The probability of a one-wavelength slip there is about \(5.5\times10^{-6}\), only 0.04 expected in 8000 trials, so most likely none is drawn, and the sample mean square falls below its expectation; by the mixture approximation the expectation is about −70.0 dB, above the bound. That the bound leaves the carrier line before the trials do is partly a matter of the number of trials. The bound constrains the mean square; it does not mark the turning point of the trial curve.

At large noise the bound tends to the variance of the half-window prior, \(L^2/12\), or 12.7 dB; at \(\sigma = 2.34\,\mathrm{rad}\) the equal-grid bound is 8.1 dB, still below that value. When the estimate carries no information and the error is uniform over the whole search window (29.979 m), the mean square is 18.74 dB, the search-window floor of Figure 1. The equal-grid coherent trials are still on the bandwidth line at \(\sigma = 0.50\,\mathrm{rad}\) (\(6\,\sigma_{\mathrm{jump}}\)), about 12 dB above it at \(\sigma = 0.84\,\mathrm{rad}\) (\(10\,\sigma_{\mathrm{jump}}\)), and at the floor by \(\sigma = 2.34\,\mathrm{rad}\). This is the second threshold; it decides when the error reaches the floor and does not change \(\sigma_{\mathrm{jump}}\).

At the threshold the bound and the trials are about 5 dB apart. One is a lower bound averaged over the prior and the other a conditional mean square, so the two cannot simply be subtracted; whether these 5 dB belong to the bound or to the estimator is one of the open questions in Section 8.

## 6 Unequal spacing

What can unequal spacing change? A frequency set has three design quantities, and none of them substitutes for another.

1. \(S/\bar f^2\) sets the phase scale of the first threshold and the gap between the two lines (\(R = 1 + N\pi^2/\sigma_{\mathrm{jump}}^2\)). Spreading the tones away from the mean frequency reduces the height of the jump and moves it to larger phase noise.
2. \(J(\delta)\) at the bottom of each deep valley, together with the offset \(\delta\), decides whether the coherent mean square stays on the bandwidth line, and whether a slip costs one wavelength squared or several metres squared. When \(J(\delta)\) is comparable to \(J(\lambda)\), slipping into that valley is as likely as a one-wavelength error.
3. The greatest common divisor of the tones sets the exact period on which all phases return to whole cycles together, that is, the unambiguous range. A coprime or Chinese-remainder layout changes only this period and does not postpone the first threshold. The approximate periods left by nearly equal gaps do not belong to this item; they appear as deep valleys under the second item (Section 6.2).

### 6.1 Holding \(S\) fixed: moment matching

The second set has the same number of tones, mean frequency, and \(S\) as the equal grid (the moment-matched set below): 2390, 2418, 2425, 2432, 2443, 2449, 2454, and 2469 MHz. The tones are whole megahertz, the greatest common divisor of the gaps is 1 MHz, and the exact period is about 300 m. \(R\) and \(\sigma_{\mathrm{jump}}\) are therefore those of the equal grid.

Within 2 m the two \(J(\delta)\) curves coincide wherever \(J \le 14\): near the valleys in this stretch no residual wraps, so \(J\) depends only on \(N\), \(\bar f\), and \(S\), which the two sets share. At \(\sigma = 0.0836\,\mathrm{rad}\) the half-window integral also splits into the same segments as on the equal grid. Making the gaps unequal while holding \(S\) fixed does not remove the jump; the two sets part only at larger noise (Section 6.3).

### 6.2 Geometric spacing: a far valley as deep as one wavelength

The third set is a geometric progression with ratio 1.02, from 2272 to 2610 MHz, with geometric centre 2435 MHz and arithmetic mean 2437.5 MHz. Its span, about 338 MHz, is 4.8 times that of the equal grid; \(S\) is about 23 times that of the equal grid, almost entirely because of the wider span. \(R\) falls to 26.87 dB, \(\sigma_{\mathrm{jump}}\) rises to 0.403 rad, and the mean wavelength is still 0.1230 m. This is a change in the first design quantity: the threshold moves to larger phase noise, and the two lines move closer together.

The change also falls on the second design quantity. Figure 4 shows \(J(\delta)\) for this set: at the bottom of the valley at 6.211 m, \(J = 0.650\), the same as \(J(\lambda) = 0.650\), and none of the eight phase residuals exceeds 0.44 rad in magnitude. The valley comes from the gaps. The seven gaps lie between 45.4 and 51.2 MHz, with mean 48.26 MHz, and \(c\) divided by 48.26 MHz is 6.212 m; at 6.211 m the phase difference between neighbouring tones is between 0.94 and 1.06 cycles. With steps of 2%, the gaps stay close enough to equal to keep this approximate period. The equal grid's gap is 10 MHz, and its period of this kind lies at 29.98 m, outside the half-window. At 6.211 m, about 50 mean wavelengths, a slip is as likely as a one-wavelength error, but its squared error is about 34 dB larger.

![Figure 4: J(δ) for geometric spacing](https://franklinzhang.github.io/multifrequency-phase-ranging/figures/fig05_geometric_cost.png)

Figure 4. \(J(\delta)\) for the geometric set, drawn only where \(J \le 10\). The orange vertical line marks the mean wavelength and the green vertical line 6.21 m; both valley bottoms are at 0.65. Near 12.4 m there is another cluster of valleys, the deepest at about 2.6.

In the half-window integral this far valley has a smaller \(L-\delta\), so its weighted error probability is lower than at the mean wavelength and does not fill back to the main lobe; the segment it sets runs from the mean wavelength out to 6.211 m and carries 96% of the integral. The contribution is spread along the range: half of the integral comes from beyond 4.3 m and nine tenths from within 5.9 m, rather than from the single point at 6.211 m. Figure 5 compares the cumulative fractions of the three sets at the noise of the equal-grid threshold. At that noise the half-window bound of the geometric set is about −51 dB, and its carrier line −64.8 dB.

![Figure 5: cumulative fraction of the half-window bound](https://franklinzhang.github.io/multifrequency-phase-ranging/figures/fig08_bound_cumulative.png)

Figure 5. Cumulative fraction of the half-window bound along range offset, at the equal-grid threshold noise 0.0836 rad. The equal and moment-matched curves coincide, reaching 96% at 0.25 m and 100% within 0.4 m; the geometric set (green) completes only at 6.21 m.

The trials use the same absolute noise levels as on the equal grid, so \(\sigma = 0.0836\,\mathrm{rad}\) is only 0.21 of the geometric set's own threshold. There the coherent sample mean square is −64.5 dB, still on the carrier line, while the bound is already −51 dB. The reason is the one given in Section 5.2: the probability of slipping to \(\pm 6.211\,\mathrm{m}\) is about \(1.4\times10^{-6}\), and 1600 trials do not draw it. At 0.31 of its own threshold, the coherent and incoherent sample mean squares both reach about −13 dB (−13.2 and −13.1), about 21 dB above this set's bandwidth line (−34.4 dB), while the slip rate is only 0.0031. The same slip rate charged at one wavelength would give about −43 dB; \(0.0031\times(6.21\,\mathrm{m})^2\) is about −9 dB, the same order as −13 dB. After leaving the carrier line, the coherent mean square does not stay on the bandwidth line. The mixture leaving point computed with the mean wavelength, 0.0967 rad, gives the noise at which slips begin, not how far they go. Figure 6 shows these curves.

![Figure 6: mean-square error for geometric spacing](https://franklinzhang.github.io/multifrequency-phase-ranging/figures/fig06_geometric_mse.png)

Figure 6. Mean-square error for the geometric set, against \(\sigma\) normalized by the set's own threshold, 0.403 rad. Black dashed: carrier line; grey dash-dot: bandwidth line; green: half-window bound. The coherent estimator (blue circles) leaves the carrier line and passes straight over the bandwidth line; the incoherent estimator (orange squares) starts on the bandwidth line and leaves it at the same noise.

### 6.3 A cluster of valleys at larger noise

The moment-matched set and the equal grid part only when the phase noise is an order of magnitude larger. At \(\sigma = 0.836\,\mathrm{rad}\) their half-window bounds are −2.4 and −3.8 dB. The moment-matched set has a cluster of valleys from 10.6 to 11.3 m whose bottoms have \(J\) about 14; on the equal grid no valley bottom in the same stretch is below 22.9. Within the cluster the weighted error probability is highest at 10.65 m, and the segment it sets carries 29% of the moment-matched integral on its own; without that segment the bound drops from −2.4 dB to about −3.9 dB, in line with the equal grid's −3.8 dB. The deepest valley of the cluster lies slightly farther out, at 11.0 m; \(L-\delta\) favours the nearer valley, so the new running maximum from the right falls at 10.65 m. The 1.4 dB comes from this far cluster, not from the near valleys.

The bottoms of these valleys have \(J\) about 14, not 0, so they are not exact whole-cycle coincidences either. Their phase residuals reach about 2.5 rad, close to half a cycle, where the Gaussian vector test of Section 3.1 is least reliable. The trials at the same noise are also higher on the moment-matched set: a coherent sample mean square of 10.3 dB, against 7.8 dB on the equal grid.

## 7 Engineering rules

The rules below apply only where the phase model holds, that is, for \(\sigma\) well below 1 rad; at the order of 1 rad the budget moves to complex observations (the second threshold). To convert to SNR, use the high-SNR relation of the complex observations, \(\sigma^2 \approx 1/(2\,\mathrm{SNR})\), where SNR is that of a single tone's complex phasor.

With a common phase reference:

- Budget on the carrier line only for operating points below about \(0.2\,\sigma_{\mathrm{jump}}\). On the equal grid that is a phase standard deviation of about 1°, or 0.017 rad, a per-tone SNR of about 32 dB.
- From about \(0.5\,\sigma_{\mathrm{jump}}\) until the phase model fails, budget on the bandwidth line; on the equal grid it is 40.53 dB worse than the carrier line.
- In between, the mean square is close to the slip rate times \(\lambda^2\) (Section 3.3).

Without a common phase, the high-SNR budget is the bandwidth line from the start, and the jump off the carrier line has no operating region. If a far valley comparable to \(J(\lambda)\) also exists, removing the common phase removes only the carrier line, not the rise that the valley causes; the incoherent estimator on the geometric set shows this (Section 6.2).

When placing tones, solve, for each pair of local maxima \(\pm\delta\) of \(P_e\) (that is, each pair of valleys),

\[
2P_e(\delta)\,\delta^2 = \mathrm{HMSE},
\]

and take the earliest \(\sigma\). The left side is the expected squared error from slipping into that pair of valleys, the factor 2 counting both \(+\delta\) and \(-\delta\); once it exceeds HMSE, the mean square leaves the carrier line. Quoting \(\sigma_{\mathrm{jump}}\) alone is not enough. On the equal grid and the moment-matched set the earliest solution comes from the mean-wavelength valley, \(\sigma = 0.0170\,\mathrm{rad}\) (\(0.204\,\sigma_{\mathrm{jump}}\)), the same as the mixture leaving point of Section 3.2. On the geometric set it comes from the valley at 6.211 m, at about 0.069 rad (0.17 of the set's own threshold), earlier than the 0.0967 rad computed from the mean wavelength; beyond that point, an undetected slip is as likely as a one-wavelength error but costs about 34 dB more in mean square.

When the gaps are nearly equal, check \(J(\delta)\) at \(c\) divided by the typical gap: a valley as deep as \(J(\lambda)\) may remain there, and neither the greatest common divisor nor the Chinese remainder theorem sees it.

| Tone set | \(\bar f\) | \(S\) | \(R\) | \(\sigma_{\mathrm{jump}}\) | Earliest leaving point | Deciding valley |
| --- | --- | --- | --- | --- | --- | --- |
| Equal spacing | 2435 MHz | 4200 MHz² | 40.53 dB | 0.0836 rad | 0.0170 rad | 0.123 m |
| Moment-matched | 2435 MHz | 4200 MHz² | 40.53 dB | 0.0836 rad | 0.0170 rad | 0.123 m |
| Geometric | 2437.5 MHz | \(9.78\times 10^4\) MHz² | 26.87 dB | 0.403 rad | 0.069 rad | 6.211 m |

Table 1. The three tone sets. The earliest leaving point is the earliest solution of \(2P_e(\delta)\,\delta^2 = \mathrm{HMSE}\) over all valleys in the half-window.

Assess the three design quantities separately: widening \(S/\bar f^2\) changes the first threshold and the gap between the two lines; the deep valleys' \(J(\delta)\) and \(\delta\) decide whether the coherent mean square stays on the bandwidth line and how far a slip goes; the greatest common divisor changes the unambiguous range. A longer exact period moves only the third quantity and leaves the first threshold where it was.

## 8 What the existing tools settle, and what remains open

No new estimation theory is proposed here. All the tools are existing ones from parameter estimation, each applied until it has nothing further to say about this problem:

1. Linear least squares inside the main lobe, which is the perturbation analysis of Li and Qi, gives the heights of the two lines; \(R = F/S = 1 + N\pi^2/\sigma_{\mathrm{jump}}^2\) does not depend on \(\sigma\), so the lines are parallel.
2. The Gaussian test between two hypotheses gives the phase scale of the first threshold, \(\sigma_{\mathrm{jump}} = \pi\sqrt{S}/\bar f\), the noise scale for telling apart two ranges one mean wavelength apart.
3. The mixture approximation for slips and the archived trials give the mean square near the threshold: the mixture leaves the carrier line at \(0.204\,\sigma_{\mathrm{jump}}\); the first slip in 8000 trials appears at \(0.25\,\sigma_{\mathrm{jump}}\), consistent with the \(0.26\,\sigma_{\mathrm{jump}}\) set by the number of trials; from there to the threshold the sample mean square equals the slip rate times \(\lambda^2\) within 0.1 dB.
4. The Fisher information of the wrapped normal and the Cramér–Rao bound describe only the curvature of the main lobe and do not jump; the jump is outside what local tools describe.
5. The incoherent estimator is a counterexample on the same observations: it has no first threshold, so the threshold belongs to estimators that use the absolute carrier phase, not to the observations.
6. The Ziv–Zakai bound with the valley filling of Bellini and Tartara gives the half-window bound and its valley-by-valley decomposition. It is a lower bound on the mean square under a uniform prior; it leaves the carrier line before the conditional trials do and does not mark the turning point of the trial curve.
7. Bounds and trials on three tone sets separate the design quantities: unequal spacing at fixed \(S\) does not remove the first threshold; nearly equal gaps leave a far valley as deep as \(J(\lambda)\), which keeps the coherent mean square off the bandwidth line and stretches a slip to metres; tone placement takes the earliest solution of \(2P_e(\delta)\,\delta^2 = \mathrm{HMSE}\).
8. The unambiguous range is set by the greatest common divisor of the frequencies, and a Chinese-remainder layout lengthens it. Number theory has already answered this completely; \(\sigma_{\mathrm{jump}}\) depends only on \(S/\bar f^2\), and lengthening the unambiguous range does not postpone the first threshold.

The mean square of global wrapped least squares at finite \(N\) has no closed form, and none is sought here. The five questions below remain open; each can be answered with the same tools, by computation and trials.

### 8.1 Placing tones against the worst valley

Under constraints on bandwidth, on \(N\), and on an integer-megahertz synthesizer step, can the earliest solution of \(2P_e(\delta)\,\delta^2 = \mathrm{HMSE}\) be made later than the equal grid's 0.0170 rad without leaving a valley with \(J(\delta)\approx J(\lambda)\) at \(\delta \gg \lambda\)? Moment matching shows that shuffling the gaps at fixed \(S\) does nothing; the geometric set pushes the earliest solution to about 0.069 rad through a 4.8-fold span and leaves exactly such a far valley at 6.211 m. \(S\) and the positions of the deep valleys have to be adjusted together. The question is answered by a tone set that meets both conditions under the same constraints, or by a complete search that confirms none exists.

### 8.2 Bayes mean square under one prior

At the threshold the half-window bound is −27.6 dB and the coherent trials −22.8 dB, about 5 dB apart. The bound is for the half-window uniform prior and the trials are a conditional mean square at the centre of the search window, so the two cannot simply be subtracted. Computing the mean square of the posterior mean under the same prior, and comparing it with coherent maximum likelihood on the same noise grid, decides whether the 5 dB is looseness in the bound or suboptimality in the estimator. The same comparison also shows whether integer least squares for the cycle ambiguity, or integration over the ambiguity edges, is worth pursuing.

### 8.3 A wrapped likelihood ratio only on valleys with residuals near half a cycle

The Gaussian vector test is a good approximation when the residuals are far below \(\pi\). At the first threshold \(\sigma = 0.08\,\mathrm{rad}\), and the residuals at the geometric set's 6.211 m valley stay below 0.44 rad; the approximation holds at both. It is not reliable on the cluster from 10.6 to 11.3 m at \(\sigma = 0.836\,\mathrm{rad}\), where residuals reach about 2.5 rad. Replacing the Gaussian vector test by the wrapped likelihood ratio on that cluster alone, recomputing the 1.4 dB gap between the moment-matched and equal-grid bounds, and seeing whether it moves answers the question; error probabilities at other offsets need not be recomputed first.

### 8.4 The first threshold on complex observations

The phase model is the high-SNR limit of the complex observations, and that limit holds at the first threshold. Running coherent maximum likelihood directly on the complex observations of the same equal grid checks whether the first threshold stays near \(\sigma_{\mathrm{jump}}\): the first jump is expected to stay, and the floor and the second threshold to move. After that check, amplitude fading, unequal SNR across tones, and correlated phase noise from a shared local oscillator enter the budget as corrections.

### 8.5 The second moment once slips exceed one wavelength

From \(1.5\,\sigma_{\mathrm{jump}}\) on, the equal-grid sample mean square exceeds the slip rate times \(\lambda^2\), because some slips are longer than one wavelength. The first few equal-grid valleys lie near \(k\lambda\), with \(J \approx k^2 J(\lambda)\) at the bottom and pairwise error probability \(Q(k\sigma_{\mathrm{jump}}/\sigma)\). The valleys at the right-to-left running maxima of the weighted error probability in Section 5.1 can be weighted by \(\delta^2\) and \(Q(k\sigma_{\mathrm{jump}}/\sigma)\) in place of the single \(\lambda^2\); for unequal spacing, use the same running-maximum points rather than integer multiples of the wavelength. The question is answered if, beyond \(1.5\,\sigma_{\mathrm{jump}}\), this weighting tracks the sample mean square clearly better than the \(\lambda^2\) product and locates the stretch where it still does not match.

## References

Wei Li and Wangdong Qi, Optimization in Multi-Frequency Interferometry Ranging: Theory and Experiment, arXiv:1202.1424, 2012. Equations (30) and (42) give HMSE and MMSE; the analysis is by perturbation.

A. J. Weiss and E. Weinstein, Fundamental Limitations in Passive Time-Delay Estimation—Part I: Narrow-Band Systems, IEEE Transactions on Acoustics, Speech, and Signal Processing, vol. 31, no. 2, pp. 472–486, 1983. E. Weinstein and A. J. Weiss, Part II: Wide-Band Systems, same journal, vol. 32, no. 5, pp. 1064–1078, 1984. Part II uses an improved Ziv–Zakai bound and finds high-, moderate-, and low-SNR regions in continuous-time delay estimation, the moderate region a factor \(12(\omega_0/W)^2\) above the Cramér–Rao bound; the observation model differs from the one here.

J. Ziv and M. Zakai, Some Lower Bounds on Signal Parameter Estimation, IEEE Transactions on Information Theory, vol. 15, no. 3, pp. 386–391, 1969. Lower bounds on estimation error from the error probability of binary hypothesis tests.

S. Bellini and G. Tartara, Bounds on Error in Signal Parameter Estimation, IEEE Transactions on Communications, vol. 22, no. 3, pp. 340–342, 1974. The valley-filled form of the Ziv–Zakai bound.

D. C. Rife and R. R. Boorstyn, Single-Tone Parameter Estimation from Discrete-Time Observations, IEEE Transactions on Information Theory, vol. 20, no. 5, pp. 591–598, 1974. The same structure in single-tone frequency estimation: the local bound is tight, and the global mean square leaves it at a threshold; the parameter is frequency, not range.

## License and citation

Text, figures, notebooks, and archived computations are licensed under [CC BY 4.0](LICENSE). Cite as (also in [`CITATION.cff`](CITATION.cff)):

```bibtex
@misc{zhang2026multifrequency,
  author       = {Zhang, Franklin},
  title        = {Range Estimation from Wrapped Multi-Frequency Phase Observations},
  year         = {2026},
  howpublished = {\url{https://franklinzhang.github.io/multifrequency-phase-ranging/}},
  note         = {In Chinese and English; source https://github.com/franklinzhang/multifrequency-phase-ranging}
}
```

A value taken from `computation/data/` should be cited to this repository as an archived trial result.

## Reproduce

```bash
pip install -r requirements.txt
python3 site/build_site.py               # HTML pages and explorer
jupyter notebook notes/threshold.ipynb   # Chinese: notes/多频相位测距的门限.ipynb
```

The numbers above are read from [`computation/data/`](computation/data/) (seed 7); [`computation/README.md`](computation/README.md) describes the scripts that wrote them. A re-run that changes a number should state the difference.
