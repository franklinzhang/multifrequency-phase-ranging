# Range Estimation from Wrapped Multi-Frequency Phase Observations

[**Read online**](https://franklinzhang.github.io/multifrequency-phase-ranging/)

October 2026

In interferometric ranging, stepped-frequency radar, and similar applications, range is estimated from phase readings at a set of discrete frequencies. Monte Carlo trials of the maximum-likelihood estimator show the mean-square range error following a straight line proportional to the square of the phase noise, then jumping onto a higher line parallel to the first. In an ordinary estimation problem the Cramér–Rao bound is usually a single line of that kind.

Four questions follow. Is the jump a lower bound for this problem? What is the closed form of the estimator's mean square before the jump, through the jump, and after it? At what noise level does the jump occur, and by what mechanism? How do the choice of tones and the frequency at which the band starts enter those expressions?

The observation model is a wrapped Gaussian phase. Fisher information inside the main lobe produces two lines: the carrier line when the carrier phase is known, and the bandwidth line when it is not. Each is the Cramér–Rao bound of its model. The sample mean square of the coherent maximum-likelihood estimator follows the carrier line, then leaves it because of an alias one mean wavelength away. The rising stretch has a closed form in the probability of a one-wavelength slip, and after the estimator has landed the mean square coincides with the bandwidth line. The local bound, the carrier line, does not jump. The bound that does rise with the noise is the half-window Ziv–Zakai bound. It leaves the carrier line earlier than the trials, and its graph is not the broken line formed by the two straight lines. The second moment of the tones about their mean sets the threshold scale and the gap between the lines; the locations of the deep valleys set how far a slip goes; with the spacings held fixed, the start of the band enters only through the mean frequency.

Wei Li and Wangdong Qi pointed out this double-threshold behaviour in multi-frequency interferometry and obtained the heights of the two lines by perturbation analysis. Those heights are the carrier line and the bandwidth line below. The definitions here are written directly from the phase model. The relationship between the two accounts is in the references.

## 1 Observation model

A phase is measured at each of \(N\) frequencies \(f_1,\dots,f_N\), and the range \(d\) is estimated from these phases. The reading at the \(i\)-th frequency is

\[
\varphi_i = \mathrm{wrap}(2\pi f_i d/c + e_i), \qquad e_i \sim \mathcal{N}(0,\sigma^2),
\]

with independent \(e_i\) and \(c = 299\,792\,458\,\mathrm{m/s}\). The map \(\mathrm{wrap}\) reduces an angle modulo \(2\pi\) to \((-\pi,\pi]\). Write \(\bar f = N^{-1}\sum_i f_i\) for the mean frequency and \(\lambda = c/\bar f\) for the mean wavelength.

The true phase \(2\pi f_i d/c\) is a number of whole cycles plus a fraction of a cycle, and the reading keeps only the fraction. When the range changes by \(\delta\), the \(i\)-th reading turns by \(2\pi f_i\delta/c\). If every frequency turns by nearly a whole number of cycles, the readings stay nearly the same, and \(d\) and \(d+\delta\) are hard to tell apart. At \(\delta = \lambda\) the \(i\)-th reading turns by \(2\pi + 2\pi(f_i-\bar f)/\bar f\), only a small angle past a whole cycle. When every \(f_i\delta/c\) is an integer, the readings do not change at all. A finite set of wrapped phases therefore determines the range only up to a set of aliases: approximate aliases near integer multiples of the mean wavelength, and an exact period set by the greatest common divisor of the frequencies.

The coherent estimator uses the absolute phase of each tone relative to one local oscillator and minimizes the sum of squared wrapped residuals:

\[
\hat d = \arg\min_{d} \sum_i \mathrm{wrap}\bigl(\varphi_i - 2\pi f_i d/c\bigr)^2 .
\]

The wrapped normal density is the normal density on the principal interval, summed over translates by multiples of \(2\pi\). When \(\sigma\) is far below \(\pi\), the translated terms are negligible on that interval, and this minimization is the maximum-likelihood estimator.

If the tones also share an unknown phase \(\theta\), as when transmitter and receiver have no common phase reference, the reading is \(\mathrm{wrap}(2\pi f_i d/c + \theta + e_i)\). The incoherent estimator removes \(\theta\) by maximizing the modulus of the phasor sum:

\[
\hat d_{\mathrm{inc}} = \arg\max_{d} \Bigl|\sum_i e^{\mathrm{j}(\varphi_i - 2\pi f_i d/c)}\Bigr| .
\]

Multiplying every phasor by \(e^{\mathrm{j}\theta}\) does not change the modulus. The carrier term \(2\pi\bar f d/c\) is the same constant at every tone, cannot be separated from \(\theta\), and is removed with it. The range information that remains is in \(f_i-\bar f\).

For a round trip the phase is \(4\pi f d/c\), which replaces \(c\) by \(c/2\) in the expressions above. The wavelength and the exact period halve, and both mean-square lines below divide by 4. The phase scale of the threshold depends only on frequency ratios and keeps its value in radians.

## 2 Two lines inside the main lobe

The main lobe is the stretch around the true range that has not yet crossed into a neighbouring alias. While the error stays there and the residuals are far below \(\pi\), the wrap can be dropped. Plotted against frequency, the phase is a linear Gaussian model. Both mean squares in this section are Cramér–Rao bounds of that linear model, and both are the variances of the corresponding maximum-likelihood estimators. Each is proportional to \(\sigma^2\). Neither jumps.

### 2.1 The carrier line

When the absolute phase is available,

\[
\varphi_i = \frac{2\pi f_i}{c}\, d + e_i .
\]

Write \(F = \sum_i f_i^2\). Least squares fits a line through the origin. The estimate of the slope \(\alpha = 2\pi d/c\) is \(\hat\alpha = \sum_i f_i\varphi_i/F\), with variance \(\sigma^2/F\). In range,

\[
\mathrm{HMSE}
= \mathrm{Var}(\hat d)
= \frac{c^2\sigma^2}{4\pi^2 F}.
\]

HMSE means high-SNR mean-square error. The line is called the carrier line: most of the information is in the carrier.

The Fisher information of the same model gives the same bound. The derivative of the mean with respect to \(d\) is \(2\pi f_i/c\),

\[
\mathcal{I}(d) = \sum_i \frac{(2\pi f_i/c)^2}{\sigma^2} = \frac{4\pi^2 F}{c^2\sigma^2},
\]

and the Cramér–Rao bound \(\mathcal{I}(d)^{-1}\) is HMSE. Least squares through the origin attains it.

### 2.2 The bandwidth line

When a common phase \(\theta\) is unknown,

\[
\varphi_i = \frac{2\pi f_i}{c}\, d + \theta + e_i .
\]

Write \(S = \sum_i (f_i-\bar f)^2\). After the mean phase is subtracted, the slope depends only on \(f_i-\bar f\) and has variance \(\sigma^2/S\). In range,

\[
\mathrm{MMSE}
= \frac{c^2\sigma^2}{4\pi^2 S}.
\]

MMSE means moderate-SNR mean-square error. The line is called the bandwidth line. It is not a minimum mean-square error. \(S\) describes only how the tones are spread about their mean frequency, not the carrier itself.

The Fisher information is now the \(2\times 2\) matrix of \((d,\theta)\). The two derivatives are \(\partial\mu_i/\partial d = 2\pi f_i/c\) and \(\partial\mu_i/\partial\theta = 1\),

\begin{align*}
\mathcal{I}_{dd} &= \frac{4\pi^2 F}{c^2\sigma^2}, &
\mathcal{I}_{d\theta} &= \frac{2\pi N\bar f}{c\sigma^2}, &
\mathcal{I}_{\theta\theta} &= \frac{N}{\sigma^2}.
\end{align*}

After \(\theta\) is eliminated, the information for \(d\) is

\[
\mathcal{I}_{dd} - \frac{\mathcal{I}_{d\theta}^2}{\mathcal{I}_{\theta\theta}}
= \frac{4\pi^2(F - N\bar f^2)}{c^2\sigma^2}
= \frac{4\pi^2 S}{c^2\sigma^2}.
\]

The Cramér–Rao bound is MMSE. For residuals far below \(\pi\), the incoherent estimator is this same least-squares fit with a free intercept: the modulus of the phasor sum is invariant to a common phase, and maximizing it fits that intercept out. The high-SNR variance of the incoherent estimator is therefore MMSE. There is no carrier line beneath it.

### 2.3 Why the lines are parallel

Since \(F = S + N\bar f^2\),

\[
R = \frac{\mathrm{MMSE}}{\mathrm{HMSE}} = \frac{F}{S} = 1 + \frac{N\bar f^2}{S}.
\]

\(R\) does not depend on \(\sigma\). On a log-log plot of mean square against \(\sigma\), or against SNR, the lines have the same slope and are \(10\log_{10} R\) decibels apart. For the complex phasor of a single tone, \(\sigma^2 \approx 1/(2\,\mathrm{SNR})\) at high SNR, so both lines are also proportional to \(1/\mathrm{SNR}\).

For \(N\) equally spaced tones spanning a bandwidth \(B = (N-1)\Delta f\), \(S = N(N^2-1)\Delta f^2/12\), and

\[
R = 1 + \frac{12(N-1)}{N+1}\left(\frac{\bar f}{B}\right)^2.
\]

The coefficient tends to 12 as \(N\) grows. In continuous-time passive time-delay estimation, Weinstein and Weiss found a moderate-SNR region a factor \(12(\omega_0/W)^2\) above the Cramér–Rao bound, the same limit of a different observation model. Discrete tones are not a flat continuous spectrum: at \(N = 8\) the coefficient is 9.33, about 1.1 dB below the continuous-spectrum value at the same \(\bar f/B\).

## 3 The jump of the coherent estimator

The two lines are local bounds of two different models. The coherent estimator uses the first model, but its global maximum leaves the main lobe at some noise level. Past that point the sample mean square no longer follows the carrier line.

Take eight tones from 2400 to 2470 MHz in steps of 10 MHz, the equal grid below. Here \(\bar f = 2435\,\mathrm{MHz}\), \(S = 4200\,\mathrm{MHz}^2\), \(F = S + N\bar f^2\), \(R = 11295\), or 40.53 dB, and \(\lambda = 0.12312\,\mathrm{m}\). Section 4 gives the phase scale \(\sigma_{\mathrm{jump}} = 0.0836\,\mathrm{rad}\). The horizontal axis of Figure 1 is \(\sigma/\sigma_{\mathrm{jump}}\).

![Figure 1: mean-square error on the equal grid](https://franklinzhang.github.io/multifrequency-phase-ranging/figures/fig02_equal_mse.png)

Figure 1. Mean-square error of the coherent estimator on eight equally spaced tones. Blue circles: trials, 1600 or 8000 per point, random seed 7. Black dashed: carrier line (HMSE). Grey dash-dot: bandwidth line (MMSE). Green: half-window Ziv–Zakai bound (Section 5). Grey dotted: search-window floor. Orange vertical line: \(\sigma = \sigma_{\mathrm{jump}}\).

Read the figure from left to right. Up to \(0.22\,\sigma_{\mathrm{jump}}\) the circles lie on the carrier line, within 0.35 dB: every whole cycle is resolved, and the estimator is using the model of Section 2.1. As \(\sigma\) grows from there to about \(0.65\,\sigma_{\mathrm{jump}}\), by less than a factor of three, the circles rise onto the bandwidth line. From about \(0.65\,\sigma_{\mathrm{jump}}\) to \(6\,\sigma_{\mathrm{jump}}\) they stay at most 1.4 dB above it. Beyond that, the mean square leaves the bandwidth line and rises to the search-window floor, the mean square when the estimate carries no information and the error is uniform over the whole search window.

The jump is the sample mean square of the coherent maximum-likelihood estimator, not a lower bound. The Cramér–Rao bound of Section 2.1 follows the carrier line, and the black dashed line remains that local bound after the jump. The green curve is a lower bound on a different functional. Section 5 says what it constrains.

Counted in the direction of growing noise, the jump is the first threshold and the later departure from the bandwidth line is the second. The second appears at phase noise of order \(1\,\mathrm{rad}\). The phase model of this note is the high-SNR limit of the complex observations \(y_i = \exp(\mathrm{j}\,2\pi f_i d/c) + n_i\). At that order the limit itself no longer holds, and the second threshold has to be computed from the complex observations. Only the first threshold is fixed below.

## 4 The scale, the rise, and the line after the landing

### 4.1 Alias cost and \(\sigma_{\mathrm{jump}}\)

Two ranges \(\delta\) apart differ in phase by \(2\pi f_i\delta/c\) at the \(i\)-th frequency. The squared length of that difference after wrapping is

\[
J(\delta) = \sum_i \mathrm{wrap}(2\pi f_i\delta/c)^2 .
\]

\(J(0) = 0\). A local minimum of \(J\) at \(\delta \ne 0\) is a valley. The lower the bottom, the harder it is to tell the two ranges apart.

![Figure 2: J(δ) for the equal grid](https://franklinzhang.github.io/multifrequency-phase-ranging/figures/fig01_equal_cost.png)

Figure 2. \(J(\delta)\) for the equal grid, in rad², drawn only where \(J \le 2\). The orange vertical line marks the mean wavelength \(\lambda\). The valleys lie near integer multiples of \(\lambda\). At the mean wavelength \(J = 0.028\), almost on zero at this scale.

Take \(\delta = \lambda\). The \(i\)-th component is \(2\pi f_i/\bar f = 2\pi + 2\pi(f_i-\bar f)/\bar f\). If \(|f_i-\bar f| < \bar f/2\), the wrap removes exactly that one cycle and leaves \(2\pi(f_i-\bar f)/\bar f\), so

\[
J(\lambda) = \frac{4\pi^2 S}{\bar f^2}.
\]

Define

\[
\sigma_{\mathrm{jump}} = \frac{\pi\sqrt{S}}{\bar f}.
\]

Then \(\sqrt{J(\lambda)} = 2\sigma_{\mathrm{jump}}\). The equal grid meets the condition above: \(\sigma_{\mathrm{jump}} = 0.0836\,\mathrm{rad}\), about 4.8°, and \(J(\lambda) = 0.0280\).

When the residuals are far below \(\pi\), the maximum-likelihood test between two equally likely hypotheses treats the phase-difference vector as a Gaussian observation. The error probability is

\[
P_e(\delta) = Q\bigl(\sqrt{J(\delta)}/(2\sigma)\bigr),
\]

where \(Q\) is the right tail of the standard normal. At \(\delta = \lambda\) the argument is exactly \(\sigma_{\mathrm{jump}}/\sigma\):

\[
P_e(\lambda) = Q\bigl(\sigma_{\mathrm{jump}}/\sigma\bigr).
\]

\(\sigma_{\mathrm{jump}}\) is the noise scale for telling apart two ranges one mean wavelength apart. At \(\sigma = \sigma_{\mathrm{jump}}\) the probability of a one-wavelength error to either side is \(2Q(1) = 0.317\) in total. At that noise no single reading is near \(\pm\pi\). The wrap matters because the range hypothesis has crossed a carrier cycle, not because a reading has been folded.

This \(P_e\) is the Gaussian vector test, not the likelihood ratio after wrapping. The two are close when the residuals are far below \(\pi\). Valleys whose residuals approach half a cycle are separated out in Section 6.3.

Under the same condition \(|k(f_i-\bar f)/\bar f| < 1/2\), the \(k\)-th mean wavelength has \(k\) whole cycles removed, and the residual is \(k\) times the residual at one wavelength:

\[
J(k\lambda) = k^2 J(\lambda), \qquad P_e(k\lambda) = Q\bigl(k\,\sigma_{\mathrm{jump}}/\sigma\bigr).
\]

The first few valleys of the equal grid agree: \(J\) at the bottom grows as \(k^2\).

From \(R = 1 + N\bar f^2/S\) and the definition of \(\sigma_{\mathrm{jump}}\), the gap between the lines and the scale of the threshold are one quantity,

\[
R = 1 + \frac{N\pi^2}{\sigma_{\mathrm{jump}}^2}.
\]

Substituting \(\sigma = \sigma_{\mathrm{jump}}\) into the bandwidth line cancels \(S\):

\[
\mathrm{MMSE}(\sigma_{\mathrm{jump}}) = \frac{c^2}{4\bar f^2} = \frac{\lambda^2}{4}.
\]

On the equal grid this is \(-24.2\) dB relative to \(1\,\mathrm{m}^2\). Decibel values below are relative to \(1\,\mathrm{m}^2\). The probability of a one-wavelength slip, times \(\lambda^2\), is then the same order as the bandwidth line.

### 4.2 The valley bottom is not exactly one wavelength

While \(\delta\) is near \(\lambda\) and every tone still loses one whole cycle,

\[
J(\delta) = 4\pi^2\sum_i\bigl(f_i\delta/c - 1\bigr)^2.
\]

Set \(\delta = \lambda(1-\varepsilon)\). Differentiating in \(\varepsilon\) puts the minimum at \(\varepsilon = S/F = 1/R\), so

\[
\delta_\star = \lambda\Bigl(1 - \frac{1}{R}\Bigr), \qquad
J(\delta_\star) = \Bigl(1 - \frac{1}{R}\Bigr) J(\lambda).
\]

The cost differs from \(J(\lambda)\) by the fraction \(1/R\). On the equal grid \(R = 11295\), the bottom is at \(0.12310\,\mathrm{m}\) and \(\lambda = 0.12312\,\mathrm{m}\). The difference is not visible in \(P_e\).

### 4.3 The rising stretch

Suppose a trial slips by one mean wavelength with probability

\[
p = 2Q\bigl(\sigma_{\mathrm{jump}}/\sigma\bigr)
\]

and otherwise stays in the main lobe, with variance HMSE. The mixture approximation of the mean square is

\[
\mathrm{MSE} \approx (1-p)\,\mathrm{HMSE} + p\lambda^2.
\]

While \(p\) is small, the second term is the one that first exceeds the first. Setting \(p\lambda^2 = \mathrm{HMSE}\) gives the noise at which the sample mean square leaves the carrier line. On the equal grid the solution is \(\sigma = 0.0170\,\mathrm{rad} = 0.204\,\sigma_{\mathrm{jump}}\). As a per-tone SNR that is about 32 dB, a phase standard deviation of about 1°.

At \(\sigma = \sigma_{\mathrm{jump}}\), \(p = 0.317\) and \(\mathrm{MMSE} = \lambda^2/4\), so \(p\lambda^2\) is \(10\log_{10}(4p) = 1.0\) dB above the bandwidth line. The rise meets the bandwidth line at the threshold scale because of the identity in Section 4.1. There is no third asymptote.

Any other valley is held to the same test. The expected squared error from slipping into the pair of valleys at \(\pm\delta\) is \(2P_e(\delta)\,\delta^2\). Solve

\[
2P_e(\delta)\,\delta^2 = \mathrm{HMSE}(\sigma)
\]

at every valley. The earliest \(\sigma\) is where the mean square leaves the carrier line. Quoting \(\sigma_{\mathrm{jump}}\) alone is not enough. \(\sigma_{\mathrm{jump}}\) is the scale on which one mean wavelength is still distinguishable, and the orange line already sits on the bandwidth line. It is not the point at which the curve leaves the carrier line.

### 4.4 Why the line after the landing is the bandwidth line

In the mixture, \(p\lambda^2\) tends to the constant \(\lambda^2\) as \(\sigma\) grows. On the equal grid \(\lambda^2\) is \(-18.2\) dB. The bandwidth line keeps rising in proportion to \(\sigma^2\). From about \(1.5\,\sigma_{\mathrm{jump}}\) on, \(p\lambda^2\) lies below the bandwidth line and no longer describes the trials.

The model after the landing is Section 2.2. The absolute cycle is no longer reliable, so the carrier phase common to every tone plays the role of an unknown common phase, and the range information that remains is \(S\). The mean square is therefore MMSE, still proportional to \(\sigma^2\), and larger than the carrier line by the factor \(R\). The incoherent estimator removes that common phase from the start, so it has no rising stretch and lies on the bandwidth line at high SNR. The coherent estimator joins it once the absolute cycle is lost.

This description assumes that the valley in the half-window as deep as \(J(\lambda)\) lies near one mean wavelength. If a valley with \(J(\delta)\approx J(\lambda)\) also exists at \(\delta \gg \lambda\), a slip into it is as likely as a one-wavelength error and costs \(\delta^2\). The coherent mean square then passes over the bandwidth line. Section 6.2 is that case.

### 4.5 The equal-grid trials

The coherent estimator takes the global minimum in a search window of half-length 14.990 m centred on the true range 12.345 m: the minimum on a grid of step \(\lambda/8\), then iterative refinement about that point. The window length, 29.979 m, is the equal grid's exact period \(c/\Delta f\). The six noise levels from 0.22 to 0.50 times the threshold have 8000 trials each, and the others 1600. The random seed is 7. A slip is an error larger than half a mean wavelength in magnitude.

From 0.10 to 0.22 times the threshold there are no slips, and the sample mean square is within 0.35 dB of the carrier line. At 0.25 times the threshold, one slip occurs in 8000 trials and the sample mean square is \(-57.2\) dB. One copy of \(\lambda^2\), averaged over 8000 trials, is exactly \(-57.2\) dB. One expected slip in 8000 trials needs \(p = 1/8000\), which is \(0.26\,\sigma_{\mathrm{jump}}\). The rise here is a single large error.

From 0.25 to 1.0 times the threshold, the sample mean square equals the slip rate times \(\lambda^2\), within 0.1 dB. At the threshold the slip rate is 0.344 and the sample mean square is \(-22.8\) dB, 1.4 dB above the bandwidth line at \(-24.2\) dB. The mixture with \(p = 0.317\) predicts 1.0 dB above. From 1.5 times the threshold on, some slips are longer than one wavelength, the sample mean square exceeds the slip rate times \(\lambda^2\), and it moves onto the bandwidth line. At 6 times the threshold (\(\sigma = 0.50\,\mathrm{rad}\)) it is only 0.1 dB above the bandwidth line, while \(p\lambda^2\) is already about 10 dB lower. At 10 times the threshold it is about 12 dB above the bandwidth line, and at \(\sigma = 2.34\,\mathrm{rad}\) it reaches the search-window floor, 18.7 dB.

## 5 The lower bound does not follow the broken line

### 5.1 The local bound does not jump

The Fisher information of a wrapped-normal phase about its mean is a smooth function of \(\sigma\). The Cramér–Rao bound for range built from it follows the carrier line. On the equal grid it is still within 0.03 dB of the carrier line at \(\sigma = 0.84\,\mathrm{rad}\): there \(\pi/\sigma \approx 3.7\), so little probability mass lies outside \(\pm\pi\), and the wrap barely changes the local curvature.

A local bound describes the curvature of the likelihood near the true range. It does not contain the valleys outside the main lobe. The jump of the coherent estimator comes from the global maximum leaving the main lobe. There is no theorem that the Cramér–Rao bound breaks from HMSE to MMSE at \(\sigma_{\mathrm{jump}}\). Maximum-likelihood estimation of a single tone's frequency has the same structure. Rife and Boorstyn showed that at high SNR the variance coincides with the Cramér–Rao bound, and that below a threshold the estimate falls on a wrong spectral peak. There the parameter is frequency.

### 5.2 The half-window Ziv–Zakai bound

The bound that rises with the noise is a Bayesian one. Let the true range be uniform on an interval of length \(L\). Ziv and Zakai bounded the estimation error by the error probability of a binary test. Bellini and Tartara added valley filling:

\[
\mathrm{ZZB}
= \frac{1}{L}\int_0^L \delta\, G\bigl[(L-\delta)\,P_e(\delta)\bigr]\,\mathrm{d}\delta,
\qquad
G[h](\delta) = \max_{\eta\ge\delta} h(\eta).
\]

\(G\) replaces the weighted error probability by its running maximum from the right. \(P_e\) is still the Gaussian test of Section 4.1. Averaged over this prior, the mean square of any estimator is at least ZZB. Here \(L = 14.990\,\mathrm{m}\), half the search window, hence the name. No trial error can exceed half the window. If \(L\) were the whole period, the replica valley at the far end of the period would raise the bound.

At high SNR only the main lobe contributes. \(\sqrt{J(\delta)}\approx\beta\delta\) with \(\beta = 2\pi\sqrt{F}/c\), \(L\) is far larger than any \(\delta\) that contributes, and valley filling does not yet reach an alias, so

\[
\mathrm{ZZB} \approx \int_0^\infty \delta\, Q\bigl(\beta\delta/(2\sigma)\bigr)\,\mathrm{d}\delta.
\]

The substitution \(u = \beta\delta/(2\sigma)\) turns the integral into \((2\sigma/\beta)^2\int_0^\infty u\,Q(u)\,\mathrm{d}u\). Changing the order of integration,

\[
\int_0^\infty u\,Q(u)\,\mathrm{d}u
= \int_0^\infty \frac{t^2}{2}\,\phi(t)\,\mathrm{d}t
= \frac14,
\]

where \(\phi\) is the standard normal density, and the last step uses \(\int_0^\infty t^2\phi(t)\,\mathrm{d}t = 1/2\). Multiplying by \(1/4\) gives \(\mathrm{ZZB} = \sigma^2/\beta^2 = \mathrm{HMSE}\). On the equal grid at \(\sigma = 0.008\,\mathrm{rad}\) the two differ by 0.00 dB.

Read \(h(\delta) = (L-\delta)P_e(\delta)\) from right to left. Each new running maximum at \(\delta_a\) sets the stretch to its left down to the next higher point \(\delta_b\), and that stretch contributes

\[
\frac{L-\delta_a}{L}\, P_e(\delta_a)\, \frac{\delta_a^2-\delta_b^2}{2}.
\]

On the equal grid at \(\sigma = \sigma_{\mathrm{jump}}\), the three largest contributions come from one, two, and three mean wavelengths: \(-29.2\), \(-32.9\), and \(-43.0\) dB. Together they give \(-27.6\) dB, the value of the numerical integral. Split by the numerical integral, one wavelength carries 67% and two wavelengths 29%. Valleys farther out barely enter. \(J\) grows as \(k^2\), and \(P_e\) falls faster than \(\delta^2\) grows.

### 5.3 It is not the broken line in Figure 1

The half-window bound constrains the Bayes mean square under a uniform prior. The trials fix the true range at the centre of the search window and report a sample of the conditional mean square. On the equal-grid archive the two sit as follows.

| \(\sigma/\sigma_{\mathrm{jump}}\) | Half-window ZZB | Carrier line | Bandwidth line | Coherent trials |
| --- | --- | --- | --- | --- |
| 0.10 | −84.7 | −84.7 | −44.2 | −84.4 |
| 0.22 | −74.8 | −77.9 | −37.4 | −77.8 |
| 1 | −27.6 | −64.7 | −24.2 | −22.8 |
| 6 | −9.4 | −49.2 | −8.7 | −8.5 |
| 10 | −3.8 | −44.7 | −4.2 | +7.8 |

At 0.22 times the threshold the bound is already about 3 dB above the carrier line, while 8000 trials are still on the carrier line. There \(p \approx 5.5\times 10^{-6}\), only 0.04 expected slips in 8000 trials, so the sample usually draws none. The mixture puts the expectation near \(-70.0\) dB, still above the bound. At the threshold the bound is 3.3 dB below the bandwidth line and 4.8 dB below the trials. At 6 times the threshold it is still 0.8 dB below the bandwidth line. At 10 times the threshold it leaves the bandwidth line itself and rises toward the variance \(L^2/12 = 12.7\) dB of the half-window prior. The mean square of an error uniform on the whole search window is 18.7 dB.

The bandwidth line is therefore not the middle piece of this bound, and the green curve in Figure 1 is not a lower bound on the blue points. Before the jump the blue points equal the carrier line, through the rise they equal \(p\lambda^2\), and after the landing they equal the bandwidth line. The green curve leaves the carrier line earlier and more gradually, tracks a few decibels below the bandwidth line, and then rises toward the prior floor. How much of the 5 dB at the threshold belongs to looseness of the bound, and how much to the estimator, is still open in Section 8. The Gaussian test ceases to be the likelihood ratio when residuals approach \(\pi\), and the guarantee on this bound stops at valleys whose residuals are far below \(\pi\).

### 5.4 The incoherent estimator has no such jump

The incoherent estimator uses the same observations at the same tones as the coherent estimator, without the absolute carrier phase. In Figure 3 it lies on the bandwidth line from the lowest noise on: \(-44.4\) dB at 0.10 times the threshold, against \(-44.2\) dB for the bandwidth line. From about 1.5 times the threshold the two curves coincide, and at 10 times the threshold they leave the bandwidth line together.

![Figure 3: coherent and incoherent estimators](https://franklinzhang.github.io/multifrequency-phase-ranging/figures/fig07_coherent_incoherent.png)

Figure 3. Mean-square error of the coherent estimator (blue circles) and the incoherent estimator (orange squares) on the equal grid, against \(\sigma/\sigma_{\mathrm{jump}}\).

\(\sigma_{\mathrm{jump}}\) is therefore a threshold of estimators that use the absolute carrier phase. The observations do not force every estimator to jump there.

## 6 Tones and the start of the band

A frequency set has three quantities, and none of them substitutes for another.

1. \(S/\bar f^2\) sets \(\sigma_{\mathrm{jump}}\) and \(R\). It contains both how widely the tones are spread and where the band is placed.
2. \(J(\delta)\) at each deep valley, together with \(\delta\), decides whether the coherent mean square stays on the bandwidth line, and whether a slip costs \(\lambda^2\) or \(\delta^2\). When \(J(\delta)\) is comparable to \(J(\lambda)\), slipping into that valley is as likely as a one-wavelength error.
3. The greatest common divisor of the tones sets the exact period on which every phase returns to a whole cycle, the unambiguous range. A coprime or Chinese-remainder layout changes only this period. The approximate periods left by nearly equal gaps do not belong to this item. They are deep valleys under the second item.

Three ways of changing the set separate cleanly.

Shift the whole band by \(\Delta\) with the spacings fixed. Then \(S\) is unchanged and \(\bar f\) becomes \(\bar f+\Delta\), so

\[
\sigma_{\mathrm{jump}}' = \frac{\pi\sqrt{S}}{\bar f+\Delta}, \qquad
R' = 1 + \frac{N(\bar f+\Delta)^2}{S}.
\]

The bandwidth line does not move. The carrier line falls through \(F' = S + N(\bar f+\Delta)^2\). Raising the frequency at which the band starts moves the jump to smaller phase noise and widens the gap between the lines.

Hold the mean frequency fixed and spread the tones. Then \(S\) grows, \(\sigma_{\mathrm{jump}}\) grows, \(R\) falls, and the bandwidth line drops. When \(N\bar f^2 \gg S\) the carrier line barely moves. The jump becomes shorter and moves to larger phase noise.

Multiply every frequency by the same factor \(\alpha\). Then \(S\) and \(F\) both scale by \(\alpha^2\), \(\sigma_{\mathrm{jump}}\) and \(R\) are unchanged, and both lines divide by \(\alpha^2\).

### 6.1 Shuffling the gaps only

The second set has the same number of tones, the same mean frequency, and the same \(S\) as the equal grid. Call it moment-matched: 2390, 2418, 2425, 2432, 2443, 2449, 2454, and 2469 MHz. The tones are whole megahertz, the greatest common divisor of the gaps is 1 MHz, and the exact period is about 300 m. \(R\) and \(\sigma_{\mathrm{jump}}\) are those of the equal grid.

Within 2 m, the two curves \(J(\delta)\) coincide wherever \(J \le 14\). Near those valleys no residual wraps, so \(J\) depends only on \(N\), \(\bar f\), and \(S\). At \(\sigma = \sigma_{\mathrm{jump}}\) the half-window integral splits into the same segments as on the equal grid. Holding \(S\) fixed and only making the gaps unequal does not remove the first threshold.

### 6.2 A wider span that leaves a far valley

The third set is a geometric progression with ratio 1.02, from 2272 to 2610 MHz, geometric centre 2435 MHz, arithmetic mean 2437.5 MHz. Its span, about 338 MHz, is 4.8 times the equal grid, and \(S\) is about 23 times larger. \(R\) falls to 26.87 dB, \(\sigma_{\mathrm{jump}}\) rises to 0.403 rad, and the mean wavelength is still 0.1230 m. This is a change in the first quantity: the threshold moves to larger phase noise, and the two lines move closer.

The second quantity changes at the same time. In Figure 4 the valley at 6.211 m has \(J = 0.650\), the same as \(J(\lambda) = 0.650\), and none of the eight phase residuals exceeds 0.44 rad in magnitude. The seven gaps lie between 45.4 and 51.2 MHz, with mean 48.26 MHz, and \(c\) divided by 48.26 MHz is 6.212 m. Steps of 2% do not break this approximate period. At 6.211 m, about 50 mean wavelengths, a slip is as likely as a one-wavelength error and costs about 34 dB more in squared error. The equal grid's 10 MHz gap puts the period of the same kind at 29.98 m, outside the half-window.

![Figure 4: J(δ) for geometric spacing](https://franklinzhang.github.io/multifrequency-phase-ranging/figures/fig05_geometric_cost.png)

Figure 4. \(J(\delta)\) for the geometric set, drawn only where \(J \le 10\). The orange line marks the mean wavelength and the green line 6.21 m. Both valley bottoms are at 0.65.

Solving \(2P_e(\delta)\,\delta^2 = \mathrm{HMSE}\) at every valley in the half-window, the earliest solution comes from 6.211 m, at about 0.069 rad, 0.17 of the set's own threshold. The solution computed from the mean wavelength, 0.0967 rad, is later. In the half-window integral the segment set by this far valley runs from the mean wavelength out to 6.211 m and carries 96% of the integral.

The trials use the same absolute noise levels as the equal grid, so \(\sigma = 0.0836\,\mathrm{rad}\) is only 0.21 of the geometric set's own threshold. There the coherent sample mean square is \(-64.5\) dB, still on the carrier line (\(-64.8\) dB), while the bound is already \(-51\) dB. The probability of slipping to \(\pm 6.211\,\mathrm{m}\) is about \(1.4\times 10^{-6}\), and 1600 trials do not draw it. At 0.31 of its own threshold the coherent and incoherent sample mean squares are \(-13.2\) and \(-13.1\) dB, about 21 dB above this set's bandwidth line (\(-34.4\) dB), while the slip rate is only 0.0031. The same slip rate charged at one wavelength would give about \(-43\) dB. The value \(0.0031\times(6.21\,\mathrm{m})^2\) is about \(-9\) dB, the same order as \(-13\) dB. After leaving the carrier line, the coherent mean square does not stay on the bandwidth line.

![Figure 5: cumulative fraction of the half-window bound](https://franklinzhang.github.io/multifrequency-phase-ranging/figures/fig08_bound_cumulative.png)

Figure 5. Cumulative fraction of the half-window bound along range offset, at the equal-grid threshold noise 0.0836 rad. The equal and moment-matched curves coincide and reach 96% by 0.25 m. The geometric set (green) completes only at 6.21 m.

![Figure 6: mean-square error for geometric spacing](https://franklinzhang.github.io/multifrequency-phase-ranging/figures/fig06_geometric_mse.png)

Figure 6. Mean-square error for the geometric set, against \(\sigma\) normalized by the set's own threshold, 0.403 rad. After leaving the carrier line, the coherent estimator passes straight over the bandwidth line. The incoherent estimator starts on the bandwidth line and leaves it at the same noise.

### 6.3 Valleys that appear only at larger noise

The moment-matched set and the equal grid part only when the phase noise is an order of magnitude larger. At \(\sigma = 0.836\,\mathrm{rad}\) their half-window bounds are \(-2.4\) and \(-3.8\) dB. The moment-matched set has a cluster of valleys from 10.6 to 11.3 m whose bottoms have \(J\) about 14. On the equal grid no valley bottom in the same stretch is below 22.9. The weighted error probability is highest at 10.65 m, and the segment it sets carries 29% of the moment-matched integral. Without that segment the bound drops from \(-2.4\) dB to about \(-3.9\) dB. Phase residuals in these valleys reach about 2.5 rad, where the Gaussian test of Section 4.1 is least reliable. At the same noise the coherent sample mean square is 10.3 dB on the moment-matched set and 7.8 dB on the equal grid.

## 7 The closed forms, and how to use them

Where the phase model holds, that is for \(\sigma\) well below 1 rad, the mean square of the coherent estimator has three stretches.

\[
\begin{aligned}
\sigma &\le \sigma_{\mathrm{leave}}: & \mathrm{MSE} &= \mathrm{HMSE} = \dfrac{c^2\sigma^2}{4\pi^2 F}, \\[4pt]
\sigma_{\mathrm{leave}} \lesssim \sigma \lesssim \sigma_{\mathrm{jump}}: &
\mathrm{MSE} &= p\lambda^2, \quad p = 2Q(\sigma_{\mathrm{jump}}/\sigma), \\[4pt]
\sigma_{\mathrm{jump}} \lesssim \sigma \ll 1: &
\mathrm{MSE} &= \mathrm{MMSE} = \dfrac{c^2\sigma^2}{4\pi^2 S}.
\end{aligned}
\]

Here

\[
\sigma_{\mathrm{jump}} = \frac{\pi\sqrt{S}}{\bar f}, \qquad
R = 1 + \frac{N\pi^2}{\sigma_{\mathrm{jump}}^2}, \qquad
\mathrm{MMSE}(\sigma_{\mathrm{jump}}) = \frac{\lambda^2}{4},
\]

and \(\sigma_{\mathrm{leave}}\) is the earliest solution of \(2P_e(\delta)\,\delta^2 = \mathrm{HMSE}\) over the valleys in the search range. When there is no far valley, that solution comes from the mean wavelength, \(p\lambda^2 = \mathrm{HMSE}\). The third stretch assumes that the nearest deep valley lies near \(\lambda\). Otherwise the third stretch is the mean square lifted by that valley's \(\delta^2\), and the bandwidth line is not where the estimator lands.

The incoherent estimator has neither the first stretch nor the rise. While the phase model holds, \(\mathrm{MSE} = \mathrm{MMSE}\).

The local bounds do not take part in the jump. For the unwrapped Gaussian model the Cramér–Rao bound equals HMSE when the absolute phase is known and MMSE when a common phase is unknown. The range bound of the wrapped normal moves smoothly along HMSE. The half-window Ziv–Zakai bound equals HMSE at high SNR, then leaves the carrier line earlier, and does not coincide with the broken line above.

With a common phase reference and no far valley, budget the carrier line only below about \(0.2\,\sigma_{\mathrm{jump}}\), the bandwidth line from about \(0.5\,\sigma_{\mathrm{jump}}\) until the phase model fails, and \(p\lambda^2\) in between. Without a common phase reference, the budget starts on the bandwidth line. Convert to SNR by \(\sigma^2 \approx 1/(2\,\mathrm{SNR})\).

| Tone set | \(\bar f\) | \(S\) | \(R\) | \(\sigma_{\mathrm{jump}}\) | Earliest leaving point | Deciding valley |
| --- | --- | --- | --- | --- | --- | --- |
| Equal spacing | 2435 MHz | 4200 MHz² | 40.53 dB | 0.0836 rad | 0.0170 rad | 0.123 m |
| Moment-matched | 2435 MHz | 4200 MHz² | 40.53 dB | 0.0836 rad | 0.0170 rad | 0.123 m |
| Geometric | 2437.5 MHz | \(9.78\times 10^4\) MHz² | 26.87 dB | 0.403 rad | 0.069 rad | 6.211 m |

Table 1. The three tone sets. The earliest leaving point is the earliest solution of \(2P_e(\delta)\,\delta^2 = \mathrm{HMSE}\) over all valleys in the half-window.

When the gaps are nearly equal, check \(J(\delta)\) at \(c\) divided by the typical gap. Neither the greatest common divisor nor the Chinese remainder theorem sees that valley. A longer exact period changes the unambiguous range and does not change \(\sigma_{\mathrm{jump}}\).

## 8 Calculations still open

At finite \(N\), the mean square of global wrapped least squares has no closed form tighter than the three stretches of Section 7. The items below are numerical or experimental. The tools are \(J\), \(P_e\), and the half-window integral already written down.

Under constraints on bandwidth, on \(N\), and on an integer-megahertz step, can the earliest leaving point be made later than the equal grid's 0.0170 rad without leaving a valley with \(J(\delta)\approx J(\lambda)\) at \(\delta\gg\lambda\)? Moment matching shows that shuffling the gaps at fixed \(S\) does nothing. The geometric set pushes the earliest solution to about 0.069 rad by a 4.8-fold span and leaves exactly the far valley at 6.211 m. \(S\) and the positions of the deep valleys have to be adjusted together.

At the threshold the half-window bound is \(-27.6\) dB and the coherent trials are \(-22.8\) dB, under different priors. The mean square of the posterior mean under the same half-window uniform prior, compared with coherent maximum likelihood, decides whether the 5 dB belongs to the bound or to the estimator.

At \(\sigma = 0.836\,\mathrm{rad}\), residuals in the cluster from 10.6 to 11.3 m are about 2.5 rad, and the Gaussian test is not reliable. Replacing it by the wrapped likelihood ratio on that cluster alone shows whether the 1.4 dB gap between the moment-matched and equal-grid bounds moves.

At the first threshold the phase model, as the high-SNR limit of the complex observations, still holds. Coherent maximum likelihood on the complex observations of the same equal grid checks whether the first jump stays near \(\sigma_{\mathrm{jump}}\). The first jump is expected to stay. The floor and the second threshold are expected to move.

From \(1.5\,\sigma_{\mathrm{jump}}\) on, the equal-grid sample mean square exceeds the slip rate times \(\lambda^2\), because some slips are longer than one wavelength. The right-to-left running maxima of Section 5.2 carry weights \(Q(k\sigma_{\mathrm{jump}}/\sigma)\) and \(\delta^2\). Replacing the single \(\lambda^2\) by those valleys tests whether, past 1.5 times the threshold, the description tracks the sample mean square better than \(p\lambda^2\). The stretch after the landing is already described by MMSE. This correction is for the tail of the rise, not for the bandwidth line.

## References

Wei Li and Wangdong Qi, Optimization in Multi-Frequency Interferometry Ranging: Theory and Experiment, arXiv:1202.1424, 2012. They pointed out the double-threshold behaviour of least squares in multi-frequency interferometry and, by perturbation analysis, obtained closed forms for the mean square at high and moderate SNR, the high-SNR form coinciding with the Cramér–Rao bound of the complex observations. The carrier line and the bandwidth line of this note are those two expressions. Numbering thresholds in the direction of increasing SNR, they call the transition from the bandwidth line down to the carrier line the second threshold, which is the first threshold here. They do not write the phase noise of that transition as a function of \(S\) and \(\bar f\).

A. J. Weiss and E. Weinstein, Fundamental Limitations in Passive Time-Delay Estimation—Part I: Narrow-Band Systems, IEEE Transactions on Acoustics, Speech, and Signal Processing, vol. 31, no. 2, pp. 472–486, 1983. E. Weinstein and A. J. Weiss, Part II: Wide-Band Systems, same journal, vol. 32, no. 5, pp. 1064–1078, 1984. Part II uses an improved Ziv–Zakai bound and finds high-, moderate-, and low-SNR regions in continuous-time delay estimation, the moderate region a factor \(12(\omega_0/W)^2\) above the Cramér–Rao bound. The observation model differs from the one here. For large \(N\), the equal-spacing ratio in Section 2.3 meets that factor.

J. Ziv and M. Zakai, Some Lower Bounds on Signal Parameter Estimation, IEEE Transactions on Information Theory, vol. 15, no. 3, pp. 386–391, 1969.

S. Bellini and G. Tartara, Bounds on Error in Signal Parameter Estimation, IEEE Transactions on Communications, vol. 22, no. 3, pp. 340–342, 1974. Section 5.2 uses this valley-filled form.

D. C. Rife and R. R. Boorstyn, Single-Tone Parameter Estimation from Discrete-Time Observations, IEEE Transactions on Information Theory, vol. 20, no. 5, pp. 591–598, 1974. The same structure, a tight local bound and a global mean square that leaves it at a threshold. The parameter is frequency, not range.

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
