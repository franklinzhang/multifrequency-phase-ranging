# Range Estimation from Wrapped Multi-Frequency Phase Observations

[**Read online**](https://franklinzhang.github.io/multifrequency-phase-ranging/) (中文 / English)

## Abstract

Phases are measured at a finite set of frequencies, each known only modulo \(2\pi\), and range is estimated from these wrapped phases. As the phase noise \(\sigma\) grows, the mean-square error of the coherent maximum-likelihood estimator first follows the carrier line HMSE (the accuracy when every whole cycle is resolved) and then rises almost vertically onto the bandwidth line MMSE (the accuracy that the tone spacings alone provide); on a log-log plot the two lines have the same slope. The phase scale of this jump is fixed by the frequency set:

\[
\sigma_{\mathrm{jump}} = \pi\sqrt{S}/\bar f,
\]

where \(\bar f\) is the mean frequency and \(S = \sum_i (f_i-\bar f)^2\) is the spread of the tones about it. For eight equally spaced tones from 2400 to 2470 MHz in 10 MHz steps, \(\sigma_{\mathrm{jump}} = 0.0836\,\mathrm{rad}\) (4.8°). This is the first threshold met as the noise grows; a second appears at phase noise of order 1 rad, where the phase model itself fails. With mean wavelength \(\lambda = c/\bar f\) and \(p\) the probability of slipping by one wavelength, the mixture approximation \((1-p)\,\mathrm{HMSE} + p\lambda^2\) leaves the carrier line at \(0.204\,\sigma_{\mathrm{jump}}\). With 8000 trials per noise level between 0.22 and 0.50 \(\sigma_{\mathrm{jump}}\), the first slip (an error beyond half a mean wavelength) appears at \(0.25\,\sigma_{\mathrm{jump}}\); from there to \(1.0\,\sigma_{\mathrm{jump}}\) the sample mean square equals the slip rate times \(\lambda^2\), within 0.1 dB. At the threshold the slip rate is 0.344, the sample mean square −22.8 dB re 1 m², and the bandwidth line −24.2 dB.

With the true range uniform on an interval 14.99 m long, half the length of the search window, the Ziv–Zakai bound (the half-window bound) at the same threshold is −27.6 dB. It is a Bayes risk under that prior, whereas the trials are conditional on a fixed range at the centre of the search window; the bound leaves the carrier line before the trials do. The Cramér–Rao bound built from local Fisher information does not jump. An incoherent estimator, which does not use the absolute carrier phase, lies on the bandwidth line from the lowest noise tried, so \(\sigma_{\mathrm{jump}}\) is a threshold only for estimators that use the absolute carrier phase.

Unequal spacing does not remove the jump. Holding \(N\), \(\bar f\), and \(S\) fixed (moment matching) and only breaking the equal spacing leaves the first threshold where it was. A geometric progression with ratio 1.02 enlarges \(S\) about 23-fold and raises \(\sigma_{\mathrm{jump}}\) to 0.403 rad; but on that set two ranges 6.211 m apart are as hard to tell apart as two ranges one mean wavelength apart, so after leaving the carrier line the mean square does not stay on the bandwidth line, and a slip spans about 50 mean wavelengths. At ten times the equal-grid threshold, ranges 10.6 to 11.3 m apart are harder to tell apart on the moment-matched set than on the equal grid, which raises its half-window bound about 1.4 dB above the equal-grid one.

For budgets: with a common phase reference, use the carrier line only below about \(0.2\,\sigma_{\mathrm{jump}}\), and the bandwidth line from about \(0.5\,\sigma_{\mathrm{jump}}\) until the phase model fails; without one, use the bandwidth line from the start. When placing tones, let \(P_e(\delta)\) be the probability of confusing two ranges \(\delta\) apart; at every pair of local maxima \(\pm\delta\) of \(P_e\), solve \(2P_e(\delta)\,\delta^2 = \mathrm{HMSE}\) for \(\sigma\), and keep the earliest. A longer unambiguous range does not postpone the first threshold.

Five questions remain open: placing tones on an integer-megahertz grid so that the earliest of these solutions comes later than on the equal grid; the posterior mean against coherent maximum likelihood under the bound's prior, to assign the 5 dB gap at the threshold to the bound or to the estimator; a wrapped likelihood ratio on the 10.6–11.3 m offsets only; a check of the first threshold on complex observations; and \(\delta^2\) weighting once slips exceed one wavelength.

## 1 Observation and two estimators

The observation at frequency \(f_i\) is

\[
\varphi_i = \mathrm{wrap}(2\pi f_i d/c + e_i), \qquad e_i \sim \mathcal{N}(0,\sigma^2),
\]

with independent \(e_i\) and \(c = 299\,792\,458\,\mathrm{m/s}\). The map \(\mathrm{wrap}\) reduces an angle modulo \(2\pi\) to \((-\pi,\pi]\): the true phase may differ from the reading by any multiple of \(2\pi\), and only the principal value is kept. A finite set of such phases therefore fixes the range only up to a set of aliases.

The coherent estimator minimizes the sum of squared wrapped residuals:

\[
\hat d = \arg\min_{d} \sum_i \mathrm{wrap}\bigl(\varphi_i - 2\pi f_i d/c\bigr)^2 .
\]

It uses absolute phase relative to a single local oscillator, so the carrier enters the range information. If all tones also share an unknown phase offset, that offset must be eliminated, and the estimator maximizes the modulus of the phasor sum instead:

\[
\hat d_{\mathrm{inc}} = \arg\max_{d} \Bigl|\sum_i e^{\mathrm{j}(\varphi_i - 2\pi f_i d/c)}\Bigr| ;
\]

the common phase takes the carrier out of the range information with it, and only the tone spacings remain. This is the incoherent estimator.

The phase model is the high-SNR limit of complex observations \(y_i = \exp(\mathrm{j}\,2\pi f_i d/c) + n_i\). In the equal-spacing example the first threshold lies near \(0.08\,\mathrm{rad}\), where that limit holds. At phase noise of order \(1\,\mathrm{rad}\) the phase model itself fails and the budget must return to the complex observations; that is the second threshold, treated separately from the first.

For a round trip the phase is \(4\pi f d/c\), which amounts to replacing \(c\) by \(c/2\): wavelengths and the ambiguity period halve, both mean-square lines below divide by 4, and \(\sigma_{\mathrm{jump}}\), which depends only on frequency ratios, keeps its value in radians.

## 2 Two parallel lines

When the noise is small and the error has not crossed into a neighbouring alias, what is the mean square of the coherent estimator? What is it once the common phase has been removed?

While the error stays in the main lobe the wrap can be dropped, and both estimators become linear least squares: with phase plotted against frequency, the coherent estimator fits a line through the origin, and the incoherent estimator also fits an intercept, the common phase. With \(F = \sum_i f_i^2\), the two mean squares are

\[
\mathrm{HMSE} = \frac{c^2\sigma^2}{4\pi^2 F}, \qquad \mathrm{MMSE} = \frac{c^2\sigma^2}{4\pi^2 S}.
\]

\(F\) comes mostly from the carrier, while \(S\) depends only on how the tones are spread about their mean; hence the names carrier line and bandwidth line. Both are proportional to \(\sigma^2\), and the ratio \(R = F/S\) does not depend on \(\sigma\), so the lines are parallel on a log-log plot. Fitting an intercept amounts to subtracting the mean frequency first: \(F = S + N\bar f^2\), so \(R = 1 + N\bar f^2/S\).

These heights are equations (30) and (42) of Wei Li and Wangdong Qi (arXiv:1202.1424), derived by perturbation. They name the mean squares in the high- and moderate-SNR regions HMSE and MMSE, and the names are kept here (MMSE is not the minimum mean-square error). They number thresholds in the direction of increasing SNR and call the transition from MMSE down to HMSE the second threshold; numbered in the direction of increasing phase noise, as here, the same transition is the first threshold. They do not write the noise level of this transition as a function of the frequency set.

For \(N\) equally spaced tones spanning \(B = (N-1)\Delta f\),

\[
R = 1 + \frac{12(N-1)}{N+1}\left(\frac{\bar f}{B}\right)^2.
\]

The coefficient tends to 12 for large \(N\), and the leading term becomes \(12(\bar f/B)^2\). Weiss and Weinstein obtained the same factor, \(12(\omega_0/W)^2\), for continuous-time delay estimation with a flat spectrum; that is the same limit reached from a different observation model, and a discrete set of wrapped phases is not a flat continuous spectrum. At \(N = 8\) the coefficient is 9.33, about 1.1 dB below the continuous-spectrum value at the same \(\bar f/B\).

The equal-spacing example (the equal grid) is eight tones from 2400 to 2470 MHz in 10 MHz steps: \(\bar f = 2435\,\mathrm{MHz}\), \(B = 70\,\mathrm{MHz}\), \(S = 4200\,\mathrm{MHz}^2\), \(R = 11295\) (40.53 dB), and mean wavelength \(\lambda = c/\bar f = 0.12312\,\mathrm{m}\).

## 3 Where the first threshold lies

At what phase noise does the coherent estimator jump from the carrier line to the bandwidth line?

### 3.1 The phase scale \(\sigma_{\mathrm{jump}}\)

Two range hypotheses \(\delta\) apart differ in phase by \(2\pi f_i\delta/c\) at tone \(i\). Write the squared length of the wrapped difference vector as

\[
J(\delta) = \sum_i \mathrm{wrap}(2\pi f_i\delta/c)^2 .
\]

Take \(\delta = \lambda\). The \(i\)-th component is \(2\pi f_i/\bar f = 2\pi + 2\pi(f_i-\bar f)/\bar f\); once the whole cycle is removed, \(2\pi(f_i-\bar f)/\bar f\) remains. Hence \(\sqrt{J(\lambda)} = 2\pi\sqrt{S}/\bar f = 2\sigma_{\mathrm{jump}}\), with

\[
\sigma_{\mathrm{jump}} = \pi\sqrt{S}/\bar f .
\]

Since \(F = S + N\bar f^2\), also \(R = 1 + N\pi^2/\sigma_{\mathrm{jump}}^2\). The phase scale of the threshold and the gap between the two lines are set by the same quantity, \(S/\bar f^2\).

In the equal-spacing example \(\sigma_{\mathrm{jump}} = 0.0836\,\mathrm{rad}\). At that noise a single phase reading is nowhere near \(\pm\pi\); the wrap matters because the range hypothesis has crossed one carrier cycle.

Treating the phase-difference vector as a Gaussian observation, the error probability of deciding between the two hypotheses is

\[
P_e(\delta) = Q\bigl(\sqrt{J(\delta)}/(2\sigma)\bigr),
\]

where \(Q\) is the standard normal right-tail probability. At \(\delta = \lambda\) the argument is exactly \(\sigma_{\mathrm{jump}}/\sigma\). At \(\sigma = \sigma_{\mathrm{jump}}\) the probability of a one-wavelength error to either side is \(2Q(1) = 0.317\) in total.

A local minimum of \(J\) at \(\delta \neq 0\) is called a valley. The smaller \(J\) at the bottom of a valley, the harder the two ranges are to tell apart and the larger \(P_e\). The valley near one mean wavelength has its bottom at \(\delta = \lambda(1-1/R)\), where \(J = (1-1/R)\,J(\lambda)\), within a fraction \(1/R\) of \(J(\lambda) = 4\sigma_{\mathrm{jump}}^2\).

This \(P_e\) comes from the Gaussian vector test, not from the wrapped likelihood ratio. Near the main lobe, where residuals are far below \(\pi\), the two agree; valleys with residuals near \(\pi\) appear in Section 5.4.

### 3.2 Where the mixture leaves the carrier line

Suppose a trial slips by one wavelength with probability \(p = 2Q(\sigma_{\mathrm{jump}}/\sigma)\) and otherwise stays in the main lobe. The mixture approximation of the mean square is

\[
\mathrm{MSE} \approx (1-p)\,\mathrm{HMSE} + p\lambda^2 .
\]

Setting \(p\lambda^2 = \mathrm{HMSE}\) gives \(\sigma = 0.0170\,\mathrm{rad}\) in the equal-spacing example, or \(0.204\,\sigma_{\mathrm{jump}}\); there the mixture leaves the carrier line. Where a measured curve starts to rise also depends on the number of trials: one expected slip in 8000 trials requires \(p = 1/8000\), which occurs at \(0.26\,\sigma_{\mathrm{jump}}\).

At \(\sigma = \sigma_{\mathrm{jump}}\), \(\mathrm{MMSE} = c^2/(4\bar f^2) = \lambda^2/4\). With \(p = 0.317\), \(p\lambda^2\) is 1.0 dB above the bandwidth line.

### 3.3 The equal-spacing trials

The coherent estimator takes the global minimum over a search window of half-length 14.990 m centred on the true range, 12.345 m; the window length, 29.979 m, is exactly the equal grid's ambiguity period \(c/\Delta f\). Each noise level from 0.22 to 0.50 \(\sigma_{\mathrm{jump}}\) has 8000 trials and the others 1600; the random seed is 7, and the results are in [`computation/data/`](computation/data/). A slip is an error larger than half a mean wavelength in magnitude.

At the three noise levels from 0.10 to 0.22 \(\sigma_{\mathrm{jump}}\) there are no slips, and the sample mean square is within 0.35 dB of the carrier line. At 0.25 \(\sigma_{\mathrm{jump}}\) one slip occurs in 8000 trials and the sample mean square rises to −57.2 dB, which is exactly the contribution of a single one-wavelength error averaged over 8000 trials. That rise comes from one large error, not from another asymptote between the two lines.

From 0.25 to 1.0 \(\sigma_{\mathrm{jump}}\) the sample mean square equals the slip rate times \(\lambda^2\), within 0.1 dB. At the threshold the slip rate is 0.344 and the sample mean square −22.8 dB, 1.4 dB above the bandwidth line (−24.2 dB); the mixture with \(p = 0.317\) predicts 1.0 dB, in the same direction. Beyond 1.5 \(\sigma_{\mathrm{jump}}\) some slips are longer than one wavelength, and the sample mean square rises above that product. Figure 1 shows the curve.

![Figure 1: mean-square error, equal spacing](https://franklinzhang.github.io/multifrequency-phase-ranging/figures/fig02_equal_mse.png)

Figure 1. Mean-square error on eight equally spaced tones, against \(\sigma\) normalized by \(\sigma_{\mathrm{jump}} = 0.0836\,\mathrm{rad}\). Black dashed: carrier line (HMSE); grey dash-dot: bandwidth line (MMSE); the coherent estimator leaves the first and rises onto the second. Green: half-window bound (Section 4); grey dotted: search-window floor; orange dashed: threshold.

The incoherent estimator lies on the bandwidth line from the lowest noise tried: −44.4 dB at 0.10 \(\sigma_{\mathrm{jump}}\), against a line at −44.2 dB. The common phase removes the carrier from the start, so there is no lower line for it to leave. \(\sigma_{\mathrm{jump}}\) is therefore a threshold only for estimators that use the absolute carrier phase.

### 3.4 Local information does not jump

The Fisher information of a single wrapped-normal phase about its mean is a smooth function of \(\sigma\). The Cramér–Rao bound for range built from it follows the carrier line and does not jump; in the equal-spacing example it is still within 0.03 dB of the carrier line at \(\sigma = 0.84\,\mathrm{rad}\). The jump comes from the global estimator leaving the main lobe; Fisher information describes only the curvature of the likelihood near the true range and says nothing about valleys outside the main lobe.

Single-tone frequency estimation by maximum likelihood has the same structure: Rife and Boorstyn showed that at high SNR the variance coincides with the Cramér–Rao bound, and below a threshold the estimator picks a wrong spectral peak. There the parameter is frequency, and the observations are not wrapped phases.

## 4 The half-window Ziv–Zakai bound

In Figure 1 the green curve leaves the carrier line before the trials do. What does it constrain?

Let the true range be uniform on an interval of length \(L\). Ziv and Zakai bounded the estimation error by the error probability of deciding between two hypotheses; in the valley-filled form of Bellini and Tartara the mean-square bound is

\[
\mathrm{ZZB} = \frac{1}{L}\int_0^L \delta\, G\bigl[(L-\delta)\,P_e(\delta)\bigr]\,\mathrm{d}\delta,
\qquad G[h](\delta) = \max_{\eta\ge\delta} h(\eta),
\]

so \(G\) replaces a function by its running maximum from the right. The pairwise error probability \(P_e\) is the Gaussian vector test of Section 3.1. Here \(L = 14.990\,\mathrm{m}\), half the length of the search window, hence the name half-window bound.

At high SNR only the main lobe contributes: \(\sqrt{J(\delta)} \approx \beta\delta\) with \(\beta = 2\pi\sqrt{F}/c\), and the integral becomes \((2\sigma/\beta)^2\int_0^\infty u\,Q(u)\,\mathrm{d}u\). Since \(\int_0^\infty u\,Q(u)\,\mathrm{d}u = 1/4\), the bound equals \(\sigma^2/\beta^2 = \mathrm{HMSE}\). In the equal-spacing example at \(\sigma = 0.008\,\mathrm{rad}\) the two differ by 0.00 dB, and only range differences up to about 0.5 mm contribute.

The bound is the Bayes risk when the true range is uniform on the half-window; the trials fix the true range at the centre of the search window and give a conditional mean square. At the threshold the bound is −27.6 dB, the bandwidth line −24.2 dB, and the coherent trials −22.8 dB. At 0.22 \(\sigma_{\mathrm{jump}}\), where the trials have no slips yet, the bound is already about 3 dB above the carrier line. The bound leaves the carrier line first and the trials later; it constrains the mean square and does not mark a turning point of the trial curve.

The half-window prior variance is \(L^2/12\), or 12.7 dB, the value the bound approaches at large noise. The trial mean square is capped by the variance of an error uniform over the whole search window (29.979 m), 18.74 dB: the search-window floor. At \(\sigma = 2.34\,\mathrm{rad}\) the equal-spacing bound is 8.1 dB, still below its prior variance. Taking \(L\) as the whole ambiguity period would let the replica valley at the far end of the period raise the bound, so the figures show the half-window bound.

The second threshold is a separate matter. Once \(\sigma\) reaches the order of 1 rad, the phase model fails and the budget must return to complex observations. The equal-spacing coherent trials are still on the bandwidth line at \(\sigma = 0.50\,\mathrm{rad}\) and reach the 18.7 dB search-window floor only at \(\sigma = 2.34\,\mathrm{rad}\). This threshold sets when the error reaches the floor; it does not change \(\sigma_{\mathrm{jump}}\).

## 5 What unequal spacing changes

What can unequal spacing change? A frequency set has three design quantities, and none substitutes for another.

1. \(S/\bar f^2\) sets the phase scale of the first threshold and the gap between the two lines (\(R = 1 + N\pi^2/\sigma_{\mathrm{jump}}^2\)). Spreading the tones away from the mean makes the jump shorter and moves it to larger phase noise.
2. The value of \(J(\delta)\) at the bottom of each deep valley, together with its offset \(\delta\), decides whether the coherent mean square stays on the bandwidth line at all, and whether a slip costs one wavelength squared or several metres squared. A valley with \(J(\delta)\approx J(\lambda)\) is as easy to slip into as a one-wavelength error.
3. The greatest common divisor of the tones sets the exact period on which every phase returns to a whole cycle, that is, the unambiguous range. A coprime or Chinese-remainder layout changes only this period; it does not postpone the first threshold.

### 5.1 Holding \(S\) fixed does not remove the jump

The second set has the same \(N\), \(\bar f\), and \(S\) as the equal grid (moment matching): 2390, 2418, 2425, 2432, 2443, 2449, 2454, and 2469 MHz. The tones are integer megahertz, the greatest common divisor of the gaps is 1 MHz, and the exact period is about 300 m. \(R\) and \(\sigma_{\mathrm{jump}}\) are therefore those of the equal grid.

At \(\sigma = 0.0836\,\mathrm{rad}\) the half-window integral splits into the same segments as on the equal grid, and the two \(J(\delta)\) curves coincide within 2 m. Breaking the equal spacing while holding \(S\) does not remove the jump; the two sets part only at larger noise (Section 5.4).

### 5.2 How valley segments enter the half-window integral

Valley filling does not add up the error probabilities of individual valleys. By the definition of \(G\), read the weighted error probability \(h(\delta) = (L-\delta)P_e(\delta)\) from right to left: each new running maximum at \(\delta_a\) sets the integrand over the whole stretch to its left, down to the first point \(\delta_b\) on the left where \(h\) is higher. That stretch contributes

\[
\frac{L-\delta_a}{L}\, P_e(\delta_a)\, \frac{\delta_a^2-\delta_b^2}{2}
\]

to the bound. When the left end falls in the main lobe, \(\delta_b \approx \sqrt{J(\delta_a)}/\beta\); for the equal grid at its threshold that end is only 1.2 mm.

At \(\sigma = 0.0836\,\mathrm{rad}\) the equal grid's one-wavelength segment contributes −29.2 dB, the two-wavelength segment −32.9 dB, and the three-wavelength segment −43.0 dB; together they give −27.6 dB, the value of the numerical half-window integral. Split by the numerical integral, the one-wavelength valley carries 67% and the two-wavelength valley 29%. Valleys further out at \(k\lambda\) barely enter: \(J\) at the bottom of the \(k\)-th valley grows as \(k^2\), and \(P_e\) falls much faster than \(\delta^2\) grows.

Geometric spacing behaves differently. Its far valley lies where \(L-\delta\) is smaller, so its weighted error probability is below that at one mean wavelength and does not extend back to the main lobe; the far valley sets a long stretch from one mean wavelength out to itself. Figure 2 compares the three sets at the same noise.

![Figure 2: cumulative fraction of the half-window bound](https://franklinzhang.github.io/multifrequency-phase-ranging/figures/fig08_bound_cumulative.png)

Figure 2. Cumulative fraction of the half-window bound along range offset, at the equal-grid threshold 0.0836 rad. The equal and moment-matched curves coincide, reaching 96% at 0.25 m and 100% within 0.4 m; the geometric set (green) completes only at 6.21 m.

### 5.3 Geometric spacing: a far valley as deep as one wavelength

The third set is a geometric progression with ratio 1.02, geometric centre 2435 MHz and arithmetic mean 2437.5 MHz. \(S\) is about 23 times that of the equal grid, \(R\) drops to 26.87 dB, \(\sigma_{\mathrm{jump}}\) rises to 0.403 rad, and the mean wavelength is still 0.1230 m. This is a change in the first design quantity, \(S/\bar f^2\): the threshold moves to larger phase noise and the two lines move closer.

At the bottom of the valley at 6.211 m, \(J = 0.650\), the same as \(J(\lambda) = 0.650\), and none of the eight phase residuals exceeds 0.44 rad in magnitude. It is a deep sidelobe, not a point where all tones return to a whole cycle. It lies about 50 mean wavelengths out, so its squared error is about 34 dB larger: a slip there is as likely as a one-wavelength error and costs 6.2 m.

Because \(L - 6.211\,\mathrm{m}\) is smaller, this valley's weighted error probability is below that of the one-wavelength valley and does not extend back to the main lobe; the stretch it sets runs from one mean wavelength to 6.211 m and carries 96% of the half-window integral. The contribution is spread along that stretch: half of the integral comes from beyond 4.3 m and nine tenths from inside 5.9 m; it is not concentrated at 6.211 m. The half-window bound is therefore −51 dB, while the carrier line is −64.8 dB.

The trials use the same absolute noise levels as for the equal grid, so \(\sigma = 0.0836\,\mathrm{rad}\) is only 0.21 of the geometric set's own threshold. There the coherent trials are still on the carrier line (−64.5 dB) while the bound is already −51 dB. At 0.31 of its own threshold the coherent and incoherent estimators both reach −13.2 dB, with a slip rate of only 0.0031. The same slip rate charged at one wavelength would give about −43 dB; \(0.0031\times(6.21\,\mathrm{m})^2\) is about −9 dB, the same order as −13 dB. The coherent mean square does not stay on the bandwidth line. The mixture leaving point computed with one mean wavelength, 0.0967 rad, gives the noise at which slips start; it does not give how far they go. Figure 3 shows these curves.

![Figure 3: mean-square error, geometric spacing](https://franklinzhang.github.io/multifrequency-phase-ranging/figures/fig06_geometric_mse.png)

Figure 3. Mean-square error for the geometric set, against \(\sigma\) normalized by its own threshold 0.403 rad. Black dashed: carrier line; grey dash-dot: bandwidth line. The coherent estimator (blue circles) leaves the carrier line and crosses the bandwidth line without staying on it; the incoherent estimator (orange squares) starts on the bandwidth line and crosses it at the same noise. Green: half-window bound.

| Design | \(\bar f\) | \(S\) | \(R\) | \(\sigma_{\mathrm{jump}}\) |
| --- | --- | --- | --- | --- |
| Equal spacing | 2435 MHz | 4200 MHz² | 40.53 dB | 0.0836 rad |
| Moment-matched | 2435 MHz | 4200 MHz² | 40.53 dB | 0.0836 rad |
| Geometric | 2437.5 MHz | \(9.78\times 10^4\) MHz² | 26.87 dB | 0.403 rad |

Table 1. The three tone sets. Moment matching shares \(\bar f\), \(S\), and the threshold with the equal grid; the geometric set has larger \(S\), a higher threshold, and a smaller gap between the lines.

### 5.4 A cluster of valleys at larger noise

The moment-matched and equal sets part only when the phase noise is ten times larger. At \(\sigma = 0.836\,\mathrm{rad}\) their half-window bounds are −2.4 and −3.8 dB. The moment-matched set has a cluster of valleys from 10.6 to 11.3 m with \(J\) about 14 at the bottom; the equal grid has about 25 over the same stretch. Within the cluster the weighted error probability is highest at 10.65 m, and the segment it sets carries 29% of the moment-matched integral; removing that segment drops the bound from −2.4 dB to about −3.9 dB, matching the equal grid's −3.8 dB. The deepest valley in the window is not at 10.65 m but slightly farther out, at 11.0 m; the weight \(L-\delta\) favours the nearer one, so the running maximum from the right falls at 10.65 m. The 1.4 dB comes from this far cluster, not from nearby valleys. \(J\) at the bottom of these valleys is about 14, not 0, so they are not exact whole-cycle coincidences either. Their phase residuals reach about 2.5 rad, close to half a cycle, where the Gaussian vector test of Section 3.1 is least reliable.

## 6 Engineering rules

These rules apply where the phase model holds, that is, for \(\sigma\) well below 1 rad; at the order of 1 rad the budget moves to complex observations (the second threshold).

With a common phase reference:

- Budget the carrier line only below about \(0.2\,\sigma_{\mathrm{jump}}\). In the equal-spacing example that is a phase standard deviation of about one degree, 0.017 rad.
- From about \(0.5\,\sigma_{\mathrm{jump}}\) until the phase model fails, budget the bandwidth line, 40.53 dB worse than the carrier line in the equal-spacing example.
- In between, the mean square is close to the slip rate times \(\lambda^2\) (Section 3.3).

Without a common phase, the high-SNR budget is the bandwidth line from the start; the jump off the carrier line has no operating region. If a far valley with \(J(\delta)\approx J(\lambda)\) exists, removing the common phase removes the carrier line but not the rise that the valley causes.

When placing tones, solve \(2P_e(\delta)\,\delta^2 = \mathrm{HMSE}\) at every pair of local maxima \(\pm\delta\) of \(P_e\) (every pair of valleys) and keep the earliest \(\sigma\); quoting only \(\sigma_{\mathrm{jump}}\) is not enough. On the equal grid the solution coincides with \(0.20\,\sigma_{\mathrm{jump}}\), because the deciding valley is at one mean wavelength. On the geometric set the far valley raises the bound earlier: at 0.21 of its own threshold the bound is already −51 dB while the carrier line is −64.8 dB, and an undetected slip, as likely as a one-wavelength error, costs about 34 dB more in mean square.

Assess the three design quantities separately: spreading \(S/\bar f^2\) changes the first threshold and the gap between the lines; the deep valleys' \(J(\delta)\) and \(\delta\) decide whether the mean square stays on the bandwidth line and how far a slip goes; the greatest common divisor changes the unambiguous range. A longer exact period moves only the third quantity and leaves the first threshold where it was.

## 7 Conclusions and open questions

### 7.1 Conclusions

1. The heights of the two lines are set by \(F\) and \(S\); \(R = F/S = 1 + N\pi^2/\sigma_{\mathrm{jump}}^2\) does not depend on \(\sigma\), so the lines are parallel on a log-log plot.
2. The phase scale of the first threshold is \(\sigma_{\mathrm{jump}} = \pi\sqrt{S}/\bar f\), the scale of the Gaussian test between two ranges one mean wavelength apart. It is a threshold only for estimators that use the absolute carrier phase: the incoherent estimator lies on the bandwidth line from the lowest noise, so no theorem places every estimator's jump at \(\sigma_{\mathrm{jump}}\).
3. On the equal grid the mixture leaves the carrier line at \(0.204\,\sigma_{\mathrm{jump}}\); the first slip in 8000 trials is at \(0.25\,\sigma_{\mathrm{jump}}\); from there to the threshold the sample mean square equals the slip rate times \(\lambda^2\) within 0.1 dB.
4. Local Fisher information and the Cramér–Rao bound do not jump.
5. The half-window Ziv–Zakai bound is a Bayes risk under a uniform prior; it leaves the carrier line before the conditional trials and does not mark a turning point of the trial curve.
6. Unequal spacing at fixed \(S\) does not remove the first threshold; a far valley with \(J(\delta)\approx J(\lambda)\) keeps the mean square from staying on the bandwidth line and stretches a slip to metres.
7. The unambiguous range and the first threshold are separate design quantities: a Chinese-remainder layout lengthens the former and does not postpone the latter.

The global wrapped least-squares mean square at finite \(N\) has no closed form. What is established is the phase scale of the threshold, the valley-filled bound, and the comparison of both with the archived trials.

### 7.2 Open questions

Each of the following can be settled by computation and trials, without that closed form.

#### Placing tones against the worst valley

Under a band limit, a fixed \(N\), and an integer-megahertz synthesizer step, can the earliest leaving noise over all deep valleys be pushed above the equal grid's 0.0170 rad without leaving a valley with \(J(\delta)\approx J(\lambda)\) at \(\delta\gg\lambda\)? Moment matching shows that shuffling the gaps at fixed \(S\) does nothing; \(S\) and the positions of the deep valleys have to be moved together. A tone set that meets both conditions answers the question.

#### Bayes mean square under one prior

At the threshold the half-window bound is −27.6 dB and the coherent trials −22.8 dB, about 5 dB apart. The bound is a Bayes risk under the half-window uniform prior; the trials are a conditional mean square at the centre of the search window, so the two cannot be compared directly. Computing the mean square of the posterior mean under the same prior, and comparing it with coherent maximum likelihood on the same noise grid, decides whether the 5 dB belongs to a loose bound or to a suboptimal estimator. The same comparison shows whether integer least squares for the cycle ambiguity, or integrating over the ambiguity edges, is worth doing.

#### A wrapped likelihood ratio only where residuals approach half a cycle

The Gaussian vector test is a good approximation when residuals are far below \(\pi\). At the first threshold \(\sigma = 0.08\,\mathrm{rad}\), and the residuals of the geometric set's 6.211 m valley stay below 0.44 rad; the approximation holds at both. It is least reliable on the cluster from 10.6 to 11.3 m at \(\sigma = 0.836\,\mathrm{rad}\), with residuals up to about 2.5 rad. Replacing the Gaussian vector test by the wrapped likelihood ratio on that cluster alone, and checking whether the 1.4 dB gap between the moment-matched and equal-grid bounds moves, answers the question; error probabilities at other offsets need not be recomputed first.

#### The first threshold on complex observations

The phase model is the high-SNR limit of the complex observations and holds at the first threshold. Running coherent maximum likelihood directly on the complex observations of the same equal grid checks whether the first threshold stays near \(\sigma_{\mathrm{jump}}\); the expectation is that the first jump stays and the floor and the second threshold move. Amplitude fading, unequal SNR across tones, and correlated phase noise from a shared oscillator are budget corrections to add after that check.

#### Second moment once slips exceed one wavelength

Beyond \(1.5\,\sigma_{\mathrm{jump}}\) the equal-grid sample mean square exceeds the slip rate times \(\lambda^2\), because some slips are longer than one wavelength. The first few equal-grid valleys sit near \(k\lambda\), with \(J = k^2 J(\lambda)\) at the bottom and pairwise error probability \(Q(k\sigma_{\mathrm{jump}}/\sigma)\). The valleys at right-to-left running maxima of the weighted error probability (Section 5.2) can be weighted by \(\delta^2\) and \(Q(k\sigma_{\mathrm{jump}}/\sigma)\) in place of the single \(\lambda^2\); for unequal spacing, use the same running-maximum points instead of integer multiples of \(\lambda\). The question is answered when, beyond 1.5 \(\sigma_{\mathrm{jump}}\), this weighting tracks the sample mean square clearly better than the \(\lambda^2\) product and the remaining mismatch can be located.

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

## References

Wei Li and Wangdong Qi, Optimization in Multi-Frequency Interferometry Ranging: Theory and Experiment, arXiv:1202.1424, 2012. Equations (30) and (42) give HMSE and MMSE; the transition between them is their second threshold (numbered with increasing SNR), the first threshold here; the analysis is by perturbation and gives no \(\sigma_{\mathrm{jump}}\).

A. Weiss and E. Weinstein, Fundamental Limitations in Passive Time-Delay Estimation — Part I, IEEE Transactions on Acoustics, Speech, and Signal Processing, 1983; Part II, 1984. The factor \(12(\omega_0/W)^2\) for a flat continuous spectrum, equal to the large-\(N\) limit of the equal-grid \(R\); different observation model.

J. Ziv and M. Zakai, Some Lower Bounds on Signal Parameter Estimation, IEEE Transactions on Information Theory, 1969. Lower bounds on estimation error from the error probability of binary hypothesis tests.

S. Bellini and G. Tartara, Bounds on Error in Signal Parameter Estimation, IEEE Transactions on Communications, 1974. The valley-filled form of the Ziv–Zakai bound.

D. C. Rife and R. R. Boorstyn, Single-Tone Parameter Estimation from Discrete-Time Observations, IEEE Transactions on Information Theory, vol. 20, no. 5, pp. 591–598, 1974. Tight local bound and global departure at a threshold, for frequency rather than range.
