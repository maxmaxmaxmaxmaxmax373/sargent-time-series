# Linear Time Series Analysis

This book is a modernized and extended treatment of the linear time series methods that run
through macroeconomics and dynamic econometrics. Its core is **Chapter XI: Time Series** from

> Thomas J. Sargent, *Macroeconomic Theory*, 2nd edition (1987), Academic Press,

now preceded by two chapters — **Chapter IX** and **Chapter X** — that develop the elementary
tools the time series theory rests on, and followed by several new sections that extend the
original.

## A book in three movements

The material is organized as a single arc that runs from elementary algebra and geometry to a
full theory of linear prediction and its uses.

- {doc}`Chapter IX — Difference Equations and Lag Operators <chapters/ch09_difference_equations>`
  introduces the **lag operator** $L$, the algebra of polynomials in $L$, and the use of that
  algebra to solve linear difference equations, to factor characteristic polynomials into
  stable and unstable roots, and to solve the linear-quadratic optimization (Euler equation)
  problems at the heart of dynamic economics.

- {doc}`Chapter X — Linear Least Squares Projections (Regressions) <chapters/ch10_regressions>`
  introduces the **linear least squares projection** — the regression viewed geometrically
  through the orthogonality principle — together with recursive projection, the law of
  iterated projections, and the signal-extraction problem.

These two chapters supply the two elementary methods — **lag operators** and **linear least
squares projections** — that, combined, set the stage for the linear prediction theory of
Chapter XI. The prediction theory is, in essence, what one obtains by projecting the solution
of a difference equation driven by white noise onto its own past: the lag-operator algebra of
Chapter IX organizes the dynamics, while the projection geometry of Chapter X organizes the
forecasting.

- **Chapter XI — Linear Time Series** (the bulk of this book) puts the two together. It studies
  covariance stationary stochastic processes built from white noise by linear difference
  equations, their {doc}`spectra <chapters/06_spectrum>` and
  {doc}`cross spectra <chapters/07_cross_spectrum>`, the
  {doc}`Wold decomposition <chapters/13_representation_theory>`, the
  {doc}`Wiener–Kolmogorov theory of linear prediction <chapters/14_linear_prediction>`,
  {doc}`signal extraction <chapters/19_signal_extraction>`,
  {doc}`Granger causality <chapters/27_granger_causality>`, and the cross-equation restrictions
  of {doc}`rational expectations <chapters/22_rational_expectations>`.

## What is new in this edition

Beyond modernizing the 1987 text, this edition adds several sections that extend the theory or
connect it to recent work:

- {doc}`The uncertainty principle for Fourier transform pairs <chapters/05a_uncertainty_principle>` —
  the time–frequency trade-off that limits what any filter or window can resolve.
- {doc}`The residue theorem behind partial fractions <chapters/18a_partial_fractions>` — the
  complex-analysis foundation of the partial-fraction calculus used throughout the prediction
  and signal-extraction chapters.
- **Whittle's spectral factorization** (in
  {doc}`Representation Theory <chapters/13_representation_theory>`) — a constructive, FFT-based
  way to recover the Wold moving-average kernel and the innovation variance from a spectral
  density; the computational complement to the existence theorems of
  {doc}`Chapters 16–18 <chapters/17_wold_ma>`.
- {doc}`A difficulty in interpreting vector autoregressions <chapters/36a_interpreting_vars>` —
  following Hansen and Sargent (1991), an example in which the innovations a vector
  autoregression recovers are *not* the shocks that hit agents, with dynamic supply and demand
  curves derived by the stable-roots-backward / unstable-roots-forward method of
  {doc}`Chapter IX <chapters/ch09_difference_equations>`.
- {doc}`Complex demodulation <chapters/41_comp_demod>` — a tool for estimating a *moving*
  spectrum and cross spectrum, applied to the changing seasonal in U.S. interest rates.

## What this book adds throughout

Relative to the 1987 original, this version:

- **Corrects typographical and mathematical errors** present in the LaTeX source.
- **Adds Python code** that generates modern versions of every figure, using current U.S. and
  international data; each figure links to the script that produced it.
- **Extends several sections** with updated empirical examples and commentary.
- **Cross-references QuantEcon lectures** where related code already exists
  (see [ARMA](https://python-advanced.quantecon.org/arma.html) and
  [Spectral Estimation](https://python-advanced.quantecon.org/estspec.html)).

## How to read this book

A reader new to the material can proceed linearly: Chapters IX and X first, then Chapter XI. A
reader already comfortable with lag operators and regressions can begin directly at
{doc}`chapters/01_introduction` and refer back to
{doc}`Chapter IX <chapters/ch09_difference_equations>` and
{doc}`Chapter X <chapters/ch10_regressions>` as needed. Two mathematical digressions —
{doc}`Fourier and z-transforms <chapters/04_fourier_z_transforms>` and
{doc}`the uncertainty principle <chapters/05a_uncertainty_principle>` — are marked *optional on
first reading* and may be skipped without loss of continuity.

## Notation

Throughout, $L$ denotes the **lag operator**, $L x_t = x_{t-1}$, introduced in
{doc}`Section 1 of Chapter IX <chapters/ch09_difference_equations>`. $E$ denotes the
mathematical expectation operator, and $\hat E$ or $P[\,\cdot \mid \cdot\,]$ the **linear least
squares projection** operator of {doc}`Chapter X <chapters/ch10_regressions>`. All stochastic
processes are discrete time unless otherwise stated.

## References

- Sargent, T.J. (1987). *Macroeconomic Theory*, 2nd ed. Academic Press.
- Hansen, L.P., and T.J. Sargent (1991). Two difficulties in interpreting vector
  autoregressions. In *Rational Expectations Econometrics*, ch. 4. Westview Press.
- Whittle, P. (1983). *Prediction and Regulation*, 2nd ed. University of Minnesota Press.
- Wold, H. (1938). *A Study in the Analysis of Stationary Time Series*. Almqvist & Wiksell.
- Anderson, T.W. (1971). *The Statistical Analysis of Time Series*. Wiley.
