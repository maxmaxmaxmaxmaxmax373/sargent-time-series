# Linear Time Series Analysis

This is a modernized and extended version of **Chapter XI: Time Series** from

> Thomas J. Sargent, *Macroeconomic Theory*, 2nd edition (1987), Academic Press.

The original chapter develops the theory of linear stochastic difference equations,
spectral analysis, linear prediction, Wold decomposition, Granger causality, and
rational expectations—tools that remain central to modern macroeconometrics and
quantitative economics.

## What this book adds

Relative to the 1987 original, this version:

- **Corrects typographical and mathematical errors** present in the LaTeX source.
- **Adds Python code** that generates modern versions of all figures, using
  current U.S. and international data.
- **Extends several sections** with updated empirical examples and commentary.
- **Cross-references QuantEcon lectures** where related code already exists
  (see [ARMA](https://python-advanced.quantecon.org/arma.html) and
  [Spectral Estimation](https://python-advanced.quantecon.org/estspec.html)).

## How to read this book

The chapter can be read linearly. Section 4 (Fourier and $z$-transforms) is marked
*optional on first reading* and may be skipped without loss of continuity.

## Notation

Throughout, $L$ denotes the **lag operator**: $L x_t = x_{t-1}$.
All stochastic processes are assumed to be discrete-time unless otherwise stated.
$E$ denotes the mathematical expectation operator.

## References

- Sargent, T.J. (1987). *Macroeconomic Theory*, 2nd ed. Academic Press.
- Whittle, P. (1983). *Prediction and Regulation*, 2nd ed. University of Minnesota Press.
- Wold, H. (1938). *A Study in the Analysis of Stationary Time Series*. Almqvist & Wiksell.
- Anderson, T.W. (1971). *The Statistical Analysis of Time Series*. Wiley.
