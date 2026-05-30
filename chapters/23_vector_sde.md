# Vector Stochastic Difference Equations

Let $x_t$ be an $(n \times 1)$-vector wide-sense stationary stochastic process that is governed by the matrix difference equation

```{math}
:label: eq-100
C(L)x_t = \epsilon_t
```

where $\epsilon_t$ is now an $n \times 1$ vector of white noises with means of zero and contemporaneous covariance matrix $E \epsilon_t \epsilon_t' = V$, an $n \times n$ matrix. We assume $E \epsilon_t \epsilon_{t-s}' = 0_{n \times n}$ for all $s \neq 0$. In {eq}`eq-100`, $C(L)$ is an $n \times n$ matrix of (finite order) polynomials in the lag operator $L$:

$$
C(L) = \begin{bmatrix}
    C_{1 1}(L) & C_{1 2}(L) & \cdots & C_{1 n}(L) \\
    \vdots & & & \\
    C_{n 1}(L) & \cdots & & C_{n n}(L)
\end{bmatrix}
$$

where each $C_{i j}(L)$ is a finite order polynomial in the lag operator.

We assume that the matrix $C(L)$ has an inverse under convolution $C(L)^{-1} \equiv B(L)$; $C(L)^{-1}$ is defined as the matrix that satisfies

$$
C(L)^{-1}C(L) = I_{n \times n}
$$

where $I_{n \times n}$ is the $n \times n$ identity matrix. If it exists, $C(L)^{-1}$ can be found as follows. Evaluate the matrix $z$ transform $C(z)$ at $z = e^{-i\omega}$ to get $C(e^{-i\omega})$. Then invert $C(e^{-i\omega})$, frequency by frequency, to get $C(e^{-i\omega})^{-1}$. Finally, the matrix coefficients $C(L)^{-1} = B(L) = \sum_{j=0}^\infty B_j L^j$, $B_j$ being an $n \times n$ matrix, can be found from the inversion formula

$$
B_j = \frac{1}{2\pi}\int_{-\pi}^{\pi} C(e^{-i\omega})^{-1}e^{i \omega j} d\omega,
$$

where by integrating a matrix we mean to denote element-by-element integration.

A solution of {eq}`eq-100` is found by premultiplying by $B(L)$ to obtain

```{math}
:label: eq-101
x_t = B(L)\epsilon_t.
```

The vector stochastic difference equation $C(L)x_t = \epsilon_t$ is said to be an *autoregressive representation* for the vector process $x_t$. The solution $x_t = B(L)\epsilon_t$ is said to be a vector *moving average representation* for the process $x_t$. The cross-spectral density matrix of the $n\times 1$ $x_t$ process (which has the cross spectrum between the $i$th and $j$th components of $x$ in the $(i,j)$th position) is given by

```{math}
:label: eq-102
g_{xx}(e^{-i\omega}) = B(e^{-i\omega})V B(e^{+i\omega})'
```

where the prime denotes transposition. Formula {eq}`eq-102` is analogous to the univariate equation (5), and can be derived by comparable methods.

Equation {eq}`eq-102` is a very compact formula for calculating the cross spectra of the $n\times 1$ $x_t$ process as a function of the fundamental parameters, the covariance matrix $V$ and the coefficients in $C(L)$ (or $B(L)$). Equation {eq}`eq-100` is quite a general representation and is flexible enough to incorporate exogenous variables and serially correlated noises.

In Equation {eq}`eq-100` a variable $x_{it}$ is said to be exogenous if $C_{ij}(L) = 0$ for all $j$ not equal to $i$. This means that the row of Equation {eq}`eq-100` corresponding to $x_{it}$ becomes $C_{ii}(L)x_{it}$ so that $x_{it}$ is governed by only its own past interacting with the random shock $\epsilon_{it}$. In this sense the evolution of $x_{it}$ is not affected by interactions with other variables in $x_t$. This is not to say however that $x_{it}$ is uncorrelated with other components of $x_t$, since $\epsilon_{it}$ can be correlated contemporaneously with other $\epsilon$'s (i.e., $V$ need not be diagonal). The definition of exogeneity given here turns out to be precisely the one used by econometricians in a time series context (see Section 24).

Serially correlated errors can be incorporated by suitably redefining the errors as components of $x_t$ and then modeling them as exogenous processes that affect but are not affected by other components of $x_t$.
