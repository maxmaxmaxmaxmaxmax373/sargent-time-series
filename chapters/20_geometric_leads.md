# Predicting Geometric Distributed Leads

It is important to know the solution of the following problem in order to use a variety of linear rational expectations models. Let $x_t$ be a covariance stationary stochastic process with Wold moving average representation

```{math}
:label: eq-82
x_t = c(L)\epsilon_t, \quad c_0 = 1
```

where $\epsilon_t$ is a fundamental white noise for $x$ and $c(L) = \sum_{j=0}^\infty c_j L^j$ is square summable. We further assume that $c(L)$ has an inverse $a(L) = c(L)^{-1}$ which is one-sided in nonnegative powers of $L$ and square summable. Thus, $x_t$ has the autoregressive representation

```{math}
:label: eq-83
a(L)x_t = \epsilon_t
```

where $a(L) = 1 - a_1 L - a_2 L^2 - \cdots$.

We want to calculate the following linear projection

```{math}
:label: eq-84
y_t = P\left[\sum_{j=0}^\infty \lambda^j x_{t+j}\Big|x_t, x_{t-1},\ldots\right] \equiv P_t\sum_{j=0}^\infty \lambda^j x_{t+j}
```

where $|\lambda|<1$. Projections of such geometric distributed leads occur in a variety of linear rational expectations models. We begin by noting that $y_t$ defined by {eq}`eq-84` satisfies the stochastic difference equation

```{math}
:label: eq-85
y_t = \lambda P_t y_{t+1} + x_t.
```

That is, $y_t$ is the stationary solution of the difference equation {eq}`eq-84` as can be verified by repeated substitution in {eq}`eq-85`. We seek expressions for $y_t$ of the forms

```{math}
:label: eq-86
y_t = g(L)x_t
```

and

```{math}
:label: eq-87
y_t = d(L)\epsilon_t
```

where

$$
d(L) = \sum_{j=0}^\infty d_j L^j,\quad g(L) = \sum_{j=0}^\infty g_j L^j,\quad \sum_{j=0}^\infty g_j^2 < +\infty,\quad \sum_{j=0}^\infty d_j^2 < + \infty
$$

We know that representation {eq}`eq-86` exists by definition, and therefore that $d(L) = g(L) c(L)$ also exists. That is, a representation of the form {eq}`eq-87` exists because $\{x_t, x_{t-1}, \ldots\}$ and $\{\epsilon_t, \epsilon_{t-1},\ldots\}$ span the same space.

We shall solve for $d(L)$ using {eq}`eq-85` and prediction theory. Using {eq}`eq-87`, we have that

$$
P_t y_{t+1} = \left[\frac{d(L)}{L}\right]_{+}\epsilon_t
$$

or

$$
P_t y_{t+1} = \left[\frac{d(L)}{L} - \frac{d_0}{L}\right]\epsilon_t.
$$

Substituting this and {eq}`eq-82` and {eq}`eq-87` into {eq}`eq-85` gives

$$
d(L)\epsilon_t = \lambda\left[\frac{d(L)}{L} - \frac{d_0}{L}\right]\epsilon_t + c(L)\epsilon_t.
$$

Since this equation holds for all $\epsilon_t$ realizations, it implies, after rearranging, that

$$
(1 - \lambda L^{-1})d(L) = c(L) - \lambda d_0 L^{-1}
$$

an equation that we desire to solve for $d(L)$ as a function of $c(L)$. We determine $d_0$ by evaluating the above equation at $L = \lambda$, to get $c(\lambda) = d_0$. Using this value for $d_0$ gives

```{math}
:label: eq-88
d(L) = \frac{c(L) - \lambda c(\lambda) L^{-1}}{1 - \lambda L^{-1}}
```

using $g(L) = d(L)c(L)^{-1}$ and $c(L)^{-1} = a(L)$, we get

```{math}
:label: eq-89
g(L) = \frac{1 - \lambda a(\lambda)^{-1} a(L) L^{-1}}{1 - \lambda L^{-1}}.
```

For the case in which $a(L)$ is an $r$th order polynomial $a(L) = 1 - \sum_{j=1}^r a_j L^j$, Hansen and Sargent (1980) show using polynomial long division that {eq}`eq-89` can be expressed

```{math}
:label: eq-90
g(L) = a(\lambda)^{-1}\left[1 + \sum_{j=1}^{r-1}\left(\sum_{k=j+1}^r \lambda^{k-j} a_k\right)L^j\right]
```

so that

$$
g(L) = \sum_{j=0}^{r-1} g_j L^j,
$$

with

$$
g_0 = a(\lambda)^{-1},\quad g_j = a(\lambda)^{-1} \sum_{k=j+1}^r \lambda^{k-j} a_k
$$

for $j = 1,\ldots, r-1$. Evidently, the coefficients $g_j$ can be computed recursively via the following formulas:

```{math}
:label: eq-91
\begin{aligned}
g_0 &= a(\lambda)^{-1} \\
g_r &= 0 \\
g_{j-1} &= \lambda g_j + \lambda g_0 a_j \quad j=r,r-1, \ldots, 2.
\end{aligned}
```

Various versions of formulas {eq}`eq-88`, {eq}`eq-89`, and {eq}`eq-90` were originally derived in papers by Saracoglu and Sargent (1978), Hansen and Sargent (1980), and Futia (1981).
