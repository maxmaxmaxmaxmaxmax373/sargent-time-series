# Linear Least Squares Prediction[^fn-pred-1]

It is common in economics to assume that $x_t$ is purely (linearly) indeterministic, which means that $\eta_t = 0$ for all $t$ or else that $\eta_t$ has been removed.[^fn-pred-2] Wold's theorem says that *any* indeterministic covariance stationary stochastic process $x_t$ has *the moving average representation*

$$
x_t = \sum_{j=0}^\infty d_j \epsilon_{t-j}
$$

or

```{math}
:label: eq-59
x_t = d(L)\epsilon_t, \qquad d(L) = \sum_{j=0}^\infty d_j L^j
```

where $\{\epsilon_t \}$ is the sequence of one-step-ahead linear least squares forecasting errors (innovations) in predicting $x_t$ as a linear function of $\{x_{t-1}, x_{t-2}, \ldots\}$, i.e., $\epsilon_t = x_t - P[x_t | x_{t-1}, x_{t-2}, \ldots]$. (As we have seen, it is natural to normalize $d(L)$ so that $d_0 = 1$, in which case $\sigma^2 = E \epsilon_t^2$ is the variance of the one-step-ahead prediction error.)

Now suppose that $d(L)$ has an inverse that is one-sided in nonnegative powers of $L$. Where $d(L) = \sum_{j=0}^{\infty} d_j L^j$, a necessary condition for $d(L)$ to have such a one-sided inverse is that the roots $\mu$ of $\sum_{j=0}^n d_j\mu^j = 0$ all like outside the unit circle, i.e., all have absolute values exceeding unity. An inverse $a(L) \equiv d(L)^{-1}$ of $d(L)$ satisfies $a(L)d(L) = d(L) a(L) = I$ where $I$ is the identity lag operator $I = 1 + 0L + 0L^2 + \ldots$. Operating on both sides of {eq}`eq-59` with $a(L) = d(L)^{-1}$ gives

```{math}
:label: eq-60
a(L)x_t = \epsilon_t, \qquad a(L) = a_0 - \sum_{j=1}^\infty a_j L^j
```

or

$$
a_0 x_t = a_1 x_{t-1} + a_2 x_{t-2} + \ldots + \epsilon_t.
$$

Since $d_0$ is unity, it turns out that $a_0$ is unity also. Equation {eq}`eq-60` is termed the autoregressive representation for $x_t$. While every linearly indeterministic covariance stationary process has a moving average representation, not all of them have an autoregressive representation. Still, those that do have both a moving average and an autoregressive representation constitute a very wide class, and we shall henceforth assume that we are dealing with a member of this class.[^fn-pred-3]

We now derive some formulas due to Wiener and Kolmogorov for linear least squares predictors. Let $P_{t-j}x_t$ be the linear least squares projection of $x_t$ on the space spanned by $\{ x_{t-j}, x_{t-j-1}, \ldots\}$ i.e.,

$$
P_{t-1}x_t \equiv P[x_t | x_{t-j}, x_{t-j-1}, \ldots ].
$$

Now project both sides of {eq}`eq-59` against $\{x_{t-1}, x_{t-2}, \ldots\}$ to get

$$
P_{t-j}x_t = \sum_{j=0}^\infty d_j P_{t-1} \epsilon_{t-j} = \sum_{j=1}^{\infty} d_j \epsilon_{t-j},
$$

which follows since $P_{t-1}\epsilon_t = 0$, because $\epsilon_t$ is orthogonal to lagged $x$'s; and since $P_{t-1} \epsilon_{t-j} = \epsilon_{t-j}$ for all $j \geq 1$ because $\epsilon_{t-j}$ is in the space spanned by $\{x_{t-1}, x_{t-2}, \ldots\}$. We write the above equation as

$$
P_{t-1}x_t = \left(\frac{d(L)}{L}\right)_{+} \epsilon_{t-1}
$$

where $(\,)_{+}$ means "ignore negative powers of $L$," i.e., $\left(\sum_{j = -\infty}^{\infty} h_j L^j\right)_{+} \equiv \sum_{j=0}^\infty h_j L^j$. Now assuming that $x_t$ has an autoregressive representation, we can write $\epsilon_{t-1} = a(L)x_{t-1} = d(L)^{-1} x_{t-1}$. Substituting this into the above equation gives

```{math}
:label: eq-61
P_{t-1}x_t = \left(\frac{d(L)}{L}\right)_{+} \frac{1}{d(L)} x_{t-1}
```

which is a compact formula for the one-step-ahead linear least squares forecast of $x_t$ based on it own past.

To get a formula for the general $k$-step-ahead linear least squares forecast, project both sides of {eq}`eq-50` against $\{ x_{t-k}, x_{t-k-1}, \ldots\}$ to get

$$
P_{t-k}x_t = \sum_{j=k}^\infty d_j \epsilon_{t-j} = \left(\frac{d(L)}{L^k}\right)_+ \epsilon_{t-k}
$$

```{math}
:label: eq-62
P_{t-k}x_t = \left(\frac{d(L)}{L^k}\right)_{+} \frac{1}{d(L)}x_{t-k}
```

which generalizes formula {eq}`eq-61`. Equation {eq}`eq-62` is the Wiener-Kolmogorov formula for $k$-step-ahead linear least squares predictions.

[^fn-pred-1]: A key reference on the subject of this section is Whittle (1983).

[^fn-pred-2]: For example, by suitable detrending and seasonal adjustment.

[^fn-pred-3]: We remarked earlier in general the sequence of the $a_j^n$ in

    $$
    P[x_t | x_{t-1}, x_{t-2}, \ldots, x_{t-n}] = \sum_{j=1}^n a_j^n x_{t-j}
    $$

    does not converge as $n \to \infty$. However, under the roots condition given in the text, the $a_j^n$ do converge. In particular, they converge to $a_j$ of Equation {eq}`eq-60`, so that $\lim_{n\to \infty} a_j^n = a_j$ for all $j = 1, 2, \ldots$
