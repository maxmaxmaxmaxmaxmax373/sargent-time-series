# Optimal Filtering Formula

It is convenient to have a formula for the projection of a random variable $y_t$ against current and past values of a covariance stationary, indeterministic random process $x_t$. We assume that $y_t$ and $x_t$ have means of zero and are jointly covariance stationary, indeterministic processes. That is, we seek the $h_j$ that characterize the one-sided projection

```{math}
:label: eq-111
y_t = \sum_{j=0}^\infty h_j x_{t-j} + u_t
```

where $E x_{t-j} u_t = 0$ for all $j \geq 0$. First, suppose that $x_t$ has the moving average representation

$$
x_t = d(L)\epsilon_t,\qquad d(L) = \sum_{j=0}^\infty d_j L^j
$$

where $\{\epsilon_t\}$ is a serially uncorrelated process of innovations in $x$, i.e., $\epsilon_t$ is fundamental for $x$. As an intermediate step,[^fn-filt-1] think of projecting $y_t$ on current and past $\epsilon$'s:

```{math}
:label: eq-112
y_t = \sum_{j=0}^\infty \phi_j\epsilon_{t - j} + u_t
```

where $E u_t \epsilon_{t-j} = 0$ for all $j \geq 0$. We assume that $x_t$ has both a moving average and an autoregressive representation, so that it is easy to see that $\{ x_t, x_{t-1},\ldots\}$ and $\{\epsilon_t, \epsilon_{t-1},\ldots\}$ span the same space. For this reason, $u_t$ in {eq}`eq-111` equals $u_t$ in {eq}`eq-112`. Since the $\epsilon$'s form an orthogonal process, we have that the $\phi_j$ are the simple least squares coefficients:

$$
\phi_j = E y_t \epsilon_{t-j}/E \epsilon_t^2 = E y_t \epsilon_{t-j}/\sigma^2
$$

where $\sigma^2 = E \epsilon_t^2$. Thus we can write

```{math}
:label: eq-113
\phi(L) = \sum_{j=0}^\infty \phi_j L^j, \qquad \phi(L) = \sigma^{-2}[g_{y \epsilon}(L)]_{+}
```

where $[\,]_{+}$ again means "ignore negative powers of $L$" and $g_{y \epsilon}(L)$ is the cross-covariance generating function

$$
g_{y \epsilon}(L) = \sum_{k = - \infty}^{\infty} E (y_t \epsilon_{t-k})L^k.
$$

We can relate $g_{y \epsilon}(L)$ to the cross-covariance generating function $g_{y x}(L)$ as follows:

$$
\begin{aligned}
g_{y x}(z) &= \sum_{k}(E y_t x_{t-k})z^k \\
&= \sum_k \left(E y_t d(L) \epsilon_{t-k}\right)z^k \\
&= \sum_k \left(E y_t(d_0 \epsilon_{t-k} + d_1 \epsilon_{t-k-1} + \cdots)\right)z^k \\
&= d_0 \sum_k (E y_t \epsilon_{t-k})z^k + d_1 \sum_k (E y_t \epsilon_{t-k-1})z^k + d_2 \sum_k (E y_t \epsilon_{t-k-2})z^k + \cdots \\
&= d_0 g_{y \epsilon}(z) + d_1 z^{-1} g_{y \epsilon}(z) + d_2 z^{-2} g_{y \epsilon}(z) + \cdots \\
&= d(z^{-1})g_{y \epsilon}(z).
\end{aligned}
$$

Thus we have $g_{y \epsilon}(z) = g_{y x}(z)/d(z^{-1})$. Substituting this into {eq}`eq-113`, we obtain

$$
\phi(L) = \frac{1}{\sigma^2}\left(\frac{g_{y x}(L)}{d(L^{-1})}\right)_{+}.
$$

So we have

```{math}
:label: eq-114
\begin{aligned}
y_t &= \frac{1}{\sigma^2}\left(\frac{g_{y x}(L)}{d(L^{-1})}\right)_{+}\epsilon_t + u_t \\
&= \frac{1}{\sigma^2}\left(\frac{g_{y x}(L)}{d(L^{-1})}\right)_{+} \frac{1}{d(L)}x_t + u_t,
\end{aligned}
```

so that in {eq}`eq-111` we have

```{math}
:label: eq-115
h(L) = \frac{1}{\sigma^2}\left(\frac{g_{y x}(L)}{d(L^{-1})}\right)_{+} \frac{1}{d(L)}.
```

A classic application of this formula is due to Muth (1960). Suppose that income evolves according to $x_t = y_t + \epsilon_t$, where $y_t = \rho y_{t-1} + u_t$, $|\rho| < 1$, and where $u_t$ and $\epsilon_t$ are mutually orthogonal at all lags and serially uncorrelated. Here $x_t$ is measured income, while $y_t$ is "systematic" or permanent income. The consumer only "sees" $x_t, x_{t-1}, \ldots$ and desires to estimate systematic income $y_t$ by a linear function of $x_t, x_{t-1}, \ldots$. The consumer is assumed to know all the relevant moments. This problem can be solved quickly using formula {eq}`eq-115`, and the reader is invited to do so.

[^fn-filt-1]: This is the method that Kolmogorov used to derive the formula we are after. See Whittle (1983, p. 42).
