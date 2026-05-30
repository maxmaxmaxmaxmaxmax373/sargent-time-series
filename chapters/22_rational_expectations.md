# Some Applications to Rational Expectations Models

Let us return to the example of Cagan's portfolio balance schedule, only now we assume that $m_t$ is a covariance stationary stochastic process and the log of the price level now expected for next period is the linear least squares projection of $p_{t+1}$ on information available at time $t$. We then have the difference equation

```{math}
:label: eq-93
m_t - p_t = \alpha P_t p_{t+1} - \alpha p_t,\quad \alpha < 0
```

where $P_t p_{t+1}$ is the linear least squares forecast of $p_{t+1}$ given information available at time $t$. This difference equation can be rewritten as

$$
p_t = \left(\frac{-\alpha}{1-\alpha}\right)P_t p_{t+1} - \left(\frac{1}{1-\alpha}\right)m_t
$$

or

$$
p_t = \lambda P_t p_{t+1} + (1-\lambda)m_t
$$

where $\lambda = -\alpha/(1-\alpha)$, which implies that $0 < \lambda < 1$ since $\alpha < 0$. The stationary solution of the above difference equation obeys[^fn-1]

```{math}
:label: eq-94
p_t = (1 - \lambda)\sum_{j=0}^\infty \lambda^j P_t m_{t+j}
```

Let us assume that $m_t$ has the autoregressive representation

$$
a(L)m_t = \epsilon_t
$$

where $\epsilon_t$ is fundamental for $m$, and $a(L) = 1 - a_1 L - \ldots - a_r L^r$. Then from formula {eq}`eq-90` of the preceding section we have that {eq}`eq-94` implies

```{math}
:label: eq-95
p_t = (1-\lambda)a(\lambda)^{-1}\left[1 + \sum_{j=1}^{r-1}\left(\sum_{k=j+1}^r \lambda^{k-j} a_k \right) L^j \right]m_t
```

```{math}
:label: eq-96
a(L)m_t = \epsilon_t.
```

These two equations express how the stochastic process for $p_t$ depends on $m_t, m_{t-1}, \ldots, m_{t-r+1}$ via coefficients that partly reflect the stochastic process {eq}`eq-96` that governs $m_t$. As an examples, we set $a(L) = 1 - a_1 L - a_2 L^2 - a_3 L^3$. Then {eq}`eq-95` and {eq}`eq-96` become

$$
p_t = (1-\lambda)(1- a_1 \lambda - a_2 \lambda^2 - a_3 \lambda^3)^{-1}[1 + (a_2 \lambda + a_3 \lambda^2)L + (a_3\lambda)L^2]m_t
$$

$$
m_t = a_1 m_{t-1} + a_2 m_{t-2} + a_3 m_{t-3} + \epsilon_t
$$

Let us reconsider the supply-demand example of Chapter IX where $x_t$ is now a covariance stationary, indeterministic random process with mean zero and autoregressive representation $a(L)x_t = \epsilon_t$, where $\epsilon_t$ is a fundamental white noise for $x_t$. Our system is naturally modified to become

$$
\begin{aligned}
C_t &= -\beta p_t, & \beta > 0 \\
Y_t &= \gamma P_{t-1} p_t + x_t & \gamma > 0 \\
I_t &= \alpha(P_t p_{t+1} - p_t), & \alpha > 0 \\
Y_t &= C_t + I_t - I_{t-1},
\end{aligned}
$$

where $Y_t$ is production, $C_t$ demand for consumption, and $I_t$ holdings of inventories. Substituting the first three equations into the fourth gives

```{math}
:label: eq-97
(\gamma + \alpha)P_{t-1}p_t + (\alpha + \beta)p_t = \alpha p_{t-1} - x_t.
```

Taking projections of both sides against information available at time $t - 1$ gives

$$
\alpha P_{t-1}p_{t+1} - (\gamma + \beta + 2\alpha)P_{t-1}p_t + \alpha P_{t-1}p_{t-1} = P_{t-1}x_t
$$

or

$$
(B^{-1} - \phi + B)P_{t-1}p_t = \alpha^{-1} P_{t-1} x_t
$$

where

$$
B^{-1}P_{t-1}z_t \equiv P_{t-1} z_{t+1}, \quad B P_{t-1}z_t \equiv P_{t-1} z_{t-1},
$$

and where

$$
\phi = ((\beta + \gamma)/\alpha) + 2 > 0.
$$

Multiplying by $B$ gives

```{math}
:label: eq-98
\begin{aligned}
(1 - \phi B + B^2)P_{t-1}p_t &= \alpha^{-1}P_{t-1}x_{t-1} \\
(1 - \lambda^{-1}B)(1 - \lambda B)P_{t-1}p_t &= \alpha^{-1}P_{t-1}x_{t-1}
\end{aligned}
```

where $|\lambda| < 1$ satisfies $\lambda + \lambda^{-1} = \phi$. To ensure covariance stationarity of the solution, we shall insist that all lag distributions be square summable. Operating on both sides of {eq}`eq-98` with the forward inverse of $(1-\lambda^{-1} B)$ gives

$$
(1 - \lambda B)P_{t-1}p_t = \frac{- \lambda \alpha^{-1}}{1-\lambda B^{-1}}P_{t-1}x_t
$$

or

$$
P_{t-1}p_t - \lambda p_{t-1} = \frac{-\lambda}{\alpha}\sum_{i=0}^{\infty}\lambda^i P_{t-1} x_{t+i}.
$$

Substituting this solution for $P_{t-1}p_t$ into {eq}`eq-97` gives

```{math}
:label: eq-99
p_t = \lambda p_{t-1} + \frac{1}{\alpha + \beta - \alpha \lambda}\left[\alpha^{-1} \lambda(\gamma + \alpha)\sum_{i=0}^{\infty} \lambda^i P_{t-1}x_{t+i} - \sum_{i=0}^{\infty} \lambda^i P_t x_{t+i}\right]
```

We have assumed that $x_t$ has the autoregressive representation $a(L)x_t = \epsilon_t$. Now by using methods similar to those used to derive {eq}`eq-90`, it can be established that

$$
P_{t-1}\sum_{j=0}^\infty \lambda^j x_{t+j} = \left(\frac{L^{-1}I - L^{-1} a(\lambda)^{-1}a(L)}{1-\lambda L^{-1}}\right)x_{t-1}
$$

Substituting this and {eq}`eq-90` into {eq}`eq-99` we have the following formula for the equilibrium stochastic process for price $p_t$ as a function of the $x_t$ process

$$
p_t = \lambda p_{t-1} + \frac{1}{\alpha + \beta - \alpha \lambda}\left\{\alpha^{-1} \lambda(\gamma + \alpha)\left[\frac{L^{-1}I - L^{-1}a(\lambda)^{-1}a(L)}{1-\lambda L^{-1}}\right]x_{t-1} - \left(\frac{1 - \lambda a(\lambda)^{-1}a(L) L^{-1}}{1-\lambda L^{-1}}\right)x_t\right\}
$$

$$
a(L)x_t = \epsilon_t
$$

This is the solution to the stochastic difference equation {eq}`eq-97` which expresses $p_t$, as a function of current and lagged $x$'s and $p$'s, and which gives a covariance stationary process for $p_t$.

[^fn-1]: in this and the following example we set "transient terms" of the form $c\lambda^t$ to zero because we are interested in obtaining solutions that are covariance stationary processes.
