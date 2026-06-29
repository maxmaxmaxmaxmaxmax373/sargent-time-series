# Signal Extraction Problems[^fn-1]

This section extends to dynamic, serially correlated settings the static
{ref}`signal-extraction problem <sec-10-4>` of {doc}`Chapter X <ch10_regressions>`, where an
agent estimates an unobserved variable from a noisy observation by linear least squares
projection. Here the signal and the noise are themselves stochastic processes, and the
projection runs over their entire histories.

Let $y_t$ be a covariance stationary stochastic process with $m$th order moving average representation

$$
y_t = a(L) u_t,
$$

where

$$
a(L) = \sum_{j=0}^m a_j L^j,
$$

where $u_t$ is a white noise with variance $\sigma_u^2$ that is not necessarily fundamental for $y_t$. Suppose that $x_t$ is the sum of $y_t$ and an orthogonal serially uncorrelated white noise $\eta_t$ with variance $\sigma_\eta^2$, where

$$
x_t = y_t + \eta_t
$$

where

$$
E \eta_t u_{t-s} = 0 \qquad \text{for all } s.
$$

Suppose that an agent observes $\{x_t, x_{t-1}, \ldots\}$ at $t$, and wishes to construct linear least squares forecasts of $x$'s on the basis of this information set. To construct the linear least squares forecast for $x_{t+k}$ given $\{x_t, x_{t-1}, \ldots\}$, one uses the Wiener-Kolmogorov formula {eq}`eq-62`, which requires that a Wold moving average representation $x_t = d(L)\epsilon_t$ be obtained for $x_t$.

To obtain the Wold representation for $x_t$, we simply use the method of Section 16. In particular, the covariance generating function of $x_t$ is

$$
g_x(z) = a(z) a(z^{-1})\sigma_u^2 + \sigma_\eta^2.
$$

We find the zeros of $g_x(z)$, which come in reciprocal pairs, and prepare the factorization

$$
g_x(z) = d(z)d(z^{-1})\sigma_\epsilon^2
$$

where the zeros of $d(z) = (1 - \lambda_1 z)\cdots(1-\lambda_m z)$ do not lie inside the unit circle, where $\sigma_\epsilon^2$ solves

$$
\sigma_\epsilon^2 = \frac{g_x(1)}{d(1)^2}
$$

The Wiener-Kolmogorov formula {eq}`eq-62` can then be used to calculate $P[x_{t+k}|x_t, x_{t-1},\ldots]$.

Moving into a richer class of examples, we now let $y_t$ be a process with mixed moving average, autoregressive representation

$$
y_t = \frac{a(L)}{b(L)}u_t
$$

where $u_t$ is a white noise with variance $\sigma_u^2$, and

$$
\begin{aligned}
a(L) &= (1 - \alpha_1 L)\cdots (1 - \alpha_m L) \\
b(L) &= (1 - \mu_1 L)\cdots(1 - \mu_m L),\quad |\mu_j| < 1
\end{aligned}
$$

where the $\alpha_j$'s can be either side of the unit circle. Suppose that $x_t$ is the sum of $y_t$ and a serially uncorrelated white noise $\eta_t$ with variance $\sigma_\eta^2$,

$$
x_t = y_t + \eta_t
$$

where $E \eta_t u_{t-s} = 0$ for all $s$. Again we desire to find $P[x_{t+k}|x_t, x_{t-1}, \ldots]$, so we need to find a Wold representation for $x_t$. We use the method of Section 17.

The covariance generating function of $x$ is

$$
g_x(z) = \frac{a(z) a(z^{-1})}{b(z) b(z^{-1})}\sigma_u^2 + \sigma_\eta^2.
$$

Taking the right-hand side to a common denominator gives

```{math}
:label: eq-77
g_x(z) = \frac{\sigma_u^2 a(z) a(z^{-1}) + \sigma_\eta^2 b(z) b(z^{-1})}{b(z) b(z^{-1})}.
```

The numerator polynomial is of order $p = \max(n,m)$, and can be factored to be of the form

```{math}
:label: eq-78
\sigma_u^2 a(z) a(z^{-1}) + \sigma_\eta^2 b(z) b(z^{-1}) = \sigma_\epsilon^2 d(z) d(z^{-1})
```

where

$$
d(z) = (1 - \lambda_1 L)\cdots(1 - \lambda_p L),\quad |\lambda_j| \leq 1, \quad j=1,\ldots,p
$$

and where $\sigma_\epsilon^2$ solves

$$
\sigma_\epsilon^2 = \frac{\sigma_u^2 a(1)^2 + \sigma_\eta^2 b(1)^2}{d(1)^2}
$$

The Wold moving average representation for $x_t$ is then

```{math}
:label: eq-79
x_t = \frac{d(L)}{b(L)}\epsilon_t
```

The Wiener-Kolmogorov formula can be applied to {eq}`eq-79`.

A famous application of the preceding analysis is due to Muth (1960). Muth assumed that income $x_t$ is the sum of a first order Markov process $[1/(1-\rho L)]u_t$, $|\rho| < 1$, and an uncorrelated white noise $\eta_t$. The agent's problem was to predict his future income. Setting $a(L)=1$, $b(L) = (1 - \rho L)$, we find that equation {eq}`eq-78` becomes

$$
\sigma_u^2 + \sigma_\eta^2(1 - \rho z)(1 - \rho z^{-1}) = \sigma_\epsilon^2(1 - \lambda_1 z)(1 - \lambda_1 z^{-1}).
$$

The expression on the left can be written

$$
\rho z^{-1} \sigma_\eta^2 \left[-z^2 + \left(\frac{\sigma_u^2}{\sigma_\eta^2 \rho} + \left(\frac{1}{\rho} + \rho\right)\right)z - 1\right].
$$

Applying the quadratic formula, and setting $\lambda_1$ equal to the root that is smaller in absolute value, we have

```{math}
:label: eq-80
\lambda_1 = \frac{1}{2}\left[\left(\frac{\sigma_u^2}{\sigma_\eta^2 \rho}\right) + \left(\frac{1}{\rho} + \rho\right) - \left\{\left[\left(\frac{\sigma_u^2}{\sigma_\eta^2 \rho}\right) + \left(\frac{1}{\rho} + \rho\right)\right]^2 - 4 \right\}^{1/2}\right].
```

The limiting value of $\lambda_1$ as $\rho$ approaches 1 from below is

```{math}
:label: eq-81
\lambda_1 = 1 + \frac{1}{2}\left(\frac{\sigma_u^2}{\sigma_\eta^2}\right) - \left\{\frac{\sigma_u^2}{\sigma_\eta^2}\left(1 + \frac{1}{4}\frac{\sigma_u^2}{\sigma_\eta^2}\right)\right\}^{1/2},
```

which is the expression obtained by Muth (1960). Thus we have that $x_t$ has the first-order moving average, first-order autoregressive representation

$$
x_t = \frac{1 - \lambda_1 L}{1 - \rho L}\epsilon_t,
$$

where $\epsilon_t$ is a fundamental white noise for $x_t$ with variance $\sigma_\epsilon^2$ that solves

$$
\sigma_\epsilon^2 = \frac{\sigma_u^2 + \sigma_\eta^2(1-\rho)^2}{(1 - \lambda_1)^2}.
$$

The result from page 294 now applies with $\beta \equiv \rho$ and $\lambda_1 \equiv -a$. Thus we have

$$
P_t x_{t + k} = [\rho^{k-1}(\rho - \lambda_1)/(1 - \lambda_1 L)]x_t
$$

so that projections of future $x$'s are a geometric average of past $x$'s.

[^fn-1]: Multivariate versions of signal extraction problems require the factorization of matrix polynomials in the lag operator, of the same kind encountered in {ref}`Chapter IX, Section 10 <sec-9-10>`. As in the univariate case, there is an equivalence between the mathematics involved in solving linear quadratic control problems and linear signal extraction problems. In the control literature, this equivalence is called "duality". Hansen and Sargent (1981a) describe practical methods for factoring matrix polynomials in the lag operator which arise from an optimal control or signal extraction problem. Blanchard and Kahn (1980), Whiteman (1983), and Dagli and Taylor (1984) discuss problems in which the matrix polynomial in the lag operator does not come from an optimum control or prediction problem (i.e., the matrix polynomial lacks $b$ symmetry in $L$). These three works contain results on the existence, uniqueness, and computation of solutions of linear systems generated by such nonsymmetric matrix polynomials.
