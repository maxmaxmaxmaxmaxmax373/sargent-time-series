# Finding a Wold Representation: An $m$th Order Moving Average, $n$th Order Autoregression

The following problem is a useful input into solving an interesting class of "signal extraction" problems.

Consider a covariance stationary process $x_t$ with representation

$$
x_t = \frac{a(L)}{b(L)}u_t
$$

where $u_t$ is a (not necessarily fundamental) white noise and

$$
\begin{aligned}
b(L) &= (1 - \mu_1 L)\cdots(1-\mu_n L),\quad |\mu_j| < 1\\
a(L) &= (1 - \alpha_1 L)\cdots(1 - \alpha_m L).
\end{aligned}
$$

Note that we assume that the zeros of $b(z)$ are outside the unit circle, but those of $a(z)$ are unrestricted. Our problem is to find a Wold moving average representation for $x_t$.

The solution of this problem is simply

```{math}
:label: eq-76
x_t = \frac{d(L)}{b(L)}\epsilon_t
```

where $d(L) = (1 - \lambda_1 L)\cdots(1 - \lambda_m L)$, where $\lambda_1, \ldots, \lambda_m$ are the zeros of $a(z) a(z^{-1})$ that do not lie outside the unit circle, and where $\epsilon_t$ is the fundamental white noise for $x_t$, with variance $\sigma_\epsilon^2$ given by

$$
\sigma_\epsilon^2 = \frac{a(1)^2\sigma_u^2}{d(1)^2}.
$$

In other words, the denominator polynomial $b(L)$ is left unaltered while the methods of the preceding section are applied to factor the numerator polynomial. The reader should convince himself that this method delivers an $\epsilon_t$ process that is a white noise, and that lies in the linear space spanned by $\{x_t, x_{t-1}, \ldots\}$. This can be done by constructing an argument along the lines of the one in the preceding section, by assuming $|\lambda_j| < 1$ for $j = 1, \ldots, m$ and by premultiplying {eq}`eq-76` by $b(L)/d(L)$.
