# Finding a Wold Representation: $m$th Order Moving Average

More generally, suppose that we are given a process $x_t$ with finite order moving average representation

$$
x_t = a(L)u_t
$$

where

$$
a(L) = \sum_{j=0}^m a_j L^j,
$$

and where $u_t$ is a white noise that is not fundamental for $x_t$, i.e., $u_t \neq x_t - P[x_t|x_{t-1}, x_{t-2},\ldots]$. For convenience, we assume that the zeros of the polynomial $a(z)$ are distinct. The covariance generating function of $x_t$ is given by

```{math}
:label: eq-67
g_x(z) = a(z) a(z^{-1})\sigma_u^2.
```

We know that this process also possesses a Wold moving average representation

```{math}
:label: eq-68
x_t = d(L)\epsilon_t, \qquad d_0 = 1,
```

where $\epsilon_t = x_t - P[x_t|x_{t-1}, x_{t-2},\ldots]$. The condition that $\epsilon_t$ lie in the linear space spanned by $x_{t-1}$ is equivalent with the condition that the zeros of $d(z)$ not be less than unity in absolute value. To see this heuristically, represent $d(L)$ as

```{math}
:label: eq-69
d(L) = (1 - \lambda_1 L)\cdots(1-\lambda_m L)
```

where $\lambda_j$ is the reciprocal of the $j$th zero of $d(z)$. Represent $d(L)^{-1}$ as

```{math}
:label: eq-70
d(L)^{-1} = \sum_{j=1}^m \frac{A_j}{(1 - \lambda_j L)}
```

Then {eq}`eq-68` and {eq}`eq-70` imply

```{math}
:label: eq-71
\epsilon_t = \sum_{j=1}^m A_j \sum_{k=0}^{\infty} \lambda_j^k x_{t-k}.
```

The geometric sums on the right side fail to converge if $|\lambda_j| > 1$. Imposing that $|\lambda_j| \leq 1$ for $j = 1, \ldots, m$ is necessary and sufficient for $\epsilon_t$ to lie in the space spanned by current and lagged $x$'s. (This argument fails to reveal why when $|\lambda_j| = 1$ for some $j$, $\epsilon_t$ lies in the space spanned by current and lagged $x_t$'s. When $|\lambda_j| = 1$ for some $j$, then although $x_t$ possesses a moving average representation, it possesses no autoregressive representation. In this case, although $\epsilon_t$ is in the closure of the linear space spanned by $\{x_t, x_{t-1}, \ldots\}$, it can be expressed in the form $\sum_{j=0}^\infty w_j x_{t-j}$) for no square summable sequence of $w_j$'s, but only as the limit of a sequence of such expressions.)

To find a fundamental moving average representation, we note that {eq}`eq-68` and {eq}`eq-69` imply that $g_x(z) = \sigma_\epsilon^2(1 - \lambda_1 z)\cdots(1-\lambda_m z)(1 - \lambda_1 z^{-1})\cdots(1 - \lambda_m z^{-1})$. Then we equate this to $g_x(z)$ given by {eq}`eq-67` to get

$$
\sigma_u^2 a(z) a(z^{-1}) = \sigma_\epsilon^2 d(z) d(z^{-1})
$$

or

```{math}
:label: eq-72
\sigma_u^2 a(z) a(z^{-1}) = \sigma_\epsilon^2 (1 - \lambda_1 z)\cdots(1 - \lambda_m z)(1 - \lambda_1 z^{-1})\cdots (1 - \lambda_m z^{-1}).
```

Equation {eq}`eq-72` asserts that $d(z)d(z^{-1})\sigma_\epsilon^2$ is a *symmetric factorization* of $\sigma_u^2 a(z)a(z^{-1})$ in which the zeros of $\sigma_u^2 a(z)a(z^{-1})$ that are not inside the unit circle are placed in $d(z)$, while those that are not outside are put into $d(z^{-1})$, (Note that since $\sigma_u^2 a(z)a(z^{-1})$ evaluated at $z_0$ equals its value at the zeros of $\sigma_u^2 a(z)a(z^{-1})$ come in reciprocal pairs. Thus, it is reminiscent of the characteristic polynomial of an Euler equation in the undiscounted case, when $b = 1$.)

The preceding tells us how to achieve a Wold moving average representation for a finite-order moving average process. First, find the zeros of $g_x(z)$. By the reciprocal pairs properties of these roots, half will not be outside the unit circle, while half will not be inside the unit circle. (Excuse the cumbersome wording, which is designed to cover the case of roots on the unit circle.) Let $\lambda_1, \ldots \lambda_m$ be the roots that are not outside the unit circle, Set

```{math}
:label: eq-73
d(L) = (1-\lambda_1 L)\cdots(1-\lambda_m L).
```

Then to find $\sigma_\epsilon^2$, solve equation {eq}`eq-72` at $z=1$ to get

```{math}
:label: eq-74
\sigma_\epsilon^2 = \frac{g_x(1)}{d(1)^2}
```

Hansen and Sargent (1990, p. 102) give a quick but equivalent method of finding $d(L)$. Given $\sigma_u^2 a(z) a(z^{-1})$, the problem is to find $\sigma_\epsilon^2 d(z) d(z^{-1})$, where the zeros of $d(z)$ are not inside the unit circle. Let $z_1, \ldots z_k$ be the zeros of $a(z)$ that are inside the unit circle, where $0 \leq k \leq m$. Then $d(z)$ satisfies

```{math}
:label: eq-75
\frac{\sigma_{\epsilon}}{\sigma_u}d(z) = a(z)\prod_{j=1}^k \frac{1 - z_j z}{z - z_j}
```

For example, suppose $x_t = (1 + 2 L)u_t$, where $\sigma_u^2 = 1$ and where $u_t$ is a white noise. Then application of {eq}`eq-75` gives the Wold moving average representation $x_t = (1 + (1/2)L)\epsilon_t$, with $\sigma_{\epsilon} = 2$, where $\epsilon_t$ is a fundamental white noise for $\{x_t\}$.

The factorizations of this section locate the zeros of $g_x(z)$ and sort them by the unit circle by hand. When the process is summarized instead by its *spectral density*, the same fundamental factor $d(L)$ and innovation variance $\sigma_\epsilon^2$ can be recovered numerically, in a handful of FFTs, by **Whittle's spectral factorization** of {doc}`Representation Theory <13_representation_theory>`. Applied to the spectral density $g_x(e^{-i\omega}) = |1 + 2e^{-i\omega}|^2$ of the example just worked, Whittle's method returns exactly $d = (1, \tfrac12, 0, \dots)$ and $\sigma_\epsilon^2 = 4$ — selecting the fundamental factor automatically.
