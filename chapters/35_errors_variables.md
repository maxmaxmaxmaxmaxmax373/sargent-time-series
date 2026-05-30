# Effects of Errors in Variables

We have already considered some aspects of the ways that measurement errors distort the autoregressive or Wold moving average representation of an error ridden variable relative to that of the true variable. For example, the solution of the dynamic signal extraction problem described in Section 18 hinged on deducing a Wold representation for the error-ridden data from knowledge of the covariance generating functions of the true process and the measurement error process. Such a Wold representation was derived by obtaining the spectral density of the error-ridden data, then "factoring" it so that the associated moving average kernel generated a white noise contained in the space of square summable linear combinations of current and lagged values of the error-ridden variable. In interpreting autoregressions and vector autoregressions, one often traces out the "impulse response function," i.e., the coefficients of a Wold moving average representation. (For example, see the "innovation accounting" methods introduced by Sims (1980).) In interpreting such representations where some or all of the data are error-ridden, it is important to remember the signal extraction or spectral factorization problem which those autoregressions are implicitly solving. This section briefly considers the effects of errors in variables on the nature of the projections in Sims's theorem 2. We shall show how measurement errors make error-ridden observations $x$ on a variable $X$ which is not Granger caused by $Y$ appear to be Granger caused by either $Y$ or an error-ridden measure of it, $y$.

Let $(Y_t, X_t)$ be a jointly covariance stationary process of two "true" (accurately measured) variables with covariance generating functions $g_{YX}(z)$, $g_X(z)$, $g_Y(z)$. The two-sided projection of $Y$ on the $X$ process is given by

$$
Y_t=B(L)X_t+U_t
$$

where $EU_tX_{t-s} = 0$ for all $s$, and

$$
B(z) = \frac{g_{YX}(z)}{g_X(z)}.
$$

Now let $(y_t,x_t)$ be a jointly covariance stationary process of error-ridden measured values corresponding to $(Y_t, X_t)$. Assume that $(y_t,x_t)$ is related to $(Y_t,X_t)$ by

$$
y_t=Y_t + \eta_t, \qquad x_t = X_t + \epsilon_t
$$

where $\eta_t$ and $\epsilon_t$ are jointly covariance stationary stochastic processes that satisfy the orthogonality conditions

```{math}
:label: eq-199
E\epsilon_t\eta_{t-s}=EY_t\eta_{t-s} = EX_t\eta_{t-s} = EY_t\epsilon_{t-s} = EX_t\epsilon_{t-s} = 0
```

for all $t$ and $s$. Let the covariance generating functions of $\epsilon$ and $\eta$ be $g_{\epsilon}(z)$ and $g_\eta(z)$. The two-sided projection of $y$ on the $x$ process is given by

$$
y_t=b(L)x_t+u_t
$$

where $Eu_tx_{t-s} = 0$ for all $s$, and

$$
b(z) = \frac{g_{yx}(z)}{g_x(z)}.
$$

But the orthogonality conditions {eq}`eq-199` imply that $g_{yx}(z) = g_{YX}(z)$ and $g_x(z) = g_X(z) + g_\epsilon(z)$. Therefore, we have that

```{math}
:label: eq-200
b(z) = \frac{g_{YX}(z)}{g_X(z) + g_\epsilon(z)} =B(z) \frac{g_{X}(z)}{g_X(z) + g_\epsilon(z)}
```

Defining $w(z) = g_X(z)/(g_X(z) + g_\epsilon(z))$, we have

```{math}
:label: eq-201
b(z) = B(z)w(z)
```

or

```{math}
:label: eq-202
b_k = \sum_{j=-\infty}^\infty w_jB_{k-j}
```

where $w_j$ is given by the inversion formula

$$
w_j = \frac{1}{2 \pi i} \int w(z)z^{-j-1}dz.
$$

Equation {eq}`eq-200`, {eq}`eq-201`, or {eq}`eq-202` shows that in general $b(z) \neq B(z)$, and that $\{b_k\}$ is the convolution of $\{B_k\}$ with a weighting function $\{w_k\}$ that is two-sided and symmetric ($w_k = w_{-k}$) and that in general is of infinite order. Therefore, $b_k$ for a given $k$ in general is a function of $B_h$ for all $h$'s. It follows from these observations, or directly from {eq}`eq-200`, that even if $B(z)$ is one-sided, $b(z)$ will in general be two-sided. By Sims's theorem, this means that even if $Y$ fails to Granger cause $X$, $y$ will Granger cause the error-ridden series $x$. To take a simple example that illustrates this phenomenon, let $g_X(z) =[(1 - \rho z)(1 - \rho z^{-1})]^{-1}$, $|\rho|<1$, $g_{YX}(z)=(1 + az)/(1 - \rho z^{-1})$, $g_\eta(z) = 0$, $g_\epsilon(z) = \sigma^2_\epsilon$. We have that $B(z)=(1 + az)(1 - \rho z)$, so that the projection of $Y$ on $X$ is one-sided, and $X$ fails to Granger cause $Y$. We have that $w(z) = 1/[1 + \sigma^2_\epsilon(1 - \rho z)(1 - \rho z^{-1})] = 1/\lambda_0(1 - \lambda_1z)(1 - \lambda_1z^{-1})$ where $|\lambda_1|<1$ and $\lambda_0(1 - \lambda_1z)(1 - \lambda_1z^{-1}) = 1 + \sigma^2_\epsilon(1 - \rho z)(1 - \rho z^{-1})$. We see that $b(z)$ is two-sided even though $B(z)$ is one-sided.
