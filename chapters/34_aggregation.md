# Aggregation over Time

This section briefly introduces Sims's (1971) framework for studying the effects of aggregation over time. The idea is to assume that there is an underlying stochastic process that is evolving in continuous time, on which the economist has observations at discrete points in time. Sims studied the question of how closely the projection of $y$ on $x$ in discrete time approximates the corresponding projection in continuous time. Discussing Sims's work permits us briefly to introduce the concept of a continuous time stochastic process.

We now suppose that $y(t)$, $x(t)$ is a continuous time jointly covariance stationary stochastic process. The linear structure of such a process is characterized by the covariance functions

```{math}
:label: eq-184
\begin{aligned}
c_x(\tau) &= Ex_tx_{t-\tau} \\
c_y(\tau) &= Ey_ty_{t-\tau} \\
c_{yx}(\tau) &= Ey_tx_{t-\tau}
\end{aligned}
```

where now $\tau \in \mathbf{R}$, where $\mathbf{R}$ is the set of real numbers. To be the matrix covariance function of a real-valued covariance stationary process it is necessary and sufficient that the matrix function

$$
\begin{bmatrix}
    c_y(\tau) & c_{yx}(\tau) \\
    c_{xy}(\tau) & c_x(\tau) \\
\end{bmatrix}
$$

must be nonnegative definite. This is equivalent with the requirement that for any finite $n$, any finite index set of times $t_1,\ldots,t_n$, and any sequence of $(2\times 1)$ real-valued vectors $h_1, \ldots ,h_n$ the covariance function {eq}`eq-184` must imply that

$$
Ez^2 \geq 0
$$

where

$$
z \equiv \sum_{j=1}^n h_j^T \begin{bmatrix} y_{t_j} \\ x_{t_j} \end{bmatrix}.
$$

The projection of $y_t$ on the entire $x$ process is defined by

```{math}
:label: eq-185
y_t = \int_{-\infty}^\infty b(s)x(t-s)ds + u(t)
```

where $u(t)$ is a least-squares residual that satisfies the orthogonality conditions $Eu(t)x(t - s) = 0$ for all $s \in \mathbf{R}$. Multiplying {eq}`eq-185` by $x(t - k)$ and taking expectations gives the normal equation

```{math}
:label: eq-186
c_{yx}(k) = \int_{-\infty}^\infty b(s)c_x(k-s)ds
```

for all $k \in \mathbf{R}$. The normal equation {eq}`eq-186` determines the function $b(s)$. Now suppose that one has observations on $(y_t, x_t)$ only at the discrete points of time $t \in I$, where $I$ is the set of integers $\{0, \pm 1, \pm 2, \ldots\}$. Thus, we view the discrete time process as representing a sampled version of the underlying continuous time process. Sims (1971) posed and answered the following question. Consider the discrete time projections of $y_t$ on all past, present, and future $x_t$'s:

```{math}
:label: eq-187
y_t = \sum_{j=-\infty}^\infty B_jx_{t-j} +\epsilon_t
```

where $E\epsilon_t x_{t-j}=0$ for $j\in I$. What is the relationship between the sequence $\{B_j\}_{-\infty}^\infty$ in {eq}`eq-187` and the function $(b(s), s \in \mathbf{R})$ in {eq}`eq-185`? Under what circumstances does the $B_j$ sequence "resemble $b(s)$," i.e., under what circumstances is $B_j$ close to being a sampled version of $b(s)$, so that $B_j \simeq b(j)$, $j \in I$? To answer this question, Sims derived a formula expressing $B_j$ as a particular kind of average of $b(s)$.

Multiplying {eq}`eq-187` by $x_{t-j}$ for integer $\tau$ and taking expectations of both sides gives the normal equations

```{math}
:label: eq-188
c_{yx}(\tau) = \sum_{j=-\infty}^\infty B_jc_x(\tau -j),\quad \tau = 0, \pm 1, \pm2,\ldots
```

Since {eq}`eq-186` holds for all $k \in \mathbf{R}$, {eq}`eq-188` and {eq}`eq-186` together imply that

```{math}
:label: eq-189
\sum_{j=-\infty}^\infty B_jc_x(\tau -j) = \int_{-\infty}^\infty b(s)c_x(\tau-s)ds,\quad \tau = 0, \pm 1, \pm2,\ldots
```

Multiplying both sides by $z^\tau$ and summing over $\tau$ gives

```{math}
:label: eq-190
B(z)g_x(z)=m(z)
```

where

$$
B(z) = \sum_{\tau=-\infty}^\infty B_\tau z^\tau, \qquad g_x(z) = \sum_{\tau=-\infty}^\infty c_x(\tau) z^\tau,
$$

and

$$
m(z)=\sum_{\tau=-\infty}^\infty z^\tau \int b(s)c_x(\tau-s)ds.
$$

Now let $g_x(z)^{-1} = \sum_{j=-\infty}^\infty \gamma_jz^j \equiv \gamma_x(z)$.[^fn-agg-1] Operating on both sides of {eq}`eq-190` gives

$$
B(z)=\gamma_x(z)m(z)
$$

or

$$
B_j=\sum_{k=-\infty}^\infty\gamma_x(k)\int_{-\infty}^\infty b(s)c_x(j-k-s)ds
$$

or

```{math}
:label: eq-191
B_j=\int_{-\infty}^\infty b(s)\sum_{k=-\infty}^\infty\gamma_x(k) c_x(j-k-s)ds
```

Equation {eq}`eq-191` is a version of Sims's formula. Define

```{math}
:label: eq-192
r_{x}(j-s) = \sum_{k=-\infty}^\infty \gamma_x(k)c_x(j-k-s).
```

Then from {eq}`eq-191`, we have

```{math}
:label: eq-193
B_j=\int_{-\infty}^\infty b(s)r_x(j-s)ds,\quad j=0, \pm 1, \pm2, \ldots
```

Formula {eq}`eq-193` has an interpretation as a version of Theil's omitted variables theorem. Since

```{math}
:label: eq-194
P[y_t|\{x_{t-k}\}_{k=-\infty}^\infty] = \sum_{k=-\infty}^\infty B_kx_{t-k},
```

use {eq}`eq-185` and the orthogonality of $u(t)$ to $x(t - s)$ for all $s$ to write

$$
P[y_t|\{x_{t-k}\}_{k=-\infty}^\infty] = \int_{-\infty}^\infty b(s)P[x(t-s)|\{x_{t-k}\}_{k=-\infty}^\infty] ds
$$

or from {eq}`eq-194`

```{math}
:label: eq-195
\sum_{k=-\infty}^\infty B_kx_{t-k} = \int_{-\infty}^\infty b(s) \sum_{k=-\infty}^\infty \delta_k^s x_{t-k} ds
```

where $P[x(t-s)|\{x_{t-k}\}_{k=-\infty}^\infty] = \sum_{k=-\infty}^\infty \delta_k^s x_{t-k}$, $s \in \mathbf{R}$. The $\{\delta_k^s\}_{k=-\infty}^\infty$ sequence for each $s \in \mathbf{R}$ is determined from the orthogonality conditions

$$
Ex(t - s)x(t - j) = \sum_{k=-\infty}^\infty \delta_k^s Ex_{t-j}x_{t-j}
$$

or

```{math}
:label: eq-196
c_x(s-j)=\sum_{k=-\infty}^\infty \delta_k^s c_x(j-k) \quad \text{for all integer } j.
```

Recalling that we have defined $\gamma_x(k)$ to be the inverse under convolution of the sequence $c_x(k)$, convoluting each side of {eq}`eq-196` with the sequence $\gamma_x(k)$ gives

```{math}
:label: eq-197
\delta_h^s=\sum_{k=-\infty}^\infty \gamma_x(k)c_x(h-k-s),
```

which shows that $\delta_h^s = r_x(h-s)$. Now exchanging order of integration and summation in {eq}`eq-195` gives

$$
\sum_{k=-\infty}^\infty B_kx_{t-k} = \sum_{k=-\infty}^\infty\left[\int_{-\infty}^\infty b(s)\delta_k^sds\right]x_{t-k}
$$

This implies that

```{math}
:label: eq-198
B_k = \int_{-\infty}^\infty b(s) \delta_k^sds
```

which is equivalent with {eq}`eq-193` because $\delta_k^s=r_x(k - s)$. Equation {eq}`eq-198` expresses $B_k$ as the convolution of the function $b(s)$ with the function $\{\delta_k^s\}$, where $\delta_k^s$ is the coefficient on $x_k$ in the projection of $x(t - s)$ on $\{x_j\}_{j=-\infty}^\infty$.

Equation {eq}`eq-193` or {eq}`eq-198` shows how $B_k$ depends on $b(s)$ for all $s$, not just $s$ close to $k$. Sims (1971) describes examples in which the $r_x(k - s)$ function is such that $B_k$ is very different from $b(k)$. In general, {eq}`eq-193` reveals that even if $b(s)$ is one-sided on the past ($b(s) = 0$ for $s < 0$), $B_k \neq 0$ for $k < 0$. This means that if $y$ fails to Granger cause $x$ in continuous time, in general $y$ in discrete time will Granger cause $x$; how strongly $y$ in discrete time Granger causes $x$ depends on the function $r_x(s)$. Sims's work was extended to a multivariate setting by Geweke (1978). Hansen and Sargent (1982) apply an approach related to Sims's in order to study the connection between vector autoregressions in continuous time and their discrete time counterparts. Methods for estimating continuous time models from discrete time observations are described by Bergstrom (1983) and Hansen and Sargent (1981b).

[^fn-agg-1]: The sequence $\{\gamma_j\}$ is called the "inverse under convolution" of the sequence $\{c_x(\tau)\}$.
