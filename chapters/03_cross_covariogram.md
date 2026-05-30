# The Cross Covariogram

Suppose we have two wide-sense stationary stochastic processes $y_t$ and $x_t$.
The processes are said to be *jointly wide-sense stationary* if the cross covariance
$E(y_t - Ey_t)(x_{t-k} - Ex_{t-k})$ depends only on $k$ and not on $t$.
The *cross covariogram* is this list of covariances as a function of $k$:

$$
c_{yx}(k) = E(y_t - Ey_t)(x_{t-k} - Ex_{t-k}).
$$

Now suppose $y_t$ and $x_t$ can be expressed as (possibly two-sided) distributed lags
of a single white-noise process $\epsilon_t$:

$$
y_t = B(L)\epsilon_t, \qquad x_t = D(L)\epsilon_t,
$$

where

$$
B(L) = \sum_{j=-\infty}^{\infty} b_j L^j, \quad
D(L) = \sum_{j=-\infty}^{\infty} d_j L^j, \quad
\sum_{j=-\infty}^{\infty} b_j^2 < \infty, \quad
\sum_{j=-\infty}^{\infty} d_j^2 < \infty.
$$

Since $E\epsilon_t = 0$,

$$
c_{yx}(k) = Ey_t x_{t-k} = \sigma_\epsilon^2 \sum_{j=-\infty}^{\infty} b_j d_{j-k}.
$$

The **cross-covariance generating function** $g_{yx}(z)$ is defined by

$$
g_{yx}(z) = \sum_{k=-\infty}^{\infty} c_{yx}(k)\, z^k.
$$

In the present case, we have

$$
g_{yx}(z) = \sigma_\epsilon^2 \sum_{k=-\infty}^{\infty} \sum_{j=-\infty}^{\infty} b_j d_{j-k}\, z^k.
$$

Letting $h = j-k$ so that $k = j-h$, we have

$$
g_{yx}(z) = \sigma_\epsilon^2 \sum_{j=-\infty}^{\infty}\sum_{h=-\infty}^{\infty} b_j d_h\, z^{j-h}
= \sigma_\epsilon^2 \sum_{j=-\infty}^{\infty} b_j z^j \sum_{h=-\infty}^{\infty} d_h z^{-h}.
$$

The last equation gives

```{math}
:label: eq-16
g_{yx}(z) = \sigma_\epsilon^2\, B(z)\, D(z^{-1}).
```

This is a counterpart to {eq}`eq-5` and includes it as a special case.

## A More General Bivariate System

Suppose instead we have the system

```{math}
:label: eq-17
y_t = A(L)\epsilon_t + B(L)u_t, \qquad x_t = C(L)\epsilon_t + D(L)u_t,
```

where $\epsilon_t$ and $u_t$ are two mutually uncorrelated white-noise processes with
variances $\sigma_\epsilon^2$ and $\sigma_u^2$, and $Eu_t\epsilon_{t-k} = 0$ for all $k$.
By calculations analogous to those above, the cross-covariance generating function is

```{math}
:label: eq-18
g_{yx}(z) = \sigma_\epsilon^2\, A(z)\, C(z^{-1}) + \sigma_u^2\, B(z)\, D(z^{-1}).
```

The representation {eq}`eq-17` is in fact very general—it includes all jointly
wide-sense stationary, indeterministic bivariate processes.[^fn-cross-1]

[^fn-cross-1]: Namely, all jointly wide-sense stationary, indeterministic processes.

## Symmetry Relations

We also define the cross-covariance going the other way:

$$
c_{xy}(k) = E(x_t - Ex_t)(y_{t-k} - Ey_{t-k}) = c_{yx}(-k).
$$

Correspondingly, we have

$$
g_{yx}(z) = \sum_{k=-\infty}^{\infty} c_{yx}(k)\, z^k
          = \sum_{h=-\infty}^{\infty} c_{xy}(h)\, z^{-h},
$$

where the second equality follows by substituting $k = -h$ and using
$c_{xy}(h) = c_{yx}(-h)$.

The particular system {eq}`eq-17` implies the symmetric counterpart:

$$
g_{xy}(z) = \sigma_\epsilon^2\, A(z^{-1})\, C(z) + \sigma_u^2\, B(z^{-1})\, D(z).
$$
