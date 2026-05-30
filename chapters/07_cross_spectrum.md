# The Cross Spectrum

An alternative representation of the cross covariogram is provided by the cross
spectrum. Recall that the cross-covariance generating function between the jointly
stationary process $y$ and $x$ is defined by

$$
g_{yx}(z) = \sum_{k=-\infty}^{\infty} c_{yx}(k)\, z^k.
$$

If we evaluate $g_{yx}(z)$ at the value $e^{-i\omega}$, we have the cross spectrum

$$
g_{yx}(e^{-i\omega}) = \sum_{k=-\infty}^{\infty} c_{yx}(k)\, e^{-i\omega k}.
$$

Viewed as a function of angular frequency $\omega$, $g_{yx}(e^{-i\omega})$ is called
the cross spectrum between $y$ and $x$.

The cross spectrum is of course a cross-covariance generating function. Given an
expression for $g_{yx}(e^{-i\omega})$, it is possible to recover the cross covariances
from the inversion formula

$$
c_{yx}(k) = \frac{1}{2\pi}\int_{-\pi}^{\pi} g_{yx}(e^{-i\omega})\, e^{i\omega k}\, d\omega.
$$

The validity of this inversion formula can be checked by following calculations
analogous to those used to verify the inversion formula for the spectrum.

Unlike the spectrum, the cross spectrum is in general a complex number at each
frequency, this being a consequence of the fact that $c_{yx}(k)$ is in general
*not* symmetric ($c_{yx}(k)$ does *not* in general equal $c_{yx}(-k)$). In place of
the symmetry property we have the readily verified property

```{math}
:label: eq-43
g_{xy}(e^{-i\omega}) = \overline{g_{yx}(e^{-i\omega})} = g_{yx}(e^{i\omega}),
```

where the bar denotes complex conjugation and

$$
g_{xy}(e^{-i\omega}) = \sum_{k=-\infty}^{\infty} c_{xy}(k)\, e^{-i\omega k}
$$

and $c_{xy}(k) = Ex_t y_{t-k}$. Notice that $c_{xy}(k) = c_{yx}(-k)$.

## Cross Spectrum and Linear Projections

Suppose that the stationary stochastic process $y_t$ is related to the stochastic
processes $x_t$ and $\epsilon_t$ by

```{math}
:label: eq-44
y_t = \sum_{j=-\infty}^{\infty} h_j x_{t-j} + \epsilon_t,
```

where $E\epsilon_t = Ex_t = 0$, and $E\epsilon_t x_{t-s} = 0$ for all $s$, an
orthogonality condition that characterizes $\sum h_j x_{t-j}$ as the projection of
$y_t$ on the space spanned by $\{x_{t-\infty}, \ldots, x_0, \ldots, x_{t+\infty}\}$.
Then we have already seen that the spectrum of $y$ satisfies

$$
g_y(e^{-i\omega}) = \left|h(e^{-i\omega})\right|^2 g_x(e^{-i\omega}) + g_\epsilon(e^{-i\omega}),
$$

where

$$
h(e^{-i\omega}) = \sum_{j=-\infty}^{\infty} h_j\, e^{-i\omega j}.
$$

To find the cross spectrum between $y$ and $x$, first use {eq}`eq-44` to calculate
the $k$th lagged covariance as

$$
Ey_t x_{t-k} = \sum_{j=-\infty}^{\infty} h_j E(x_{t-j} x_{t-k}), \qquad
c_{yx}(k) = \sum_{j=-\infty}^{\infty} h_j\, c_x(k-j).
$$

Thus the cross covariogram between $y$ and $x$ is the convolution of the sequence
$\{h_j\}$ with the sequence $\{c_x(j)\}$. From the convolution property we immediately
have

$$
g_{yx}(e^{-i\omega}) = h(e^{-i\omega})\, g_x(e^{-i\omega})
$$

since the Fourier transform of a convolution of two sequences is the product of the
Fourier transforms of the two sequences. That is, taking the Fourier transforms on
each side (i.e., multiplying by $e^{-i\omega k}$ and summing over $k$) gives

$$
\sum_{k=-\infty}^{\infty} c_{yx}(k)\, e^{-i\omega k}
= \sum_{j=-\infty}^{\infty}\sum_{k=-\infty}^{\infty} h_j\, c_x(k-j)\, e^{-i\omega k}.
$$

Noting that $e^{-i\omega k} = e^{-i\omega(k-j)}\, e^{-i\omega j}$, the above can be
written as

$$
g_{yx}(e^{-i\omega})
= \sum_{j=-\infty}^{\infty} h_j\, e^{-i\omega j}
  \sum_{k=-\infty}^{\infty} c_x(k-j)\, e^{-i\omega(k-j)},
$$

or

```{math}
:label: eq-45
g_{yx}(e^{-i\omega}) = h(e^{-i\omega})\, g_x(e^{-i\omega}).
```

Notice that the covariance between $y$ and $x$ can be recovered from the inversion
formula

$$
c_{yx}(k) = \frac{1}{2\pi}\int_{-\pi}^{\pi}
  h(e^{-i\omega})\, g_x(e^{-i\omega})\, e^{i\omega k}\, d\omega.
$$

Further, notice that given $g_{yx}(e^{-i\omega})$ and $g_x(e^{-i\omega})$ are used in
the above equation, the resulting estimator of the $h_k$ is known as Hannan's
inefficient estimator.

## An Example

As an example, suppose that the jointly covariance stationary process $(y, x)$ has
covariance generating functions

$$
g_x(z) = \sigma_\epsilon^2 \left(\frac{1}{1-0.9z}\right)\!\left(\frac{1}{1-0.9z^{-1}}\right),
\qquad
g_y(z) = \sigma_u^2(1-0.8z)(1-0.8z^{-1}),
$$
$$
g_{yx}(z) = \sigma_{u\epsilon}(1-0.8z)(1+0.5z^{-1}).
$$

Notice that this is equivalent with

$$
g_x(e^{-i\omega}) = \sigma_\epsilon^2 \frac{1}{(1-0.9e^{-i\omega})(1-0.9e^{+i\omega})}
= \frac{\sigma_\epsilon^2}{1.81 - 1.8\cos\omega},
$$
$$
g_y(e^{-i\omega}) = \sigma_u^2(1.64 - 1.6\cos\omega),
$$
$$
g_{yx}(e^{-i\omega}) = \sigma_{u\epsilon}(0.6 - 0.8e^{-i\omega} + 0.5e^{+i\omega}).
$$

Let us now use formula {eq}`eq-45` to calculate the coefficient generating function
$h(z)$ in the projection of $y_t$ on the entire $x$ process. Using $z$ instead of
$e^{-i\omega}$, formula {eq}`eq-45` becomes

$$
h(z) = g_{yx}(z)/g_x(z).
$$

For our example this gives

$$
g_{yx}(z)/g_x(z)
= (\sigma_{u\epsilon}/\sigma_\epsilon^2)(1-0.8z)(1+0.5z^{-1})(1-0.9z^{-1})(1-0.9z).
$$

The reader can easily multiply this polynomial in $z$ and verify that $h_j = 0$ for
$|j| > 2$ and that $h_j \neq 0$ for $j = -2,-1,0,1,2$. Notice that, as in general,
$h(z)$ is "two-sided," having nonzero coefficients on negative powers of $z$.

Now let us calculate the coefficients of $x_t$ on the entire $y$ process:

$$
x_t = \sum_{j=-\infty}^{\infty} f_j\, y_{t-j} + \eta_t,
$$

$E\eta_t\, y_{t-j} = 0$ for all $j$. Applying formula {eq}`eq-45`, exchanging the
roles of $y$ and $x$, gives

$$
f(z) = g_{xy}(z)/g_y(z) = g_{yx}(z^{-1})/g_y(z).
$$

In our example this gives

$$
f(z) = \frac{\sigma_{u\epsilon}(1-0.8z^{-1})(1+0.5z)}
            {\sigma_u^2(1-0.8z^{-1})(1-0.8z)}
= \frac{\sigma_{u\epsilon}}{\sigma_u^2}\frac{(1+0.5z)}{(1-0.8z)}.
$$

It is readily verified that $f_j = 0$ for $j < 0$, so that $f(z)$ is "one-sided on
the past and present."

We shall shortly study the conditions under which the projection of $y$ on $x$ or of
$x$ on $y$ are "one-sided" or "two-sided."

## A Cross-Spectral Formula

Equation {eq}`eq-16` can be generalized as follows. Let

```{math}
:label: eq-46
y_{1t} = B_1(L)x_t, \qquad y_{2t} = B_2(L)x_t,
```

where $x_t$ is a covariance stationary process and $B_1(L)$ and $B_2(L)$ are the lag
generating functions for square summable lag distributions. Then

```{math}
:label: eq-47
g_{y_1 y_2}(e^{-i\omega}) = B_1(e^{-i\omega})\, B_2(e^{+i\omega})\, g_x(e^{-i\omega}).
```

We invite the reader to verify this formula by using {eq}`eq-46` to calculate
$Ey_{1t}y_{2t-k}$, multiplying by $e^{-i\omega k}$, and summing over $k$. The
derivation mimics the derivation of Equation {eq}`eq-16` above.[^fn-cs-1]

[^fn-cs-1]: Alternatively, write $x_t$ in terms of its moving average representation
$x_t = c(L)\eta_t$, where $g_x(z) = \sigma_\eta^2 c(z)c(z^{-1})$. Then apply
{eq}`eq-16` to the system $y_{1t} = B_1(L)c(L)\eta_t$, $y_{2t} = B_2(L)c(L)\eta_t$.

We now use formula {eq}`eq-47` to show that the spectrum reflects a decomposition of
$x_t$ into processes that are orthogonal across frequencies. Thus let
$y_{1t} = B_1(L)x_t$ and $y_{2t} = B_2(L)x_t$, where $B_1(L)$ and $B_2(L)$ are
chosen to satisfy

$$
B_1(e^{-i\omega}) =
\begin{cases}
1, & \omega \in [-b,-a] \cup [a,b], \\
0, & \omega \notin [-b,-a] \cup [a,b];
\end{cases}
\qquad
B_2(e^{-i\omega}) =
\begin{cases}
1, & \omega \in [-d,-c] \cup [c,d], \\
0, & \omega \notin [-d,-c] \cup [c,d].
\end{cases}
$$

To find the individual distributed lag coefficients, Equation {eq}`eq-36` can be used.
From {eq}`eq-47`, we note that if $[-b,-a] \cup [a,b]$ does not intersect the set of
frequencies $[-d,-c] \cup [c,d]$, then $B_1(e^{-i\omega})B_2(e^{i\omega}) = 0$ for all
$\omega$, so that $g_{y_1 y_2}(e^{-i\omega}) = 0$. This in turn implies that $y_1$ and
$y_2$ are processes that are orthogonal (uncorrelated) at all lags, as can be verified
directly from the inversion formula. In this sense, the spectrum $g_x(e^{-i\omega})$
decomposes the variance of $x$ into a set of mutually orthogonal processes across
frequencies.

## Polar Representation: Phase, Gain, and Coherence

The cross spectrum is a complex quantity that is usually characterized by real numbers
in various ways. One characterization is in terms of its real and imaginary parts

$$
g_{yx}(e^{-i\omega}) = co(\omega) + i\, qu(\omega),
$$

where $co(\omega)$ is called the **cospectrum** and $qu(\omega)$ is called the
**quadrature spectrum**. A more usual representation is the polar one

```{math}
:label: eq-48
g_{yx}(e^{-i\omega}) = r(\omega)\, e^{i\theta(\omega)},
```

where

$$
r(\omega) = \sqrt{co(\omega)^2 + qu(\omega)^2}, \qquad
\theta(\omega) = \tan^{-1}\!\left[\frac{qu(\omega)}{co(\omega)}\right].
$$

The phase statistic gives the lead of $y$ over $x$ at frequency $\omega$, while the
"gain" $r(\omega)$ tells how the amplitude in $x$ is multiplied in contributing to the
amplitude of $y$ at frequency $\omega$. Another interesting number is the coherence

$$
coh(\omega) = \frac{|g_{yx}(e^{-i\omega})|^2}{g_x(e^{-i\omega})\, g_y(e^{-i\omega})},
$$

which, being essentially the ratio of a covariance squared to the product of two
variances, is analogous to an $R^2$ statistic. It indicates the proportion of the
variance in one series at frequency $\omega$ that is accounted for by variation in the
other series.

## Interpretation of the Phase

Notice that from {eq}`eq-45` and from the fact that the spectrum $g_x(e^{-i\omega})$
is real, the phase of the cross spectrum equals the phase of
$h(e^{-i\omega}) = \sum h_j e^{-i\omega j}$, which is the Fourier transform of the
$h_j$. That is, writing {eq}`eq-45` and {eq}`eq-48`, we have

$$
r(\omega)\, e^{i\theta(\omega)} = g_{yx}(e^{-i\omega}) = h(e^{-i\omega})\, g_x(e^{-i\omega}),
$$

or

$$
h(e^{-i\omega}) = \frac{r(\omega)}{g_x(e^{-i\omega})}\, e^{i\theta(\omega)},
$$

which shows that the phase of $g_{yx}(e^{-i\omega})$ equals the phase of
$h(e^{-i\omega})$. For convenience, represent $h(e^{-i\omega})$ in polar form

$$
h(e^{-i\omega}) = s(\omega)\, e^{i\theta(\omega)},
$$

where $s(\omega) = r(\omega)/g_x(e^{-i\omega})$.

The following provides a heuristic device for interpreting $\theta(\omega)$. Suppose
we consider as an input into the system {eq}`eq-44` an $x$ series consisting of a pure
cosine wave of frequency $\omega$:

$$
x_t = 2\cos\omega t = e^{i\omega t} + e^{-i\omega t}.
$$

For this input path, suppressing the disturbances $\epsilon_t$, {eq}`eq-44` becomes

$$
y_t = \sum h_j[e^{i\omega(t-j)} + e^{-i\omega(t-j)}]
= e^{i\omega t}\sum h_j e^{-i\omega j} + e^{-i\omega t}\sum h_j e^{+i\omega j}.
$$

But $\sum h_j e^{-i\omega j} = s(\omega)e^{i\theta(\omega)}$ and
$\sum h_j e^{+i\omega j}$, being the complex conjugate of $\sum h_j e^{-i\omega j}$,
equals $s(\omega)e^{-i\theta(\omega)}$. Therefore, we have

$$
y_t = e^{i\omega t} s(\omega) e^{i\theta(\omega)} + e^{-i\omega t} s(\omega) e^{-i\theta(\omega)}
= s(\omega)\left[e^{i\omega t} e^{i\theta(\omega)} + e^{-i\omega t} e^{-i\theta(\omega)}\right]
= s(\omega)\cdot 2\cos(\omega t + \theta(\omega)).
$$

Therefore, the response of {eq}`eq-44` to an input in the form of a cosine wave of
frequency $\omega$ is a cosine wave at the same frequency with amplitude multiplied by
$s(\omega)$ and phase shifted by $\theta(\omega)$. The input cosine wave is at its peak
at $t = 0$, while the output is at its peak at $\omega t + \theta(\omega) = 0$ or
$t = -\theta(\omega)/\omega$ units of time. Thus, for $\theta(\omega) > 0$, the output
*leads* the input by $-\theta(\omega)/\omega$ units of time (where we adopt the usual
convention that $\theta(\omega)$ is constrained to be between $-\pi$ and $+\pi$, a
convention needed to make the arctangent function single-valued).

While useful, the preceding interpretation of the phase has to be used cautiously. The
reason is that the stochastic difference equations that we have been studying generate
random processes with spectral power distributed across a continuum of frequencies
between $-\pi$ and $+\pi$. It is really only over a non-negligible *band* of frequencies
that there occurs a positive contribution to variance. Thus, for such processes, there
really do not occur input processes that are pure cosines, though this situation could be
approached if the spectral density did display a very sharp peak at a given frequency.
Processes with positive spectral power at a single given frequency do exist, and
realizations of these processes do consist of (sums of) sine and cosine waves. But such
processes are not generated by the stochastic equations that we are studying.

## Two Facts about $h(e^{-i\omega})$

It is interesting to note the following two facts about $h(e^{-i\omega})$.

**First**, from the definition of $h(e^{-i\omega})$

$$
h(e^{-i\omega}) = \sum_j h_j\, e^{-i\omega j},
$$

we note that $h(e^{-i\omega})$ evaluated at $\omega = 0$ is the sum of the lag weights,
i.e.,

$$
h(e^{-i\cdot 0}) = \sum_j h_j.
$$

Notice that since

$$
\sum h_j\, e^{-i\omega j} = \sum h_j \cos\omega j - i\sum h_j \sin\omega j
$$

and that since $\sin 0 = 0$, we have that

$$
h(e^{-i\cdot 0}) = s(0) = \sum h_j.
$$

Since $h(e^{-i\omega})$ is real at zero frequency, the phase statistic $\theta(\omega)$
is zero at zero frequency, provided $\sum h_j \neq 0$:

```{math}
:label: eq-49
\theta(\omega) = \tan^{-1}\!\left[\frac{-\sum h_j \sin\omega j}{\sum h_j \cos\omega j}\right],
\qquad \theta(0) = \tan^{-1}[0] = 0.
```

**Second**, it is possible to show that the derivative of the phase statistic with
respect to $\omega$ evaluated at $\omega = 0$ equals minus the mean lag. Recall that

$$
\frac{d}{dx}\tan^{-1} u = \frac{1}{1+u^2}\frac{du}{dx}.
$$

Applying this to {eq}`eq-49` gives

$$
\theta'(\omega) = \frac{1}{1 + \left[-\dfrac{\sum h_j \sin\omega j}{\sum h_j \cos\omega j}\right]^2}
\times \frac{-\sum h_j \cos\omega j \sum h_j j \cos\omega j
             - \sum h_j j \sin\omega j \sum h_j \sin\omega j}
            {(\sum h_j \cos\omega j)^2}.
$$

Evaluating $\theta'(\omega)$ at $\omega = 0$ gives

$$
\theta'(0) = -\frac{\sum h_j\, j}{\sum h_j}.
$$

(Here we have used the facts that $\cos 0 = 1$, $\sin 0 = 0$.) The right-hand side of
this equation is minus the "mean lag" of the lag distribution formed by the $h$'s, a
statistic often reported in econometric studies involving the estimates of distributed
lags.
