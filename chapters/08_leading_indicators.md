# A Digression on Leading Indicators

For years, the National Bureau of Economics Research (NBER) has employed a number of
heuristic techniques designed to isolate "leading indicators" of business cycle
movements, presumably as an aid in the early recognition and prediction of cyclical
movements.[^fn-1] To translate into our vocabulary, essentially a good leading indicator
displays a sizable phase lead at the low business cycle frequencies over some important
"coincident" measures of the cycle, such as unemployment or GNP (as well as a large
coherence with those coincident measures—so that the phase lead is not only large on
average but is regular in its occurrence). While searching for the leading indicators is
perhaps an important thing to do in terms of categorizing data, it is important to
recognize that a series $y_t$ that displays a sizable phase lead over another series
$x_t$, at the most important business cycle frequencies does *not* necessarily help in
predicting $x_t$ any better that can be done by using past $x$'s alone to predict $x$.
We illustrate this fact with two examples.

[^fn-1]: Leading indicators are published in *Business Conditions Digest*, published by
the Department of Commerce.

First suppose we have the system governed by

```{math}
:label: eq-50
x_t = \lambda x_{t-1} + u_t, \qquad |\lambda| < 1, \qquad y_t = h_0 x_t + h_1 x_{t-1} + \epsilon_t
```

where $E\epsilon_t = Eu_t\epsilon_{t-s} = 0$ for all $t$ and $s$, and where both $u$
and $\epsilon$ are serially uncorrelated. The cross spectrum between $y$ and $x$ is
given by

$$
\theta_{yx}(e^{-i\omega}) = (h_0 + h_1 e^{-i\omega})\, g_x(e^{-i\omega})
= (h_0 + h_1\cos\omega - ih_1\sin\omega)\, g_x(e^{-i\omega}) \\
= r(\omega)\, e^{i\theta(\omega)}\, g_x(e^{-i\omega})
$$

where

$$
r(\omega) = \sqrt{(h_0 + h_1\cos\omega)^2 + (h_1\sin\omega)^2}, \qquad
\theta(\omega) = \tan^{-1}\!\left[\frac{-h_1\sin\omega}{h_0 + h_1\cos\omega}\right].
$$

Now by suitably choosing $h_0$ and $h_1$, at a given frequency $\theta(\omega)$ can be
set arbitrarily in the interval $(-\pi, \pi)$. This is in spite of the fact that the
model {eq}`eq-50` implies that $y_t$ is of no use in terms of predicting $x_t$, for
$x_t$ is governed by a pure "autoregression," and depends only on itself lagged and the
unpredictable random term $u_t$. Thus, even if $y_t$ leads $x_t$ at the low business
cycle frequencies, it is of no use in predicting $x_t$.

To specialize this example somewhat, suppose we have

$$
x_t = \lambda x_{t-1} + u_t, \qquad y_t = (x_t - x_{t-1}) + \epsilon_t,
$$

where as before $u$ and $\epsilon$ are mutually orthogonal (at all lags) white-noise
processes. Calculating $h(e^{-i\omega})$, we have

$$
h(e^{-i\omega}) = 1 - e^{-i\omega} = e^{-i\omega/2}(e^{i\omega/2} - e^{-i\omega/2}) \\
= e^{-i\omega/2}\, 2i\sin(\omega/2) = 2e^{-i\omega/2}\, e^{i\pi/2}\sin(\omega/2) \\
= 2e^{i(\pi/2 - \omega/2)}\sin(\omega/2).
$$

For $0 < \omega < \pi$, the phase angle is positive, implying that the output $y$ leads
$x$ at all frequencies between zero and $\pi$. In spite of the fact that $y$ leads at
all of these frequency components, $y$ is of no use in predicting $x$ once lagged
$x$'s are taken into account.

As our second example, consider the system

$$
y_t = \sum_{j=-\infty}^{\infty} h_j x_{t-j} + \epsilon_t, \qquad y_t = \lambda y_{t-1} + u_t
$$

where we assume $E\epsilon_t x_s = 0$ for all $t, s$, $Eu_t = 0$ and $u_t$ is a
white-noise stationary process. We further assume that

$$
h_j = h_{-j} \qquad \text{for all } j \geq 1.
$$

The cross spectrum between $y$ and $x$ is calculated to be

$$
g_{yx}(e^{-i\omega}) = \{h_0 + h_1(e^{i\omega} + e^{-i\omega}) + h_2(e^{2i\omega} + e^{-2i\omega}) + \cdots\}\, g_x(e^{-i\omega}) \\
= \left(h_0 + 2\sum_{j=1}^{\infty} h_j\cos\omega j\right) g_x(e^{-i\omega})
$$

which is *real* for all $\omega$. Therefore, the phase shift $\theta(\omega) = 0$, so
that $y$ and $x$ are perfectly in phase at all frequencies. Despite this, by using a
theorem due to Sims (see pp. 314–324 in original check now) it is possible to show that
even given the past of $x$, past $y$ does help to predict present and future $x$'s.
This is a consequence of the lag distribution of the $h_j$ being two-sided and of
Sim's theorem 2, which we will describe in detail presently.

Taken altogether, these two examples illustrate the fact that displaying a phase lead is
neither a necessary nor a sufficient condition for one series to be of use in predicting
another.
