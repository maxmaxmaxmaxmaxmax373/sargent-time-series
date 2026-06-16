# Seasonal Adjustment and Forecasting Geometric Distributed Leads

A stochastic process $x_t$ is said to have an important seasonal component if its spectral density has peaks in the vicinity of the seasonal frequencies, meaning that the process has a substantial portion of its variance occurring over the seasonal frequencies. It is a common practice to "seasonally adjust" the series in order to remove the seasonal components. Unless the seasonal components are thought mainly to reflect measurement errors that are unrelated to economic activity, seasonal adjustment is typically misleading. We briefly indicate how using seasonally adjusted data can distort analysis. We have often encountered, and shall encounter even more often later, models in which economic agents decide to make an endogenous variable $y_t$ a geometric sum of expected future values of a variable $x_t$:

```{math}
:label: eq-170
y_t=P_t\sum_{j=0}^\infty \lambda^j x_{t+j},
```

where $P_t$ is the linear least squares projection operator conditioned on $[x_t,x_{t-1},\ldots]$. Now let $x_t$ be a covariance stationary stochastic process with Wold representation

```{math}
:label: eq-171
x_t=c(L)\epsilon_t,
```

Let $x_t$ have a substantial seasonal component, say reflecting weather or seasonal movements in tastes (Christmas and Easter). These seasonal movements are assumed not to be due to measurement error, but are truly occurring in the process that $x_t$ is measuring. We have seen above that {eq}`eq-170` and {eq}`eq-171` imply that $y_t$ will obey

```{math}
:label: eq-172
y_t=\left[\frac{c(L)-\lambda c(\lambda)L^{-1}}{1-\lambda L^{-1}}\right]\epsilon_t,
```

or

```{math}
:label: eq-173
y_t=\left[\frac{1-\lambda c(\lambda)c(L)^{-1}L^{-1}}{1-\lambda L^{-1}}\right]x_t,
```

Now suppose that an analyst has seasonally adjusted $x_t$ to obtain a process $x_t^a$. Typically this is accomplished by taking

```{math}
:label: eq-174
x_t^a=h(L)x_t
```

where $h(L)$ is a long two-sided filter that attenuates (but does not completely eradicate) power at the seasonal frequencies. We indicated how to construct such a filter in Section 6 above. Typically, these seasonal adjustment filters are symmetric with $h_j=h_{-j}$. (See Exercise 2 and Section 6.)

If the economic agent is choosing $y_t$ to obey {eq}`eq-173`, then $y_t$ is related to the seasonally adjusted series $x_t^a$ by

$$
h(L)y_t=\left[\frac{1-\lambda c(\lambda)c(L)^{-1}L^{-1}}{1-\lambda L^{-1}}\right]h(L)x_t,
$$

or

```{math}
:label: eq-175
h(L)y_t=\left[\frac{1-\lambda c(\lambda)c(L)^{-1}L^{-1}}{1-\lambda L^{-1}}\right]x_t^a.
```

Suppose that seasonally adjusted data on $y_t$ denoted $y_t^a$ are created via the *same* filter $h(L)$ used to create $x_t^a$ so that

```{math}
:label: eq-176
y_t^a=h(L)y_t.
```

In view of {eq}`eq-176`, {eq}`eq-175` becomes

```{math}
:label: eq-177
y_t^a=\left[\frac{1-\lambda c(\lambda)c(L)^{-1}L^{-1}}{1-\lambda L^{-1}}\right]x_t^a.
```

where $c(L)$ is the Wold representation for the *seasonally unadjusted* series $x_t$, given in {eq}`eq-171`. (In Sims (1974) and Exercise 10, the effects are studied of adjusting $y_t$ and $x_t$ by using *different* seasonal filters.)

If agents are behaving according to {eq}`eq-173` or {eq}`eq-177`, the restrictions across the $y_t^a$ and $x_t^a$ processes cannot be computed from knowledge of the autoregression for $x_t^a$ alone, without knowledge of the seasonal adjustment filter $h(L)$. If one routinely applies the Hansen-Sargent formula {eq}`eq-90` based on using the moving average or autoregressive representation for the seasonally adjusted series, an error is committed, whose nature we now briefly explore. We continue to suppose that economic agents are actually choosing $y_t$ to satisfy {eq}`eq-170` and therefore {eq}`eq-173`. Suppose that the economist mistakenly uses the seasonally adjusted data to create the model

```{math}
:label: eq-178
y_t=P_t \sum_{j=0}^\infty \lambda^j x^a_{t+j}
```

where now $P_t$ denotes the least squares projection on the space $[x_t^a, x^a_{t-1},\ldots]$. To deduce the implications of model {eq}`eq-178`, we need the Wold representation of $x_t^a$,

```{math}
:label: eq-179
x_t^a=d(L)a_t,\quad d(L)=\sum_{j=0}^\infty d_jL^j
```

$$
a_t = x^a_t - P[x_t^a| x^a_{t-1},x^a_{t-2},\ldots]
$$

Now $d(L)$ of {eq}`eq-179` is related to $c(L)$ of {eq}`eq-171` by the spectral factorization equation

```{math}
:label: eq-180
\sigma^2_ad(L)d(L^{-1}) =\sigma_\epsilon^2h(L)h(L^{-1})c(L)c(L^{-1})
```

where $\sigma^2_a=Ea_t^2$, $\sigma_\epsilon^2=E\epsilon_t^2$, and $d(L)$ has all of its zeros outside the unit circle. Note that in general $h(L)$ is two-sided in powers of $L$, and that $d(L) \neq h(L)c(L)$. Equations {eq}`eq-178` and {eq}`eq-179` imply that

```{math}
:label: eq-181
y_t=\left[\frac{1-\lambda d(\lambda)d(L)^{-1}L^{-1}}{1-\lambda L^{-1}}\right]x_t^a.
```

Using $x_t^a=h(L)x_t$, we have that

```{math}
:label: eq-182
y_t=\left[\frac{1-\lambda d(\lambda)d(L)^{-1}L^{-1}}{1-\lambda L^{-1}}\right]h(L)x_t.
```

A comparison of equations {eq}`eq-173`, {eq}`eq-180`, and {eq}`eq-182` shows how {eq}`eq-182` misspecifies the relationship between $y_t$ and $\{x_t\}$ if agents are responding to forecasts of the seasonally unadjusted series while the economist attributes to agents their responding to forecasts of the seasonally adjusted series. Alternatively, a comparison of {eq}`eq-175` with {eq}`eq-181` shows how deducing the cross-equation restrictions on the basis of the seasonally adjusted data misspecifies the relationship between $y_t^a$ and $x_t^a$. Continuing to suppose that $y_t^a=h(L)y_t$, applying $h(L)$ to both sides of {eq}`eq-181` implies that

```{math}
:label: eq-183
y_t^a=h(L)\left[\frac{1-\lambda d(\lambda)d(L)^{-1}L^{-1}}{1-\lambda L^{-1}}\right]x_t^a.
```

whereas the actual value of $y_t^a$ is determined by the right side of {eq}`eq-175`.

We illustrate these calculations as follows. For monthly U.S. money supply (M2) for the period 1959 to the present (Sargent's original used 1959:1–1986:2) we constructed an 18th order autoregression for both seasonally adjusted and seasonally unadjusted M2. The spectral densities of seasonally adjusted and unadjusted M2 are reported in Figure 8. Notice the dips in the spectral density of the adjusted series at the seasonal frequencies (see Exercise 2). The autoregressive coefficients are plotted for seasonally adjusted and unadjusted series in Figure 9. In Figure 10, we report the values of the two filters in $x$ associated with {eq}`eq-177` and {eq}`eq-181`; namely,

$$
b^{nsa}(L)=h(L)\left[\frac{1-\lambda c(\lambda)c(L)^{-1}L^{-1}}{1-\lambda L^{-1}}\right]x_t^a, \qquad b^{sa}(L)=h(L)\left[\frac{1-\lambda d(\lambda)d(L)^{-1}L^{-1}}{1-\lambda L^{-1}}\right]x_t^a
$$

where we set $\lambda=0.9$ and where we have used the estimated autoregressive representations reported in Figure 9 for $c(L)$ and $d(L)$. Note the difference in $b^{nsa}(L)$ and $b^{sa}(L)$ in the vicinity of 12 lags. Figures 11 and 12 summarize the magnitude and phase of the transfer function from innovations in M2 to "$y_t$" that are associated with $b^{sa}(L)$. That is, let the filters associated with the autoregressive representations be denoted as $a^{sa}(L)$ and $a^{nsa}(L)$, respectively. Then Figures 11–12 report the magnitude and phase of

$$
h(e^{-i\omega})h(e^{+i\omega})
$$

where $h(L)$ is given by $b^{nsa}(L)/a^{nsa}(L)$ or $b^{sa}(L)/a^{sa}(L)$, respectively. Note the substantially different responses revealed in these figures.

Some of these ideas are explored further in Exercises 22 and 23.

```{figure} ../figures/fig8_m2_spectra.png
:name: fig-m2-spectra
:align: center
:width: 100%

**Figure 8.** AR(18) spectral densities of seasonally adjusted (`M2SL`) and
unadjusted (`M2NS`) monthly M2, **1959–present** (Sargent's original used
1959:1–1986:2). An 18th-order autoregression is fit by OLS to each log series;
the curves are $\log g(\omega) = \log\!\big(\sigma^2/|A(e^{-i\omega})|^2\big)$.
The unadjusted series shows sharp peaks at the seasonal periods (12, 6, 4, 3,
2.4, 2 months); seasonal adjustment carves matching dips there — the shaded gap
is exactly the seasonal power removed. (Sargent's caption reads "M1"; the text
describes the construction with M2, which we follow.) Generated by
[`code/fig8_m2_spectra.py`](https://github.com/maxmaxmaxmaxmaxmax373/sargent-time-series/blob/main/code/fig8_m2_spectra.py).
```

```{figure} ../figures/fig9_m2_ar_coeffs.png
:name: fig-ar-coefficients
:align: center
:width: 100%

**Figure 9.** OLS AR(18) coefficients $a_1,\dots,a_{18}$ for unadjusted (`M2NS`)
and adjusted (`M2SL`) log M2, 1959–present. Both share the dominant near-unit-root
coefficients at lags 1–2, but the unadjusted series carries pronounced **seasonal**
coefficients at lags 12–14 that the adjusted series lacks (e.g. $a_{12}\approx 0.76$
unadjusted versus $0.17$ adjusted). Generated by
[`code/fig9_m2_ar_coeffs.py`](https://github.com/maxmaxmaxmaxmaxmax373/sargent-time-series/blob/main/code/fig9_m2_ar_coeffs.py).
```

```{figure} ../figures/fig10_m2_lead_filter.png
:name: fig-filter-comparison
:align: center
:width: 100%

**Figure 10.** Geometric-lead forecast filters $b^{nsa}(L)$ and $b^{sa}(L)$
(Hansen–Sargent formula, $\lambda = 0.9$) implied by the two AR(18)
representations, for the projection $y_t = P_t\sum_j \lambda^j x_{t+j}$. The
filters agree at short lags but diverge sharply near lag 12: forecasting future
*unadjusted* money inherits the seasonal autoregressive structure (the large
spike at lags 12–13), which seasonal adjustment removes. Generated by
[`code/fig10_m2_lead_filter.py`](https://github.com/maxmaxmaxmaxmaxmax373/sargent-time-series/blob/main/code/fig10_m2_lead_filter.py).
```

```{figure} ../figures/fig11_m2_transfer_nsa.png
:name: fig-transfer-magnitude
:align: center
:width: 90%

**Figure 11.** Amplitude (log scale) and phase of the transfer function
$h(e^{-i\omega}) = b^{nsa}(e^{-i\omega})/a^{nsa}(e^{-i\omega})$ — from a money
innovation to the geometric-lead forecast $y_t$ — for the **unadjusted** series
(the adjusted case is overlaid faintly for comparison). The amplitude is
dominated by the near-unit-root pole at low frequency; the unadjusted response
carries extra bumps at the seasonal frequencies (12, 6, 4, 3, 2.4, 2 months).
Generated by
[`code/fig11_12_m2_transfer.py`](https://github.com/maxmaxmaxmaxmaxmax373/sargent-time-series/blob/main/code/fig11_12_m2_transfer.py).
```

```{figure} ../figures/fig12_m2_transfer_sa.png
:name: fig-transfer-phase
:align: center
:width: 90%

**Figure 12.** The same transfer function for the **adjusted** series,
$h(e^{-i\omega}) = b^{sa}(e^{-i\omega})/a^{sa}(e^{-i\omega})$ (the unadjusted case
overlaid faintly). Comparing Figures 11 and 12, the responses coincide away from
the seasonal frequencies and differ there: seasonal adjustment smooths the
filter's response to innovations at the seasonal frequencies. (In Sargent's
1959–1986 sample the two responses looked more dramatically different; over the
longer, more persistent modern M2 the difference is concentrated at the seasonal
bands.) Generated by
[`code/fig11_12_m2_transfer.py`](https://github.com/maxmaxmaxmaxmaxmax373/sargent-time-series/blob/main/code/fig11_12_m2_transfer.py).
```
