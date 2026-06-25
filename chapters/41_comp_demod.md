# Complex Demodulation

The spectrum $g_x(e^{-i\omega})$ and the cross spectrum $g_{yx}(e^{-i\omega})$ of the
earlier chapters are properties of a *covariance-stationary* process: the autocovariances
$c_x(k)=Ex_t x_{t-k}$ that generate them do not depend on calendar time $t$. Yet many
economic series carry seasonal or cyclical components whose strength and timing visibly
drift across a sample. The seasonal in a short-term interest rate, for example, may be
large in one decade and small in another, and its peak may slide from one month to
another as institutions change. For such series we would like a *moving* spectrum and a
*moving* cross spectrum — an estimate, local in time, of how much variance sits in a
narrow band of frequencies and of how two series are related there.

**Complex demodulation** is a tool for exactly this. It isolates a narrow band around a
chosen frequency $\omega_0$ and returns a slowly varying *amplitude* and *phase* for that
band as functions of time. Time-averaging its squared modulus recovers the ordinary
spectrum; letting the average run only over a moving window produces an evolutionary
spectrum. The technique was put to work in economics by Godfrey (1965), as a step in
estimating bispectra, and by Sargent (1968), to study the changing seasonal behavior of
U.S. interest rates during the 1950s. This chapter develops the machinery in the notation
of the preceding chapters and reproduces its logic with code.

## The complex demodulate

Fix a center frequency $\omega_0 \in (0,\pi)$. Complex demodulation of a real series $x_t$
at $\omega_0$ consists of two operations: a **frequency shift** followed by **low-pass
filtering**.

Earlier chapters reserve $L$ for the lag operator, while Godfrey and Sargent write $L$ for
the low-pass filter. To avoid the collision we denote the low-pass filter by $W(L)$, a
polynomial (or rational function) in the lag operator whose frequency response
$W(e^{-i\omega})$ is concentrated near $\omega = 0$. The complex demodulate of $x_t$ at
$\omega_0$ is

```{math}
:label: eq-cd-1
d_t(\omega_0) = W(L)\!\left[x_t\, e^{-i\omega_0 t}\right].
```

Multiplying by $e^{-i\omega_0 t}$ slides the whole frequency axis down by $\omega_0$: power
that sat in the band $[\omega_0-\Delta,\ \omega_0+\Delta]$ is moved to $[-\Delta,\Delta]$,
centered on the origin. The low-pass filter $W(L)$ then keeps that band and discards
everything else. What survives is a complex-valued series that varies *slowly* in $t$.

Writing $e^{-i\omega_0 t} = \cos\omega_0 t - i\sin\omega_0 t$, the real and imaginary parts
of {eq}`eq-cd-1` are built from two real demodulates,

```{math}
:label: eq-cd-2
a_t = W(L)\!\left[x_t\cos\omega_0 t\right], \qquad
b_t = W(L)\!\left[x_t\sin\omega_0 t\right], \qquad
d_t(\omega_0) = a_t - i\,b_t .
```

From the demodulate we read off a time-varying **amplitude** and **phase**,

```{math}
:label: eq-cd-3
R_t(\omega_0) = 2\left|d_t(\omega_0)\right| = 2\sqrt{a_t^2 + b_t^2}, \qquad
\theta_t(\omega_0) = \arg d_t(\omega_0).
```

The factor of $2$ in $R_t$ is conventional and is explained by the calculation below.

```{note}
Sargent (1968, fn. 9) and Godfrey (1965, eq. 2.1) define the demodulate with the opposite
sign, $d_t = W(L)[x_t e^{+i\omega_0 t}] = a_t + i\,b_t$. This merely conjugates
{eq}`eq-cd-1`; it flips the sign of the phase $\theta_t$ but leaves the amplitude $R_t$,
the coherence, the gain, and the remodulated component (below) unchanged. We use
$e^{-i\omega_0 t}$ so that phase and lead carry the same sign convention as the cross
spectrum $g_{yx}(e^{-i\omega})$ of {doc}`07_cross_spectrum`.
```

## Why it works: a single moving component

To see what the demodulate measures, feed in a pure cosine whose amplitude and phase drift
slowly,

```{math}
:label: eq-cd-4
x_t = A_t \cos(\omega_0 t + \phi_t),
```

and suppose $A_t,\phi_t$ change little over the span of the filter $W(L)$. Then

$$
x_t\, e^{-i\omega_0 t}
= \tfrac{1}{2} A_t\, e^{i\phi_t}
+ \tfrac{1}{2} A_t\, e^{-i(2\omega_0 t + \phi_t)} .
$$

The first term is slowly varying (it sits near frequency $0$); the second oscillates near
$2\omega_0$. The low-pass filter $W(L)$ passes the first and annihilates the second, so

```{math}
:label: eq-cd-5
d_t(\omega_0) \approx \tfrac{1}{2} A_t\, e^{i\phi_t}, \qquad
R_t(\omega_0) \approx A_t, \qquad \theta_t(\omega_0) \approx \phi_t .
```

The demodulate is a slowly turning phasor whose length tracks the local amplitude of the
$\omega_0$ component and whose angle tracks its local phase. This is the sense in which
complex demodulation estimates a *moving* spectrum: $R_t(\omega_0)$ is a local measure of
how much of the variance of $x_t$ is concentrated at frequency $\omega_0$ at date $t$. If
$x_t$ were exactly covariance stationary, $R_t$ and $\theta_t$ would fluctuate only because
of sampling noise around constants; drift in $R_t$ or $\theta_t$ is the visible signature
of nonstationarity within the band.

Connecting back to the spectrum, if $x_t$ is stationary with spectrum
$g_x(e^{-i\omega})$, then the shifted-and-filtered series has

```{math}
:label: eq-cd-6
E\left|d_t(\omega_0)\right|^2
= \frac{1}{2\pi}\int_{-\pi}^{\pi}
   \left|W(e^{-i\omega})\right|^2 g_x\!\left(e^{-i(\omega+\omega_0)}\right) d\omega ,
```

which, because $|W|^2$ is a narrow bump at the origin, is approximately $g_x(e^{-i\omega_0})$
times the filter's bandwidth. Averaging $|d_t|^2$ over the whole sample therefore estimates
the band power at $\omega_0$ — the ordinary spectrum — while averaging over a moving window
estimates the evolutionary spectrum.

## The low-pass filter

Following Godfrey and Sargent we build $W(L)$ by iterating a centered moving average. One
pass of the equal-weight average of length $2m+1$,

```{math}
:label: eq-cd-7
y_t = \frac{1}{2m+1}\sum_{k=-m}^{m} x_{t-k},
```

has the Dirichlet frequency response
$D_m(\omega) = \dfrac{\sin\!\big((2m+1)\omega/2\big)}{(2m+1)\sin(\omega/2)}$. Iterating it
$r$ times convolves the rectangular data window with itself and gives the response
$\big[D_m(\omega)\big]^{r}$. The case $r=2$ is the Bartlett (triangular) window and $r=4$
the Parzen window (Godfrey 1965, eqs. 2.7–2.8). Sargent (1968) used $r=2$ with $m=18$ on
monthly data. Iteration is cheap — each pass is a running sum — and it gives a smooth,
non-negative response with small side lobes, which keeps neighboring bands from leaking in.

The half-width $m$ sets a **time–frequency trade-off** that pervades this subject. A large
$m$ makes $W(e^{-i\omega})$ narrow, sharpening frequency resolution and cleanly separating
$\omega_0$ from its neighbors; but it also smooths the demodulate heavily in time, blurring
any genuine movement in $A_t,\phi_t$ and enlarging the unreliable region at the ends of the
sample (about $rm$ observations at each end). A small $m$ tracks fast changes in time at the
cost of letting adjacent frequencies contaminate the band. The choice is the same one faced
in choosing the lag window for a spectral estimate in {doc}`06_spectrum`.

## Remodulation

Demodulation is invertible on the band it keeps. **Remodulation** shifts the demodulate
back up to $\omega_0$ and takes the real part,

```{math}
:label: eq-cd-8
r_t(\omega_0) = 2\,\mathrm{Re}\!\left[d_t(\omega_0)\, e^{i\omega_0 t}\right]
             = 2\big(a_t \cos\omega_0 t + b_t \sin\omega_0 t\big)
             = R_t(\omega_0)\cos\!\big(\omega_0 t + \theta_t(\omega_0)\big).
```

For the moving component {eq}`eq-cd-4` this returns $A_t\cos(\omega_0 t + \phi_t)$ — the
original band component, with its drifting amplitude and phase intact. To extract a seasonal
that lives at a fundamental and its harmonics, one demodulates at each frequency in the
seasonal set $S = \{\omega_0, 2\omega_0, 3\omega_0,\dots\}$, remodulates each, and adds:

```{math}
:label: eq-cd-9
s_t = \sum_{\omega_j \in S} r_t(\omega_j).
```

Sargent (1968) demodulated monthly interest rates at the twelve-month cycle
$\omega_0 = 2\pi/12 = \pi/6$ and its first three harmonics ($\pi/3$, $\pi/2$, $2\pi/3$),
remodulated, and summed to obtain a *moving seasonal* $s_t$; the deseasonalized series is
then $x_t - s_t$. Because $A_t$ and $\phi_t$ are allowed to evolve, this is a seasonal
adjustment with a seasonal pattern that itself changes through time — exactly what a fixed
symmetric filter of {doc}`33_seasonal_adjustment` cannot deliver.

Figure 1 puts {eq}`eq-cd-1`–{eq}`eq-cd-8` to work on a synthetic monthly
series whose seasonal amplitude swells and fades and whose phase drifts, buried under a
low-frequency cycle, a harmonic, and noise. Demodulation at $\omega_0=\pi/6$ recovers the
moving amplitude and phase, and remodulation reconstructs the seasonal band.

```{figure} ../figures/ch41_demod_seasonal.png
:name: fig-cd-seasonal
:align: center
:width: 100%

**Figure 1.** Complex demodulation at $\omega_0 = 2\pi/12$ of a synthetic monthly series
$x_t = A_t\cos(\omega_0 t + \phi_t) + (\text{10-yr cycle}) + (\text{6-mo harmonic}) +
\text{noise}$. (a) the series with its moving seasonal and the envelope $\pm A_t$;
(b) the recovered amplitude $R_t = 2|d_t|$ tracks the true $A_t$ — a moving spectrum at
the seasonal band (shaded ends are corrupted by the filter); (c) the recovered phase
$\theta_t$ tracks the true $\phi_t$; (d) the remodulated band $r_t$ recovers the seasonal,
and $x_t - r_t$ is the deseasonalized series. The low-pass filter is a length-$25$ moving
average iterated twice (Bartlett window). Generated by
[`code/ch41_demod_seasonal.py`](https://github.com/maxmaxmaxmaxmaxmax373/sargent-time-series/blob/main/code/ch41_demod_seasonal.py).
```

## Locating a spectral peak

The phase of a demodulate carries information about *where* the power actually sits, a fact
Godfrey (1965, eqs. 2.9–2.12) exploited to locate peaks. Suppose the series contains a pure
component at an unknown frequency $\omega_1$, $x_t = a\cos\omega_1 t$, and we demodulate at a
trial frequency $\omega_0$. Then $x_t e^{-i\omega_0 t}$ has a slow part proportional to
$e^{i(\omega_1-\omega_0)t}$, so after low-pass filtering

```{math}
:label: eq-cd-10
d_t(\omega_0) \approx \tfrac{1}{2} a\, e^{i(\omega_1 - \omega_0)t}, \qquad
\theta_t(\omega_0) = (\omega_1 - \omega_0)\, t .
```

The demodulate phase is **linear in time**, with slope equal to the offset between the trial
frequency and the true peak. Demodulating just below the peak gives a phase that climbs at a
constant rate; just above, it falls; exactly on the peak it is flat. This converts
peak-finding into reading the slope of a straight line, and it warns that a demodulate
centered slightly off a sharp spectral peak will *appear* to drift in phase even though the
underlying process is perfectly stationary — a trap to keep in mind when interpreting moving
phases. Figure 2 illustrates the three cases.

```{figure} ../figures/ch41_demod_peak.png
:name: fig-cd-peak
:align: center
:width: 100%

**Figure 2.** A stationary sinusoid of period $12$ demodulated at three trial frequencies.
(a) the unwrapped demodulate phase is linear in time with slope $\omega_1-\omega_0$ (dotted
theory lines lie under the data): demodulating below the peak gives a rising line, above the
peak a falling line, on the peak a flat one. (b) the demodulate amplitude $2|d_t|$ is largest
on the peak. Generated by
[`code/ch41_demod_peak.py`](https://github.com/maxmaxmaxmaxmaxmax373/sargent-time-series/blob/main/code/ch41_demod_peak.py).
```

## A moving cross spectrum

Demodulation extends to two series, producing a cross spectrum that is local in time. Demodulate
$y_t$ and $x_t$ at the same $\omega_0$ to get $d^y_t$ and $d^x_t$, and form the local cross
product, smoothed over a moving window $\langle\cdot\rangle$ in time:

```{math}
:label: eq-cd-11
C^{yx}_t(\omega_0) = \big\langle\, d^y_t(\omega_0)\,\overline{d^x_t(\omega_0)} \,\big\rangle,
\qquad
P^{x}_t(\omega_0) = \big\langle |d^x_t(\omega_0)|^2 \big\rangle,
\qquad
P^{y}_t(\omega_0) = \big\langle |d^y_t(\omega_0)|^2 \big\rangle .
```

This is the time-varying analog of the cross spectrum $g_{yx}(e^{-i\omega_0})$, and from it
come the same three real summaries introduced in {doc}`07_cross_spectrum` — now functions of
$t$:

```{math}
:label: eq-cd-12
\mathrm{coh}_t(\omega_0) = \frac{|C^{yx}_t(\omega_0)|^2}{P^{y}_t(\omega_0)\,P^{x}_t(\omega_0)},
\qquad
\mathrm{gain}_t(\omega_0) = \frac{|C^{yx}_t(\omega_0)|}{P^{x}_t(\omega_0)},
\qquad
\theta^{yx}_t(\omega_0) = \arg C^{yx}_t(\omega_0).
```

The phase converts to a **lead of $y$ over $x$**, measured in time units, by
$\tau_t = \theta^{yx}_t(\omega_0)/\omega_0$; at the twelve-month frequency this is a lead in
months, the conversion Sargent reports in his Table 2. The extra time-smoothing
$\langle\cdot\rangle$ in {eq}`eq-cd-11` is essential: without it, the squared coherence of a
single demodulate pair is identically $1$, exactly as the raw squared coherence of one
realization is (the same point made for the spectrum in {doc}`06_spectrum`). Smoothing over a
window short enough to preserve genuine movement but long enough to average out noise is the
cross-spectral counterpart of the bandwidth choice for $W(L)$.

Figure 3 demodulates two seasonal series whose lead and whose strength of
comovement both drift through time.

```{figure} ../figures/ch41_demod_cross.png
:name: fig-cd-cross
:align: center
:width: 100%

**Figure 3.** A moving cross spectrum at the seasonal band. Two monthly series share a
twelve-month seasonal, but the lead of $y$ over $x$ drifts and their comovement is degraded
by incoherent seasonal noise that is strongest in the middle of the sample. (a) an excerpt of
the two series; (b) the moving coherence {eq}`eq-cd-12` dips where the comovement is noisiest;
(c) the moving phase lead $\tau_t = \theta^{yx}_t/\omega_0$ in months tracks the true drifting
lead (and wobbles where coherence is low). Generated by
[`code/ch41_demod_cross.py`](https://github.com/maxmaxmaxmaxmaxmax373/sargent-time-series/blob/main/code/ch41_demod_cross.py).
```

## Application: the changing seasonal in interest rates, 1951–1960

Sargent (1968) studied monthly rates on eleven money- and capital-market instruments —
three-month Treasury bills, governments of one to twenty years' maturity, commercial paper,
and Moody's Aaa and Baa corporates — over 1951–1960, working with first differences to remove
the trend. The ordinary spectral and cross-spectral statistics (coherence, phase, gain of
{doc}`07_cross_spectrum`) revealed strong seasonal peaks at twelve months and its harmonics,
generally declining in amplitude as maturity lengthens, and a tendency for the longer rates to
*lead* the bill rate, the lead increasing as maturity increases — a pattern reminiscent of
Macaulay's (1938) finding for 1890–1913.

The seasonal peaks raised a question that stationary cross-spectra could not settle: were the
seasonal in the long–short rate *spread* and the seasonal in the average maturity of the
publicly held debt moving *together*, as a segmented-markets or hedging theory of the term
structure would predict? To answer it Sargent turned to complex demodulation, "the relatively
new tool" of his title's program. He demodulated the rate spread and the average maturity at
the twelve-month cycle and its first three harmonics, using the iterated moving-average
low-pass filter of {eq}`eq-cd-7` with $m=18$, and then **remodulated** ({eq}`eq-cd-8`–{eq}`eq-cd-9`)
to construct the estimated seasonal of each series as a function of calendar time. The
resulting moving seasonals — his Figure 3 — show the seasonal pattern of the bill rate, of the
ten-year-minus-bill spread, and of the average maturity evolving month by month from 1954 to
1963, with amplitudes that grow and shrink across the decade rather than repeating identically.

Reading the demodulates jointly delivered the economic conclusion. The seasonals in the rate
spread and in the average maturity were close in phase over the first several years,
consistent with the two being directly related at the seasonal frequencies; the moving cross
spectrum ({eq}`eq-cd-11`–{eq}`eq-cd-12`, his Table 2) showed high coherence at the seasonal and
a small phase difference, with the rate spread leading the average maturity by one to five
months at the lower frequencies. Sargent read this as the Treasury adjusting the maturity
structure of the debt in response to the spread — lengthening maturity when long rates were
low — while the tight seasonal comovement was at best transitory. None of this is visible in a
single stationary spectrum; it requires the local-in-time amplitude and phase that
demodulation provides.

### A replication on FRED data

We can re-run Sargent's exercise on data anyone can download. The variables he used are
*interest rates*, which the St. Louis Fed publishes **not seasonally adjusted** — exactly
what a study of the seasonal requires. The closest monthly NSA series reaching back to the
1950s are the three-month Treasury bill rate (`TB3MS`, from 1934), the one-, five-, and
twenty-year Treasury constant-maturity rates (`GS1`, `GS5`, `GS20`, from April 1953), the
ten-year rate (`GS10`), and Moody's Aaa and Baa corporate yields (`AAA`, `BAA`, from 1919).
Following Sargent we first-difference each rate to remove the trend, then demodulate at the
twelve-month frequency $\omega_0 = 2\pi/12$ and its first three harmonics with the iterated
moving-average filter of {eq}`eq-cd-7` ($m=18$, $r=2$).

Figure 4 reproduces his three seasonal findings. Panel (a) — a smoothed
periodogram of the differenced bill rate over 1951–1972 — shows the sharp peaks at periods
$12,6,4,3$ months that establish the seasonal (his Figure 1a). Panel (b) is the moving
seasonal component $s_t$ of {eq}`eq-cd-9`, the running analog of his Figure 3a. Panel (c)
measures, for each maturity, the local *seasonal share* — the variance of the remodulated
seasonal band as a fraction of the variance of the differenced rate — over 1957–1967, and
finds it declining with maturity (bill $0.27$, one-year $0.19$, five-year $0.16$, twenty-year
$0.16$), Sargent's finding (b). The volatility-normalized share is the honest way to compare
across maturities and eras: the *absolute* seasonal swing in basis points is largest in the
high-rate 1970s and early 1980s simply because every movement was larger then, but relative to
total variation the bill-rate seasonal was strongest in the 1950s–60s and has largely faded
from short rates since the mid-1980s.

```{figure} ../figures/ch41_fred_seasonal.png
:name: fig-cd-fred-seasonal
:align: center
:width: 100%

**Figure 4.** Sargent (1968) replicated on FRED. (a) The spectral density of the
first-differenced three-month bill rate (`TB3MS`, NSA) over 1951–1972 has peaks at the
seasonal periods $12,6,4,3,2$ months. (b) The moving seasonal $s_t = \sum_{k=1}^{4}
r_t(k\omega_0)$ of the differenced bill rate, 1951–1975, with the moving amplitude $\pm R_t$
of the twelve-month fundamental (basis points; his Figure 3a). (c) The seasonal share by
maturity over 1957–1967 declines as maturity lengthens (his finding b). Generated by
[`code/ch41_fred_seasonal.py`](https://github.com/maxmaxmaxmaxmaxmax373/sargent-time-series/blob/main/code/ch41_fred_seasonal.py).
```

The moving cross spectrum reproduces his cross-spectral results just as closely.
Figure 5 demodulates the differenced bill and ten-year rates at the seasonal
and forms {eq}`eq-cd-11`–{eq}`eq-cd-12` over 1953–1975. The coherence is very high throughout
(median $0.97$ — Sargent: "almost all very high"), and the phase lead of the bill over the
ten-year is negative, hovering around $-1.5$ months: **the ten-year rate leads the bill rate**
at the seasonal by about a month and a half. Sargent's Table 1 reports a phase of $-0.7978$
radians for this pair at period 12, which is $-0.7978 \times 12/(2\pi) = -1.52$ months; our
demodulation gives a median of $-1.49$. That the longer rate should lead the shorter at the
seasonal is the pattern Sargent traced back to Macaulay (1938) and read as evidence that the
seasonal was, to a degree, forecast and impounded in the longer rate.

```{figure} ../figures/ch41_fred_cross.png
:name: fig-cd-fred-cross
:align: center
:width: 90%

**Figure 5.** Moving cross spectrum of the first-differenced bill (`TB3MS`) and ten-year
(`GS10`) rates at the twelve-month band, 1953–1975, both NSA. (a) the moving coherence
{eq}`eq-cd-12` stays near $0.97$; (b) the moving phase lead $\theta^{yx}_t/\omega_0$ of the bill
over the ten-year is about $-1.5$ months — the ten-year leads — matching Sargent's Table 1 value
of $-1.52$ months (red dashed). Generated by
[`code/ch41_fred_cross.py`](https://github.com/maxmaxmaxmaxmaxmax373/sargent-time-series/blob/main/code/ch41_fred_cross.py).
```

The one piece of Sargent's study we cannot reproduce from FRED is the average maturity of the
publicly held debt, which is not available as a monthly series back to the 1950s; the
bill-versus-long comparison above stands in for the cross-spectral half of his analysis.

## Exercise

```{admonition} Exercise: detecting a seasonal break
:class: tip

Many seasonal patterns shift abruptly when an institution changes — a payment rule, a market
convention, a policy regime. Complex demodulation should reveal such a break as a jump in the
demodulate amplitude or a kink in its phase.

1. Simulate a monthly series with a twelve-month seasonal that breaks at $t^\star$: let the
   amplitude jump (say from $A=1$ to $A=2$) and the phase shift (say by $\phi=\pi/2$) at
   $t=t^\star$, and add a low-frequency cycle and white noise so the seasonal is not obvious in
   the raw plot:
   :::{math}
   x_t = A(t)\cos\!\big(\tfrac{\pi}{6}t + \phi(t)\big)
          + 0.5\cos\!\big(\tfrac{2\pi}{120}t\big) + \varepsilon_t .
   :::
2. Demodulate at $\omega_0=\pi/6$ with the iterated moving-average filter {eq}`eq-cd-7`
   ($r=2$; try $m=6$ and $m=18$). Plot $R_t = 2|d_t|$ and $\theta_t = \arg d_t$ and locate the
   break. How does the apparent *sharpness* of the break in $R_t$ depend on $m$? Relate your
   answer to the time–frequency trade-off and to the $\sim rm$ observations lost at each end.
3. Remodulate ({eq}`eq-cd-8`) to extract the seasonal component on each side of the break and
   confirm that $x_t - r_t$ is free of the twelve-month oscillation.
4. *Harder.* Add a second series $z_t$ that shares the seasonal but whose break is *lagged* by
   three months, and use the moving cross spectrum {eq}`eq-cd-11`–{eq}`eq-cd-12` to estimate the
   coherence and the phase lead of $x$ over $z$ before and after $t^\star$. Verify that the
   estimated lead $\theta^{xz}_t/\omega_0$ recovers the three-month lag, and watch the coherence
   fall transiently around the break where the two seasonals are momentarily out of step.

The scripts behind Figure 1 and Figure 3 provide a `lowpass`,
`demodulate`, and cross-spectral template you can adapt.
```

## References

- Godfrey, M. D. (1965). An exploratory study of the bi-spectrum of economic time series.
  *Journal of the Royal Statistical Society, Series C (Applied Statistics)* 14(1), 48–69.
- Macaulay, F. R. (1938). *The Movements of Interest Rates, Bond Yields and Stock Prices in the
  United States Since 1856.* National Bureau of Economic Research.
- Sargent, T. J. (1968). Interest rates in the nineteen-fifties. *The Review of Economics and
  Statistics* 50(2), 164–172.
- Blackman, R. B., and J. W. Tukey (1959). *The Measurement of Power Spectra.* Dover.
