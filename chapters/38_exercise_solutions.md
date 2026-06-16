# Solutions to Exercises

Worked solutions to the exercises of the previous chapter. Each heading links
back to the corresponding problem. Equation numbers in `{eq}` refer to the
chapter text; e.g. {eq}`eq-45` is the cross-spectral projection formula
$g_{yx}(z) = h(z)\,g_x(z)$.

---

## Exercise 1

[→ {ref}`Exercise 1 <ex-1>`]  *(Sims's approximation error formula.)*

Write the true (population) projection of $y_t$ on the whole $x$ process as

$$
y_t = \sum_{j=-\infty}^{\infty} b_j^0 x_{t-j} + \epsilon_t,
\qquad E\epsilon_t x_{t-j} = 0 \ \text{ for all } j,
$$

so that $\epsilon_t$ is orthogonal to every $x_{t-j}$. This orthogonality is the
content of the projection.

**(a)** For any constrained coefficient sequence $\{b_j^1\}$,

$$
y_t - \sum_j b_j^1 x_{t-j}
= \Big(\sum_j b_j^0 x_{t-j} + \epsilon_t\Big) - \sum_j b_j^1 x_{t-j}
= \sum_j (b_j^0 - b_j^1)\,x_{t-j} + \epsilon_t \equiv z_t .
$$

Because $\epsilon_t \perp x_{t-j}$ for all $j$, the cross term vanishes when we
square and take expectations:

$$
E\Big(y_t - \sum_j b_j^1 x_{t-j}\Big)^2
= E(z_t^2)
= E\Big(\sum_j (b_j^0 - b_j^1) x_{t-j}\Big)^2 + E\epsilon_t^2 .
$$

**(b)** The first term is the variance of the filtered process
$w_t \equiv \big(b^0(L) - b^1(L)\big)x_t$. By the filtering formula {eq}`eq-33`
(with no extra noise), the spectrum of $w_t$ is

$$
g_w(e^{-i\omega}) = \big|b^0(e^{-i\omega}) - b^1(e^{-i\omega})\big|^2\, g_x(e^{-i\omega}),
$$

and since $\epsilon_t \perp x$ at all lags, $z_t = w_t + \epsilon_t$ has spectrum

$$
g_z(e^{-i\omega}) = \big|b^0(e^{-i\omega}) - b^1(e^{-i\omega})\big|^2\, g_x(e^{-i\omega})
+ g_\epsilon(e^{-i\omega}).
$$

**(c)** By the inversion formula {eq}`eq-20` evaluated at lag $0$, the variance is
the integral of the spectrum:

$$
E z_t^2 = \frac{1}{2\pi}\int_{-\pi}^{\pi} g_z(e^{-i\omega})\,d\omega
= \frac{1}{2\pi}\int_{-\pi}^{\pi}
   \big|b^0(e^{-i\omega}) - b^1(e^{-i\omega})\big|^2\, g_x(e^{-i\omega})\,d\omega
 + E\epsilon_t^2 .
$$

The term $E\epsilon_t^2$ does not involve $b^1$. Hence least squares — which
minimizes $E z_t^2$ over the admissible class $\{b_j^1\}$ — chooses $b^1$ to
minimize

$$
\int_{-\pi}^{\pi}
   \big|b^0(e^{-i\omega}) - b^1(e^{-i\omega})\big|^2\, g_x(e^{-i\omega})\,d\omega,
$$

which is **Sims's formula**. The fitted $b^1$ approximates the true $b^0$ in a
mean-square sense *weighted by the spectral density $g_x$*: the fit is most
accurate at the frequencies where $x$ has most of its power. $\blacksquare$

---

## Exercise 2

[→ {ref}`Exercise 2 <ex-2>`]  *("Optimal" seasonal adjustment via signal extraction.)*

Here $X_t = x_t + u_t$ with $x \perp u$ at all leads and lags, so the cross- and
auto-covariance generating functions satisfy $g_{xX}(z) = g_x(z)$ and
$g_X(z) = g_x(z) + g_u(z)$.

**A.** The coefficients of the two-sided projection of $x_t$ on the $X$ process
follow from {eq}`eq-45` (with the roles "$y\to x$, $x\to X$"):

$$
h(z) = \frac{g_{xX}(z)}{g_X(z)} = \frac{g_x(z)}{g_x(z) + g_u(z)},
\qquad
h(e^{-i\omega}) = \frac{g_x(e^{-i\omega})}{g_x(e^{-i\omega}) + g_u(e^{-i\omega})} .
$$

This is the **Wiener filter**; the $h_j$ are its Fourier coefficients,
$h_j = \tfrac{1}{2\pi}\int_{-\pi}^{\pi} h(e^{-i\omega}) e^{i\omega j}\,d\omega$.

**B.** Since $\hat{x}_t = h(L)X_t$, the filtering formula {eq}`eq-33` gives
$g_{\hat{x}} = |h(e^{-i\omega})|^2 g_X$. Both $g_x$ and $g_u$ are real and
nonnegative, so $h$ is real and $|h|^2 = g_x^2/(g_x+g_u)^2$. Using $g_X = g_x+g_u$,

$$
g_{\hat{x}}(e^{-i\omega})
= \frac{g_x^2}{(g_x+g_u)^2}\,(g_x+g_u)
= \frac{g_x(e^{-i\omega})^2}{g_x(e^{-i\omega}) + g_u(e^{-i\omega})}
= g_x(e^{-i\omega})\cdot\frac{g_x}{g_x+g_u}.
$$

Therefore

$$
\frac{g_{\hat{x}}(e^{-i\omega})}{g_x(e^{-i\omega})}
= \frac{g_x(e^{-i\omega})}{g_x(e^{-i\omega}) + g_u(e^{-i\omega})} < 1
\quad\text{for all }\omega,
$$

because $g_u(e^{-i\omega}) > 0$ everywhere. Hence
$g_{\hat{x}}(e^{-i\omega}) < g_x(e^{-i\omega})$ for all $\omega$: the optimal
estimate is always *less* variable, frequency by frequency, than the target.
$\blacksquare$

**C.** Write $g_{\hat{x}} = g_x \cdot h$ with the gain $h = g_x/(g_x+g_u)$.
At the seasonal frequencies $g_u$ is large (its power is concentrated there), so
$h = g_x/(g_x+g_u)$ is small — the gain has deep dips at the seasonal
frequencies. Away from the seasonal frequencies $g_u \approx 0$, so $h \approx 1$
and $g_{\hat{x}} \approx g_x$. If $g_x$ is itself relatively smooth across both
bands, then $g_{\hat{x}} = g_x h$ inherits the shape of $h$: it nearly equals
$g_x$ at the nonseasonal frequencies and dips sharply at the seasonal ones.
This is precisely the pattern of spectral dips produced by seasonal adjustment
seen in Figure 8. $\blacksquare$

---

## Exercise 3

[→ {ref}`Exercise 3 <ex-3>`]  *(The backward Wold representation.)*

```{note}
As printed, the representation reads $x_t = \sum_{j\ge0} c_j u_{t-j} + \theta_t$;
since $u_t = x_t - P[x_t\mid x_{t+1},x_{t+2},\dots]$ is the *backward* innovation,
the sum must run over **future** innovations, $x_t = \sum_{j\ge0} c_j u_{t+j} + \theta_t$.
(For $j\ge1$, $u_{t-j}$ is orthogonal to $x_t$, so past backward-innovations
cannot build $x_t$.) We prove the corrected statement.
```

**A.** Apply Wold's theorem (chapter §13) to the **time-reversed** process
$\check{x}_t := x_{-t}$. It is covariance stationary with the same covariogram
($c(\tau)=c(-\tau)$), so Wold gives

$$
\check{x}_t = \sum_{j=0}^{\infty} c_j \check{\epsilon}_{t-j} + \check{\eta}_t,
\quad c_0 = 1,\ \textstyle\sum_j c_j^2 < \infty,\quad
\check{\epsilon}_t = \check{x}_t - P[\check{x}_t\mid \check{x}_{t-1},\dots],
$$

with $\check{\epsilon}$ serially uncorrelated, $\check{\eta}$ linearly
deterministic, and $E\check{\eta}_t\check{\epsilon}_s = 0$. Reversing time
($t \mapsto -t$) and setting

$$
u_t := \check{\epsilon}_{-t} = x_t - P[x_t\mid x_{t+1}, x_{t+2},\dots],
\qquad
\theta_t := \check{\eta}_{-t},
$$

(the conditioning set $\{\check{x}_{t-1},\dots\}$ becomes $\{x_{t+1},\dots\}$),
gives

$$
x_t = \sum_{j=0}^{\infty} c_j\, u_{t+j} + \theta_t ,
$$

with $u$ white noise, $E\theta_t u_s = 0$ for all $t,s$, and $\theta_t$
predictable arbitrarily well from *future* $x$'s. $\blacksquare$

**B.** The forward Wold theorem writes $x_t = d(L)\epsilon_t$ with spectrum
$g_x(e^{-i\omega}) = \sigma_\epsilon^2\, d(e^{-i\omega})d(e^{i\omega})$, where
$d(z)=\sum_{j\ge0}d_j z^j$ is *fundamental* ($d_0=1$, zeros of $d(z)$ on or
outside the unit circle). The backward representation gives the **same**
spectrum,

$$
g_x(e^{-i\omega}) = \sigma_u^2\, c(e^{i\omega})c(e^{-i\omega}),
$$

with $c(z)=\sum_{j\ge0}c_j z^j$, $c_0=1$, fundamental. The factorization of a
spectral density into $\sigma^2\,\phi(z)\phi(z^{-1})$ with $\phi$ one-sided,
$\phi(0)=1$, and zeros outside the unit circle is **unique**. Comparing the two
factorizations forces

$$
c_j = d_j \ \text{ for all } j, \qquad \sigma_u^2 = \sigma_\epsilon^2 . \qquad\blacksquare
$$

**C.** **No**, $u_t \neq \epsilon_t$ in general: $\epsilon_t$ is the innovation
relative to the *past* ($x_t - P[x_t\mid x_{t-1},\dots]$) while $u_t$ is the
innovation relative to the *future*; they are built from different information
sets. But **yes**, $Eu_t^2 = E\epsilon_t^2$: by part B, $\sigma_u^2 = \sigma_\epsilon^2$.
(Equivalently, the one-step prediction-error variance is the same forward and
backward — it equals $\exp\!\big[(2\pi)^{-1}\int_{-\pi}^{\pi}\ln g_x(e^{-i\omega})\,d\omega\big]$,
the Kolmogorov formula, which is symmetric in time.)

**D.** **Yes**, $\theta_t = \eta_t$. The linearly deterministic component is a
process (e.g. a sum of fixed-frequency sinusoids with random amplitudes) that
can be predicted arbitrarily well from *either* its past or its future; it is
the unique "remote" component common to both decompositions. Hence the
deterministic part obtained looking backward coincides with the one obtained
looking forward. (For a purely indeterministic $x$, both are zero.) $\blacksquare$

---

## Exercise 4

[→ {ref}`Exercise 4 <ex-4>`]  *(The explosive first-order Markov process.)* Write $\lambda \equiv A > 1$.

**A.** From $y_t = \lambda y_{t-1} + \epsilon_t$ and the initial condition $y_0$,
iterating gives the causal solution

$$
y_t = \lambda^t y_0 + \sum_{k=1}^{t}\lambda^{t-k}\epsilon_k
    = \lambda^t\Big(y_0 + \sum_{k=1}^{t}\lambda^{-k}\epsilon_k\Big).
$$

Define the (path-dependent) random variable

$$
\eta_0 := y_0 + \sum_{k=1}^{\infty}\lambda^{-k}\epsilon_k ,
$$

which converges in mean square because $\lambda > 1$. Then

$$
y_t = \lambda^t \eta_0 - \lambda^t\!\!\sum_{k=t+1}^{\infty}\lambda^{-k}\epsilon_k
    = \lambda^t \eta_0 \;-\; \sum_{i=1}^{\infty}\lambda^{-i}\epsilon_{t+i}
    \;\equiv\; \lambda^t \eta_0 + s_t ,
$$

so the deviation from the explosive trend, $s_t := y_t - \lambda^t\eta_0 =
-\sum_{i\ge1}\lambda^{-i}\epsilon_{t+i}$, is a covariance-stationary process. Its
covariogram is

$$
c_s(\tau) = E s_t s_{t-\tau}
= \sigma_\epsilon^2 \sum_{i\ge1}\lambda^{-i}\lambda^{-(i+|\tau|)}
= \frac{\sigma_\epsilon^2}{\lambda^2-1}\,\lambda^{-|\tau|},
$$

which is exactly the covariogram of a **stable** first-order Markov process with
root $\lambda^{-1}$. Hence $s_t$ has the fundamental (one-sided, past) Wold
representation

$$
s_t = \frac{1}{1-\lambda^{-1}L}\,u_t, \qquad u_t = s_t - \lambda^{-1}s_{t-1},
$$

and therefore

$$
y_t = \lambda^t\eta_0 + \frac{1}{1-\lambda^{-1}L}\,u_t .
$$

Matching the covariogram $c_s(\tau)=\dfrac{\sigma_u^2}{1-\lambda^{-2}}\lambda^{-|\tau|}$
fixes the innovation variance,

$$
\sigma_u^2 = (1-\lambda^{-2})\,\frac{\sigma_\epsilon^2}{\lambda^2-1}
          = \frac{\sigma_\epsilon^2}{\lambda^2}.
$$

(A simulation confirms $c_s(\tau)=\sigma_\epsilon^2\lambda^{-|\tau|}/(\lambda^2-1)$,
that $u_t$ is white, and that $\sigma_u^2=\sigma_\epsilon^2/\lambda^2$.)

**B.** The original $\epsilon_t$ is the innovation in the *causal* (explosive)
representation, with the unstable root $\lambda$; it is **not** the fundamental
white noise. The white noise $u_t$ *is* fundamental for the stationary part
$\{s_t\}$ (equivalently, for $\{y_t\}$ once the deterministic explosive term
$\lambda^t\eta_0$ is removed): its representation
$s_t = (1-\lambda^{-1}L)^{-1}u_t$ uses the stable root $\lambda^{-1}$, so $u_t$
lies in the space spanned by current and past $s$'s and equals
$s_t - P[s_t\mid s_{t-1},s_{t-2},\dots]$. Note $u_t \neq \epsilon_t$ and
$\sigma_u^2 = \sigma_\epsilon^2/\lambda^2 \neq \sigma_\epsilon^2$: "solving the
unstable root forward" trades the original innovation for the fundamental one and
shrinks the innovation variance by $\lambda^{-2}$. $\blacksquare$

---

## Exercise 5

[→ {ref}`Exercise 5 <ex-5>`]  *(Companion form and prediction of an ARMA(1,1).)*

**A.** With state $x_t = (z_t,\ a_t)'$ and $\epsilon_t = (a_t,\ a_t)'$, the
process $z_{t+1} = \lambda z_t - \beta a_t + a_{t+1}$ together with the identity
$a_{t+1}=a_{t+1}$ is the first-order system {eq}`eq-104`,
$x_{t+1} = A x_t + \epsilon_{t+1}$, with

$$
A = \begin{bmatrix} \lambda & -\beta \\ 0 & 0 \end{bmatrix},
\qquad
\begin{bmatrix} z_{t+1} \\ a_{t+1} \end{bmatrix}
= \begin{bmatrix} \lambda & -\beta \\ 0 & 0 \end{bmatrix}
  \begin{bmatrix} z_{t} \\ a_{t} \end{bmatrix}
+ \begin{bmatrix} a_{t+1} \\ a_{t+1} \end{bmatrix}.
$$

**B.** Because $a_t$ is fundamental for $z$, it is recoverable from current and
past $z$'s, $a_t = z_t - P[z_t\mid z_{t-1},\dots]$; conditioning on
$\{z_t,z_{t-1},\dots\}$ is the same as conditioning on $x_t=(z_t,a_t)'$.
The prediction formula {eq}`eq-108`, $P[x_{t+\tau}\mid x_t]=A^\tau x_t$, gives for
$\tau=2$

$$
A^2 = \begin{bmatrix} \lambda & -\beta \\ 0 & 0 \end{bmatrix}^2
    = \begin{bmatrix} \lambda^2 & -\lambda\beta \\ 0 & 0 \end{bmatrix},
\qquad
P[z_{t+2}\mid z_t,z_{t-1},\dots] = e_1' A^2 x_t = \lambda^2 z_t - \lambda\beta\, a_t,
$$

where $e_1'=(1,0)$ and $a_t = \dfrac{1-\lambda L}{1-\beta L}\,z_t$.

*Verification via Wiener–Kolmogorov {eq}`eq-62`.* The Wold representation is
$z_t = c(L)a_t$ with $c(L) = \dfrac{1-\beta L}{1-\lambda L}$, whose coefficients
are $c_0=1$ and $c_j = \lambda^{j-1}(\lambda-\beta)$ for $j\ge1$. Then

$$
\Big[\tfrac{c(L)}{L^2}\Big]_+ = \sum_{k\ge0} c_{k+2}L^k
= \lambda(\lambda-\beta)\sum_{k\ge0}\lambda^k L^k
= \frac{\lambda(\lambda-\beta)}{1-\lambda L},
$$

so {eq}`eq-62` gives

$$
P[z_{t+2}\mid z_t,\dots]
= \Big[\tfrac{c(L)}{L^2}\Big]_+\frac{1}{c(L)}\,z_t
= \frac{\lambda(\lambda-\beta)}{1-\lambda L}\cdot\frac{1-\lambda L}{1-\beta L}\,z_t
= \frac{\lambda(\lambda-\beta)}{1-\beta L}\,z_t .
$$

The companion-form answer agrees: substituting $a_t = \frac{1-\lambda L}{1-\beta L}z_t$,

$$
\lambda^2 z_t - \lambda\beta\,a_t
= \frac{\lambda^2(1-\beta L) - \lambda\beta(1-\lambda L)}{1-\beta L}\,z_t
= \frac{\lambda(\lambda-\beta)}{1-\beta L}\,z_t . \qquad\blacksquare
$$
