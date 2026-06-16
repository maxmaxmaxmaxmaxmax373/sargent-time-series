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

---

```{note}
**Granger-causality criterion used in Exercises 6–10.** By Sims's theorem
(chapter §27), "$P$ Granger-causes $Q$" if and only if the coefficient
generating function of the projection of $P$ on the *entire* $Q$ process,
$g_{PQ}(z)/g_Q(z)$, is **two-sided** (has nonzero coefficients on negative powers
of $z$). If it is one-sided in nonnegative powers, $P$ fails to Granger-cause $Q$.
Here $g_{PQ}(z)=\sum_k E[P_t Q_{t-k}]z^k$ and $g_{QP}(z)=g_{PQ}(z^{-1})$.
```

## Exercise 6

[→ {ref}`Exercise 6 <ex-6>`]

For each pair, form $b(z)=g_{yx}(z)/g_x(z)$ (projection of $y$ on the $x$ process;
its two-sidedness tests "$y$ causes $x$") and $f(z)=g_{yx}(z^{-1})/g_y(z)$
(projection of $x$ on the $y$ process; tests "$x$ causes $y$").

**A.**

$$
b(z)=\frac{g_{yx}(z)}{g_x(z)}
=\frac{\sigma_{\epsilon u}}{\sigma_\epsilon^2}(1-0.8z)(1+0.5z^{-1})(1-0.9z)(1-0.9z^{-1}),
$$

which contains $z^{-1}$ and $z^{-2}$ — **two-sided**, so **$y$ Granger-causes $x$**.

$$
f(z)=\frac{g_{yx}(z^{-1})}{g_y(z)}
=\frac{\sigma_{\epsilon u}(1-0.8z^{-1})(1+0.5z)}{\sigma_u^2(1-0.8z)(1-0.8z^{-1})}
=\frac{\sigma_{\epsilon u}}{\sigma_u^2}\frac{1+0.5z}{1-0.8z},
$$

one-sided in nonnegative powers — **$x$ does not Granger-cause $y$**. (This is the
example of chapter §7.)

**B.**

$$
b(z)=\frac{\sigma_{\epsilon u}(1+0.2z)(1+0.99z^{-1})}{\sigma_\epsilon^2(1+0.99z)(1+0.99z^{-1})}
=\frac{\sigma_{\epsilon u}}{\sigma_\epsilon^2}\frac{1+0.2z}{1+0.99z},
$$

one-sided — **$y$ does not Granger-cause $x$**. Whereas

$$
f(z)=\frac{\sigma_{\epsilon u}(1+0.2z^{-1})(1+0.99z)}
          {\sigma_u^2}(1-0.7z+0.3z^2)(1-0.7z^{-1}+0.3z^{-2})
$$

contains negative powers (from $1+0.2z^{-1}$ and $1-0.7z^{-1}+0.3z^{-2}$) —
**two-sided, so $x$ Granger-causes $y$**.

**C.**

$$
b(z)=\frac{\sigma_{\epsilon u}}{\sigma_\epsilon^2}\frac{1-0.7z}{1-0.8z},
\qquad
f(z)=\frac{\sigma_{\epsilon u}}{\sigma_u^2}\frac{1-0.8z}{1-0.7z},
$$

both one-sided in nonnegative powers (the $z^{-1}$ factors cancel against $g_x$,
$g_y$). Hence **neither variable Granger-causes the other**. (All three parts
confirmed numerically by inverse-FFT of $b(z)$, $f(z)$ and inspection of the
negative-lag coefficients.) $\blacksquare$

## Exercise 7

[→ {ref}`Exercise 7 <ex-7>`]

Write the consumption multiplier $A(L)\equiv(1-b(L))^{-1}$, one-sided in
nonnegative powers by hypothesis. Substituting $c_t=b(L)Y_t+\epsilon_t$ into
$Y_t=c_t+I_t$ gives the reduced forms

$$
Y_t = A(L)\,(\epsilon_t + I_t),
\qquad
c_t = Y_t - I_t = A(L)\epsilon_t + \big(A(L)-1\big)I_t ,
$$

with $\epsilon\perp I$ at all lags.

**A.** The cross generating function of $Y$ with $I$ is $g_{YI}(z)=A(z)g_I(z)$
(the $\epsilon$ part drops out because $\epsilon\perp I$). Hence the projection of
$Y$ on the entire $I$ process has coefficient function

$$
\frac{g_{YI}(z)}{g_I(z)} = A(z),
$$

which is one-sided in nonnegative powers. **$Y$ does not Granger-cause $I$** — as
expected, since $I$ is exogenous and the only extra information in past $Y$ is
past $\epsilon$, which is orthogonal to $I$.

**B.** Using $g_{cY}(z)=A(z)A(z^{-1})g_\epsilon+(A(z)-1)A(z^{-1})g_I$ and
$g_Y(z)=A(z)A(z^{-1})(g_\epsilon+g_I)$, the projection of $c$ on the whole $Y$
process is

$$
\frac{g_{cY}(z)}{g_Y(z)} = 1 - \frac{g_I(z)}{A(z)\,(g_\epsilon(z)+g_I(z))} ,
$$

which is **two-sided** (the factor $g_I/(g_\epsilon+g_I)$ is a two-sided spectral
ratio), so **$c$ Granger-causes $Y$**. Symmetrically, the projection of $Y$ on the
whole $c$ process, $g_{Yc}(z)/g_c(z)$, is also two-sided, so **$Y$ Granger-causes
$c$**. Thus $c$ and $Y$ Granger-cause each other — the feedback inherent in the
simultaneous consumption–income system. (Confirmed numerically for
$b(L)=0.5L$, white $\epsilon$, $I$ an AR(1).)

**C.** **No.** A projection (regression) equation requires the disturbance to be
orthogonal to every regressor. Here $\epsilon_t$ is orthogonal to *lagged*
$Y$'s, but the consumption function regresses on *current* income $Y_t$, and
$Y_t=A(L)(\epsilon_t+I_t)$ contains $\epsilon_t$, so
$E[\epsilon_t Y_t]=\sigma_\epsilon^2\neq 0$. The contemporaneous regressor is
correlated with the disturbance, so $c_t=b(L)Y_t+\epsilon_t$ is a **structural
(behavioral) equation**, not a least-squares projection; OLS on it would be
inconsistent (simultaneity bias). $\blacksquare$

## Exercise 8

[→ {ref}`Exercise 8 <ex-8>`]

The representation $y_t=a(L)\epsilon_t+ka(L)u_t$, $x_t=c(L)\epsilon_t$ has the
generating functions ($\epsilon\perp u$, variances $\sigma_\epsilon^2,\sigma_u^2$):

$$
g_x = \sigma_\epsilon^2\, c(z)c(z^{-1}),\quad
g_y = (\sigma_\epsilon^2+k^2\sigma_u^2)\,a(z)a(z^{-1}),\quad
g_{yx} = \sigma_\epsilon^2\, a(z)c(z^{-1}).
$$

**A.** Projection of $y$ on $x$:

$$
\frac{g_{yx}(z)}{g_x(z)} = \frac{a(z)}{c(z)} ,
$$

one-sided (both $a$ and $c^{-1}$ are one-sided in nonnegative powers), so **$y$
does not Granger-cause $x$**. Projection of $x$ on $y$:

$$
\frac{g_{yx}(z^{-1})}{g_y(z)}
= \frac{\sigma_\epsilon^2}{\sigma_\epsilon^2+k^2\sigma_u^2}\,\frac{c(z)}{a(z)} ,
$$

also one-sided, so **$x$ does not Granger-cause $y$**. **Neither causes the
other.** (The original representation is already lower-triangular — $x$ loads only
on $\epsilon$ — which by Sims's theorem 1 means $y\not\Rightarrow x$; part D
exhibits the rotated triangular form giving $x\not\Rightarrow y$.)

**B.** The coefficient generating function of the projection of $y$ on the entire
$x$ process is $\;b(z)=g_{yx}(z)/g_x(z)=a(z)/c(z)$, one-sided.

**C.** The projection of $x$ on the entire $y$ process is
$\;f(z)=g_{yx}(z^{-1})/g_y(z)=\dfrac{\sigma_\epsilon^2}{\sigma_\epsilon^2+k^2\sigma_u^2}\dfrac{c(z)}{a(z)}$,
also one-sided.

**D.** Following the hint, set $\eta_{1t}=\epsilon_t+ku_t$ and let $\eta_{2t}$ be
the residual in projecting $\epsilon_t$ on $\eta_{1t}$,
$\epsilon_t=\rho\,\eta_{1t}+\eta_{2t}$ with
$\rho=\dfrac{E\epsilon_t\eta_{1t}}{E\eta_{1t}^2}=\dfrac{\sigma_\epsilon^2}{\sigma_\epsilon^2+k^2\sigma_u^2}$,
so $\eta_1\perp\eta_2$. Then $y_t=a(L)(\epsilon_t+ku_t)=a(L)\eta_{1t}$ and
$x_t=c(L)\epsilon_t=\rho\,c(L)\eta_{1t}+c(L)\eta_{2t}$, i.e.

$$
\begin{bmatrix} y_t \\ x_t \end{bmatrix}
= \begin{bmatrix} a(L) & 0 \\ \rho\,c(L) & c(L) \end{bmatrix}
  \begin{bmatrix} \eta_{1t} \\ \eta_{2t} \end{bmatrix},
\qquad
\mathrm{var}(\eta_1)=\sigma_\epsilon^2+k^2\sigma_u^2,\ \
\mathrm{var}(\eta_2)=\frac{\sigma_\epsilon^2 k^2\sigma_u^2}{\sigma_\epsilon^2+k^2\sigma_u^2}.
$$

This second Wold representation is lower-triangular with $y$ on top ($y$ loads
only on $\eta_1$), confirming $x\not\Rightarrow y$. Together with the original
representation (lower-triangular with $x$ on top), it shows the process admits
*both* triangular factorizations — neither variable Granger-causes the other.
$\blacksquare$

## Exercise 9

[→ {ref}`Exercise 9 <ex-9>`]

Since $p_t=\sum_i w_i p_{t-i}+\epsilon_t$ with $P[\epsilon_t\mid\Omega_{t-1}]=0$,
the price *surprise* is the innovation, $p_t-P[p_t\mid\Omega_{t-1}]=\epsilon_t$,
and the supply curve becomes

$$
y_t = \lambda y_{t-1} + \gamma\epsilon_t + u_t .
$$

**A.** If $u_t$ is serially uncorrelated with $P[u_t\mid\Omega_{t-1}]=0$, then both
$\epsilon_t$ and $u_t$ are orthogonal to $\Omega_{t-1}\supseteq\{y_{t-1},\dots,p_{t-1},\dots\}$.
Projecting on past $y$'s and past $p$'s,

$$
P[y_t\mid y_{t-1},\dots,p_{t-1},\dots]
= \lambda y_{t-1} + \gamma\,P[\epsilon_t\mid\cdot] + P[u_t\mid\cdot]
= \lambda y_{t-1},
$$

which involves **only** $y_{t-1}$ — no past $p$'s. Hence past $p$ does not help
predict $y$ given past $y$: **$p$ fails to Granger-cause $y$**. (Only the
contemporaneous surprise $\epsilon_t$ enters, and it is unforecastable; this holds
for *any* process $p$, not just the Markov one — the Lucas neutrality result.)

**B.** Now $u_t=\rho u_{t-1}+\xi_t$ with $P[\xi_t\mid\Omega_{t-1}]=0$. Then

$$
y_t = \lambda y_{t-1} + \gamma\epsilon_t + \rho u_{t-1} + \xi_t .
$$

From the supply curve one period back, $u_{t-1}=y_{t-1}-\lambda y_{t-2}-\gamma\epsilon_{t-1}$,
and $\epsilon_{t-1}=p_{t-1}-\sum_i w_i p_{t-1-i}$ is a function of past prices.
Both lie in the past information set, while $\epsilon_t,\xi_t\perp$ past. Therefore

$$
P[y_t\mid y_{t-1},\dots,p_{t-1},\dots]
= (\lambda+\rho)\,y_{t-1} - \rho\lambda\,y_{t-2}
  - \rho\gamma\Big(p_{t-1}-\textstyle\sum_i w_i p_{t-1-i}\Big).
$$

The coefficients on past prices are nonzero whenever $\rho\gamma\neq0$, so **past
$p$ helps predict $y$ given past $y$: $p$ Granger-causes $y$**. Serially
correlated supply shocks destroy the neutrality of part A, because past prices now
reveal the persistent component $u_{t-1}$ of the current shock. $\blacksquare$

## Exercise 10

[→ {ref}`Exercise 10 <ex-10>`]

Because $y$ fails to Granger-cause $x$, by Sims's theorem the projection of $y$ on
the whole $x$ process is one-sided: $b(z)=g_{yx}(z)/g_x(z)$ has nonnegative powers
only. Filtering by $y_t^a=f(L)y_t$, $x_t^a=g(L)x_t$ transforms the generating
functions to (formula {eq}`eq-47` style)

$$
g_{y^a x^a}(z)=f(z)g(z^{-1})\,g_{yx}(z),\qquad
g_{x^a}(z)=g(z)g(z^{-1})\,g_x(z),
$$

so the projection of $y^a$ on the whole $x^a$ process has coefficient function

$$
b^a(z)=\frac{g_{y^a x^a}(z)}{g_{x^a}(z)}
=\frac{f(z)}{g(z)}\,\frac{g_{yx}(z)}{g_x(z)}
=\frac{f(z)}{g(z)}\,b(z).
$$

The filters are finite-order, two-sided and symmetric ($f_j=f_{-j}$,
$g_j=g_{-j}$), so $g(z)$ has its zeros in reciprocal pairs and $1/g(z)$ — hence
$f(z)/g(z)$ — is **two-sided** unless $f\equiv g$. Multiplying the one-sided
$b(z)$ by the two-sided $f(z)/g(z)$ generally produces a $b^a(z)$ with nonzero
negative-power coefficients. By Sims's theorem this means **$y^a$ Granger-causes
$x^a$**. Only when the *same* filter is used on both series ($f=g$, so
$f/g\equiv1$ and $b^a=b$) is the one-sidedness — and the absence of
Granger-causality — preserved. This is the chapter's warning (§30) that using
*different* seasonal-adjustment filters on the two series manufactures spurious
Granger-causality. $\blacksquare$
