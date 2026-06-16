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

---

## Exercise 11

[→ {ref}`Exercise 11 <ex-11>`]

```{note}
This problem shows that a long distributed lag of prices on money is *fully
consistent* with a frictionless classical model: the lag distribution reflects
the public's *forecasting of future money*, not sticky prices.
```

**The reduced form for $p_t$.** Substitute $y_t=\bar y$ (constant) into the
portfolio-balance schedule and collect the $p_t$ terms:

$$
m_t-p_t=\alpha\bigl(P_tp_{t+1}-p_t\bigr)+\bar y+u_t
\;\Longrightarrow\;
(1-\alpha)\,p_t+\alpha\,P_tp_{t+1}=m_t-\bar y-u_t .
$$

Divide by $(1-\alpha)$ (which is positive, since $\alpha<0$) and write
$\beta\equiv\dfrac{-\alpha}{1-\alpha}=\dfrac{\alpha}{\alpha-1}$. Because
$\alpha<0$ we have $0<\beta<1$, so the equation can be solved *forward*:

$$
p_t-\beta\,P_tp_{t+1}=\frac{1}{1-\alpha}\bigl(m_t-\bar y-u_t\bigr)
\;\Longrightarrow\;
p_t=\frac{1}{1-\alpha}\sum_{j=0}^\infty\beta^{\,j}\,P_t\bigl(m_{t+j}-\bar y-u_{t+j}\bigr).
$$

The constant $\bar y$ contributes only a constant $-\bar y\sum_j\beta^j/(1-\alpha)
=-\bar y/[(1-\alpha)(1-\beta)]$, which is absorbed into the intercept $a$. The
shocks $u_{t+j}$ satisfy $Eu_tm_{t-s}=0$ for all $s$, so the entire
$u$-contribution is orthogonal to the money process and is absorbed into the
regression residual $\epsilon_t$. The **money component** of $p_t$ is therefore

$$
p_t^{\,m}=\frac{1}{1-\alpha}\sum_{j=0}^\infty\beta^{\,j}\,P_t\,m_{t+j}.
$$

**Applying the geometric-lead (Hansen–Sargent) formula.** With
$m_t=d(L)e_t$, $e_t$ the fundamental innovation, the one-step forecasts are
$P_tm_{t+j}=\sum_{i\ge0}d_{i+j}e_{t-i}$, so the geometrically weighted sum is the
one-sided filter (formula {eq}`eq-90`)

$$
\sum_{j=0}^\infty\beta^{\,j}P_tm_{t+j}=h^\ast(L)\,e_t,\qquad
h^\ast(L)=\frac{L\,d(L)-\beta\,d(\beta)}{L-\beta},
$$

derived as follows. Writing $d(L)=\sum_kd_kL^k$ and changing the order of summation,

$$
\sum_{i\ge0}\Bigl(\sum_{j\ge0}\beta^{\,j}d_{i+j}\Bigr)L^i
=\sum_{k\ge0}d_k\sum_{i=0}^{k}\beta^{\,k-i}L^i
=\sum_{k\ge0}d_k\,\frac{L^{k+1}-\beta^{\,k+1}}{L-\beta}
=\frac{L\,d(L)-\beta\,d(\beta)}{L-\beta},
$$

which is one-sided because the numerator vanishes at $L=\beta$, so $(L-\beta)$
divides it.

**The regression coefficients.** Since $\{m_{t-s}\}$ and $\{e_{t-s}\}$ span the
same Hilbert space (the MA is fundamental), the projection of $p_t$ on current
and past money is $h(L)\,m_t=p_t^{\,m}$, i.e. $h(L)\,d(L)\,e_t=\tfrac1{1-\alpha}
h^\ast(L)e_t$. Hence

$$
\boxed{\;h(L)=\frac{1}{1-\alpha}\;\frac{L\,d(L)-\beta\,d(\beta)}{(L-\beta)\,d(L)},
\qquad \beta=\frac{-\alpha}{1-\alpha}.\;}
$$

**Is the economist right?** No. Classical theory with rational expectations does
*not* imply $h_0=1,\;h_j=0\;(j\neq0)$. It implies the rich distributed lag above,
whose shape is governed jointly by the interest semi-elasticity $\alpha$ *and* by
the money-supply dynamics $d(L)$. The $h_j$ are nonzero for many $j$ precisely
because rational agents forecast future money from its serial correlation. Only
in the knife-edge case of **serially uncorrelated money**, $d(L)=1$, does the
filter collapse: then $h^\ast(L)=\frac{L-\beta}{L-\beta}=1$ and

$$
h(L)=\frac{1}{1-\alpha},\qquad h_0=\frac{1}{1-\alpha}<1,\quad h_j=0\;(j\ge1),
$$

and even here $h_0=1/(1-\alpha)\neq1$. The observed long lag is an artefact of
the money process, not evidence of price stickiness; the macroeconomist's
inference is unwarranted. $\blacksquare$

---

## Exercise 12

[→ {ref}`Exercise 12 <ex-12>`]

**A. The forward solution.** Project the Cagan schedule, written at date
$t+j$ ($j\ge1$), onto information available at $t$. Using
$P_t\bigl(P_{t+j-1}x_{t+j}\bigr)=P_tx_{t+j}$ for $j\ge1$ and
$P_t\eta_{t+j}=0$ for $j\ge1$ (since $P_{t-1}\eta_t=0$ makes $\eta$ a one-step
innovation):

$$
P_t\mu_{t+j}=(1-\alpha)\,P_tx_{t+j}+\alpha\,P_tx_{t+j+1},\qquad j\ge1 .
$$

Let $g_j\equiv P_tx_{t+j}$ and $\theta\equiv\dfrac{-\alpha}{1-\alpha}\in(0,1)$
(again because $\alpha<0$). Dividing by $(1-\alpha)$,

$$
g_j-\theta\,g_{j+1}=\frac{1}{1-\alpha}P_t\mu_{t+j}
\;\Longrightarrow\;
g_j=\frac{1}{1-\alpha}\sum_{k=0}^\infty\theta^{\,k}P_t\mu_{t+j+k}.
$$

Setting $j=1$ and re-indexing $m=k+1$ gives exactly

$$
P_tx_{t+1}=\frac{1}{1-\alpha}\sum_{j=1}^\infty
\Bigl(\frac{-\alpha}{1-\alpha}\Bigr)^{j-1}P_t\mu_{t+j}. \qquad\blacksquare
$$

**B. Rational expectations reproduce Cagan's adaptive formula.** The first row
of the VAR is $x_t=x_{t-1}+a_{1t}-\lambda a_{1,t-1}$, i.e.
$(1-L)x_t=(1-\lambda L)a_{1t}$, an IMA(1,1) with innovation
$a_{1t}=x_t-P_{t-1}x_t=\dfrac{1-L}{1-\lambda L}x_t$. Leading one period,
$x_{t+1}=x_t+a_{1,t+1}-\lambda a_{1t}$, and $P_ta_{1,t+1}=0$, so

$$
\pi_t\equiv P_tx_{t+1}=x_t-\lambda a_{1t}
=x_t-\lambda\frac{1-L}{1-\lambda L}x_t
=\frac{(1-\lambda L)-\lambda(1-L)}{1-\lambda L}x_t
=\frac{1-\lambda}{1-\lambda L}x_t .
$$

This is exactly Cagan's geometric-distributed-lag ("adaptive expectations")
formula — here *derived* from rational expectations, because the IMA(1,1)
structure of $x$ is the one process for which the optimal forecast happens to
coincide with an exponentially weighted moving average. $\blacksquare$

**C. $\mu$ fails to Granger-cause $x$.** The $x$-equation
$x_t=x_{t-1}+a_{1t}-\lambda a_{1,t-1}$ contains *no* $\mu$ terms: $x$ is a
univariate IMA(1,1) driven only by its own innovation $a_{1t}$. The bivariate
system is therefore block-triangular — $\mu$ depends on lagged $x$, but $x$ does
not depend on lagged $\mu$ — so the projection of $x_t$ on the *entire* $\mu$
process adds nothing beyond $x$'s own past. By Sims's criterion, $\mu$ does not
Granger-cause $x$ (while $x$ *does* cause $\mu$). $\blacksquare$

**D. The projection of $\mu_t-x_t$ on the $x$ process versus Cagan's equation.**
First reduce both rows to innovations. From the second row,
$(1-\lambda L)\mu_t=(1-\lambda)x_{t-1}+(1-\lambda L)a_{2t}$, so
$\mu_t=\dfrac{(1-\lambda)L}{1-\lambda L}x_t+a_{2t}
=\dfrac{(1-\lambda)L}{1-L}a_{1t}+a_{2t}$ (using
$x_t=\tfrac{1-\lambda L}{1-L}a_{1t}$). Therefore

$$
\mu_t-x_t=\frac{(1-\lambda)L-(1-\lambda L)}{1-L}a_{1t}+a_{2t}
=\frac{-(1-L)}{1-L}a_{1t}+a_{2t}=a_{2t}-a_{1t},
$$

a *white noise* (verified numerically to machine precision). Now the projection
of $w_t\equiv\mu_t-x_t$ on the whole $x$ process has coefficient generating
function $b(z)=g_{wx}(z)/g_x(z)$. With $w_t=(-1,1)\,a_t$ and
$x_t=\bigl(\tfrac{1-\lambda L}{1-L},0\bigr)a_t$, and
$\Sigma=\operatorname{Var}(a_t)=\left(\begin{smallmatrix}\sigma_1^2&\sigma_{12}\\
\sigma_{12}&\sigma_2^2\end{smallmatrix}\right)$,

$$
g_{wx}(z)=(\sigma_{12}-\sigma_1^2)\frac{1-\lambda z^{-1}}{1-z^{-1}},\qquad
g_x(z)=\sigma_1^2\frac{(1-\lambda z)(1-\lambda z^{-1})}{(1-z)(1-z^{-1})},
$$

so

$$
b(z)=\frac{g_{wx}(z)}{g_x(z)}
=\frac{\sigma_{12}-\sigma_1^2}{\sigma_1^2}\,\frac{1-z}{1-\lambda z}
\;\Longrightarrow\;
P\bigl[\mu_t-x_t\mid x\bigr]
=\frac{\sigma_{12}-\sigma_1^2}{\sigma_1^2}\,\frac{1-L}{1-\lambda L}\,x_t .
$$

Compare Cagan's equation, whose coefficient filter (substituting
$P_tx_{t+1}-P_{t-1}x_t=\tfrac{(1-\lambda)(1-L)}{1-\lambda L}x_t$) is

$$
\mu_t-x_t=\alpha(1-\lambda)\,\frac{1-L}{1-\lambda L}\,x_t+\eta_t .
$$

The two share the *same lag shape* $\dfrac{1-L}{1-\lambda L}$ but **different
scalars**: the projection slope is $(\sigma_{12}-\sigma_1^2)/\sigma_1^2$, whereas
Cagan's structural slope is $\alpha(1-\lambda)$. They coincide only if
$\sigma_{12}-\sigma_1^2=\alpha(1-\lambda)\sigma_1^2$. In general
$\eta_t=a_{2t}-[1+\alpha(1-\lambda)]a_{1t}$ is correlated with $x_t$'s innovation
$a_{1t}$ ($E\eta_ta_{1t}=\sigma_{12}-[1+\alpha(1-\lambda)]\sigma_1^2\neq0$), so
**Cagan's equation is not a valid projection**. Estimating it as a regression
recovers the projection slope, so the implied semi-elasticity is biased to

$$
\hat\alpha=\frac{\sigma_{12}-\sigma_1^2}{(1-\lambda)\,\sigma_1^2}
\neq\alpha ,
$$

a simultaneity (errors-in-the-equation) bias driven by the covariance
$\sigma_{12}$ between the inflation and money-growth innovations. $\blacksquare$

---

## Exercise 13

[→ {ref}`Exercise 13 <ex-13>`]

```{note}
A cost-of-adjustment investment problem with a *general* gestation/depreciation
lag $g(L)$. The Euler equation is a two-sided (symmetric) deterministic
difference equation, solved by the "stable-roots-backward, unstable-roots-
forward" factorization; its characteristic polynomial is the same one Muth
factored in his signal-extraction reading of the permanent-income hypothesis.
```

**A. The Euler equation.** Capital is $K(t)=g(L)I(t)=\sum_{k\ge0}g_kI(t-k)$, so
$\partial K(t)/\partial I(s)=g_{t-s}$ for $t\ge s$ and $0$ otherwise.
Differentiating the present value with respect to $I(s)$ ($\rho\equiv p$, the
fixed output price):

$$
\frac{\partial}{\partial I(s)}
=\sum_{t\ge s}\rho\bigl[A_1-A_2K(t)\bigr]g_{t-s}
-\bigl[B_0+B_1I(s)+\epsilon(s)\bigr]=0 .
$$

Now $\sum_{t\ge s}\rho A_1g_{t-s}=\rho A_1g(1)$, and
$\sum_{t\ge s}\rho A_2K(t)g_{t-s}=\rho A_2\sum_{k\ge0}g_kK(s+k)
=\rho A_2\,g(L^{-1})K(s)=\rho A_2\,g(L^{-1})g(L)I(s)$. Collecting terms,

$$
\bigl\{-B_0-\epsilon(t)+\rho A_1g(1)\bigr\}
=\bigl\{\rho A_2\,g(L^{-1})g(L)+B_1\bigr\}I(t). \qquad\blacksquare
$$

**B. The factorization.** Let
$\Psi(L)\equiv\rho A_2\,g(L^{-1})g(L)+B_1$. It is *symmetric*
($\Psi(L)=\Psi(L^{-1})$), and on the unit circle
$\Psi(e^{-i\omega})=\rho A_2\,|g(e^{-i\omega})|^2+B_1\ge B_1>0$, so $\Psi$ has no
zeros on $|z|=1$ and its zeros occur in reciprocal pairs $(z_i,z_i^{-1})$ with
$|z_i|<1$. Factor

$$
\Psi(L)=\kappa\,\pi(L)\,\pi(L^{-1}),\qquad
\pi(L)=\prod_i(1-z_iL),\;|z_i|<1 ,
$$

with $\pi$ collecting the *stable* roots. The bounded (square-summable) solution
of $\Psi(L)I(t)=\{\cdots\}$ is

$$
I(t)=\frac{1}{\kappa}\,\pi(L)^{-1}\pi(L^{-1})^{-1}
\bigl\{-B_0-\epsilon(t)+\rho A_1g(1)\bigr\},
$$

where $\pi(L)^{-1}=\prod_i(1-z_iL)^{-1}$ is expanded in nonnegative powers of $L$
("stable roots backward") and $\pi(L^{-1})^{-1}$ in nonnegative powers of
$L^{-1}$ ("unstable roots forward"). The forward expansion is legitimate because
$\{\epsilon(t)\}$ is bounded and known (no uncertainty). $\blacksquare$

**C. The geometric case $g(L)=L/(1-\mu L)$.** Here $g(1)=1/(1-\mu)$ and
$g(L)g(L^{-1})=\dfrac{1}{(1-\mu L)(1-\mu L^{-1})}$. Multiplying the Euler
equation by $(1-\mu L)(1-\mu L^{-1})$ gives

$$
\Bigl[\rho A_2+B_1(1-\mu L)(1-\mu L^{-1})\Bigr]I(t)
=(1-\mu L)(1-\mu L^{-1})\Bigl[\frac{\rho A_1}{1-\mu}-B_0-\epsilon(t)\Bigr].
$$

The bracket on the left,
$B_1(1+\mu^2)+\rho A_2-B_1\mu(L+L^{-1})$, is symmetric and factors as
$\dfrac{B_1\mu}{\zeta}(1-\zeta L)(1-\zeta L^{-1})$, where $\zeta\in(0,1)$ is the
stable root of $B_1\mu z^2-[B_1(1+\mu^2)+\rho A_2]z+B_1\mu=0$:

$$
\zeta=\frac{[B_1(1+\mu^2)+\rho A_2]
-\sqrt{[B_1(1+\mu^2)+\rho A_2]^2-4B_1^2\mu^2}}{2B_1\mu},
\qquad \zeta\cdot\zeta^{-1}=1 .
$$

Solving stable-root backward and unstable-root forward,

$$
I(t)=\frac{\zeta}{B_1\mu}\,
\frac{(1-\mu L)(1-\mu L^{-1})}{(1-\zeta L)(1-\zeta L^{-1})}
\Bigl[\frac{\rho A_1}{1-\mu}-B_0-\epsilon(t)\Bigr],
$$

with $(1-\zeta L^{-1})^{-1}=\sum_{k\ge0}\zeta^kL^{-k}$ acting forward on the
(known) shocks. The constant part is
$\bar I=\dfrac{\zeta}{B_1\mu}\dfrac{(1-\mu)^2}{(1-\zeta)^2}
\bigl(\tfrac{\rho A_1}{1-\mu}-B_0\bigr)$, and investment responds to the supply
shocks through the two-sided geometric filter
$-\tfrac{\zeta}{B_1\mu}\tfrac{(1-\mu L)(1-\mu L^{-1})}{(1-\zeta L)(1-\zeta
L^{-1})}\epsilon(t)$. $\blacksquare$

**D. Equivalence with Muth's signal-extraction problem.** The polynomial
factored in C,

$$
\frac{\rho A_2}{(1-\mu L)(1-\mu L^{-1})}+B_1 ,
$$

is identical in form to the spectral density that Muth (1960) factors when
extracting *permanent income* from observed income. Consider

$$
y_t=P_t+u_t,\qquad P_t=\frac{1}{1-\mu L}\,v_t,\qquad
v_t,u_t\ \text{mutually orthogonal white noises},
$$

with $u_t$ transitory and $P_t$ a "permanent" component with geometric
persistence $\mu$. The spectral density of observed income is

$$
g_y(z)=\frac{\sigma_v^2}{(1-\mu z)(1-\mu z^{-1})}+\sigma_u^2 ,
$$

which is exactly the characteristic polynomial of part C under the
correspondence $\rho A_2\leftrightarrow\sigma_v^2$ (the curvature/benefit of
adjusting) and $B_1\leftrightarrow\sigma_u^2$ (the adjustment cost / transitory
noise). The Wiener factorization yields the same stable root $\zeta$ and hence
the same exponentially weighted distributed lag: the optimal-investment problem
is the *dual* of Muth's optimal-prediction (signal-extraction) problem, with the
adjustment-cost ratio $B_1/\rho A_2$ playing the role of the noise-to-signal
ratio $\sigma_u^2/\sigma_v^2$. $\blacksquare$

---

## Exercise 14

[→ {ref}`Exercise 14 <ex-14>`]

```{note}
Coherence is the frequency-by-frequency $R^2$ of the two-sided projection. The
statements printed in the problem have their denominators swapped; the
internally consistent identities (each residual normalized by the spectrum of
the variable being *explained*) are
$\mathrm{coh}=1-g_u/g_x=1-g_\epsilon/g_y$, used throughout below.
```

**A. Coherence as one minus the residual share.** Define coherence
$\mathrm{coh}(\omega)=\dfrac{|g_{yx}(\omega)|^2}{g_x(\omega)g_y(\omega)}$. For the
two-sided projection $x_t=b(L)y_t+u_t$ with $u\perp y$, the coefficient filter is
$B(\omega)=g_{xy}(\omega)/g_y(\omega)$ and spectra add:

$$
g_x=|B(\omega)|^2g_y+g_u
=\frac{|g_{xy}(\omega)|^2}{g_y(\omega)}+g_u
\;\Longrightarrow\;
g_u=g_x\Bigl(1-\frac{|g_{yx}|^2}{g_xg_y}\Bigr)=g_x\bigl(1-\mathrm{coh}\bigr).
$$

Hence $\mathrm{coh}(\omega)=1-\dfrac{g_u(\omega)}{g_x(\omega)}$. Symmetrically,
the projection $y_t=h(L)x_t+\epsilon_t$ with $\epsilon\perp x$,
$H(\omega)=g_{yx}/g_x$, gives
$g_\epsilon=g_y\bigl(1-\mathrm{coh}\bigr)$, i.e.
$\mathrm{coh}(\omega)=1-\dfrac{g_\epsilon(\omega)}{g_y(\omega)}$. $\blacksquare$

**B. $R^2$ of the $y$-on-$x$ equation.** Its $R^2$ is
$1-E\epsilon^2/Ey^2$. By the spectral representation of variances,
$E\epsilon^2=\frac1{2\pi}\int_{-\pi}^\pi g_\epsilon\,d\omega
=\frac1{2\pi}\int_{-\pi}^\pi g_y(1-\mathrm{coh})\,d\omega$ and
$Ey^2=\frac1{2\pi}\int_{-\pi}^\pi g_y\,d\omega$. Therefore

$$
R^2_{xy}=1-\frac{\frac1{2\pi}\int g_y(1-\mathrm{coh})\,d\omega}
{\frac1{2\pi}\int g_y\,d\omega}
=\frac{\frac1{2\pi}\int_{-\pi}^\pi\mathrm{coh}(\omega)\,g_y(\omega)\,d\omega}
{\frac1{2\pi}\int_{-\pi}^\pi g_y(\omega)\,d\omega},
$$

a $g_y$-weighted average of coherence. $\blacksquare$

**C. $R^2$ of the $x$-on-$y$ equation.** Identically, using
$g_u=g_x(1-\mathrm{coh})$,

$$
R^2_{yx}=1-\frac{\frac1{2\pi}\int g_x(1-\mathrm{coh})\,d\omega}
{\frac1{2\pi}\int g_x\,d\omega}
=\frac{\frac1{2\pi}\int_{-\pi}^\pi\mathrm{coh}(\omega)\,g_x(\omega)\,d\omega}
{\frac1{2\pi}\int_{-\pi}^\pi g_x(\omega)\,d\omega},
$$

a $g_x$-weighted average of coherence. In both cases the overall $R^2$ is a
spectrum-weighted average of the frequency-wise coherence, weighted by the power
of the variable being explained. $\blacksquare$

---

## Exercise 15

[→ {ref}`Exercise 15 <ex-15>`]

```{note}
To avoid the notation clash in the printed formula (the same symbols are used for
both the AR poles and the MA roots), write the **AR roots** as $p_1,\dots,p_n$
($A(L)=\prod_k(1-p_kL)$, $|p_k|<1$) and the **MA roots** as $q_1,\dots,q_m$
($B(L)=\prod_j(1-q_jL)$, $m\le n$). The printed $\lambda_s$ are the AR poles
$p_s$ and the printed $\mu_j$ are the MA roots $q_j$ (except in the denominator
product $\prod_{j=1}^n$, where they are again the AR poles).
```

The autocovariance generating function is
$g_y(z)=\dfrac{B(z)B(z^{-1})}{A(z)A(z^{-1})}$ (unit-variance innovations). Using
the residue formula {eq}`eq-25` with the symmetric branch that reproduces a
*decaying* covariance, $c_y(\tau)=\sum$ of residues of $g_y(z)\,z^{|\tau|-1}$ at
the poles inside the unit circle. Substituting
$A(z^{-1})=z^{-n}\prod_k(z-p_k)$ and $B(z^{-1})=z^{-m}\prod_j(z-q_j)$,

$$
g_y(z)\,z^{|\tau|-1}
=\frac{B(z)\,\prod_{j}(z-q_j)}{A(z)\,\prod_{k}(z-p_k)}\;
z^{\,n-m+|\tau|-1}.
$$

The poles inside $|z|=1$ are the simple poles at $z=p_s$ (the zeros of
$\prod_k(z-p_k)$; the zeros of $A(z)$ lie at $1/p_k$ outside). For
$|\tau|\ge m-n+1$ — automatic since $m\le n$ — the factor
$z^{\,n-m+|\tau|-1}$ contributes no pole at the origin, so summing the residues
at $z=p_s$,

$$
c_y(\tau)=\sum_{s=1}^n
\frac{p_s^{\,n+|\tau|-m-1}\,
\prod_{j=1}^m(1-q_jp_s)(p_s-q_j)}
{\prod_{k=1}^n(1-p_kp_s)\;\prod_{k=1,\,k\neq s}^n(p_s-p_k)} ,
$$

where $\prod_k(1-p_kp_s)=A(z)\big|_{z=p_s}\cdot$(grouping) comes from $B(p_s)$ and
$A(p_s)$ and $\prod_{k\neq s}(p_s-p_k)$ from differentiating the pole factor.
Relabelling $p_s\to\lambda_s$, $q_j\to\mu_j$ recovers the printed formula

$$
c_y(\tau)=\sum_{s=1}^n
\frac{\lambda_s^{\,n+|\tau|-m-1}\prod_{j=1}^m(1-\mu_j\lambda_s)(\lambda_s-\mu_j)}
{\prod_{j=1}^n(1-\mu_j\lambda_s)\prod_{j=1,\,j\neq s}^n(\lambda_s-\lambda_j)} .
$$

I verified this numerically against the autocovariances obtained by convolving
the MA($\infty$) expansion of $B(L)/A(L)$: for an ARMA(2,1)
($p=\{0.5,0.3\},q=\{0.2\}$) the formula matches $c_y(\tau)$ to machine precision
for all $\tau\ge0$. (When $m=n$ the origin factor reappears at $\tau=0$, so the
$\tau=0$ term then needs the extra residue at $z=0$; for $m<n$ the formula holds
for every $\tau$.) $\blacksquare$

---

## Exercise 16

[→ {ref}`Exercise 16 <ex-16>`]

The coefficient $b_j$ is, by the inverse-$z$-transform formula {eq}`eq-25`, the
residue at the origin of $b(z)\,z^{-j-1}$ — equivalently the coefficient of $z^j$
in the Laurent expansion of $b(z)=\dfrac{1+\mu z}{1-\lambda z}$. The only finite
pole of $b(z)$ is at $z=1/\lambda$, which lies **outside** the unit circle since
$|\lambda|<1$; hence the expansion about $z=0$ converges on a disk reaching the
unit circle and $b$ is one-sided.

* **$j<0$.** In the branch $b(z)\,z^{-j-1}$ the factor $z^{-j-1}=z^{|j|-1}$ has a
  nonnegative power (for $j\le-1$), so there is no pole at the origin, and the
  only pole $z=1/\lambda$ is outside $\Gamma$. The sum of residues inside is
  therefore $0$, so $b_j=0$.

* **$j=0$.** $b(z)z^{-1}=\dfrac{1+\mu z}{z(1-\lambda z)}$ has a simple pole at
  $z=0$ with residue $\dfrac{1+\mu\cdot0}{1-\lambda\cdot0}=1$, so $b_0=1$.

* **$j\ge1$.** The residue at the origin of $b(z)z^{-j-1}$ equals the $z^j$
  Taylor coefficient of $b(z)$. Expanding,

  $$
  \frac{1+\mu z}{1-\lambda z}=(1+\mu z)\sum_{k\ge0}\lambda^kz^k
  =\sum_{k\ge0}\lambda^kz^k+\mu\sum_{k\ge0}\lambda^kz^{k+1},
  $$

  the coefficient of $z^j$ ($j\ge1$) is $\lambda^j+\mu\lambda^{j-1}$.

Collecting,

$$
b_j=\begin{cases}0,&j<0,\\[2pt]1,&j=0,\\[2pt]
\lambda^j+\mu\lambda^{j-1},&j\ge1.\end{cases}\qquad\blacksquare
$$

---

## Exercise 17

[→ {ref}`Exercise 17 <ex-17>`]

The second-order Solow–Pascal generating function $w(z)=(1-Az)^{-2}$, $|A|<1$,
has a single pole of order two at $z=1/A$ (outside $\Gamma$), so its coefficients
$w_j$ are again one-sided and equal to the residue at the origin of
$w(z)\,z^{-j-1}=\dfrac{z^{-j-1}}{(1-Az)^2}$. For $j\ge0$ this is a pole of order
$j+1$; applying the residue formula {eq}`eq-23` with
$\phi(z)=z^{j+1}\,w(z)z^{-j-1}=(1-Az)^{-2}$,

$$
w_j=\operatorname{res}_{z=0}=\frac{\phi^{(j)}(0)}{j!},\qquad
\phi^{(j)}(z)=\frac{d^j}{dz^j}(1-Az)^{-2}=(j+1)!\,A^j(1-Az)^{-(j+2)} .
$$

Evaluating at $z=0$ and dividing by $j!$,

$$
w_j=\frac{(j+1)!\,A^j}{j!}=(j+1)A^j,\qquad j\ge0,\quad w_j=0\ (j<0).
$$

This is exactly the *second-order Pascal (Solow) lag distribution*: weights
proportional to $(j+1)A^j$, which rise then decline — the characteristic
"humped" shape. Comparing with equation (31) of Chapter IX, the lag weights
coincide up to the normalizing constant $w(1)^{-1}=(1-A)^2$ that makes the
weights sum to one: the normalized Solow lag is $(1-A)^2(j+1)A^j$. $\blacksquare$

---

## Exercise 18

[→ {ref}`Exercise 18 <ex-18>`]

```{note}
The covariogram $c(0)=1.25,\;c(\pm1)=-0.5,\;c(|\tau|\ge2)=0$ truncates at lag
one, so $x_t$ is a first-order moving average — exactly the case treated at the
start of §16. The finite-order projection coefficients converge to the
infinite-order autoregressive coefficients.
```

**B (done first, it organizes A). Wold representation by the method of §16.**
A covariogram vanishing for $|\tau|\ge2$ is that of an MA(1), $x_t=(1+\theta
L)\epsilon_t$, with $c(0)=\sigma^2(1+\theta^2)$, $c(1)=\sigma^2\theta$. Hence
$\dfrac{c(1)}{c(0)}=\dfrac{\theta}{1+\theta^2}=-0.4$, i.e.
$0.4\,\theta^2+\theta+0.4=0$, with roots $\theta=-0.5$ and $\theta=-2$. The
**fundamental (invertible)** root is $\theta=-\tfrac12$, giving

$$
x_t=\Bigl(1-\tfrac12L\Bigr)\epsilon_t,\qquad
\sigma^2=\frac{c(0)}{1+\theta^2}=\frac{1.25}{1.25}=1 .
$$

Inverting, $\epsilon_t=\dfrac{1}{1-\frac12L}\,x_t=\sum_{k\ge0}\bigl(\tfrac12\bigr)^k
x_{t-k}$, so the autoregressive representation is

$$
x_t=-\sum_{j=1}^\infty\Bigl(\tfrac12\Bigr)^j x_{t-j}+\epsilon_t,
\qquad A_j=-\Bigl(\tfrac12\Bigr)^j .
$$

**A. Finite-order projections.** Solving the Yule–Walker system
$\Gamma_n A^n=\gamma_n$ (with $\Gamma_n$ the $n\times n$ Toeplitz covariance
matrix and $\gamma_n=(c(1),\dots,c(n))'$) for $n=1,\dots,5$ gives:

| $n$ | $A_1^n$ | $A_2^n$ | $A_3^n$ | $A_4^n$ | $A_5^n$ |
|---|---|---|---|---|---|
| 1 | $-0.40000$ | | | | |
| 2 | $-0.47619$ | $-0.19048$ | | | |
| 3 | $-0.49412$ | $-0.23529$ | $-0.09412$ | | |
| 4 | $-0.49853$ | $-0.24633$ | $-0.11730$ | $-0.04692$ | |
| 5 | $-0.49963$ | $-0.24908$ | $-0.12308$ | $-0.05861$ | $-0.02344$ |

For each fixed $j$, $A_j^n$ increases (in magnitude) toward the limit
$A_j=-(\tfrac12)^j$: $A_1^n\to-0.5$, $A_2^n\to-0.25$, $A_3^n\to-0.125$, …. So the
finite-order coefficients **do** approach the infinite-order autoregressive
coefficients of part B as $n\to\infty$, confirming
$A_j^n\to A_j=-(\tfrac12)^j$.

```python
import numpy as np
cov = lambda k: {0: 1.25, 1: -0.5}.get(abs(k), 0.0)
for n in range(1, 6):
    G = np.array([[cov(i - j) for j in range(n)] for i in range(n)])
    g = np.array([cov(k) for k in range(1, n + 1)])
    print(n, np.linalg.solve(G, g).round(5))
# limit (Wold/AR):  A_j = -(1/2)**j
```

$\blacksquare$

---

## Exercise 19

[→ {ref}`Exercise 19 <ex-19>`]

By the Wiener–Kolmogorov one-step prediction formula (the $k=1$ case of
{eq}`eq-61`–{eq}`eq-62`), with $x_t=c(L)\epsilon_t$ the Wold representation:
projecting $x_{t+1}=\sum_{j\ge0}c_j\epsilon_{t+1-j}$ on
$\{x_t,x_{t-1},\dots\}$ and using $P_t\epsilon_{t+1}=0$,

$$
P[x_{t+1}\mid x_t,\dots]=\sum_{j\ge0}c_{j+1}\epsilon_{t-j}
=\frac{c(L)-c_0}{L}\,\epsilon_t .
$$

The hypothesis $P[x_{t+1}\mid x_t,\dots]=\rho x_t=\rho\,c(L)\epsilon_t$, matched
coefficient-by-coefficient in the (linearly independent) innovations $\epsilon$,
requires the generating-function identity

$$
\frac{c(L)-c_0}{L}=\rho\,c(L)
\;\Longrightarrow\;
c(L)-c_0=\rho L\,c(L)
\;\Longrightarrow\;
c(L)\,(1-\rho L)=c_0 .
$$

With the Wold normalization $c_0=1$,

$$
c(L)=\frac{1}{1-\rho L}. \qquad\blacksquare
$$

That is, the only Wold representation consistent with a geometrically decaying
one-step forecast is the AR(1).

---

## Exercise 20

[→ {ref}`Exercise 20 <ex-20>`]

Project $x_{t+1}$ and $x_{t+2}$ on $\{x_t,x_{t-1},\dots\}$ using the Wold
representation $x_t=c(L)\epsilon_t$ and $P_t\epsilon_{t+i}=0$ for $i\ge1$:

$$
P[x_{t+1}\mid x_t,\dots]=\frac{c(L)-c_0}{L}\,\epsilon_t,\qquad
P[x_{t+2}\mid x_t,\dots]=\sum_{j\ge0}c_{j+2}\epsilon_{t-j}
=\frac{c(L)-c_0-c_1L}{L^2}\,\epsilon_t .
$$

The hypothesis $P[x_{t+2}\mid x_t,\dots]=\rho\,P[x_{t+1}\mid x_t,\dots]$ gives, in
generating-function form,

$$
\frac{c(L)-c_0-c_1L}{L^2}=\rho\,\frac{c(L)-c_0}{L} .
$$

Multiply by $L^2$ and collect the $c(L)$ terms:

$$
c(L)-c_0-c_1L=\rho L\bigl(c(L)-c_0\bigr)
\;\Longrightarrow\;
c(L)\,(1-\rho L)=c_0+(c_1-\rho c_0)L ,
$$

so

$$
c(L)=\frac{c_0+(c_1-\rho c_0)L}{1-\rho L}. \qquad\blacksquare
$$

(Exercise 21 extends the same argument: imposing
$P[x_{t+k}\mid x_t,\dots]=\rho^k\,P[x_{t+1}\mid x_t,\dots]$ for *every* $k\ge1$
adds no restriction beyond the $k=2$ case — the recursion $c_{j+1}=\rho c_j$ for
$j\ge1$ is already implied — so $c(L)$ is the *same* ARMA(1,1) above.)
