# Exercises

```{note}
Worked solutions are collected in the next chapter,
[Solutions to Exercises](38_exercise_solutions.md).
```

(ex-1)=

**1.** (Sims's approximation error formula.) Let $(y_t,x_t)$ be jointly covariance stationary with means of zero. Let the projection of $y_t$ on the $x$ process be

$$
\sum_{j=-\infty}^\infty b_j^0x_{t-j}.
$$

Suppose a researcher fits by least squares

$$
y_t = \sum_{j=-\infty}^\infty b_j^1x_{t-j} + u_t
$$

where $u_t$ is a disturbance and $\{b_j^1\}$ is a constrained parameterization so that $b_j^1$ cannot equal $b_j^0$ for all $j$. Some examples of *commonly* encountered constrained parameterizations are:

(i) truncation: $b_j^1=0$ for $|j| \geq m$, $m$ a fixed positive integer;

(ii) polynomial approximation: $b_j^1 = \alpha_0 + \alpha_1j + \cdots + \alpha_mj^m$, $m$ a fixed positive integer, $\alpha_j$ free;

(iii) Pascal lag distributions (Solow):

$$
b^1(L)=\frac{1}{(1-\lambda L)^r}
$$

where $r$ is a fixed positive integer and $|\lambda|<1$.

(a) Derive Sims's formula, which asserts that in population, least squares picks $b_j^1$ to minimize

$$
\int_{-\pi}^\pi |b^0(e^{-i\omega})-b^1(e^{-i\omega})|^2g_x(e^{-i\omega})d\omega
$$

*Hints*: Write $y_t$ as

$$
y_t = \sum_{j=-\infty}^\infty b_j^0x_{t-j} + \epsilon_t
$$

$E\epsilon_t x_{t-j}=0$ for all $j$. Then show that

$$
E\left(y_t - \sum_{j=-\infty}^\infty b_j^1x_{t-j}\right)^2 = E(z_t^2)
$$

where

$$
z_t = \sum_{j=-\infty}^\infty (b_j^0 - b_j^1)x_{t-j} + \epsilon_t.
$$

(b) Apply formula (33) to calculate the spectrum of $z_t$. (c) Apply formula (20) to calculate the variance of $z_t$ (see Sims, 1972b).

---

(ex-2)=

**2.** ("Optimal" seasonal adjustment via signal extraction.) Suppose that an analyst is interested in estimating $x_t$ but only observes $X_t = x_t + u_t$, where $Ex_tu_t = 0$ for all $t$ and $s$, and where $x_t$ and $u_t$ are both covariance stationary stochastic processes with means of zero and known (to the analyst) covariance generating functions $g_x(z)$ and $g_u(z)$ respectively; $g_u(e^{-i\omega}) > 0$ for all $\omega$, but has most of its power concentrated at seasonal frequencies. The analyst estimates $x_t$ by the projection

$$
\hat{x}_t = \sum_{j=-\infty}^\infty h_jX_{t-j},
$$

the projection of the unknown $x_t$ on the $X_t$ process.

A. Derive a formula for the $h_j$ (use formula (45)).

B. Prove that $g_{\hat{x}}(e^{-i\omega}) < g_x(e^{-i\omega})$ for all $\omega$.

C. Prove that if $g_x(e^{-i\omega})$ is relatively smooth across the seasonal and nonseasonal frequencies, then since $g_u(e^{-i\omega})$ has big peaks at the seasonal frequencies, it follows that $g_{\hat{x}}(e^{-i\omega})$ will have substantial dips at the seasonal frequencies.

---

(ex-3)=

**3.** Let $x_t$ be any covariance stationary stochastic process with $Ex_t = 0$.

A. Prove that there exists a representation

$$
x_t = \sum_{j=0}^\infty c_ju_{t+j} + \theta_t
$$

where $c_0=1$, $\sum c_j^2<\infty$, $Eu_t^2\geq0$, $Eu_tu_s=0$ for $t\neq s$ and $E\theta_tu_s=0$ for all $t$ and $s$; $\theta_t$ is a process that can be predicted arbitrarily well by a linear function of only future values of $x$; and $u_t = x_t - P[x_t|x_{t+1}, x_{t+2},\ldots]$.

B. Prove that $c_j = d_j$ where $d_j$ is the object in Wold's theorem.

C. Does $u_t=\epsilon_t$ where $\epsilon_t$ is the object in Wold's theorem? Does $Eu_t^2=E\epsilon_t^2$?

D. Does $\theta_t=\eta_t$ where $\eta_t$ is the object in Wold's theorem?

---

(ex-4)=

**4.** Consider the "explosive" first-order Markov process $y_t = Ay_{t-1} + \epsilon_t$, $t= 1,2,\ldots$, $A>1$, where $\epsilon_t$ is white noise with mean zero and variance $\sigma_\epsilon^2$, and $y_0$ is given.

A. Prove that for each realization $(\epsilon_1,\epsilon_2,\ldots)$ the $y_t$ process has the representation

$$
y_t = \lambda^t\eta_0 + \frac{1}{1-\lambda^{-1}L}u_t
$$

where $u_t$ is a white noise. Find formulas for $\eta_0$ and $u_t$ in terms of the $\epsilon$ process, $\lambda$, and $y_0$. (*Hint*: solve the difference equation forward and impose the initial condition.)

B. Is the $u_t$ process "fundamental" for $y_t$?

---

(ex-5)=

**5.** Consider the univariate first-order mixed moving average, autoregressive process $z_t = \lambda z_{t-1}+ a_t - \beta a_{t-1}$ where $a_t$ is a fundamental white noise for $z$ and $0 < \beta < 1$, $0 < \lambda < 1$.

A. Write the process in the form {eq}`eq-104`. (*Hint*: try $x_t = (z_t,a_t)'$ and $\epsilon_t = (a_t, a_t)'$.)

B. Use formula {eq}`eq-108` to derive a formula for $P[z_{t+2}|z_t,z_{t-1},\ldots]$. Verify that this answer agrees with the result of applying the Wiener-Kolmogorov formula {eq}`eq-62`.

---

(ex-6)=

**6.** For the processes below, determine whether $x$ Granger causes $y$ and whether $y$ Granger causes $x$.

A.

$$
g_x(z)=\sigma^2_\epsilon\frac{1}{1-0.9z}\frac{1}{1-0.9z^{-1}},\quad g_y(z) =\sigma^2_u(1-0.8z)(1-0.8z^{-1})
$$

$$
g_{yx}(z)=\sigma_{\epsilon u}(1-0.8z)(1+0.5z^{-1})
$$

B.

$$
g_x(z)=\sigma^2_\epsilon(1+0.99z)(1+0.99z^{-1}),\quad g_y(z) =\sigma^2_u\left(\frac{1}{1-0.7z + 0.3z^2}\right)\left(\frac{1}{1-0.7z^{-1}+ 0.3z^{-2}}\right),
$$

$$
g_{yx}(z)=\sigma_{\epsilon u}(1+0.2z)(1+0.99z^{-1})
$$

C.

$$
g_x(z)=\sigma^2_\epsilon\left(\frac{1}{1-0.7z}\right)\left(\frac{1}{1-0.7z^{-1}}\right),\quad g_y(z) =\sigma^2_u\left(\frac{1}{1-0.8z}\right)\left(\frac{1}{1-0.8z^{-1}}\right),
$$

$$
g_{yx}(z)=\sigma_{\epsilon u}\left(\frac{1}{1-0.8z}\right)\left(\frac{1}{1-0.7z^{-1}}\right)
$$

---

(ex-7)=

**7.** Consider the simple Keynesian macroeconomic model

$$
c_t=\sum_{j=0}^\infty b_jY_{t-j} + \epsilon_t,\quad \sum_{j=0}^\infty b_j^2<\infty,\quad c_t + I_t = Y_t
$$

where $c_t$, $Y_t$ and $I_t$ are consumption, GNP, and investment, respectively, all measured as deviations from their means. Here $\epsilon_t$ is a stationary disturbance process that satisfies $E\epsilon_t I_s = E\epsilon_t = 0$ for all $t$ and $s$ and $I_s$ is a stationary stochastic process. Assume that $(I - b(L))$ has a one-sided, square summable inverse in nonnegative powers of $L$.

A. Determine whether $Y$ Granger causes $I$.

B. Determine whether $c$ Granger causes $Y$ and whether $Y$ Granger causes $c$. (*Hint*: solve for $c$ and $Y$, each as "reduced form" functions of $I$ and $\epsilon$, then apply formula (18) to calculate the cross spectrum and use formula (45) to investigate Granger causality.)

C. Is the consumption function a projection (regression) equation?

---

(ex-8)=

**8.** Consider a $(y,x)$ process that has a Wold moving average representation

$$
y_t = a(L)\epsilon_t + ka(L)u_t,\quad x_t = c(L)\epsilon_t,
$$

where $k$ is a constant, $a(L)$ and $c(L)$ are each one-sided on the past and present and square summable, $Eu_t = E\epsilon_t = Eu_t\epsilon_s = 0$ for all $t$ and $s$, and where $\epsilon_t$ and $u_t$ are jointly fundamental for $y$ and $x$. Finally, assume that both $a(L)$ and $c(L)$ are invertible, i.e., have square-summable inverses that are one-sided in nonnegative powers of $L$.

A. Determine whether $y$ Granger causes $x$ and whether $x$ Granger causes $y$.

B. Find the coefficient generating function for the projection of $y$ on the entire $x$ process.

C. Find the coefficient generating function for the projection of $x$ on the entire $y$ process.

D. Obtain a different Wold moving average representation for the $(y, x)$ process. (*Hint*: choose one white noise process as $\eta_{1t} \equiv \epsilon_t + k u_t$ and choose the other as $\eta_{2t}$, the error in the projection of $\epsilon_t$ on $\epsilon_t + ku_t$: $\epsilon_t = \rho(ku_t + \epsilon_t) + \eta_{2t}$ is a least squares disturbance.)

---

(ex-9)=

**9.** Consider Lucas's aggregate supply curve

$$
y_t= \gamma(p_t-P[p_t|\Omega_{t-1}]) + \lambda y_{t-1} + u_t,\quad 0<|\lambda|<1, \gamma>0
$$

where $y$ is the log of real GNP, $p$ the log of the price level, and $u_t$ a stationary random disturbance process. Suppose that $p_t$ follows the Markov process

$$
p_t = \sum_{i=1}^n w_ip_{t-i} + \epsilon_t
$$

where $P[\epsilon_t|\Omega_{t-1}]=0$. Here $\Omega_{t-1}$ is an information set including at least lagged $y$'s and lagged $p$'s.

A. Suppose that $P[u_t|\Omega_{t-1}] = 0$, so that $u_t$ is serially uncorrelated. Prove that $p$ fails to Granger cause $y$. (In fact, this can be proved where $p$ follows any arbitrary stationary stochastic process and is not dependent on $p$ following the equation above.)

B. Now assume the Markov equation for $p$ and suppose that $u_t$ is serially correlated, and in particular that

$$
u_t=\rho u_{t-1} + \xi_t,\quad 0<|\rho|<1
$$

where $P[\xi_t|\Omega_{t-1}]=0$. Prove that $p$ Granger causes $y$ by calculating $P[y_t|y_{t-1},y_{t-2},\ldots,p_{t-1},p_{t-2},\ldots]$.

---

(ex-10)=

**10.** Suppose that $y_t$ fails to Granger cause $x_t$ where both $y$ and $x$ are seasonally unadjusted processes. Suppose that an investigator studies seasonally adjusted processes $y_t^a$ and $x_t^a$ (see Sims, 1974):

$$
y_t^a=f(L)y_t,\quad x_t^a=g(L)x_t
$$

where $f(L)$ and $g(L)$ are each finite-order two-sided, symmetric ($f_j = f_{-j}$, $g_j = g_{-j}$) seasonal adjustment filters chosen so that $y_t^a$ and $x_t^a$ have less power at the seasonal frequencies than do $y_t$ and $x_t$ respectively. Assume that $y_t^a$ and $x_t^a$ are strictly linearly indeterministic, as are $y_t$ and $x_t$.

Prove that if $f(L) \neq g(L)$, then $y_t^a$ in general Granger causes $x_t^a$. (*Hint*: first calculate the coefficient generating function for the projection of $y_t$ on the $x$ process, then calculate the coefficient generating function for the projection of $y_t^a$ on the $x_t^a$ process.)

---

(ex-11)=

**11.** In a recent article, a macroeconomist reported a regression of the log of the price level ($p_t$) on current and past values of the log of the money supply ($m_t$):

$$
p_t = a + \sum_{j=0}^\infty h_jm_{t-j} + \epsilon_t,\quad E\epsilon_t m_{t-j} = 0 \quad \text{for } j\geq0, \quad E\epsilon_t=0
$$

where $\epsilon_t$ is a random disturbance. He found that the $h_j$'s were nonzero for many $j$'s. He concluded that prices are "too sticky" to be explained by an equilibrium model. According to this economist, "classical" macroeconomics implies that $h_0 = 1$ and $h_j = 0$ for $j \neq 0$.

Now consider the following classical macroeconomic model:

$$
m_t - p_t = \alpha(P_tp_{t+1}-p_{t}) + y_t + u_t \quad \text{(portfolio balance schedule)},
$$

$y_t =$ constant (extreme classical full-employment assumption). Here $\alpha < 0$, and $u_t$ is a stationary random process obeying $Eu_tm_{t-s}=0$ for all $t,s$, $Eu_t=0$.

The money supply is exogenous and has moving average representation

$$
m_t=d(L)e_t,\quad e_t=m_t - P[m_t|m_{t-1},m_{t-2},\ldots], \quad \sum_{j=0}^\infty d_j^2 < \infty.
$$

Derive a formula giving the $h(L) = \sum_{j=0}^\infty h_jL^j$ as a function of $\alpha$ and $d(L)$. Is the macroeconomist correct in his interpretation of the implications of classical theory?

---

(ex-12)=

**12.** Let the portfolio balance schedule be Cagan's

$$
\mu_t - x_t = \alpha(P_tx_{t+1} - P_{t-1}x_t) + \eta_t
$$

where $\mu_t$ is the rate of growth of the money supply, $x_t$ is the rate of inflation, and $\eta_t$ satisfies $P_{t-1}\eta_t=0$, where $P_t[y] = P[y|\mu_t,\mu_{t-1},\ldots,x_t,x_{t-1},\ldots]$ in which $y$ is any random variable.

A. Prove that a solution is

$$
P_tx_{t+1} = \frac{1}{1-\alpha}\sum_{j=1}^\infty \left(\frac{-\alpha}{1-\alpha}\right)^{j-1}P_t\mu_{t+j}.
$$

B. Suppose that $(x_t,\mu_t)$ has the bivariate vector moving average, autoregressive representation

$$
\begin{pmatrix} x_t \\ \mu_t \end{pmatrix} = \begin{pmatrix} 1 & 0 \\ 1-\lambda & \lambda \end{pmatrix} \begin{pmatrix} x_{t-1} \\ \mu_{t-1} \end{pmatrix} + \begin{pmatrix} a_{1t} - \lambda a_{1,t-1} \\ a_{2t} - \lambda a_{2,t-1} \end{pmatrix}
$$

where $a_{1t} = x_t - P_{t-1}x_t$, $a_{2t} = \mu_t - P_{t-1}\mu_t$, $|\lambda| < 1$, and $a_{1t}$ and $a_{2t}$ have finite variances and nonzero covariance. Prove that Cagan's formula for the expected rate of inflation $\pi_t$,

$$
\pi_t=\frac{1-\lambda}{1-\lambda L}x_t
$$

is implied by the hypothesis of rational expectations.

C. Prove that $\mu$ fails to Granger cause $x$.

D. Calculate the coefficients in the projection of $\mu_t-x_t$ on the $x_t$ process. Is this projection equation the same as Cagan's equation

$$
\mu_t - x_t = [\alpha(1 - \lambda)/(1 - \lambda L)](1 - L)x_t + \xi_t,
$$

where $\xi_t$ is random? If not, use your formula for the projection equation to determine the biases that would emerge from mistakenly regarding Cagan's equation as a projection equation.

---

(ex-13)=

**13.** *Depreciation, Gestation, and Delivery Lags.* Consider a firm that is a perfect competitor in the market for its one output, and a monopsonist in the market for additions to the stock of the single factor of production that it uses, capital. The firm sells its output at the price $p$, which is constant over time. Output in period $t$, $q(t)$, is produced according to the production function

$$
q(t) = A_0 + A_1K(t) - \frac{A_2}{2} K(t)^2,\quad A_0, A_1, A_2 > 0
$$

where $K(t)$ is the firm's capital stock at the beginning of period $t$. The firm's capital stock is related to its investment decisions at $t$ and earlier by

$$
K(t) = g(L)I(t)
$$

where $I(t)$ is investment in period $t$, and $g(L) = g_0 + g_1L + g_2L^2 + \cdots$. The polynomial $g(L)$ reflects depreciation, delivery, and possible gestation lags. (A common example is the geometrical model $g(L) = L/(1 - (1 - \delta)L)$ where $\delta$ is "the" depreciation rate; the formulation (2) is designed to encompass this as well as all sorts of plausible alternative models of depreciation and gestation-delivery lags.) The firm faces a supply curve for investment of

$$
J(t) = B_0 + \frac{B_1}{2} I(t) + \epsilon(t),\quad B_0, B_1 > 0
$$

where $J(t)$ is the price of investment goods, and where $\epsilon(t)$ is a bounded sequence of shocks to supply. There is no uncertainty. The firm's problem is to maximize undiscounted present value

$$
\sum_{t=0}^\infty\left\{p\left[A_0 + A_1g(L)I(t) - \frac{A_2}{2}[g(L)I(t)]^2 \right] - \left[B_0 + \frac{B_1}{2} I(t) + \epsilon(t)\right]I(t) \right\}
$$

with respect to $\{I(t);t\geq0\}$ and taking as given $\{I(-s);s\geq0\}$ and $\{\epsilon(t);t\geq0\}$.

A. Prove that the Euler equation for this problem is of the form

$$
\{-B_0 - \epsilon(t) + \rho A_1g(1)\} = \{\rho A_2g(L^{-1})g(L) + B_1\}I(t)
$$

B. Prove that if $B_1 > 0$, then the optimal investment plan is obtained by solving "the stable roots of the characteristic polynomial $\{\rho A_2g(L^{-1})g(L) + B_1\}$ backwards, the unstable roots forwards."

C. Solve the special version of this problem that results when

$$
g(L)=\frac{L}{1-\mu L}
$$

where $\mu = 1 - \delta$, $\delta$ = the depreciation rate, $0 < \mu < 1$.

D. Argue that the characteristic polynomial that must be factored in solving problem C is equivalent with the one encountered by John F. Muth in his signal extraction interpretation of the permanent income hypothesis. State a version of that signal extraction problem that leads to this same characteristic polynomial as encountered in C.

---

(ex-14)=

**14.** Let $(y_t, x_t)$ be a jointly covariance stationary process with mean of zero. Consider the projection equations

$$
x_t = \sum_{j=-\infty}^\infty b_j y_{t-j} + u_t,\quad Eu_t y_s = 0 \quad \text{for all } t,s
$$

$$
y_t = \sum_{j=-\infty}^\infty h_j x_{t-j} + \epsilon_t,\quad E\epsilon_t x_s = 0 \quad \text{for all } t,s
$$

Let $g_y(\omega)$, $g_x(\omega)$, $g_u(\omega)$, $g_\epsilon(\omega)$ be the spectral densities of $y$, $x$, $u$, and $\epsilon$, respectively.

A. Prove that the coherence satisfies

$$
\text{coh}(\omega)=1-\frac{g_u(\omega)}{g_y(\omega)}, \qquad \text{coh}(\omega)=1-\frac{g_\epsilon(\omega)}{g_x(\omega)}
$$

B. Prove that $R^2$ in equation (ii) (i.e., $1 - Eu^2/Ey^2$) is given by

$$
R^2_{xy} = \frac{(1/2\pi) \int_{-\pi}^{\pi} \text{coh}(\omega)g_y(\omega)d\omega}{(1/2\pi) \int_{-\pi}^{\pi} g_y(\omega)d\omega}
$$

C. Prove that the $R^2$ in equation (i) is given by

$$
R^2_{yx} = \frac{(1/2\pi) \int_{-\pi}^{\pi} \text{coh}(\omega)g_x(\omega)d\omega}{(1/2\pi) \int_{-\pi}^{\pi} g_x(\omega)d\omega}
$$

---

(ex-15)=

**15.** Let $y_t$ be a mixed moving average, autoregressive process $y_t=(B(L)/A(L))\epsilon_t$, where $\epsilon_t$ is a white noise with unit variance, $B(L) = \prod_{j=1}^m(1 - \lambda_jL)$, $|\lambda_j|<1$ for $j=1,\ldots,m$, $|\mu_k|<1$ for $k=1,\ldots,n$, and $m\leq n$. The autocovariance generating function for $y$ is $g_y(z) = B(z)B(z^{-1})/A(z)A(z^{-1})$. Use formula (25) to establish the formula

$$
c_y(\tau) = \sum_{s=1}^n \frac{\lambda_s^{n+|\tau|-m-1}\prod_{j=1}^m(1 - \mu_j\lambda_s)(\lambda_s-\mu_j)}{\prod_{j=1}^n(1 - \mu_j\lambda_s)\prod_{j=1,j\neq s}^n(\lambda_s-\lambda_j)}
$$

---

(ex-16)=

**16.** Let $b(L)$ be the polynomial in the lag operator $b(L)=(1 + \mu L)/(1-\lambda L)=\sum_{j=-\infty}^{\infty} b_jL^j$ where $|\lambda| < 1$. Use formula (25) to establish that

$$
b_j = \begin{cases}
  0 & j<0 \\
  1 & j=0 \\
  \lambda^j + \mu \lambda^{j-1} & j\geq 1 \\
\end{cases}
$$

---

(ex-17)=

**17.** Consider the generating function of the second-order Solow-Pascal lag distribution $w(z) = 1/(1 - Az)^2$, $|A|< 1$. Use formulas (23) and (25) to evaluate the coefficients of the lag distribution. Compare your results with equation (31) of Chapter IX.

---

(ex-18)=

**18.** Let $x_t$ be a covariance stationary stochastic process with mean zero and covariogram

$$
c(\tau) = \begin{cases}
  1.25 & \tau=0 \\
  -0.5 & \tau = \pm 1 \\
  0 & |\tau| \geq 2 \\
\end{cases}
$$

A. Use a computer to calculate the projections $P[x_t|x_{t-1},\ldots,x_{t-n}] = \sum_{j=1}^n A_j^n x_{t-j}$ for $n = 1, 2, 3, 4, 5$. As $n$ increases, is $A_j^n$ for fixed $j$ seeming to approach a limit? What value seems to be approached?

B. Use the method of Section 16 to find a Wold (fundamental) moving average representation for $x_t$. Invert this moving average representation to obtain the autoregressive representation $x_t = \sum_{j=1}^\infty A_j x_{t-j}$. Do the $A_j^n$ calculated in A seem to be approaching $A_j$ as $n \to \infty$?

---

**Three inverse optimal prediction problems.**

(ex-19)=

**19.** Suppose that $x_t$ is a stochastic process with Wold representation $x_t = c(L)\epsilon_t$, $\epsilon_t=x_t-P[x_t|x_{t-1},\ldots]$. Suppose that $x_t$ satisfies

$$
P[x_{t+1}|x_t,x_{t-1},\ldots] = \rho x_t,\quad |\rho|< 1.
$$

Use the Wiener-Kolmogorov formula to prove that $c(L)$ must be $1/(1 - \rho L)$.

---

(ex-20)=

**20.** Suppose that $x_t$ is a stochastic process with Wold representation $x_t = c(L)\epsilon_t$, where $\epsilon_t=x_t - P[x_t|x_{t-1},\ldots]$. Suppose that $x_t$ is such that

$$
P[x_{t+2}|x_t,x_{t-1},\ldots] = \rho P[x_{t+1}|x_t,x_{t-1},\ldots]
$$

where $|\rho| < 1$. Use the Wiener-Kolmogorov formula to prove that $c(L)$ must be given by

$$
c(L)=\frac{c_0 + (c_1-\rho c_0)L}{1-\rho L}
$$

---

(ex-21)=

**21.** Suppose that $x_t$ is a stochastic process with Wold representation $x_t = c(L)\epsilon_t$, where $\epsilon_t=x_t - P[x_t|x_{t-1},\ldots]$. Suppose that $x_t$ is such that

$$
P[x_{t+k}|x_t,x_{t-1},\ldots] = \rho^k P[x_{t+1}|x_t,x_{t-1},\ldots]
$$

for all $k \geq 1$ where $|\rho| < 1$. Use the Wiener-Kolmogorov formula to prove that $c(L)$ must be given by $c(L)=(c_0 + (c_1 - \rho c_0)L)/(1-\rho L)$.

---

(ex-22)=

**22.** *Seasonality.* Consider a firm that faces the following optimum problem: to maximize

$$
E_0\sum_{t=0}^\infty \left\{f_0n_t - \frac{f_1}{2}n_t^2 - \frac{d}{2}(n_t - n_{t-1})^2 - w_t n_t \right\}
$$

subject to $n_{-1}$ given, and where $f_0$, $f_1$, $d > 0$; $E_t$ is the mathematical expectation conditioned on information known at time $t$. Here $n_t$ is employment at $t$, and $w_t$ is the real wage at $t$. The firm maximizes over linear contingency plans for setting $n_t$ as a function of information available at $t$, which includes $n_{t-1}$ and $\{w_t,w_{t-1},\ldots\}$. The real wage $w_t$ is assumed to follow the Markov process

$$
(1 - \delta L^4)w_t = \epsilon_t,\quad 0 < \delta < 1
$$

where $\epsilon_t$ is a fundamental white noise for $w_t$. The data are *quarterly*.

A. Compute the *spectrum* of the $w_t$ process, and plot it. Argue that $w_t$ is characterized by a strong seasonal.

B. Compute the optimal linear decision rule of the form

$$
n_t=\lambda n_{t-1} + \sum_{j=0}^3h_jw_{t-j},
$$

giving explicit closed-form formulas for the $h_j$'s. Describe the cross-equation restrictions between the parameters of the Markov process for $w_t$ and the parameters of the decision rule.

C. Calculate the projection of $n_t$ against $\{w_t,w_{t-1},\ldots\}$. Then calculate the projection of $n_t$ against the entire $w$ process, $\{\ldots,w_{t+2},w_{t+1},w_t,w_{t-1},\ldots\}$. Argue that $n$ fails to Granger cause $w$.

D. Because of the strong seasonal in real wages, an analyst forms the "seasonally adjusted"

$$
n_t^a=(1-\delta L^4)n_t, \qquad w_t^a=(1-\delta L^4)w_t.
$$

Assuming that the economic agents look at the seasonally unadjusted series in forming their decisions, and so behave according to the decision rule, calculate the projection of $n_t^a$ on $\{w_t^a,w^a_{t-1},\ldots\}$.

E. Suppose that economic agents care about the seasonally unadjusted data, and so behave according to the decision rule. But, as often has been the practice, the econometrician uses the seasonally adjusted data. Let the econometrician construct his model by imagining that the agent is maximizing

$$
E_0\sum_{t=0}^\infty \left\{f_0n^a_t - \frac{f_1}{2}(n_t^a)^2 - \frac{d}{2}(n_t^a - n_{t-1}^a)^2 - w_t^a n_t^a \right\}
$$

subject to $n^a_{-1}$ given and the law of motion for $w_t^a$ implied by the autoregression for $w_t$ and the definition of $w_t^a$.

(i) Derive the decision rule of the form $n_t^a=\tilde{\lambda} n_{t-1}^a + \sum_{j=0}^3\tilde{h}_jw^a_{t-j}$ that the econometrician attributes to the agent. Describe the cross-equation restrictions built into the econometrician's model.

(ii) Use the projection of $n_t^a$ on $w_t^a$ that you calculated earlier to argue that the econometrician's cross-equation restrictions will not describe the agents' behavior. Describe how the econometrician could mistakenly reject the "rational expectations hypothesis" because of his procedure of using seasonally adjusted data.

---

(ex-23)=

**23.** Consider a jointly covariance stationary, linearly indeterministic stochastic process $(y_t,x_t,X_t)$ with means of zero. The variable $y_t$ is determined as

$$
y_t = \alpha\sum_{j=0}^{\infty} \lambda^j P[x_{t+j}|\Omega_t],\quad \alpha > 0, 0 < \lambda < 1
$$

where $\Omega_t = (x_t,x_{t-1},x_{t-2},\ldots)$. There is available a noisy measure of $x_t$, namely $X_t = x_t + u_t$, where $u_t$ is a mean zero stochastic process satisfying $Eu_sx_t = 0$ for all $t$ and $s$, with covariance generating function $g_u(z)$ which is known to an econometrician. The econometrician observes $(y_t, X_t)$ but not $x_t$.

A. Given the generating functions that are known, describe how to obtain a univariate Wold representation for $X_t$: $X_t=d(L)a_t$ where $a_t = X_t - P[X_t|X_{t-1},\ldots]$.

B. Given the known generating functions, give a formula for $\theta(L)$ in the one-sided projection $x_t=\theta(L)X_t + \epsilon_t$ where $E\epsilon_t X_{t-s}=0$ for $s \geq 0$. Define the "one-sided seasonally adjusted" series $\tilde{x}_t=\theta(L)X_t$. How does $\tilde{x}_t$ compare with $\hat{x}$ associated with the two-sided seasonal adjustment procedure described in Exercise 2?

C. Consider the projection of $y_t$ on $X_t,X_{t-1},\ldots$: $y_t=h(L)X_t + s_t$ where $Es_tX_{t-j} = 0$ for $j\geq0$. Prove that

$$
h(L)=\alpha\left[\frac{\theta(L) - \lambda L^{-1} \theta(\lambda)d(\lambda)d(L)^{-1}}{1-\lambda L^{-1}}\right]
$$

D. Using long division on the formula in C, derive explicit formulas for the $h_j$'s in terms of $\lambda$, the $\theta$'s and the parameters of $d(L)$.

E. Describe how to find a Wold representation for $\tilde{x}_t$. Describe special conditions on $\theta(L)$ which imply that $\tilde{x}_t=\theta(L)d(L)a_t$ is a Wold representation for the one-sided seasonally adjusted series $\tilde{x}_t$.

F. Since the econometrician does not have data on $x_t$ but has data on $\tilde{x}_t$, a friend recommends using the model

$$
y_t = \alpha \sum_{j=0}^\infty \lambda^j P[\tilde{x}_{t+j}|\tilde{x}_t,\tilde{x}_{t-1},\ldots] + \tilde{s}_t = f(L)\tilde{x}_t + \tilde{s}_t
$$

where $E\tilde{s}_t\tilde{x}_{t-j}=0$ for $j\geq 0$. Derive an explicit formula for $f(L)$ in terms of the parameters of the Wold representation for $\tilde{x}_t$.

G. Under what circumstances is $f(L)\tilde{x}_t = h(L)X_t$? Interpret these circumstances as describing conditions under which correct restrictions on the $(y_t, X_t)$ process are derived by proceeding as if the agents set $y_t$ as a function of forecasts of $\tilde{x}_{t}$.

H. Collect the restrictions on $(y_t, X_t)$ as: $X_t=d(L)a_t$; $x_t=\theta(L)X_t + \epsilon_t$; $y_t=h(L)X_t + s_t$.

(i) Argue that in general $y$ Granger causes $X$ in the third equation. Prove that $y$ *fails* to Granger cause $x$.

(ii) Argue informally (but precisely) that there is no reason to expect $s_t$ to be serially uncorrelated.

(iii) Let $s_t$ have a Wold representation $s_t = c(L)\eta_t$, where $\eta_t$ is a fundamental white noise for $s_t$. Argue that in general $E\eta_t X_{t-j} \neq 0$ for $j\geq0$. (Therefore, backward filtering designed to whiten the residual destroys the orthogonality with $X_{t-j}$.)

(iv) Prove that $s_t$ has a representation $s_t=c(L^{-1})\tilde{\eta}_t$ where $\tilde{\eta}_t=s_t - P[s_t|s_{t+1},s_{t+2},\ldots]$. Prove that $E\tilde{\eta}_t X_{t-j} = 0$ for $j\geq0$.

(v) Interpret the results in (iii) and (iv) in terms of the considerations motivating Hayashi and Sims.

I. Briefly describe an econometric strategy for estimating $\lambda$ and $\alpha$ from the orthogonality conditions associated with the system of equations.

---

(ex-24)=

**24.** Let $x_t$ have the Wold representation $x_t=c(L)\epsilon_t$, $c(z)=\sum_{j=0}^\infty c_jz^j$, where $\epsilon_t = x_t - P[x_t|x_{t-1},\ldots]$. Assume that $x_t$ has an autoregressive representation so that $a(L)x_t = \epsilon_t$, where $a(L) = c(L)^{-1}$.

Derive a formula for the $h_j$'s in

$$
P_t\left[\sum_{j=0}^\infty \lambda^j x_{t+j+1}\Bigg|x_t,x_{t-1},\ldots\right]=\sum_{j=0}^\infty h_jx_{t-j},\quad |\lambda|\leq 1.
$$

---

(ex-25)=

**25.** Assume that the first difference of a stochastic process $z_t$ is covariance stationary, purely linearly indeterministic, has mean zero and has Wold moving average representation

$$
(1-L)z_t=c(L)\epsilon_t
$$

where $c(L) = \sum_{j=0}^\infty c_jL^j$, $\sum_{j=0}^\infty c_j^2 < \infty$, and $\epsilon_t$ is a fundamental white noise for $(1-L)z_t$. Beveridge and Nelson (1981) define the "permanent component" of $z_t$ as

$$
\tilde{z}_t= z_t + \lim_{n \to \infty}E_t[\Delta z_{t+1} + \ldots + \Delta z_{t+n}]
$$

where $\Delta = (1 - L)$.

A. Prove that

$$
\lim_{n \to \infty}E_t[\Delta z_{t+1} + \ldots + \Delta z_{t+n}] = [(c(L)-c(1))/(L-1)]\epsilon_t
$$

*Hint*: use the Hansen-Sargent formula {eq}`eq-90` or the reasoning used to derive it.

B. Prove that $\tilde{z}_t$ follows the "random walk" $(1 - L)\tilde{z}_t = c(1)\epsilon_t$.

C. Prove that $\tilde{z}_t$ can be obtained from $z_t$ by "filtering" according to the formula

$$
\tilde{z}_t = c(1)c(L)^{-1}z_t.
$$

D. Beveridge and Nelson define the "cyclical" part of the series as $s_t = z_t - \tilde{z}_t$. Show that $s_t$ can be obtained from $z_t$ via the formula

$$
s_t=\left[1-\frac{c(1)}{c(L)}\right]z_t.
$$

E. Note that the preceding formula and that the Wold representation imply that

$$
s_t = \frac{c(L)-c(1)}{1-L}\epsilon_t.
$$

Is this a Wold representation for $s_t$? *Hints*: at this point the reader might want to consult Hansen and Sargent (1980). Note that this is the version of formula {eq}`eq-88` that is obtained by solving Exercise 24, then driving $\lambda$ to unity from below. Does the assumption that $c(L)$ has a square summable inverse in nonnegative powers of $L$ imply that $[c(L) - c(1)]$ has a square summable inverse in nonnegative powers of $L$?

F. In light of your answer to E, describe the effects of filtering a pair of time series $(z_{1t},z_{2t})$, each with the univariate filter applied by Beveridge and Nelson, to obtain cyclical series $(s_{1t},s_{2t})$. Suppose that $z_2$ fails to Granger cause $z_1$. Does it follow that $s_2$ fails to Granger cause $s_1$?
