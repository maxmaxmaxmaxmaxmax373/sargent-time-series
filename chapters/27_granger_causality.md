# The Relationship Between Wiener-Granger Causality and Econometric Exogeneity

We now study the conditions under which the projection of $y_t$ on the entire $x$ process equals the projection of $y_t$ on only current and all past $x$'s, i.e., the condition under which

$$
g_{yx}(z)/g_x(z) = [g_{yx}(z)/d(z^{-1})]_{+}(1/d(z)\sigma^2)
$$

where $g_x(z) = \sigma^2 d(z)d(z^{-1})$. Sims (1972a) proved the important result that these two projections are equal if and only if lagged $y$'s fail linearly to help predict $x$, given lagged $x$'s. We shall first prove Sims's result in a different way than he did. The method is in the spirit of Wiener's derivation of the Wiener-Kolmogorov prediction formula and has certain advantages when it comes to working Exercises 6–13.

We consider a jointly covariance stationary stochastic process $\{y_t,x_t\}$ with $Ex_t = Ey_t = 0$ and with covariance generating functions $g_x(z), g_y(z), g_{yx}(z)$. We assume that $x$ possesses an autoregressive representation and that both $y$ and $x$ are linearly indeterministic. We now consider the projection of $x_t$ on past values of $x$ and past values of $y$:

```{math}
:label: eq-116
x_t = \sum_{j=1}^\infty h_j x_{t-j} + \sum_{j=1}^\infty v_j y_{t-j} + u_t
```

where the least squares residual $u_t$ obeys the orthogonality conditions $E u_t x_{t-\tau} = E u_t y_{t-\tau} = 0$ for $\tau = 1, 2, \ldots$ Solving {eq}`eq-116` for $u_t$ permits the orthogonality conditions to assume the form of the normal equations

$$
E \left\{(x_t - \sum_{j=1}^\infty h_j x_{t-j} - \sum_{j=1}^\infty v_j y_{t-j} )x_{t-\tau}\right\} = 0, \qquad \tau = 1, 2, \ldots
$$

$$
E \left\{(x_t - \sum_{j=1}^\infty h_j x_{t-j} - \sum_{j=1}^\infty v_j y_{t-j} )y_{t-\tau}\right\} = 0, \qquad \tau = 1, 2, \ldots
$$

These equations can be written

```{math}
:label: eq-117
c_x(\tau) = \sum_{j=1}^\infty h_j c_x(\tau - j) + \sum_{j=1}^\infty v_j c_{y x} (\tau - j),
```

```{math}
:label: eq-118
c_{x y}(\tau) = \sum_{j=1}^\infty h_j c_{x y}(\tau - j) + \sum_{j=1}^\infty v_j c_{y} (\tau - j),
```

which are required to hold *only for positive integers* $\tau = 1, 2, \ldots$. Multiplying both sides of {eq}`eq-117` and {eq}`eq-118` by $z^{\tau}$ and summing over all $\tau$, we get the following equation in terms of the $z$ transforms

```{math}
:label: eq-119
g_x(z) - m(z) = h(z)g_x(z) + v(z)g_{yx}(z)
```

```{math}
:label: eq-120
g_{xy}(z) - n(z) = h(z)g_{xy}(z) + v(z)g_y(z)
```

where $m(z)$ and $n(z)$ are each unknown series in *nonpositive powers of $z$* only. That $m(z)$ and $n(z)$ are series in nonpositive powers of $z$ is equivalent with Equations {eq}`eq-117` and {eq}`eq-118` holding only for $\tau \geq 1$. Equations {eq}`eq-119` and {eq}`eq-120` are the normal equations for $h(z)$ and $v(z)$.

Following Wiener, Granger (1969) has proposed the terminology that "$y$ causes $x$" whenever $v(z) \neq 0$. That is, $y$ is said to cause $x$ if, given all past values of $x$, past values of $y$ help to predict $x$. The conditions under which $v(z)$ does or does not equal zero turn out to be of substantial interest to econometricians and macroeconomists, which is the reason that this concept of causality is an interesting one to study.

Consider the projection of $y_t$ on the entire $x$ process,

$$
y_t = \sum_{j=-\infty}^{\infty} b_j x_{t-j} + \epsilon_t
$$

where $E \epsilon_t x_{t-j} = 0$ for all $j$. Under what conditions will the lag distribution $\{ b_j \}$ be one-sided on the past and present, so that $b_j = 0$ for $j < 0$? From formula (45) we have that

```{math}
:label: eq-121
b(z) = g_{yx}(z)/g_x(z).
```

Suppose that $x_t$ has the Wold moving average representation

$$
x_t = d(L)\eta_t, \quad \eta_t = x_t - P[x_t | x_{t-1}, \ldots], \quad \sum_{j=0}^\infty d_j^2 < \infty
$$

Then

```{math}
:label: eq-122
g_x(z) = \sigma_{\eta}^2 d(z)d(z^{-1})
```

We have assumed that $x$ possesses an autoregressive representation so that $[d(z)]^{-1}$ is one-sided and square summable in nonnegative powers of $z$.

Combining {eq}`eq-121` and {eq}`eq-122`, we have that $g_{yx}(z)$ can be expressed as

$$
g_{yx}(z) = \sigma_{\eta}^2 b(z) d(z)d(z^{-1})
$$

Thus, if $b(z)$ is one-sided in nonnegative powers of $z$, then we have that $g_{yx}(z)$ can be factored as

```{math}
:label: eq-123
g_{yx}(z) = a(z) d(z^{-1})
```

where $a(z) = \sigma_{\eta}^2 b(z) d(z)$. Here $a(z)$ is one-sided in nonnegative powers of $z$.

We now prove the following important theorem that is due to Sims:

#### Theorem

$v(z) = 0$ if and only if $b(z)$ is one-sided in nonnegative powers of $z$.

*Proof:* First, we shall show that if $b(z)$ is one-sided in nonnegative powers of $z$, then $v(z) = 0$. We begin by noting that the solution of {eq}`eq-119` and {eq}`eq-120` for $(h(z), v(z))$ is unique. Consequently, it will suffice to show that if $b(z)$ is one-sided, then there exist $m(z)$ and $n(z)$, each one-sided in nonpositive powers of $z$, for which $v(z) = 0$ and {eq}`eq-119` and {eq}`eq-120` hold. Choose $m(z) = - d_0 \sigma_{\eta}^2 d(z^{-1})$ and $n(z) = - d_0 a(z^{-1})$. These are each one-sided in nonpositive powers of $z$, as required. Next, recall that if $b(z)$ is one-sided in nonnegative powers of $z$, then $g_{yx}(z) = a(z) d(z^{-1})$. Substituting these into {eq}`eq-119` and {eq}`eq-120` gives

```{math}
:label: eq-119p
\sigma_{\eta}^2 d(z)d(z^{-1}) - \sigma_{\eta}^2 d(z^{-1}) d_0 = h(z) \sigma_{\eta}^2 d(z)d(z^{-1})
```

```{math}
:label: eq-120p
a(z^{-1}) d(z) - d_0 a(z^{-1}) = h(z) a(z^{-1})d(z)
```

Dividing the first equation by $z\sigma_{\eta}^2 d(z^{-1})$ and the second by $z a(z^{-1})$ results in the same equation, namely

```{math}
:label: eq-124
\frac{d(z)}{z} - \frac{d_0}{z} = \frac{h(z)}{z} d(z)
```

Since $h(z) = \sum_{j=0}^\infty h_j z^j$, the right-hand side is one-sided in nonnegative powers of $z$, while $d_0/z$ involves only negative powers of $z$. Applying the annihilation operator $[\,]_{+}$ to both sides of {eq}`eq-124` and rearranging gives

```{math}
:label: eq-125
\frac{h(z)}{z} = \left[\frac{d(z)}{z}\right]_{+}\frac{1}{d(z)},
```

which is just the univariate Wiener-Kolmogorov prediction formula for the $z$-transform of the coefficients of the projection of $x_t$ on its own infinite past. With this choice of $h(z)/z$, it is directly verified that {eq}`eq-124` holds, and therefore that {eq}`eq-119p` and {eq}`eq-120p` also hold. Consequently, with $h(z)/z$ given by {eq}`eq-125` and $v(z) = 0$, {eq}`eq-119` and {eq}`eq-120` hold. This proves that if $b(z)$ is one-sided in nonnegative powers of $z$, then $v(z) = 0$.

Next, we prove the other direction of implication, namely, that if $v(z) = 0$ is the solution of {eq}`eq-119` and {eq}`eq-120`, then $b(z)$ is one-sided in nonnegative powers of $z$. Suppose that {eq}`eq-119` and {eq}`eq-120` are satisfied with $v(z) = 0$. Then equation {eq}`eq-119` implies that

$$
h(z) = z \left[\frac{d(z)}{z}\right]_{+}\frac{1}{d(z)}
$$

Substituting this into {eq}`eq-120` with $v(z) = 0$ gives

$$
g_{xy}(z) + n(z) = z\left[\frac{d(z)}{z}\right]_{+}\frac{1}{d(z)}g_{xy}(z)
$$

Dividing the above equation by $z g_{xy}(z)/d(z)$ gives

$$
\frac{d(z)}{z} + \frac{n(z) d(z)}{z g_{xy}(z)} = \left[\frac{d(z)}{z}\right]_{+}.
$$

But since $d(z)/z - d_0/z = [d(z)/z]_{+}$, the above equation implies that $-d_0 = n(z) d(z)/g_{xy}(z)$, or

```{math}
:label: eq-126
n(z) = - d_0 g_{xy}(z)/d(z)
```

Now, $n(z)$ is a polynomial in nonpositive powers of $z$. Therefore, for {eq}`eq-126` to hold, $g_{xy}(z)$ must be capable of being represented as $g_{xy}(z) = c(z^{-1})d(z)$, where $c(z^{-1})$ is one-sided in nonpositive powers of $z$. Then it follows that $g_{yx}(z) = c(z) d(z^{-1})$, and that the $z$-transform $b(z)$ of coefficients in the projection of $y$ on the $x_t$ process is given by $b(z) = c(z) d(z^{-1})/\sigma_{\eta}^2 d(z)d(z^{-1}) = c(z)/\sigma_{\eta}^2 d(z)$, which is one-sided in nonnegative powers of $z$. This completes the proof that if {eq}`eq-119` and {eq}`eq-120` are satisfied with $v(z) = 0$, then $b(z)$ is one-sided in nonnegative powers of $z$.

In words, Sims's theorem asserts that the projection of $y$ on the entire $x$ process equals the projection of $y$ on current and past $x$'s if and only if $y$ *fails* to Granger cause $x$ (i.e., $y$ fails to help predict $x$).

#### Example

Suppose that

$$
g_x(z) = (1 + 0.4 z)(1 + 0.4 z^{-1}), \quad g_y(z) = (1 - 0.2 z)(1 - 0.2 z^{-1}), \quad g_{yx}(z) = \frac{(1 - 0.2 z)}{1 - 0.9 z^{-1}}.
$$

The projection of $y$ on the entire $x$ process has coefficient generating function

$$
\frac{g_{yx}(z)}{g_x(z)} = \frac{1 - 0.2 z}{(1 - 0.9 z^{-1})(1 + 0.4 z)(1 + 0.4 z^{-1})},
$$

which has nonzero coefficients on negative powers of $z$ (expand the polynomial in $z$ by partial fractions to verify this). Therefore $y$ Granger causes (helps predict) $x$. The projection of $x$ on the entire $y$ process has the coefficient generating function

$$
\frac{g_{xy}(z)}{g_y(z)} = \frac{g_{y x}(z^{-1})}{g_y(z)} = \frac{1 - 0.2 z^{-1}}{(1 - 0.9 z)(1 - 0.2 z^{-1})(1 - 0.2 z)} = \frac{1}{(1 - 0.9 z)(1 - 0.2 z)},
$$

which involves only nonnegative powers of $z$. Therefore $x$ fails to Granger cause $y$ (i.e., given lagged $y$'s, $x$ fails to help predict $y$).

To gain additional perspective on Sims's theorem, we now undertake to indicate the proof that Sims gave. We do this not only because the theorem is very important, but because his proof provides useful practice in projection arguments and useful insights into the nature of bivariate Wold representations.

Let $\begin{bmatrix} x_t \\ y_t \end{bmatrix}$ be a bivariate, jointly covariance stationary stochastic process. Suppose that $\begin{bmatrix} x_t \\ y_t \end{bmatrix}$ is a strictly linearly indeterministic process with mean zero. Under these conditions, the bivariate version of Wold's theorem states that there exists a moving average representation of the $(x_t, y_t)$ process

$$
\begin{bmatrix} x_t \\ y_t \end{bmatrix} =
\begin{bmatrix}
c^{11}(L) & c^{12}(L) \\
c^{21}(L) & c^{22}(L)
\end{bmatrix}
\begin{bmatrix} \epsilon_t \\ u_t \end{bmatrix}
$$

where $c^{ij}(L) = \sum_{k=0}^{\infty} c_k^{ij} L^k$ are square summable polynomials in the lag operator $L$ that are one-sided in nonnegative powers of $L$; $\epsilon_t$ and $u_t$ are serially uncorrelated processes with $E u_t \epsilon_{s} = 0$ for all $t,s$; $E \epsilon_t^2 = \sigma_{\epsilon}^2$, $E u_t^2 = \sigma_u^2$; and where the one-step-ahead prediction errors are given by

$$
x_t - P[x_t |x_{t-1},\ldots, y_{t-1}, \ldots ] = c_0^{11}\epsilon_t + c_0^{12}u_t,
$$

$$
y_t - P[y_t |x_{t-1},\ldots, y_{t-1}, \ldots ] = c_0^{21}\epsilon_t + c_0^{22}u_t,
$$

i.e., $\epsilon$ and $u$ are "jointly fundamental for $x$ and $y$."[^fn-gc-1]

Wold's theorem establishes the sense in which a vector moving average is a general representation for an indeterministic covariance stationary vector process. The theorem can be proved by pursuing the same kind of projection arguments used in proving the univariate version of the theorem. Below, we shall show how to construct a Wold representation from knowledge of the covariograms of $x$ and $y$ and their cross covariogram.

We now make the further assumption that the $(x_t, y_t)$ process has an autoregressive representation. In particular think of constructing a sequence of projections

```{math}
:label: eq-127
\begin{bmatrix} x_t \\ y_t \end{bmatrix} = F_1^n \begin{bmatrix} x_{t-1} \\ y_{t-1} \end{bmatrix} + \cdots + F_n^n \begin{bmatrix} x_{t-n} \\ y_{t-n} \end{bmatrix} + \begin{bmatrix} a_{x t}^n \\ a_{y t}^n \end{bmatrix}
```

where $F_1^n, \ldots, F_n^n$ are $2 \times 2$ matrices of least squares coefficients and we have the orthogonality conditions

$$
E \begin{bmatrix} x_{t-j} \\ y_{t-j} \end{bmatrix} \begin{bmatrix} a_{x t}^n & a_{y t}^n \end{bmatrix} = \begin{bmatrix} 0 & 0 \\ 0 & 0 \end{bmatrix}
$$

for $j = 1, \ldots, n$. We assume that as $n \to \infty$ the $F_j^n$ converge to $F_j$ for each $j$. This is the assumption that $(x_t, y_t)$ possesses an autoregressive representation and is stronger than the conditions required for $(x_t, y_t)$ to have a vector moving average representation. We can write the autoregressive representation for $(x_t, y_t)$ as

$$
\begin{bmatrix} x_t \\ y_t \end{bmatrix} = \sum_{j=1}^\infty F_j \begin{bmatrix} x_{t-j} \\ y_{t-j} \end{bmatrix} + \begin{bmatrix} a_{x t} \\ a_{y t} \end{bmatrix} = F(L) \begin{bmatrix} x_{t-1} \\ y_{t-1} \end{bmatrix} + \begin{bmatrix} a_{x t} \\ a_{y t} \end{bmatrix}, \qquad F(L) = \sum_{j=1}^\infty F_j L^{j-1}
$$

where the random variables $(a_{x t}, a_{y t})$ obey the least squares orthogonality conditions

$$
E \begin{bmatrix} x_{t-j} \\ y_{t-j} \end{bmatrix} \begin{bmatrix} a_{x t} & a_{y t} \end{bmatrix} = \begin{bmatrix} 0 & 0 \\ 0 & 0 \end{bmatrix}
$$

for all $j \geq 1$. The random variables $(a_{x t}, a_{y t})$ are the one-step-ahead errors in predicting $(x_t, y_t)$ from all past values of $x$ and $y$.

Now consider obtaining the following representation for the $(x_t, y_t)$ process:

```{math}
:label: eq-128
A(L) \begin{bmatrix} x_t \\ y_t \end{bmatrix} = \begin{bmatrix} \epsilon_t \\ u_t \end{bmatrix} \quad \text{or} \quad
(A_0 - A_1 L - A_2 L^2 - \cdots) \begin{bmatrix} x_t \\ y_t \end{bmatrix} = \begin{bmatrix} \epsilon_t \\ u_t \end{bmatrix}
```

where $A_j$ is a $2 \times 2$ matrix for each $j$, where $A_0$ is chosen to be lower triangular and $\begin{bmatrix} \epsilon_t \\ u_t \end{bmatrix}$ are pairwise orthogonal processes (at all lags) that are serially uncorrelated. Can we be sure that such a representation can be arrived at, in particular one with $A_0$ being lower triangular and $\epsilon$ and $u$ being orthogonal processes?

The answer is in general yes[^fn-gc-2] as the following argument suggests. Think of projecting $x_t$ against all *lagged* $x$'s and *lagged* $y$'s. This gives the first row of $A(L)$ and gives a least squares residual process $\epsilon_t$ that is by construction orthogonal to all lagged $y$'s and all lagged $x$'s. Next project $y_t$ against *current* and lagged $x$'s and all lagged $y$'s. This gives the second row of $A(L)$ and delivers a disturbance process $u_t$ that is by construction orthogonal to current and lagged $x$'s and lagged $y$'s. This procedure produces an $A_0$ that is lower triangular as required. Further, notice that since $\epsilon_t$ is orthogonal to all lagged $x$'s and $y$'s and since the representation {eq}`eq-128` that we have achieved permits lagged $\epsilon_t$'s and $u_t$'s to be expressed as linear combinations of lagged $x$'s and $y$'s, it follows that $\epsilon_t$ is orthogonal to lagged $u$'s and $\epsilon$'s. A similar argument shows that $u_t$ is orthogonal to lagged $u$'s and $\epsilon$'s. Finally, since by construction $u_t$ is orthogonal to current and lagged $x$'s and lagged $y$'s and since $\epsilon_t$ is by definition a linear combination of current and lagged $x$'s and lagged $y$'s, it follows that $u_t$ and $\epsilon_t$ are orthogonal contemporaneously.

To check that he understands this construction, the reader is invited to verify that it would also be possible to choose $A_0$ to be upper triangular with a new and generally different error process $\begin{bmatrix}u' \\ \epsilon' \end{bmatrix}$ that satisfies the same conditions on second moments that the $\begin{bmatrix} u \\ \epsilon \end{bmatrix}$ process satisfies.

To get {eq}`eq-128` in a form that is useful for studying prediction problems, premultiply by $A_0^{-1}$ to get[^fn-gc-3]

```{math}
:label: eq-129
\begin{bmatrix} x_t \\ y_t \end{bmatrix} = A_0^{-1}[A_1 L + A_2 L^2 + \cdots]\begin{bmatrix} x_t \\ y_t \end{bmatrix} + A_0^{-1}\begin{bmatrix} \epsilon_t \\ u_t \end{bmatrix} = A_0^{-1} H(L)\begin{bmatrix} x_t \\ y_t \end{bmatrix} + A_0^{-1} \begin{bmatrix} \epsilon_t \\ u_t \end{bmatrix}
```

where $H(L) = A_1L + A_2L^2 + \cdots$. The linear least squares prediction of the $\begin{bmatrix}x_t \\ y_t \end{bmatrix}$ process based on all lagged $x$'s and all lagged $y$'s (call it $P_{t-1}\begin{bmatrix}x_t \\ y_t \end{bmatrix}$) from {eq}`eq-129` is then

$$
P_{t-1}\begin{bmatrix}x_t \\ y_t \end{bmatrix} = A_0^{-1}H(L)\begin{bmatrix} x_t \\ y_t \end{bmatrix} = F(L)L\begin{bmatrix}x_t \\ y_t \end{bmatrix}
$$

since by construction $P_{t-1}\begin{bmatrix} \epsilon_t \\ u_t \end{bmatrix}=0$. The one-step-ahead prediction errors in predicting $\begin{bmatrix} x \\ y \end{bmatrix}$ process are given by

$$
A_0^{-1} \begin{bmatrix} \epsilon_t \\ u_t \end{bmatrix}
$$

Thus $x$ prediction errors and $y$ prediction errors are contemporaneously correlated so long as $A_0$ is not diagonal. Notice that since $A_0$ is lower triangular, so is $A_0^{-1}$, so that $\epsilon_t$ is the one-step-ahead prediction error in predicting $x$ from past $x$'s and $y$'s which is what should be expected given the way the $\epsilon_t$ process was constructed above.

If $A_0^{-1}A(L)$ is lower triangular (i.e., the matrix coefficient is lower triangular for *each* power of $L$), then given lagged $x$'s, lagged $y$'s do not help predict current $x$. That is, if $A_0^{-1}A(L)$ is lower triangular, and, therefore, so is $A_0^{-1}H(L)$, then $P_{t-1}x_t$ involves only lagged $x$'s, lagged $y$'s all bearing zero regression coefficients. In the language of Wiener and Granger, $y$ is said to *cause* $x$ if given past $x$'s, past $y$'s help predict current $x$. Thus, the lower triangularity of $A_0^{-1}A(L)$ is equivalent with $y$'s *failing* to cause $x$, in the Wiener-Granger sense.

Given that $A_0^{-1}$ is lower triangular, we now claim the following: $A_0^{-1}A(L)$ is lower triangular if and only if $A(L)^{-1}$ is lower triangular. To show this, suppose first that $A_0^{-1}A(L)$ is lower triangular. Then note

$$
A(L)^{-1} = A(L)^{-1}A_0A_0^{-1}.
$$

But we know that $A(L)^{-1}A_0$, being the inverse of $A_0^{-1}A(L)$, is lower triangular, as is $A_0^{-1}$. Noting that the product of two lower triangular matrices is also lower triangular then proves that $A(L)^{-1}$ is lower triangular.[^fn-gc-4]

Now suppose that $A(L)^{-1}$ is lower triangular. Since $A_0$ is lower triangular, it follows that $A_0^{-1}A(L)$ is lower triangular. So we have proved that $A_0^{-1}A(L)$ is lower triangular if and only if $A(L)^{-1}$ is lower triangular.

This establishes that if $A_0^{-1}A(L)$ is lower triangular, then {eq}`eq-128` can be "inverted" to yield the vector moving average representation

```{math}
:label: eq-130
\begin{bmatrix} x_t \\ y_t \end{bmatrix} = C(L)\begin{bmatrix} \epsilon_t \\ u_t \end{bmatrix}
```

where $A(L)^{-1} = C(L) = C_0 + C_1L + C_2L^2 + \cdots$, $C_j$ being a $2 \times 2$ matrix, and where $C(L)$ is lower triangular. Recall the extensive orthogonality conditions satisfied by $\epsilon$ and $u$; the $\epsilon$ and $u$ process are orthogonal at all lags, even contemporaneously.[^fn-gc-5] Conversely, suppose that a moving average representation of the lower triangular form {eq}`eq-130` exists with $\epsilon_t$ and $u_t$ being serially uncorrelated processes with $E\epsilon_t u_s=0$ for all $t$ and $s$. Then assuming that $C(L)^{-1}$ exists and equals $A(L)$ gives a representation

$$
C(L)^{-1}\begin{bmatrix} x_t \\ y_t \end{bmatrix} = \begin{bmatrix} \epsilon_t \\ u_t \end{bmatrix} \quad \text{or} \quad A(L)\begin{bmatrix} x_t \\ y_t \end{bmatrix} = \begin{bmatrix} \epsilon_t \\ u_t \end{bmatrix}
$$

where $A(L)$ is lower triangular and one-sided on the present and past. It follows then that $y$ fails to Granger cause $x$.

We have now proved Sims's important theorem 1 which states:

Let $(x_t,y_t)$ be a jointly covariance stationary, strictly indeterministic process with mean zero. Then $\{y_t\}$ fails to Granger cause $\{x_t\}$ if and only if there exists a vector moving average representation

$$
\begin{bmatrix} x_t \\ y_t \end{bmatrix} =
\begin{bmatrix}
    C^{11}(L) &  0 \\
    C^{21}(L) & C^{22}(L)  \\
\end{bmatrix}
\begin{bmatrix} \epsilon_t \\ u_t \end{bmatrix}
$$

where $\epsilon_t$ and $u_t$ are serially uncorrelated processes with means zero and $E\epsilon_t u_s=0$ for all $t$ and $s$, and where the one-step-ahead prediction errors ($x_t-P[x_t|x_{t-1},\ldots,y_{t-1},\ldots]$) and ($y_t-P[y_t|x_{t-1},\ldots,y_{t-1},\ldots]$) are linear combinations of $\epsilon_t$ and $u_t$ (Sims, 1972a).

We are now in a position to state a second theorem of Sims that characterizes the relationship between the concept of strict econometric exogeneity and Granger's concept of causality. Sims's theorem is this:

$y_t$ can be expressed as a distributed lag of current and past $x$'s (with no future $x$'s) with a disturbance process that is orthogonal to past, present, and future $x$'s if and only if $y$ does not Granger cause $x$ (Sims, 1972a).

The condition that $y$ can be expressed as a one-sided distributed lag of $x$ with disturbance process that is orthogonal at all lags to the $x$ process is known as the strict econometric exogeneity of $x$ with respect to $y$. In applied work it is important to test for this condition since it is required if various estimators are to have good properties. It is interesting that engineers have long called a relationship in which $y$ is a one-sided (on the present and past) distributed lag of $x$ a "causal" relationship, and that this long-standing use of the word cause should happen to coincide with the failure of $y$ to cause $x$ in the Wiener-Granger sense.

First we prove that $y$'s not Granger causing $x$ implies that $y$ can be expressed as a one-sided distributed lag of $x$ with a disturbance process orthogonal to $x$ at all lags. The lack of Granger causality from $y$ to $x$ is equivalent with $A_0^{-1}A(L)$ being lower triangular. As we have seen, this implies that $C(L)$ in {eq}`eq-130` is lower triangular, so that

```{math}
:label: eq-131
x_t=C^{11}(L)\epsilon_t,
```

```{math}
:label: eq-132
y_t=C^{21}(L)\epsilon_t + C^{22}(L)u_t,
```

where all polynomials in $L$ involve only nonnegative powers of $L$. Inverting {eq}`eq-131` and substituting into {eq}`eq-132` gives[^fn-gc-6]

$$
y_t=C^{21}(L)C^{11}(L)^{-1}x_t + C^{22}(L)u_t,
$$

which expresses $y_t$ as a one-sided distributed lag of $x$ (no negative powers of $L$ enter) with a disturbance process $u_t$ that is orthogonal to $\epsilon_t$ and therefore to $x_t$ at all lags. This proves half of the theorem.

To prove the other half, one would start with a one-sided lag distribution and a moving average representation for $x_t$

$$
y_t = h(L)x_t + \eta_t, \qquad x_t = a(L)\epsilon_t
$$

where by hypothesis $\eta$ is orthogonal to $x$ and therefore to $\epsilon$ at all lags. Then by finding the moving average representation for $\eta_t$, say $\eta_t=m(L)u_t$ where $Eu_t\epsilon_s = 0$ for all $t$, $s$, one gets the lower triangular vector moving average representation

$$
\begin{bmatrix} x_t \\ y_t \end{bmatrix} = C(L)\begin{bmatrix} \epsilon_t \\ u_t \end{bmatrix}
$$

where $C(L)$ is lower triangular. Assuming $C(L)^{-1}$ exists then gives

$$
C(L)^{-1} \begin{bmatrix} x_t \\ y_t \end{bmatrix} = \begin{bmatrix} \epsilon_t \\ u_t \end{bmatrix}
$$

where $C(L)^{-1}$ is lower triangular and say equal to $A(L)$. Multiplying the above equation, which is in the form of {eq}`eq-130`, through by $A_0^{-1}$, which is also lower triangular, then gives

$$
A_0^{-1}A(L)\begin{bmatrix} x_t \\ y_t \end{bmatrix} = A_0^{-1}\begin{bmatrix} \epsilon_t \\ u_t \end{bmatrix}
$$

or

$$
[I - A_0^{-1}A_1L - A_0^{-1}A_2L^2 - \cdots]\begin{bmatrix} x_t \\ y_t \end{bmatrix} = A_0^{-1}\begin{bmatrix} \epsilon_t \\ u_t \end{bmatrix}
$$

The lower triangularity of the matrices on the left and the orthogonality properties of $\epsilon$ and $u$ establish that in this system $y$ does not Granger cause $x$, i.e., $y$ does not help predict $x$ given lagged $x$'s. This proves the other half of Sims's theorem 2.

[^fn-gc-1]: The statement "$\{\epsilon_t\}$ and $\{u_t\}$ are jointly fundamental for $\{x_t\}$ and $\{y_t\}$" means that the one-step-ahead errors in forecasting $(y_t, x_t)$ from past $x$'s and $y$'s are linear combinations of $\epsilon_t$ and $u_t$.

[^fn-gc-2]: We have remarked earlier that the vector moving average representation of a vector process $z_t$ in terms of the vector noise $\eta_t$, $z_t = C(L)\eta_t$, where the components of $\eta_t$ are white noises that are mutually orthogonal at all lags, is a very general representation. An autoregressive representation for $z_t$ can be obtained by inverting the preceding equation to get $A(L)z_t = \eta_t$, where $A(L) = C(L)^{-1}$ which is to say $A(e^{-i\omega}) = C(e^{-i\omega})^{-1}$ for each $\omega$ between $-\pi$ and $\pi$. The autoregressive representation exists provided that $C(e^{-i\omega})$ is invertible at each frequency between $-\pi$ and $\pi$. This condition is a restriction but is one that can usually be assumed in applied work. (For an example of a $C(e^{-i\omega})$ that violates the condition, consider the univariate example $C(e^{-i\omega}) = 1 - e^{-i\omega}$ — the transform of the first difference operator $I - L$ — which equals zero at $\omega = 0$ and so is not invertible there.)

[^fn-gc-3]: Notice that {eq}`eq-129` is identical to {eq}`eq-127` for $n = \infty$, so that we must have $F_j = A_0^{-1}A_j$ and $\begin{bmatrix} a_{x t} \\ a_{y t} \end{bmatrix} = A_0^{-1} \begin{bmatrix} \epsilon_t \\ u_t \end{bmatrix}$. Notice that $(a_{x t}, a_{y t})$ are by the orthogonality conditions serially uncorrelated and uncorrelated with one another at all nonzero lags.

[^fn-gc-4]: To make the argument in terms of ordinary matrices, write $A(e^{-i\omega})^{-1}=A(e^{-i\omega})^{-1}A_0A_0^{-1}$ and note that $A(e^{-i\omega})^{-1}A_0$ is the inverse of the lower triangular matrix $A_0^{-1}A(e^{-i\omega})$ at each frequency and so is lower triangular. It follows that $A(e^{-i\omega})^{-1}$ is lower triangular (at each frequency) being the product of two lower triangular matrices. It then follows that $A_j^{-1}=(1/2\pi)\int^\pi_{-\pi}A(e^{-i\omega})^{-1}e^{i\omega j}d\omega$ is lower triangular for $j = 0, 1, 2, \ldots$

[^fn-gc-5]: Assuming that things have been normalized so that $\epsilon$ and $u$ have unit variances, the spectral density matrix of the ($x$,$y$) process satisfying {eq}`eq-130` is $S(e^{-i\omega})=C(e^{-i\omega})IC(e^{-i\omega})'$ where the prime now denotes both complex conjugation and transposition. Now let $U$ be a $2 \times 2$ unitary matrix, i.e., a $2 \times 2$ matrix satisfying $UU'=U'U=I$ where here the prime again denotes complex conjugation and transposition. Then note that $S(e^{-i\omega})$ can also be represented $S(e^{-i\omega})=C(e^{-i\omega})UIU'C(e^{-i\omega})'=[C(e^{-i\omega})U]I[C(e^{-i\omega})U]'=D(e^{-i\omega})ID(e^{-i\omega})'$ where $D(e^{-i\omega})=C(e^{-i\omega})U$. Thus we have produced a new moving average representation, one with contemporaneously orthogonal disturbances. This proves that a moving average representation is unique only up to multiplication by a unitary matrix. Notice that multiplication of $C(e^{-i\omega})$ by $U$ will, in general, destroy the lower triangularity of $C(e^{-i\omega})$ if $C$ originally has this property.

[^fn-gc-6]: In assuming that $(x_t, y_t)$ has an autoregressive representation we have in effect assumed that $C^{11}(L)$ has an inverse that is one-sided in nonnegative powers of $L$.
