# Chapter X — Linear Least Squares Projections (Regressions)

The concept of linear regression has many important uses in macroeconomics, several of which we shall illustrate in subsequent chapters. One very important application will be its use in modeling the "signal extraction" problem faced by agents in an environment in which they have imperfect information about a variable affecting their welfare. By using a linear regression, agents can estimate that unobserved variable in a manner that is optimal, in a certain sense. Two leading applications of the signal extraction model in macroeconomics are Robert Lucas's model of the Phillips curve and Milton Friedman's theory of the consumption function. Another use we shall make of linear regression is to characterize and study the optimal control problem facing monetary and fiscal authorities.

(sec-10-1)=
## 1. Linear Least Squares Regression: The Orthogonality Condition

We consider a set of random variables $y$, $x_1$, $x_2, \ldots, x_n$. The population means of this list of random variables are denoted $Ey$, $Ex_1, \ldots, Ex_n$. We assume that these means are finite, as are the population second moments $Ey^2$, $Ex_1^2, \ldots, Ex_n^2$. By the Cauchy-Schwarz inequality the following cross second moments exist and are finite:

$$
\begin{matrix}
Eyx_1 & Ex_1^2 & Ex_1x_2 & \cdots & Ex_1x_n \\
Eyx_2 & Ex_2x_1 & Ex_2^2 & \cdots & Ex_2x_n \\
\vdots & \vdots & & & \\
Eyx_n & Ex_nx_1 & Ex_nx_2 & \cdots & Ex_n^2.
\end{matrix}
$$

Consider estimating the random variable $y$ on the basis of knowing values only for the random variables $x_1, \ldots, x_n$ as well as knowing all of the means and second moments listed above.[^fn-10-1] More specifically, suppose we restrict ourselves to estimating $y$ by the linear function[^fn-10-2] of the $x_i$,

$$
\hat y = a_0 + a_1 x_1 + \cdots + a_n x_n.
$$ (eq-10-1)

[^fn-10-1]: Luenberger (1969) and Naylor and Sell (1982) contain good treatments of inner product spaces, and of the technical results used in this chapter. The reader with some background in econometrics will note that we are *not* studying the "general linear model," (e.g., see Johnston, 1963, Chapter 4), which assumes that the right-hand side $x$ variables are nonstochastic.

[^fn-10-2]: The restriction to a linear function is in general a binding one. It is possible to show that to minimize $E\{y - g(x_1, \ldots, x_n)\}^2$ with respect to the choice of $g(x_1, \ldots, x_n)$, the optimal thing to do is to set $g(x_1, \ldots, x_n) = E[y|x_1, \ldots, x_n]$, the mathematical expectation of $y$ conditional on $x_1, \ldots, x_n$. In general, the mathematical expectation $E[y|x_1, \ldots, x_n]$ is *not* a linear function of $x_1, \ldots, x_n$. In the special case that the variates $(y, x_1, \ldots, x_n)$ follow a multivariate normal distribution, the conditional mathematical expectation $E[y|x_1, \ldots, x_n]$ is linear in $x_1, \ldots, x_n$.

We seek to choose the $a_i$ so that the random variable $\hat y$ is as "close" to $y$ as possible, in the least squares sense that $E(y - \hat y)^2$ is a minimum. Thus, our problem is to minimize

$$
E(y - (a_0 + a_1 x_1 + \cdots + a_n x_n))^2
$$ (eq-10-2)

with respect to $a_0, a_1, \ldots, a_n$. To facilitate the computations, let us define the new (trivial) random variable $x_0 \equiv 1$.

We are now in a position to state the *orthogonality principle*:

A necessary and sufficient set of conditions for $a_0, a_1, \ldots, a_n$ to minimize {eq}`eq-10-2` is

$$
E\{(y - (a_0 + a_1 x_1 + \cdots + a_n x_n))x_i\} = 0, \qquad i = 0, 1, \ldots, n.
$$ (eq-10-3)

Condition {eq}`eq-10-3` says that $E(y - \hat y)x_i = 0$ for all $i$. Two random variables $w$ and $z$ are said to be *orthogonal* if $E(wz) = 0$. Thus, {eq}`eq-10-3` asserts that $y - \hat y$ is orthogonal to each $x_i$, $i = 0, \ldots, n$. The orthogonality principle asserts that the condition $E(y - \hat y)x_i = 0$ for each $i$ uniquely determines $\hat y$. (It will also uniquely determine the $a_i$ if there is no linear dependence among the $x_i$.)

To prove the orthogonality principle, we proceed by minimizing

$$
J = E\left(y - \sum_{i=0}^n a_i x_i\right)^2.
$$

Differentiating $J$ with respect to $a_k$ gives

$$
-\frac{\partial J}{\partial a_k} = 2E\left(y - \sum_{i=0}^n a_i x_i\right)x_k = 0, \qquad k = 0, 1, \ldots, n,
$$

or

$$
Eyx_k - \sum_{i=0}^n a_i Ex_i x_k = 0, \qquad k = 0, 1, \ldots, n.
$$ (eq-10-4)

Clearly, these "normal equations" are equivalent with the orthogonality condition {eq}`eq-10-3`. Let $x$ be the $1 \times (n + 1)$ row vector $x = (x_0, x_1, \ldots, x_n)$ and $a$ the $(n + 1) \times 1$ column vector $a = (a_0, a_1, \ldots, a_n)'$. Then we can write the least squares normal equations {eq}`eq-10-4` compactly as

$$
Ex'y = (Ex'x)a.
$$

Since $Ex'x$ is a nonnegative definite matrix, we are assured by the second-order conditions that the mean squared error $J$ is minimized by choosing the $a$'s to satisfy {eq}`eq-10-4`. This completes the proof of the orthogonality principle.

The orthogonality condition {eq}`eq-10-3` in effect asserts that the "forecast error" $y - \sum_{i=0}^n a_i x_i$ is orthogonal to each of the $x_i$ and therefore is also orthogonal to any linear combination of the $x_i$. Defining the forecast error as $\varepsilon = y - \sum_{i=0}^n a_i x_i$, we therefore have

$$
y = \sum_{i=0}^n a_i x_i + \varepsilon
$$ (eq-10-5)

where $E(\varepsilon \sum_{i=0}^n a_i x_i) = 0$ and $E\varepsilon x_i = 0$ for $i = 0, 1, \ldots, n$.

Thus, {eq}`eq-10-5` decomposes $y$ into orthogonal parts. By virtue of the orthogonality of the random variables $\sum_{i=0}^n a_i x_i$ and $\varepsilon$ we have the decomposition

$$
Ey^2 = E\left(\sum_{i=0}^n a_i x_i\right)^2 + E\varepsilon^2.
$$

The random variable $\sum_{i=0}^n a_i x_i$, where the $a_i$ are chosen to satisfy the least squares orthogonality condition {eq}`eq-10-3`, is called the *projection* of $y$ on $x_0, x_1, \ldots, x_n$. We shall find it convenient to denote the projection of $y$ on $x_0, x_1, \ldots, x_n$ as

$$
\sum_{i=0}^n a_i x_i \equiv P[y | 1, x_1, x_2, \ldots, x_n],
$$

where remember that $x_0 = 1$ identically.

The orthogonality conditions {eq}`eq-10-3` can be readily rearranged in the form of the familiar least squares normal equations. Write {eq}`eq-10-3` explicitly for $i = 0, 1, \ldots, n$ to get the normal equations

$$
\begin{bmatrix} Ey \\ Eyx_1 \\ Eyx_2 \\ \vdots \\ Eyx_n \end{bmatrix} = \begin{bmatrix} 1 & Ex_1 & Ex_2 & \cdots & Ex_n \\ Ex_1 & Ex_1^2 & Ex_1x_2 & \cdots & Ex_1x_n \\ Ex_2 & Ex_1x_2 & Ex_2^2 & \cdots & Ex_2x_n \\ \vdots & & & & \\ Ex_n & Ex_1x_n & Ex_2x_n & \cdots & Ex_n^2 \end{bmatrix} \begin{bmatrix} a_0 \\ a_1 \\ a_2 \\ \vdots \\ a_n \end{bmatrix}.
$$ (eq-10-6)

Assuming that the $(n + 1) \times (n + 1)$ matrix above has an inverse, we have the following explicit equation for the $a_i$:

$$
\begin{bmatrix} a_0 \\ a_1 \\ \vdots \\ a_n \end{bmatrix} = [Ex_i x_j]^{-1} [Eyx_k]
$$ (eq-10-7)

where $[Ex_i x_j]^{-1}$ is the inverse of the matrix with $(i + 1, j + 1)$th element $Ex_i x_j$, and $[Eyx_k]$ is the $(n + 1) \times 1$ vector with $(k + 1)$th element $Eyx_k$.

As an example, consider projecting $y$ against a single variate $x_1$ (as well as the trivial variate $x_0 = 1$). Then {eq}`eq-10-6` becomes

$$
\begin{bmatrix} Ey \\ Eyx_1 \end{bmatrix} = \begin{bmatrix} 1 & Ex_1 \\ Ex_1 & Ex_1^2 \end{bmatrix} \begin{bmatrix} a_0 \\ a_1 \end{bmatrix}.
$$

The solution of these two equations turns out to be

$$
a_0 = Ey - a_1 Ex_1, \qquad a_1 = \frac{E(y - Ey)(x_1 - Ex_1)}{E(x_1 - Ex_1)^2}.
$$

Denote the covariance between $y$ and $x_1$ as $\sigma_{x_1 y} = E\{(y - Ey)(x_1 - Ex_1)\}$ and the variance of $x_1$ as $\sigma_{x_1}^2 = E\{(x_1 - Ex_1)^2\}$. Then the equations for $a_1$ and $a_0$ become the familiar

$$
a_0 = Ey - a_1 Ex_1, \qquad a_1 = \frac{\sigma_{x_1 y}}{\sigma_{x_1}^2}.
$$ (eq-10-8)

(sec-10-2)=
## 2. Recursive Projection

It happens that the simple univariate formulas {eq}`eq-10-8` can be used in a recursive way to assemble projections on many variables, e.g., $P[y | 1, x_1, \ldots, x_n]$. This often affords a computational saving, and also carries insights about sequential learning.

Write the decomposition {eq}`eq-10-5` for $n = 2$ as

$$
y = P[y | 1, x_1, x_2] + \varepsilon, \qquad y = a_0 + a_1 x_1 + a_2 x_2 + \varepsilon
$$ (eq-10-9)

where $E\varepsilon = 0$, $E\varepsilon x_1 = 0$, and $E\varepsilon x_2 = 0$. These three orthogonality conditions ensure that the $a_i$ are the least squares parameter values. Now project both sides of {eq}`eq-10-9` against 1 and $x_1$ to obtain the equation[^fn-10-3]

$$
P[y | 1, x_1] = a_0 + a_1 x_1 + a_2 P[x_2 | 1, x_1].
$$ (eq-10-10)

[^fn-10-3]: The projection operator is *linear* in the sense that for any two random variables $y_1$ and $y_2$ and any real constants $c$ and $d$,

    $$
    P[cy_1 + dy_2 | x_0, x_1, \ldots, x_n] = cP[y_1 | x_0, x_1, \ldots, x_n] + dP[y_2 | x_0, x_1, \ldots, x_n].
    $$

    To prove this equality, simply write out the normal equations for the projections on both sides of the asserted equality.

To get from {eq}`eq-10-9` to {eq}`eq-10-10` we have used the facts that

$$
P[a_0 | 1, x_1] = a_0, \qquad P[x_1 | 1, x_1] = x_1, \qquad P[\varepsilon | 1, x_1] = 0.
$$

The first two came directly from application of the orthogonality principle to the problem of computing the indicated projection. More directly, it is clear that $E\{a_0 - t_0 - t_1 x_1\}^2$ is minimized by setting $t_0 = a_0$ and $t_1 = 0$. Similarly, $E\{x_1 - t_0 - t_1 x_1\}^2$ is minimized by setting $t_0 = 0$ and $t_1 = 1$. The last of the three equalities above comes from noting that from the orthogonality conditions in {eq}`eq-10-9`, $E\varepsilon = E\varepsilon x_1 = 0$. Substituting these into the least squares normal equations for $P[\varepsilon | 1, x_1]$ shows that $P[\varepsilon | 1, x_1] = 0$.

Subtracting {eq}`eq-10-10` from {eq}`eq-10-9` gives

$$
y - P[y | 1, x_1] = a_2(x_2 - P[x_2 | 1, x_1]) + \varepsilon
$$ (eq-10-11)

where we repeat that $E\varepsilon = E\varepsilon x_1 = E\varepsilon x_2 = 0$. Let $P[x_2 | 1, x_1] = b_0 + b_1 x_1$. The orthogonality conditions imply that

$$
\begin{aligned}
E[\varepsilon(x_2 - P[x_2 | 1, x_1])] &= E[\varepsilon(x_2 - b_0 - b_1 x_1)] \\
&= E\varepsilon x_2 - b_0 E\varepsilon - b_1 E\varepsilon x_1 = 0.
\end{aligned}
$$

Thus $\varepsilon$ is orthogonal to $x_2 - P[x_2 | 1, x_1]$. The orthogonality principle therefore implies that $a_2(x_2 - P[x_2 | 1, x_1])$ must be the projection of $y - P[y | 1, x_1]$ against $(x_2 - P[x_2 | 1, x_1])$. Thus, {eq}`eq-10-11` can be rewritten

$$
y = P[y | 1, x_1] + P[(y - P[y | 1, x_1]) | (x_2 - P[x_2 | 1, x_1])] + \varepsilon.
$$ (eq-10-12)

Notice that by virtue of the orthogonality conditions on $\varepsilon$ {eq}`eq-10-12` implies

$$
P[y | 1, x_1, x_2] = P[y | 1, x_1] + P[(y - P[y | 1, x_1]) | (x_2 - P[x_2 | 1, x_1])].
$$

Let

$$
\begin{aligned}
P[y | 1, x_1] &= c_0 + c_1 x_1, \\
P[x_2 | 1, x_1] &= b_0 + b_1 x_1, \\
P[(y - P[y | 1, x_1]) | (x_2 - P[x_2 | 1, x_1])] &= d_0 + d_1(x_2 - b_0 - b_1 x_1).
\end{aligned}
$$

(Actually, from {eq}`eq-10-11`, we know that $d_0 = 0$.) Then {eq}`eq-10-12` can be written

$$
\begin{aligned}
y &= c_0 + c_1 x_1 + d_0 + d_1(x_2 - b_0 - b_1 x_1) + \varepsilon \\
&= (c_0 + d_0 - d_1 b_0) + (c_1 - b_1 d_1)x_1 + d_1 x_2 + \varepsilon.
\end{aligned}
$$ (eq-10-13)

Comparing {eq}`eq-10-9` with {eq}`eq-10-13`, we have

$$
a_0 = (c_0 + d_0 - d_1 b_0) = (c_0 - d_1 b_0), \qquad a_1 = (c_1 - b_1 d_1), \qquad a_2 = d_1.
$$ (eq-10-14)

The relations {eq}`eq-10-14` give the coefficients in the bivariate projection $P[y | 1, x_1, x_2]$ in terms of the parameters of three univariate projections.

Equation {eq}`eq-10-12` is a useful description of optimal least squares learning or sequential estimation. If at first we have data only on a variable $x_1$, the linear least squares estimates of $y$ and $x_2$ are $P[y | 1, x_1]$ and $P[x_2 | 1, x_1]$, respectively. If an observation $x_2$ subsequently becomes available, our estimate of $y$ can be improved by adding to $P[y | 1, x_1]$ the projection of the unobserved "forecast error" $y - P[y | 1, x_1]$ on the observed forecast error $x_2 - P[x_2 | 1, x_1]$. So long as these forecast errors are correlated, the new observation on $x_2$ carries information useful for estimating $y$.

By induction (or by suitably interpreting $x$ in {eq}`eq-10-12` as a vector of random variables) it is straightforward to extend {eq}`eq-10-12` to the vector form

$$
P[y | \Omega, x] = P[y | \Omega] + P[(y - P[y | \Omega]) | (x - P[x | \Omega])]
$$ (eq-10-15)

where $\Omega$ is a list of random variables and $x$ a random variable. The practical implication of {eq}`eq-10-12` and {eq}`eq-10-15` is that the multivariable regression {eq}`eq-10-15` can be built up sequentially from a set of univariate regressions.

Somewhat more generally, where $\Omega$ is one vector of random variables and $x = (x_1, \ldots, x_n)$ is another vector of random variables there obtains the recursive projection formula[^fn-10-4]

$$
P[y | \Omega, x] = P[y | \Omega] + P[y - P(y | \Omega) | x - P(x | \Omega)]
$$ (eq-10-15-prime)

where $x - P(x | \Omega)$ is the vector $x_1 - P(x_1 | \Omega), \ldots, x_n - P(x_n | \Omega)$.

[^fn-10-4]: To prove that Equation {eq}`eq-10-15-prime` holds, note that by the orthogonality principle each element of the vector $x - P(x | \Omega)$ is orthogonal to each element of $\Omega$. Further, the linear space spanned by $(\Omega, x)$ is the same linear space spanned by $(\Omega, x - P(x | \Omega))$. Therefore

    $$
    P[y | \Omega, x] = P[y | \Omega] + P[y | x - P(x | \Omega)]. \tag{$*$}
    $$

    Since $P[P(y | \Omega) | x - P(x | \Omega)] = 0$ because $x - P(x | \Omega)$ is orthogonal to $\Omega$ and therefore to $P(y | \Omega)$, formula {eq}`eq-10-15-prime` follows by subtracting $0 = P[P(y | \Omega) | x - P(x | \Omega)]$ from $(*)$.

The recursive relation {eq}`eq-10-15` is the foundation of *Kalman filtering*, a technique widely used by engineers. We shall see how the sequential learning mechanism in {eq}`eq-10-15` was exploited by Lucas in obtaining his model of the Phillips curve.

(sec-10-3)=
## 3. The Law of Iterated Projections

From {eq}`eq-10-15` it is easy to deduce a fundamental relation which has been dubbed the *law of iterated projections*. In particular, project the random variable $P[y | \Omega, x]$ against the set of information $\Omega$. Then this law states

$$
P[P[y | \Omega, x] | \Omega] = P[y | \Omega].
$$

To prove the law, write {eq}`eq-10-15` as

$$
P[y | \Omega, x] = P[y | \Omega] + a(x - P[x | \Omega])
$$

where the fixed number $a$ is given by

$$
a = \frac{E(y - P[y | \Omega])(x - P[x | \Omega])}{E(x - P[x | \Omega])^2}.
$$

Project both sides of the above equation against $\Omega$ and note that by the orthogonality principle $P[(x - P[x | \Omega]) | \Omega] = 0$, and of course $P[P[y | \Omega] | \Omega] = P[y | \Omega]$ since $P[y | \Omega]$ is a linear combination of elements in $\Omega$. This proves the law.

(sec-10-4)=
## 4. The Signal-Extraction Problem

Suppose an agent wants to estimate a random variable $s$ but only "sees" the random variable $x$ which is related to $s$ by

$$
x = s + n
$$

where $Esn = 0$; $Es^2, En^2 < \infty$; $Es = En = 0$. One interpretation is that $x$ differs from $s$ by the measurement error $n$. The linear least squares estimate of $s$ is

$$
P[s | 1, x] = a_0 + a_1 x.
$$

The least squares normal equations become

$$
a_1 = \frac{E(xs)}{Ex^2} = \frac{E((s + n)s)}{E(s + n)^2}, \qquad a_1 = \frac{Es^2}{Es^2 + En^2}, \qquad a_0 = 0.
$$

As a slightly richer example, consider a worker who wants to estimate (the log of) his real wage $w - p$. He "sees" the random variable $w$, but does not see the pertinent $p$ at the time that he makes his decision to work. Suppose that the log wage $w$ and log price $p$ obey

$$
w = z + u, \qquad p = z + v,
$$

and $Ezu = Ezv = Euv = Eu = Ez = Ev = 0$. Here $z$ represents neutral movements in the aggregate price level that leave the real wage unaltered. The variates $u$ and $v$ represent factors calling for real wage changes. The worker's linear least squares estimate of $w - p$ based on observing $w$ and knowing the first and second moments of all random variables is

$$
P[(w - p) | 1, w] = a_0 + a_1 w.
$$

We have

$$
w - p = u - v, \qquad w = z + u.
$$

So the normal equations imply

$$
a_0 = 0, \qquad a_1 = \frac{E[(u - v)(z + u)]}{E[(z + u)^2]}, \qquad a_1 = \frac{Eu^2}{Ez^2 + Eu^2}.
$$

Notice that $0 < a_1 < 1$, and that the greater is $Eu^2/Ez^2$, the closer to unity is $a_1$. That makes sense since the greater is $Eu^2/Ez^2$, the larger is the fraction of variance in $w$ that is due to variation in the real-wage determining factor $u$.

(sec-10-5)=
## 5. Signal Extraction with Dynamics

We now use the recursive projection formula {eq}`eq-10-15-prime` to solve a signal extraction problem that John Muth (1960) used to provide a rationalization for Milton Friedman's formula for permanent income. This will lead us to a version of the Kalman filter.[^fn-10-5]

We consider the structure

$$
\theta_{t+1} = \rho \theta_t + \varepsilon_{t+1}
$$ (eq-10-16)

$$
z_t = c\theta_t + u_t
$$ (eq-10-17)

where $\rho$ and $c$ are scalars, $\varepsilon_{t+1}$ is a random variable satisfying $E\varepsilon_t = 0$ for all $t$, $E\varepsilon_t^2 = \sigma_\varepsilon^2$ for all $t$, and $E\varepsilon_t \varepsilon_{t-s} = 0$ for $s \neq 0$. We assume that $u_t$ is a random variable satisfying $Eu_t = 0$ for all $t$, $Eu_t^2 = \sigma_u^2$ for all $t$, $Eu_t u_{t-s} = 0$ for $s \neq 0$ and $Eu_t \varepsilon_s = 0$ for all $t, s$. Equation {eq}`eq-10-16` states that $\theta_t$ is governed by a first-order linear stochastic difference equation, while equation {eq}`eq-10-17` states that $z_t$ is a linear combination of $\theta_t$ and a "noise" $u_t$. Given this structure, consider the following problem that is to be solved by an agent who knows the values of $(c, \rho, \sigma_u^2, \sigma_\varepsilon^2)$. At time $t$, the agent is supposed to see $(z_t, z_{t-1}, \ldots, z_0)$, but not to have seen $\theta_t$ for any $t$. Thus $\theta_t$ is a "hidden variable." At time 0, before observing $z_0$, the agent is supposed to have an initial idea about the location of $\theta_0$, which can be summarized by saying that he thinks it is distributed with mean $\hat\theta_0$ and variance about $\hat\theta_0$ of $\Sigma_0$. The agent's problem is to calculate $P[\theta_{t+1} | z_t, z_{t-1}, \ldots, z_0]$. Using {eq}`eq-10-15-prime`, we shall derive a convenient recursive formula for this projection.

[^fn-10-5]: Anderson and Moore (1979) is a good source of results on recursive methods for solving prediction and filtering problems.

We define $\hat\theta_{t+1} = P[\theta_{t+1} | z_t, z_{t-1}, \ldots, z_0]$. In {eq}`eq-10-15-prime`, at $t \geq 1$, we let $y = \theta_{t+1}$, $\Omega = (z_{t-1}, z_{t-2}, \ldots, z_0)$, $x = z_t$. Then in light of {eq}`eq-10-16`–{eq}`eq-10-17`, and using $P[\varepsilon_{t+1} | z_t, z_{t-1}, \ldots, z_0] = 0$ and the orthogonality conditions assumed for $(\varepsilon_{t+1}, u_t)$, {eq}`eq-10-15-prime` becomes

$$
\begin{aligned}
P[\theta_{t+1} | z_t, z_{t-1}, \ldots, z_0] &= P[\rho \theta_t | z_{t-1}, z_{t-2}, \ldots, z_0] \\
&\quad + P[(\rho \theta_t - P[\rho \theta_t | z_{t-1}, \ldots, z_0]) | (z_t - P[z_t | z_{t-1}, \ldots, z_0])].
\end{aligned}
$$

or

$$
\hat\theta_{t+1} = \rho \hat\theta_t + P[\rho(\theta_t - \hat\theta_t) | (c(\theta_t - \hat\theta_t) + u_t)]
$$ (eq-10-18)

where in the last line we use $z_t - P[z_t | z_{t-1}, \ldots, z_0] = c\theta_t + u_t - c\hat\theta_t$. Let us define

$$
\Sigma_t = E(\theta_t - \hat\theta_t)^2.
$$ (eq-10-19)

Then we have that

$$
P[\rho(\theta_t - \hat\theta_t) | c(\theta_t - \hat\theta_t) + u_t] = K_t[c(\theta_t - \hat\theta_t) + u_t]
$$ (eq-10-20)

where

$$
K_t = \frac{c\rho\Sigma_t}{c^2\Sigma_t + \sigma_u^2}.
$$ (eq-10-21)

In deriving {eq}`eq-10-21`, we use the orthogonality condition $Eu_t(\theta_t - \hat{\theta}_t) = 0$ which follows from {eq}`eq-10-16`–{eq}`eq-10-17` and the orthogonality conditions imposed on $(\varepsilon_{t+1}, u_t)$. Substituting {eq}`eq-10-20` into {eq}`eq-10-18` and rearranging gives

$$
\hat{\theta}_{t+1} = (\rho - K_t c)\hat{\theta}_t + K_t z_t.
$$ (eq-10-22)

Subtracting {eq}`eq-10-22` from {eq}`eq-10-16` and using {eq}`eq-10-17` gives

$$
\theta_{t+1} - \hat{\theta}_{t+1} = (\rho - K_t c)(\theta_t - \hat{\theta}_t) + \varepsilon_{t+1} - K_t u_t.
$$

Computing variances gives

$$
\Sigma_{t+1} = (\rho - K_t c)^2 \Sigma_t + \sigma_\varepsilon^2 + K_t^2 \sigma_u^2.
$$ (eq-10-23)

Equations {eq}`eq-10-21`, {eq}`eq-10-22` and {eq}`eq-10-23` are to be solved starting from the initial condition $\Sigma_0$ given. These three equations give a convenient recursive solution to our problem. The equations are a scalar version of the famous "Kalman filter."

By analyzing the pair of difference equations {eq}`eq-10-21`, {eq}`eq-10-23`, it is possible to establish the following two properties of the solution. First, for any value of $\rho$ and for any value $c \neq 0$, starting the system from any $\Sigma_0 \geq 0$ leads to a $\Sigma_t$ sequence that converges as $t \to \infty$. Second, for the same range of values of $\rho$ and $c$, the parameter $(\rho - Kc)$, where $K = \lim_{t \to \infty} K_t$, is less than unity in absolute value. This implies that for the infinite history filtering problem, in which the agent is imagined to form $\hat{\theta}_{t+1} = P[\theta_{t+1}|z_t, z_{t-1}, \ldots]$ by projecting on an infinite record of current and past $z$'s, the solution can be represented by the time invariant equation

$$
\hat{\theta}_{t+1} = (\rho - Kc)\hat{\theta}_t + K z_t
$$ (eq-10-24)

where $|(\rho - Kc)| < 1$, and $K$ is the unique stationary solution of {eq}`eq-10-21` and {eq}`eq-10-23` that is associated with a stationary solution for which $\lim \Sigma_t = \Sigma > 0$. Equation {eq}`eq-10-24` can be solved to give a version of Friedman and Cagan's formula

$$
\hat{\theta}_{t+1} = K \sum_{j=0}^{\infty} (\rho - Kc)^j z_{t-j}.
$$ (eq-10-25)

Muth considered the case in which $c = 1$, $\rho = 1$. In this case, {eq}`eq-10-25` becomes exactly the adaptive expectations mechanism that was used by Friedman and Cagan.

Note that the orthogonality conditions imposed on the $(\varepsilon, u)$ process imply that

$$
P[\theta_{t+j}|z_t, z_{t-1}, \ldots] = \rho^{j-1}\hat{\theta}_{t+1} \qquad \text{for} \quad j \geq 1.
$$ (eq-10-26)

(sec-10-6)=
## 6. The Term Structure of Interest Rates

Meiselman's[^fn-10-6] error-learning model of the term structure of interest rates can be described quite compactly and motivated elegantly by using our results on recursive regression. The term structure of interest rates refers to interest rates on assets of similar quality but varying terms to maturity viewed as a function of the yield to maturity. The *yield curve* is a graph of yields to maturity against the maturity. The *yield to maturity* on a bond is defined as the (single) yield that makes the present value of the bond's (expected) stream of payments just equal to the present market price of the bond. The yield to maturity is seen to be equivalent with Keynes's "internal rate of return."

Let $R_{nt}$ be the yield to maturity at time $t$ on a bond that will mature at time $t + n$. Irving Fisher and John Hicks suggested viewing the $n$-period yield as an average of the current one-period yield and a sequence of one-period forward rates, which we approximate by the formula

$$
R_{nt} = n^{-1}[R_{1t} + {}_{t+1}F_{1t} + {}_{t+2}F_{1t} + \cdots + {}_{t+n-1}F_{1t}], \qquad n = 1, 2, \ldots,
$$ (eq-10-27)

where ${}_{t+j}F_{1t}$ is the one-period forward rate that at time $t$ pertains to one-period loans to be made at time $t + j$ and mature at time $t + j + 1$.[^fn-10-7] Equation {eq}`eq-10-27` for $n = 1, 2, \ldots$ actually *defines* the forward rates ${}_{t+j}F_{1t}$ as functions of the observable yields $R_{1t}, R_{2t}, R_{3t}, \ldots$. Thus, using {eq}`eq-10-27`, we have

$$
2R_{2t} = R_{1t} + {}_{t+1}F_{1t} \qquad \text{or} \qquad {}_{t+1}F_{1t} = 2R_{2t} - R_{1t}.
$$

Similarly, we could calculate

$$
{}_{t+j}F_{1t} = (j + 1)R_{j+1, t} - j R_{jt}, \qquad j = 1, 2, \ldots.
$$

Now markets in forward loans (i.e., contracts executed at time $t$ for loans to extend between some times $t + j$ to $t + k, k > j > 0$) do not literally exist, as futures markets do in some commodities like wheat and corn. Fisher and Hicks' point was that it was fruitful to decompose a given long rate into the implicit one-period forward rates composing it. Thus, a loan for two periods made at time $t$ is viewed as a one-period (spot) loan made at time $t$ plus a forward commitment entered into at time $t$ to extend the loan for one additional period at time $t + 1$.

So far, all of this has been tautological since {eq}`eq-10-27` is only a *definition* of forward rates. Hicks added content to {eq}`eq-10-27` by adopting the expectations hypothesis, asserting that speculators would force forward yields into equality with the spot one-period rates that they expect to hold on the dates to which the forward rates pertain:

$$
{}_{t+j}F_{1t} = \hat{R}_{1,t+j}
$$ (eq-10-28)

where $\hat{R}_{1t+j}$ is speculators' forecast of the one-period rate which, as of time $t$, they expect to prevail at time $t + j$. Hicks's argument was that unless {eq}`eq-10-28` held speculators could always increase their expected returns by the appropriate combination of issuing and purchasing debts of various maturities. Thus, suppose that we have the following situation:

$$
\begin{aligned}
(R_{1t} = 0.05, \ R_{2t} = 0.04) \quad &\Rightarrow \quad {}_{t+1}F_{1t} = 0.03; \\
\hat{R}_{1t+1} &= 0.05.
\end{aligned}
$$

Speculators expect one-period rates to remain stable between period one and two, but the two-period rate is below the one-period rate, indicating that the one-period forward rate is below the one-period spot rate. In this situation speculators could increase their expected returns, say by borrowing for two periods (at 0.04 per period for two periods) and putting the proceeds into a one-period bond the first period (at 0.05) and another one-period bond the second period (at a return that they expect to be 0.05). By executing this transaction, the speculator expects to get a net return of 0.01 of the amount borrowed each period. (Notice that the speculator need put up no money of his own. Notice also, however, that the speculator is undertaking a risky investment since the return is not a sure thing—the investor is not in an arbitrage situation where the returns are sure.) Hicks entertained the hypothesis that speculators would dominate the market and force {eq}`eq-10-28` to hold.

To make {eq}`eq-10-28` operational suppose we adopt a version of Muth's (1961) hypothesis of "rational expectations" and assume that $\hat{R}_{1t+j}$ is formed as the projection of $R_{1t+j}$ against current and lagged spot rates $R_{1t}, R_{1t-1}, \ldots$. Then {eq}`eq-10-28` becomes

$$
{}_{t+j}F_{1t} = P[R_{1t+j}|1, R_{1t}, R_{1t-1}, \ldots].
$$ (eq-10-29)

Applying our recursive regression formula {eq}`eq-10-15`, we have

$$
\begin{aligned}
P[R_{1t+j}|1, R_{1t+1}, R_{1t}, \ldots] = {} &P[R_{1t+j}|1, R_{1t}, R_{1t-1}, \ldots] \\
&+ b_j(R_{1t+1} - P[R_{1t+1}|1, R_{1t}, R_{1t-1}, \ldots])
\end{aligned}
$$ (eq-10-30)

Substituting {eq}`eq-10-29` into {eq}`eq-10-30` gives

$$
{}_{t+j}F_{1t+1} - {}_{t+j}F_{1t} = b_j(R_{1t+1} - {}_{t+1}F_{1t}).
$$ (eq-10-31)

Equation {eq}`eq-10-31` is exactly the "error-learning" model proposed and implemented by Meiselman. Notice that {eq}`eq-10-31` has no random disturbance, so that strictly speaking {eq}`eq-10-31` should fit perfectly—the $\bar{R}^2$ statistic should be unity. In this sense high values of the $\bar{R}^2$ statistic confirm the theory. Meiselman estimated {eq}`eq-10-31` for annual U.S. data over the period 1901–1954, found positive and statistically significant $b_j$, and found zero constant terms which {eq}`eq-10-31` predicts (though as John Wood and Reuben Kessel pointed out, models other than Meiselman's might also be consistent with the zero intercept). Meiselman found moderately high values for the $\bar{R}^2$ statistics, though they were not all that close to unity (they ranged from 0.91 for $j = 1$ to 0.34 for $j = 8$).

It is straightforward to convert {eq}`eq-10-31` into a regression equation with a disturbance. To illustrate how, suppose that speculators form their expectations at time $t$ by projecting $R_{1t+j}$ on a set $R_{1t}, R_{1t-1}, \ldots, y_t, y_{t-1}, \ldots$, where $y_t$ is some random variable, distinct from $R_{1t}$, that is useful for forecasting $R_{1t+j}$. We can thus write

$$
{}_{t+j}F_{1t} = P[R_{1t+j}|\Omega_t]
$$

where $\Omega_t = R_{1t}, R_{1t-1}, \ldots, y_t, y_{t-1}, \ldots$. With little additional work, one can deduce the following bivariate version of our recursive learning formula {eq}`eq-10-15`:

$$
\begin{aligned}
P[R_{1t+j}|\Omega_{t+1}] - P[R_{1t+j}|\Omega_t] = {} &\gamma_j(R_{1t+1} - P[R_{1t+1}|\Omega_t]) \\
&+ \beta_j(y_{t+1} - P[y_{t+1}|\Omega_t]),
\end{aligned}
$$ (eq-10-32)

$\gamma_j$ and $\beta_j$ being the regression coefficients in the bivariate regression of $R_{1t+j} - P[R_{1t+j}|\Omega_t]$ regressed against $R_{1t+1} - P[R_{1t+1}|\Omega_t]$ and $y_{t+1} - P[y_{t+1}|\Omega_t]$.

Equation {eq}`eq-10-32` is obviously in the same spirit as {eq}`eq-10-15` and says that forecasts are revised according to the "surprising" information in the new observations $y_{t+1}$ and $R_{1t+1}$. Equation {eq}`eq-10-32` can be written

$$
{}_{t+j}F_{1t+1} - {}_{t+j}F_{1t} = \gamma_j(R_{1t+1} - P[R_{1t+1}|\Omega_t]) + \beta_j(y_{t+1} - P[y_{t+1}|\Omega_t]).
$$ (eq-10-33)

Now simply project the right-hand side of the above equation against $R_{1t+1} - P[R_{1t+1}|\Omega_t]$ to get the representation

$$
\begin{aligned}
\gamma_j(R_{1t+1} - P[R_{1t+1}|\Omega_t]) &+ \beta_j(y_{t+1} - P[y_{t+1}|\Omega_t]) \\
&= \phi_j(R_{1t+1} - P[R_{1t+1}|\Omega_t]) + \varepsilon_{t+1}
\end{aligned}
$$ (eq-10-34)

where

$$
\phi_j = \frac{\beta_j E(R_{1t+1} - P[R_{1t+1}|\Omega_t])(y_{t+1} - P[y_{t+1}|\Omega_t])}{E(R_{1t+1} - P[R_{1t+1}|\Omega_t])^2} + \gamma_j
$$

and where by the orthogonality principle $\varepsilon_{t+1}$ is orthogonal to $R_{1t+1} - P[R_{1t+1}|\Omega_t]$. Substituting Equation {eq}`eq-10-33` into Equation {eq}`eq-10-34` gives

$$
{}_{t+j}F_{1t+1} - {}_{t+j}F_{1t} = \phi_j(R_{1t+1} - {}_{t+1}F_{1t}) + \varepsilon_{t+1},
$$ (eq-10-35)

which is a regression equation that is in the form of Meiselman's error-learning model. The presence of the random term $\varepsilon_{t+1}$ means that there is no implication that Equation {eq}`eq-10-35` will bear a high $\bar{R}^2$ statistic.[^fn-10-8] High values of the $\bar{R}^2$ statistic would indicate that a large proportion of the information useful for forecasting interest rates is included in current and lagged one-period rates.

(sec-10-7)=
## 7. Application of the Law of Iterated Projections

Let us represent the expectations theory of the term structure in the form

$$
{}_{t+j}F_{1t} = P[R_{1t+j}|\Omega_t], \qquad j = 1, 2, \ldots, \quad \text{all integer } t,
$$ (eq-10-36)

where $\Omega_t \supset \Omega_{t-1} \supset \Omega_{t-2} \supset \cdots$, so that $\Omega_t$ is an information set that includes all information available at time $t - 1$ and maybe some additional information. Applying the law of iterated projections to {eq}`eq-10-36` gives

$$
\begin{aligned}
P[{}_{t+j}F_{1t}|\Omega_{t-1}] &= P[P[R_{1t+j}|\Omega_t]|\Omega_{t-1}] \\
&= P[R_{1t+j}|\Omega_{t-1}] = {}_{t+j}F_{t-1}
\end{aligned}
$$

or

$$
P[{}_{t+j}F_{1t}|\Omega_{t-1}] = {}_{t+j}F_{1t-1}.
$$ (eq-10-37)

A sequence of random variables

$$
{}_{t+j}F_{1t-n}, \quad {}_{t+j}F_{1t-n+1}, \quad \cdots, \quad {}_{t+j}F_{1t}
$$ (eq-10-38)

with the property {eq}`eq-10-37` is said to be a (weak) martingale with respect to $\Omega_{t-1}$. It was Samuelson who first pointed out that the rational expectations theory of future markets implies that sequences of forward prices like {eq}`eq-10-38` follow martingales (Samuelson, 1965).[^fn-10-9]

It is useful also to note that $P[({}_{t+j}F_{1t} - {}_{t+j}F_{1t-1})|\Omega_{t-1}] = 0$ as a result of the law of iterated projections. This says that the revisions in the forecast of $R_{1t+j}$ made between $t - 1$ and $t$ cannot be predicted by a linear function of the information in $\Omega_{t-1}$. This is a strong, testable implication of the theory.

(sec-10-8)=
## 8. Conclusions

The tools described in this chapter are workhorses in applied dynamic economics. We shall see in subsequent chapters how in the hands of Milton Friedman, John F. Muth, and Robert E. Lucas, Jr., the dynamic signal extraction problem became an important tool for understanding the way "distributed lags" could arise in economic models.

The next chapter combines aspects of Chapters IX and X to study linear stochastic difference equations. The main idea is to discover a setting in which a sequence of dated scalar- or vector-valued variables, $x_t$, can be thought of as consisting of random variables, for which the second moments of $x_t$ and $x_s$ are defined for all possible $t$ and $s$. In such a setting, the methods of this chapter can be used to calculate, say, the projection of $x_t$ on $[x_{t-1}, \ldots, x_{t-n}]$ for each $n > 0$. Such a projection is known as a "vector autoregression," and is widely used for macroeconomic prediction, and also to create models of agents whose behavior is partly determined by their predictions.

## EXERCISES[^fn-10-10]

**1.** Let the demand for money be governed by

$$
m_t = p_t + ky + u_t
$$

where $Eu_t = 0 = Eu_t m_t$. Here $m_t$ is the log of the money supply, $p_t$ the log of the price level, and $y$ the constant level of the log of real income. Assume that $Em_t^2$ and $Ep_t^2$ exist.

Suppose that a researcher attempts to verify the absence of money illusion in this economy by estimating

$$
m_t = \alpha p_t + \text{constant} + \text{residual}_t,
$$

by least squares, and testing whether $\alpha = 1$. In arbitrarily large samples, will this procedure lead him to conclude the truth, namely that $\alpha = 1$? If this procedure is flawed provide a better one and defend it.

**2.** Suppose that the expectations theory of the term structure is correct and that

$$
R_{2t} = \tfrac{1}{2}[R_{1t} + E_t R_{1t+1}]
$$

where $R_{nt}$ is the yield to maturity on an $n$-period bond and $E_t\{x\}$ is the mathematical expectation of $x$ conditioned on information available at time $t$, assumed to include observations on past and present $R_1$ and $R_2$. Deduce the implications that the theory delivers for the population values of $\alpha$, $\beta$, and $\lambda$ in the following least squares regression

$$
R_{2t} - \tfrac{1}{2}[R_{1t} + R_{1t+1}] = \alpha + \beta R_{2t-1} + \lambda R_{1t-1} + u_t
$$

where $u_t$ is a least squares residual obeying $Eu_t = Eu_t R_{2t-1} = Eu_t R_{1t-1} = 0$.

Use the following information to solve problems 3–8. Let $Y$ and $X \equiv (x_0, x_1, \ldots, x_n)$ be random variables with known means and variances, with $x_0 \equiv 1$.

**3.** Prove that if $Ex_1 x_2 = 0$, then $P[Y|x_1, x_2] = P[Y|x_1] + P[Y|x_2]$.

**4.** Prove that if $E(x_1, x_2, \ldots, x_n) = (0, 0, \ldots, 0)$, then

$$
P[Y|1, x_1, x_2, \ldots, x_n] = P[Y|1] + P[Y|x_1, x_2, \ldots, x_n].
$$

**5.** Prove that if $Ex_1 x_2 = 0$, then $P[x_1|x_2] = 0$.

**6.** Let $c$, $d$ be real numbers. Prove that $P[cY_1 + dY_2|X] = cP[Y_1|X] + dP[Y_2|X]$ where $Y_1$ and $Y_2$ are random variables.

**7.** Use 3, 5, and 6 to prove the Kalman filter (recursive projection) formula, {eq}`eq-10-15-prime`:

$$
P[Y|\Omega, X] = P[Y|\Omega] + P\{(Y - P[Y|\Omega])|(X - P[X|\Omega])\}
$$

where $\Omega$ and $X$ are random variables.

**8.** Interpret $Y$, $\Omega$, and $X$ as sets, and $P[Y|\Omega, X]$ as $Y (\Omega\ X)$. Derive the Kalman filter formula using Venn diagrams. (Interpret $Y - X$ as $Y \cap X^c$ where $X^c$ is the complement of $X$.)

**9.** The labor supply schedule is given by

$$
N_t = \gamma(w_t - \hat{E}p_t), \qquad \gamma > 0
$$

where $N_t$ is the logarithm of employment at time $t$, $w_t$ is the logarithm of the money wage, and $\hat{E}p_t$ is the workers' perception of the average price of goods they will buy during period $t$. At the time that they see $w_t$ and must make the decision to work or not, workers don't actually see the prices at which they will be able to buy goods. Instead, they must form their best guess about the average price at which they'll be able to buy things. Workers do know that the average price obeys

$$
p_t = k + u_t
$$

where $k$ is a constant and $u_t$ a serially uncorrelated random variable with mean zero and variance known by workers to be $\sigma_u^2$. The workers also know that $w_t$ obeys

$$
w_t = z_t + u_t
$$

where $z_t$ is a serially uncorrelated random variable with mean zero and variance $\sigma_z^2$. Assume also $Ez_t u_t = 0$. That is, $u$ and $z$ are uncorrelated. The variate $z_t$ measures changes in actual real wages, while the variate $u_t$ measures changes in $w$ and $p$ that actually leave the real wage unchanged. Workers know the value of $w_t$ in each period and the parameters $\sigma_u^2$, $\sigma_z^2$, and $k$.

**A.** Derive an operational labor supply schedule of the form

$$
N_t = \phi w_t + h
$$

where $h$ is a constant, and derive explicit formulas linking $h$ and $\phi$ to the parameters $\gamma$, $\sigma_u^2$ and $\sigma_z^2$.

**B.** A researcher estimates the labor supply schedule, and finds $\phi > 0$ and a very high $R^2$. The researcher concludes that if the monetary authority wants to make $N_t$ constant over time and equal to some high value $N^*$, it should follow a monetary policy that makes $w$ constant over time at the value $w$ that solves the equation

$$
N^* = \phi w + h.
$$

Is this a correct policy conclusion?

**10.** The consumption function is given by

$$
C_t = c(Y_t - T_t^*)
$$

where

$$
\begin{aligned}
C_t &= \text{consumption at time } t \\
Y_t &= \text{income at time } t \\
T_t^* &= \text{consumers' perception of the taxes they will pay at time } t,
\end{aligned}
$$

and $c$ is the marginal propensity to consume out of (perceived) disposable income. At the time they make consumption decisions, agents do not know the taxes ($T_t$) they must pay at time $t$. Consumers do know that income is described by

$$
Y_t = 3T_t + v_t
$$

where $v_t$ is a serially uncorrelated random variable with mean zero and variance $\sigma_v^2$. In addition, they know $c$, that $ET_t v_t = 0$ and the variance of $T_t$: $ET_t^2 = \sigma_T^2$.

**A.** Derive an operational consumption function of the form $C_t = \beta Y_t$ and give an explicit formula relating $\beta$ to the parameters $c$, $\sigma_v^2$, and $\sigma_T^2$.

**B.** Suppose the Congress increases the volatility of taxes by an extraordinary amount. What happens to the consumption function?

## References

- Anderson, B. D. O. and Moore, J. B. (1979). *Optimal Filtering*, Englewood Cliffs: Prentice-Hall.
- Cagan, P. (1956). "The Monetary Dynamics of Hyperinflation." *Studies in the Quantity Theory of Money* (M. Friedman, ed.), Chicago, Illinois: University of Chicago Press.
- Friedman, M. (1956). *A Theory of the Consumption Function*, Princeton, New Jersey: Princeton University Press.
- Hicks, J. R. (1953). *Value and Capital*, London and New York: Oxford University Press.
- Johnston, J. (1963). *Econometric Methods*, New York: McGraw-Hill.
- Luenberger, D. G. (1969). *Optimization by Vector Space Methods*, New York: Wiley.
- Meiselman, D. (1963). *The Term Structure of Interest Rates*, Englewood Cliffs, New Jersey: Prentice-Hall.
- Muth, J. F. (1960). "Optimal Properties of Exponentially Weighted Forecasts." *Journal of the American Statistical Association*, Vol. 55, No. 290, pp. 299–306.
- Muth, J. F. (1961). "Rational Expectations and the Theory of Price Movements." *Econometrica*, Vol. 29, No. 3, pp. 315–335.
- Naylor, A. W. and Sell, G. R. (1982). *Linear Operator Theory in Engineering and Science*, 2nd ed., New York: Springer-Verlag.
- Samuelson, P. A. (1965). "Proof That Properly Anticipated Prices Fluctuate Randomly." *Industrial Management Review*, Vol. 6, No. 2, pp. 41–50.
- Shiller, R. (1978). "Rational Expectations And the Dynamic Structure of Macroeconomic Models: A Critical Review." *Journal of Monetary Economics*, Vol. 4, No. 1.

[^fn-10-6]: The present description of Meiselman's model is along the lines developed after Meiselman by Mincer, Pye, Diller, Shiller, and Nelson. For a useful survey of this literature, see Shiller (1978).

[^fn-10-7]: This formula is an arithmetic approximation to Hicks's formula. Hicks's formula assumes discount bonds that have zero coupons. See Hicks (1939).

[^fn-10-8]: Further, notice that we have no restrictions on the signs of $\phi_j$; they may be either positive or negative, depending on the various covariances among "surprises" that go into composing $\phi_j$.

[^fn-10-9]: Samuelson proved that a property analogous to {eq}`eq-10-37` holds where linear projections are replaced by conditional mathematical expectations. A sequence that obeys that condition is called a martingale with respect to the conditioning information. For conditional expectations, there is a law of iterated expectations precisely paralleling the law of iterated projections. Given this law, Samuelson's theorem is proved in a similar fashion to the weak martingale theorem exhibited in the text.

[^fn-10-10]: Exercises 6–10 were written by Charles Whiteman.
