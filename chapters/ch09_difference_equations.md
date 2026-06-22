# Chapter IX — Difference Equations and Lag Operators[^fn-9-1]

Linear difference equations underlie much work in macroeconomics and applied economic dynamics. Lag operators provide a powerful tool for solving systems of linear difference equations, and for gaining insights about their structure. Lag operators constitute a language which much of macroeconomics and dynamic econometrics takes for granted. This chapter provides an introduction to the analysis of nonstochastic difference equations via lag operators. While a variety of economic examples occur in this chapter, additional examples continue to occur throughout the remainder of this book.

(sec-9-1)=
## 1. Lag Operators

The backward shift or lag operator is defined by

$$
LX_t = X_{t-1}
$$

$$
L^n X_t = X_{t-n} \qquad \text{for} \quad n = \cdots, -2, -1, 0, 1, 2, \ldots.
$$ (eq-9-1)

Multiplying a variable $X_t$ by $L^n$ thus gives the value of $X$ shifted back $n$ periods. Notice that if $n < 0$ in {eq}`eq-9-1`, the effect of multiplying $X_t$ by $L^n$ (more precisely, the effect of "operating on $X_t$ with $L^n$") is to shift $X$ *forward* in time by $-n$ periods.[^fn-9-2] This language is loose. Actually, we are starting out with a *sequence* $\{X_t\}_{t=-\infty}^{\infty}$ which associates a real number $X_t$ with each integer $t$. We are operating on the sequence $\{X_t\}$ with the operator $L^n$ to obtain the new sequence $\{y_t\}_{t=-\infty}^{\infty} = \{X_{t-n}\}_{t=-\infty}^{\infty}$. Formally, the operator $L^n$ maps one sequence into another sequence.

We shall consider polynomials in the lag operator

$$
A(L) = a_0 + a_1 L + a_2 L^2 + \cdots = \sum_{j=0}^{\infty} a_j L^j
$$

where the $a_j$'s are constants and $L^0 \equiv 1$. Operating on $X_t$ with $A(L)$ yields a moving sum of $X$'s:

$$
A(L)X_t = (a_0 + a_1 L + a_2 L^2 + \cdots)X_t
$$

$$
= a_0 X_t + a_1 X_{t-1} + a_2 X_{t-2} + \cdots = \sum_{j=0}^{\infty} a_j X_{t-j}.
$$

It is generally convenient to work with polynomials $A(L)$ that are "rational," meaning that they can be expressed as the ratio of two (finite order) polynomials in $L$:

$$
A(L) = B(L)/C(L)
$$

where

$$
B(L) = \sum_{j=0}^{m} b_j L^j, \qquad C(L) = \sum_{j=0}^{n} c_j L^j
$$

where the $b_j$ and $c_j$ are constants. Assuming that $A(L)$ is rational amounts to imposing a more economical and restrictive parametrization on the $a_j$.

To take the simplest example of a rational polynomial in $L$, consider[^fn-9-3]

$$
A(L) = \frac{1}{1 - \lambda L}.
$$ (eq-9-2)

For the scalar $|C| < 1$, we know that

$$
\frac{1}{1 - C} = 1 + C + C^2 + \cdots.
$$ (eq-9-3)

Thus suggests treating $\lambda L$ of {eq}`eq-9-2` exactly like the $C$ of {eq}`eq-9-3` to get

$$
\frac{1}{1 - \lambda L} = 1 + \lambda L + \lambda^2 L^2 + \cdots,
$$ (eq-9-4)

an expansion which is sometimes only "useful" so long as $|\lambda| < 1$. To motivate the equality {eq}`eq-9-4`, assume that $|\lambda| < 1$ and operate on both sides of {eq}`eq-9-4` with $1 - \lambda L$ to obtain

$$
\frac{1 - \lambda L}{1 - \lambda L} = 1 = (1 + \lambda L + \lambda^2 L^2 + \cdots) - \lambda L(1 + \lambda L + \lambda^2 L^2 + \cdots) = 1.
$$

The reason that we say that {eq}`eq-9-4` is sometimes "useful" only if $|\lambda| < 1$ derives from the following argument. We intend often to multiply $1/(1 - \lambda L)$ by $X_t$ to obtain the infinite moving sum

$$
\frac{1}{1 - \lambda L} X_t = (1 + \lambda L + \lambda^2 L^2 + \cdots)X_t = \sum_{i=0}^{\infty} \lambda^i X_{t-i}.
$$ (eq-9-5)

Consider this sum for a path of $X$ that is constant over time, so that $X_{t-i} = \bar{X}$ for all $i$ and all $t$. Then the sum {eq}`eq-9-5` becomes

$$
\frac{1}{1 - \lambda L} X_t = \bar{X} \sum_{i=0}^{\infty} \lambda^i.
$$

The sum $\sum_{i=0}^{\infty} \lambda^i$ equals $1/(1 - \lambda)$ if $|\lambda| < 1$. But if $|\lambda| \geq 1$ that sum is unbounded, being $+\infty$ if $\lambda \geq 1$. We shall sometimes (though not always) be applying the polynomial in the lag operator {eq}`eq-9-4` in situations in which it is appropriate to go infinitely far back in time; and we sometimes find it necessary to insist that in such cases the infinite sum in {eq}`eq-9-5` exist where $X$ has been constant through time. This is what leads to the requirement sometimes imposed that $|\lambda| < 1$ in {eq}`eq-9-4`. As we shall see, however, in standard analyses of difference equations, which take the starting point of all processes as some point only finitely far back into the past, the requirement that $|\lambda| < 1$ need not be imposed in {eq}`eq-9-4`.

It is useful to note that there is an alternative expansion for the "geometric" polynomial $1/(1 - \lambda L)$. For notice that formally

$$
\frac{1}{1 - \lambda L} = \frac{-(\lambda L)^{-1}}{1 - (\lambda L)^{-1}} = \frac{-1}{\lambda L}\left(1 + \frac{1}{\lambda}L^{-1} + \left(\frac{1}{\lambda}\right)^2 L^{-2} + \cdots\right)
$$

$$
= -\frac{1}{\lambda}L^{-1} - \left(\frac{1}{\lambda}\right)^2 L^{-2} - \left(\frac{1}{\lambda}\right)^3 L^{-3} - \cdots,
$$ (eq-9-6)

an expansion which is especially "useful" where $|\lambda| > 1$, i.e., where $|1/\lambda| < 1$. So {eq}`eq-9-6` implies that

$$
\frac{1}{1 - \lambda L} X_t = -\frac{1}{\lambda} X_{t+1} - \left(\frac{1}{\lambda}\right)^2 X_{t+2} - \cdots = -\sum_{i=1}^{\infty} \left(\frac{1}{\lambda}\right)^i X_{t+i},
$$

which shows $(1/(1 - \lambda L))X_t$ to be a geometrically declining weighted sum of *future* values of $X$. Notice that for this infinite sum to be finite for a constant time path $X_{t+i} = \bar{X}$ for all $i$ and $t$, the series $-\sum_{i=1}^{\infty}(1/\lambda)^i$ must be convergent, which requires that $|1/\lambda| < 1$.

More generally, let $\{x_t\}_{t=-\infty}^{\infty}$ be any bounded sequence of real numbers, i.e., for some $M > 0$, $|x_t| < M$ for all $t$. Then applying the operator $1/(1 - \lambda L)$ to the sequence $\{x_t\}_{t=-\infty}^{\infty}$ can be taken to give either the sequence $\{y_t\}_{t=-\infty}^{\infty}$ where

$$
y_t = \frac{1}{1 - \lambda L} x_t = \sum_{j=0}^{\infty} \lambda^j x_{t-j}
$$

or the sequence $\{z_t\}_{t=-\infty}^{\infty}$ where

$$
z_t = \frac{-(\lambda L)^{-1}}{1 - (\lambda L)^{-1}} x_t = -\sum_{j=1}^{\infty} \left(\frac{1}{\lambda}\right)^j x_{t+j}.
$$

In general, $\{y_t\}$ is a bounded sequence if $|\lambda| < 1$, while $\{z_t\}$ is a bounded sequence if $|\lambda| > 1$. In many (though not all) contexts, we want application of $(1 - \lambda L)^{-1}$ to map all bounded sequences into bounded sequences, so that we choose the "backward" expansion if $|\lambda| < 1$ and the forward expansion if $|\lambda| > 1$.

To illustrate how polynomials in the lag operator can be manipulated, consider the difference equation

$$
Y_t = \lambda Y_{t-1} + bX_t + a, \qquad t = -\infty, \ldots, 0, 1, 2, \ldots,
$$ (eq-9-7)

where $X_t$ is an exogenous variable and $Y_t$ is an endogenous variable and $\lambda \neq 1$. Here, $X_t$ is a sequence of real numbers $t = \cdots, -1, 0, 1, 2, \ldots$. Write the above equation as

$$
(1 - \lambda L)Y_t = a + bX_t.
$$

Operating on both sides of this equation by $(1 - \lambda L)^{-1}$ gives

$$
Y_t = \frac{a}{1 - \lambda L} + \frac{b}{1 - \lambda L} X_t + c\lambda^t
$$

$$
= \frac{a}{1 - \lambda} + b\sum_{i=0}^{\infty} \lambda^i X_{t-i} + c\lambda^t,
$$ (eq-9-8)

since $a/(1 - \lambda L) = a\sum_{i=0}^{\infty}\lambda^i = a/(1 - \lambda)$. Here $c$ is any constant. The reason that we must include the term $c\lambda^t$ in {eq}`eq-9-8` is that for any constant $c$, $(1 - \lambda L)c\lambda^t = c\lambda^t - c\lambda\lambda^{t-1} = 0$.[^fn-9-4] Therefore, it follows that application of $1 - \lambda L$ to both sides of {eq}`eq-9-8` gives Equation {eq}`eq-9-7` once again. Consequently, {eq}`eq-9-8` is the general "solution" of the difference equation {eq}`eq-9-7` and describes the entire time path of $Y$ associated with a given time path of $X$. In order to get a "particular solution" we must be able to tie down the constant $c$. This requires an additional bit of information in the form of a specified value of $Y_t$ at some particular time or some conditions on the path of $\{Y_t\}$ such as boundedness. Notice that for the $Y_t$ defined by {eq}`eq-9-8` to be finite, $\lambda^i X_{t-i}$ must be "small" for large $i$. More precisely, we require

$$
\lim_{n \to \infty} \sum_{i=n}^{\infty} \lambda^i X_{t-i} = 0 \qquad \text{for all } t.
$$ (eq-9-9)

For the case of $X$ constant for all time, $X_{t-i} = \bar{X}$ all $i$ and $t$, this condition requires $|\lambda| < 1$. Notice also that the infinite sum $a\sum_{i=0}^{\infty}\lambda^i$ in {eq}`eq-9-8` is finite only if $|\lambda| < 1$, in which case it equals $a/(1 - \lambda)$, or if $a = 0$, in which case it equals zero regardless of the value of $\lambda$. We tentatively assume that $|\lambda| < 1$.

For analyzing difference equations with arbitrary initial conditions given, it is convenient to rewrite {eq}`eq-9-8` for $t > 0$ as

$$
Y_t = a\sum_{i=0}^{t-1} \lambda^i + a\sum_{i=t}^{\infty} \lambda^i + b\sum_{i=0}^{t-1} \lambda^i X_{t-i} + b\sum_{i=t}^{\infty} \lambda^i X_{t-i} + c\lambda^t
$$

$$
= \frac{a(1 - \lambda^t)}{1 - \lambda} + \frac{a\lambda^t}{1 - \lambda} + b\sum_{i=0}^{t-1} \lambda^i X_{t-i} + b\lambda^t\sum_{i=0}^{\infty} \lambda^i X_{0-i} + c\lambda^t,
$$

$$
Y_t = \frac{a(1 - \lambda^t)}{1 - \lambda} + b\sum_{i=0}^{t-1} \lambda^i X_{t-i} + \lambda^t\left\{\frac{a}{1 - \lambda} + b\sum_{i=0}^{\infty} \lambda^i X_{0-i} + c\lambda^0\right\}, \qquad t \geq 1.
$$ (eq-9-10)

The term in braces equals $Y_0$, as reference to expression {eq}`eq-9-8` will confirm. So {eq}`eq-9-10` becomes

$$
Y_t = \frac{a(1 - \lambda^t)}{1 - \lambda} + b\sum_{i=0}^{t-1} \lambda^i X_{t-i} + \lambda^t Y_0
$$

or

$$
Y_t = \frac{a}{1 - \lambda} + \lambda^t\left(Y_0 - \frac{a}{1 - \lambda}\right) + b\sum_{i=0}^{t-1} \lambda^i X_{t-i}, \qquad t \geq 1.
$$ (eq-9-11)

Now textbooks on difference equations often analyze the special case in which $X_t = 0$ for all $t \geq 0$. Under this special circumstance {eq}`eq-9-11` becomes

$$
Y_t = \frac{a}{1 - \lambda} + \lambda^t\left(Y_0 - \frac{a}{1 - \lambda}\right),
$$ (eq-9-12)

which is the solution of the first-order difference equation $Y_t = a + \lambda Y_{t-1}$ subject to the initial condition that $Y$ equals the arbitrary given value $Y_0$ at time 0. Notice that if $Y_0 = a/(1 - \lambda)$, then {eq}`eq-9-12` implies $Y_t = Y_0$ for all $t \geq 0$, which shows $a/(1 - \lambda)$ to be a "stationary point" or long-run equilibrium value of $Y$. Notice also that if, as we are assuming, $|\lambda| < 1$, then {eq}`eq-9-12` implies that

$$
\lim_{t \to \infty} Y_t = \frac{a}{1 - \lambda},
$$

which shows that the system is "stable," tending to approach the stationary point as time passes.

Now consider the first-order system {eq}`eq-9-7` under the assumption that $a = 0$, so that $a\sum_{i=0}^{\infty}\lambda^i$ equals zero regardless of the value of $\lambda$. Then the appropriate counterpart to {eq}`eq-9-10` is

$$
Y_t = b\sum_{i=0}^{t-1} \lambda^i X_{t-i} + \lambda^t\left\{b\sum_{i=0}^{\infty} \lambda^i X_{0-i} + c\lambda^0\right\}.
$$

Assuming that condition {eq}`eq-9-9` is met even where $|\lambda| > 1$ (so that the second term in the equation is finite), the above equation becomes

$$
Y_t = b\sum_{i=0}^{t-1} \lambda^i X_{t-i} + \lambda^t Y_0, \qquad t \geq 1.
$$

As before we analyze the special case where $X_t = 0$ for all $t > 0$. Then the above equation becomes

$$
Y_t = \lambda^t Y_0, \qquad t \geq 1.
$$

The stationary point of this solution is zero since if $Y_0 = 0$, $Y$ will remain equal to zero forever, regardless of the value of $\lambda$. However, if $|\lambda| > 1$, the system will diverge farther and farther from this stationary point if either $Y_0 > 0$, or $Y_0 < 0$. If $\lambda > 1$, $Y_t$ will tend toward $+\infty$ as $t \to \infty$ provided $Y_0 > 0$; $Y_t$ will tend toward $-\infty$ as $t \to \infty$ if $Y_0 < 0$. If $\lambda < -1$, $Y_t$ will display explosive oscillations of periodicity two time periods.

We can also solve the difference equation {eq}`eq-9-7` by applying the "forward inverse" of $1 - \lambda L$ to get the general solution

$$
Y_t = \frac{-\lambda^{-1}L^{-1}}{1 - \lambda^{-1}L^{-1}} a + b\left(\frac{-\lambda^{-1}L^{-1}}{1 - \lambda^{-1}L^{-1}}\right)X_t + d\lambda^t,
$$

$$
Y_t = \frac{a}{1 - \lambda} - b\sum_{i=0}^{\infty} \left(\frac{1}{\lambda}\right)^{i+1} X_{t+i+1} + d\lambda^t,
$$ (eq-9-8-prime)

where $d$ is a constant to be determined from some side condition on the path of $Y_t$, such as an initial condition or terminal condition. If $a = 0$, then for *any* value of $\lambda \neq 1$, in general {eq}`eq-9-8` and {eq}`eq-9-8-prime` *both* represent solutions of the difference equation {eq}`eq-9-7`. They are simply alternative representations of the solution in the sense that for any given initial condition or other side condition, one can generally find values of $d$ and $c$ which guarantee that both {eq}`eq-9-8` and {eq}`eq-9-8-prime` satisfy {eq}`eq-9-7`. The equivalence of the solutions {eq}`eq-9-8` and {eq}`eq-9-8-prime` will hold whenever $(b/(1 - \lambda L))X_t$ and $(b\lambda^{-1}L^{-1}/(1 - \lambda^{-1}L^{-1}))X_t$ are both finite for all $t$.

It often happens, however, that either $(b/(1 - \lambda L))X_t$ or

$$
\frac{b\lambda^{-1}L^{-1}}{1 - \lambda^{-1}L^{-1}} X_t
$$

fails to be finite, i.e., the infinite sum fails to converge. In this case, one or the other of the representations {eq}`eq-9-8` or {eq}`eq-9-8-prime` breaks down, i.e., fails to give a $Y_t$ sequence that is finite for all finite $t$. For example, if the sequence $\{X_t\}$ is bounded, this is sufficient to imply that $\{(b/(1 - \lambda L))X_t\}$ is a bounded sequence if $|\lambda| < 1$, but not sufficient to imply that

$$
\frac{b\lambda^{-1}L^{-1}}{1 - \lambda^{-1}L^{-1}} X_t
$$

is a convergent sum for all $t$. Similarly, if $|\lambda| > 1$, boundedness of the sequence $\{X_t\}$ is sufficient to imply that

$$
\left\{\frac{b\lambda^{-1}L^{-1}}{1 - \lambda^{-1}L^{-1}} X_t\right\}
$$

is a bounded sequence, but fails to guarantee finiteness of $(b/(1 - \lambda L))X_t$. In instances where one of $(b/(1 - \lambda L))X_t$ or

$$
\frac{b\lambda^{-1}L^{-1}}{1 - \lambda^{-1}L^{-1}} X_t
$$

is always finite and the other is not we shall take as our solution to {eq}`eq-9-7` either {eq}`eq-9-8` where the backward sum in $X_t$ is finite, or {eq}`eq-9-8-prime` where the forward sum in $X_t$ is finite. This procedure assures us that we shall find the unique solution of {eq}`eq-9-7` that is finite for all finite $t$, provided that such a solution exists. Such a solution is guaranteed to exist where $\{X_t\}$ is a bounded sequence.

Now if we desired to impose that the $\{Y_t\}$ sequence given by {eq}`eq-9-8` or {eq}`eq-9-8-prime` is bounded—as we are free to do if no other side condition has been imposed—then it is evident that we must set $c = 0$ in {eq}`eq-9-8` or $d = 0$ in {eq}`eq-9-8-prime`. This is necessary since if $\lambda > 1$ and $c > 0$,

$$
\lim_{t \to \infty} c\lambda^t = \infty;
$$

while if $\lambda < 1$ and $c > 0$,

$$
\lim_{t \to -\infty} c\lambda^t = \infty.
$$

It follows that $Y_t$ will be bounded for all $t$ only if $c$ or $d$ is zero. For the solution in the forward direction to be finite for all finite $t$, we clearly require a condition analogous to {eq}`eq-9-9`

$$
\lim_{n \to \infty} \sum_{i=n}^{\infty} \left(\frac{1}{\lambda}\right)^i X_{t+i} = 0.
$$ (eq-9-9-prime)

The principle of solving "stable roots" ($|\lambda| < 1$) backward and "unstable roots" ($|\lambda| > 1$) forward was encountered in Chapter I. It is a device designed to ensure that the solution of the differential (or difference) equation maps bounded functions (or sequences) as driving processes into bounded functions (or sequences). Below, we shall see that a formal justification for this procedure is sometimes available in the context of difference equations that emerge from optimum problems.

(sec-9-2)=
## 2. Second-Order Difference Equations

Consider the second-order difference equation

$$
Y_t = t_1 Y_{t-1} + t_2 Y_{t-2} + a + bX_t,
$$ (eq-9-13)

where $\{X_t\}$ is again an exogenous sequence of real numbers for $t = \cdots, -1, 0, 1, \ldots$. Using lag operators, {eq}`eq-9-13` can be written as

$$
(1 - t_1 L - t_2 L^2)Y_t = a + bX_t.
$$

A solution to this difference equation is given by

$$
Y_t = \frac{a}{1 - t_1 L - t_2 L^2} + \frac{b}{1 - t_1 L - t_2 L^2} X_t
$$ (eq-9-14)

where we have temporarily ignored the terms analogous to $c\lambda^t$ which appeared in the general solution of the first-order equation. By long division it is easy to verify that

$$
\frac{b}{1 - t_1 L - t_2 L^2} = \sum_{i=0}^{\infty} w_i L^i
$$ (eq-9-15)

where $w_0 = b$, $w_1 = bt_1$, and

$$
w_j = t_1 w_{j-1} + t_2 w_{j-2} \qquad \text{for} \quad j \geq 2.
$$

That is,

$$
\begin{array}{l}
\phantom{1 - t_1 L - t_2 L^2 \big)\ } 1 + t_1 L + (t_2 + t_1^2)L^2 + (t_1(t_2 + t_1^2) + t_1 t_2)L^3 + \cdots \\[4pt]
1 - t_1 L - t_2 L^2 \,\big)\, 1 \\
\phantom{1 - t_1 L - t_2 L^2 \big)\ } \underline{1 - t_1 L - t_2 L^2} \\
\phantom{1 - t_1 L - t_2 L^2 \big)\ 1} t_1 L + t_2 L^2 \\
\phantom{1 - t_1 L - t_2 L^2 \big)\ 1} \underline{t_1 L - t_1^2 L^2 \quad\ \ - t_1 t_2 L^3} \\
\phantom{1 - t_1 L - t_2 L^2 \big)\ 1\ } (t_2 + t_1^2)L^2 + t_1 t_2 L^3 \\
\phantom{1 - t_1 L - t_2 L^2 \big)\ 1\ } \underline{(t_2 + t_1^2)L^2 - t_1(t_2 + t_1^2)L^3 - t_2(t_2 + t_1^2)L^4} \cdots.
\end{array}
$$

Notice that the weights in {eq}`eq-9-15` follow a geometric pattern if $t_2 = 0$, as we would expect, since then {eq}`eq-9-13` collapses to a first-order equation.

It is convenient to write the polynomial $1 - t_1 L - t_2 L^2$ in an alternative way, given by the *factorization*

$$
1 - t_1 L - t_2 L^2 = (1 - \lambda_1 L)(1 - \lambda_2 L)
$$

$$
= 1 - (\lambda_1 + \lambda_2)L + \lambda_1 \lambda_2 L^2,
$$ (eq-9-16)

[^fn-9-1]: It would be useful for the reader to be familiar with the material on difference equations in Allen (1960) and Baumol (1959).

[^fn-9-2]: This chapter aims to teach the reader to manipulate lag operators, while devoting little or no attention to describing their mathematical foundations. The key Riesz–Fischer theorem which justifies these methods is discussed briefly in Chapter XI, pp. 249–253. The reader interested in increasing his proficiency with these techniques is urged to consult Gabel and Roberts (1973, Chapter 4).

[^fn-9-3]: Actually, we should write $A(L) = I/(1 - \lambda L)$ where $I$ is the identity lag operator defined by $I \equiv 1 + 0L + 0L^2 + \cdots$. So $I$ satisfies $Ix_t = x_t$, and thus acts like unity.

[^fn-9-4]: Technically, we are free to add to the solution any function of time $f(t)$ for which $(1 - \lambda L)f(t) = 0$. It can be proved that $f(t) = c\lambda^t$ is the only such function.


so that $\lambda_1 + \lambda_2 = t_1$ and $-\lambda_1\lambda_2 = t_2$. To see how $\lambda_1$ and $\lambda_2$ are related to the *roots* or *zeros* of $1 - t_1 z - t_2 z^2$, notice that

$$
(1 - \lambda_1 z)(1 - \lambda_2 z) = \lambda_1\lambda_2\left(\frac{1}{\lambda_1} - z\right)\left(\frac{1}{\lambda_2} - z\right).
$$

Therefore the equation

$$
0 = (1 - \lambda_1 z)(1 - \lambda_2 z) = \lambda_1\lambda_2\left(\frac{1}{\lambda_1} - z\right)\left(\frac{1}{\lambda_2} - z\right)
$$

is satisfied at the two roots $z = 1/\lambda_1$ and $z = 1/\lambda_2$. Given the polynomial $1 - t_1 z - t_2 z^2$, the roots $1/\lambda_1$ and $1/\lambda_2$ are found from solving the *characteristic equation*

$$
1 - t_1 z - t_2 z^2 = 0 \qquad \text{or} \qquad t_2 z^2 + t_1 z - 1 = 0
$$

for two values of $z$. The roots are given by the quadratic formula

$$
z = \frac{-t_1 \pm \sqrt{t_1^2 + 4t_2}}{2t_2}.
$$ (eq-9-17)

Formula {eq}`eq-9-17` enables us to obtain the reciprocals of $\lambda_1$ and $\lambda_2$ for given values of $t_1$ and $t_2$.

We assume that $\lambda_1 \neq \lambda_2$, $\lambda_1 \neq 1$. Then without loss of generality we can write the second-order difference equation as

$$
(1 - \lambda_1 L)(1 - \lambda_2 L)Y_t = a + bX_t.
$$

The general solution to this difference equation is

$$
Y_t = \frac{a}{(1 - \lambda_1 L)(1 - \lambda_2 L)} + \frac{b}{(1 - \lambda_1 L)(1 - \lambda_2 L)}X_t + c_1\lambda_1^t + c_2\lambda_2^t
$$ (eq-9-18)

where $c_1$ and $c_2$ are any constants. To see that {eq}`eq-9-18` solves the difference equation for any values of $c_1$ and $c_2$, operate on both sides of {eq}`eq-9-18` with

$$
(1 - \lambda_1 L)(1 - \lambda_2 L)
$$

and notice that

$$
(1 - \lambda_1 L)(1 - \lambda_2 L)c_1\lambda_1^t = (1 - \lambda_1 L)(1 - \lambda_2 L)c_2\lambda_2^t = 0.
$$

In order to determine a particular solution to the difference equation, we now need two side conditions on the path of $Y_t$. For example, two initial conditions in the forms of given values for $Y$ at time 0 and 1 are sufficient to determine $c_1$ and $c_2$.

Notice that since $\lambda_1 \neq \lambda_2$,

$$
\frac{1}{(1 - \lambda_1 L)(1 - \lambda_2 L)} = \frac{1}{\lambda_1 - \lambda_2}\left(\frac{\lambda_1}{1 - \lambda_1 L} - \frac{\lambda_2}{1 - \lambda_2 L}\right),
$$

which can be verified directly. Thus {eq}`eq-9-18` can be written

$$
\begin{aligned}
Y_t &= \frac{a}{(1 - \lambda_1)(1 - \lambda_2)} + \frac{\lambda_1 b}{\lambda_1 - \lambda_2}\frac{1}{1 - \lambda_1 L}X_t \\
&\quad - \frac{\lambda_2 b}{\lambda_1 - \lambda_2}\frac{1}{1 - \lambda_2 L}X_t + c_1\lambda_1^t + c_2\lambda_2^t \\
&= a\sum_{i=0}^{\infty}\lambda_1^i\sum_{j=0}^{\infty}\lambda_2^j + \frac{\lambda_1 b}{\lambda_1 - \lambda_2}\sum_{i=0}^{\infty}\lambda_1^i X_{t-i} \\
&\quad - \frac{\lambda_2 b}{\lambda_1 - \lambda_2}\sum_{i=0}^{\infty}\lambda_2^i X_{t-i} + c_1\lambda_1^t + c_2\lambda_2^t
\end{aligned}
$$ (eq-9-19)

where we are making use of the fact that for a constant $a$

$$
H(L)a = \sum_{i=0}^{\infty}h_i L^i a = a\sum_{i=0}^{\infty}h_i = aH(1).
$$

Notice that

$$
\frac{1}{1 - \lambda_1 L}\frac{1}{1 - \lambda_2 L} = \sum_{i=0}^{\infty}\lambda_1^i L^i \sum_{j=0}^{\infty}\lambda_2^j L^j,
$$

so that the sum of the distributed lag weights $\sum_{i=0}^{\infty}\lambda_1^i\sum_{j=0}^{\infty}\lambda_2^i$ is finite and equals $1/((1 - \lambda_1)(1 - \lambda_2))$ provided that both $|\lambda_1| < 1$, $|\lambda_2| < 1$. So in writing {eq}`eq-9-19` we require either that both $|\lambda_1|$ and $|\lambda_2|$ be less than unity or that $a = 0$, so that $a\sum_{i=0}^{\infty}\lambda_1^i\sum_{j=0}^{\infty}\lambda_2^j$ is defined. Furthermore, we require that

$$
\lim_{n \to \infty}\sum_{i=n}^{\infty}\lambda_j^i X_{t-i} = 0, \qquad \text{all } t,
$$

hold for $j = 1, 2$, so that the geometric sums in {eq}`eq-9-19` are both finite.

Suppose that $a = 0$. On this assumption write {eq}`eq-9-19` for $t \geq 1$ as

$$
\begin{aligned}
Y_t &= \frac{\lambda_1 b}{\lambda_1 - \lambda_2}\sum_{i=0}^{t-1}\lambda_1^i X_{t-i} - \frac{\lambda_2 b}{\lambda_1 - \lambda_2}\sum_{i=0}^{t-1}\lambda_2^i X_{t-i} \\
&\quad + \frac{\lambda_1 b}{\lambda_1 - \lambda_2}\sum_{i=t}^{\infty}\lambda_1^i X_{t-i} - \frac{\lambda_2 b}{\lambda_1 - \lambda_2}\sum_{i=t}^{\infty}\lambda_2^i X_{t-i} + c_1\lambda_1^t + c_2\lambda_2^t
\end{aligned}
$$

or

$$
Y_t = \frac{\lambda_1 b}{\lambda_1 - \lambda_2}\sum_{i=0}^{t-1}\lambda_1^i X_{t-i} - \frac{\lambda_2 b}{\lambda_1 - \lambda_2}\sum_{i=0}^{t-1}\lambda_2^i X_{t-i} + \lambda_1^t\theta_0 + \lambda_2^t\eta_0, \quad t \geq 1,
$$ (eq-9-20)

where

$$
\theta_0 = \left\{c_1 + \frac{b\lambda_1}{\lambda_1 - \lambda_2}\sum_{i=0}^{\infty}\lambda_1^i X_{0-i}\right\} \quad \text{and} \quad \eta_0 = \left\{c_2 - \frac{b\lambda_2}{\lambda_1 - \lambda_2}\sum_{i=0}^{\infty}\lambda_2^i X_{0-i}\right\}.
$$

The case in which $X_t = 0$ for $t \geq 1$ is often analyzed, as for the first-order case. On this assumption {eq}`eq-9-20` becomes

$$
Y_t = \lambda_1^t\theta_0 + \lambda_2^t\eta_0, \qquad t \geq 1.
$$ (eq-9-21)

If $\theta_0 = \eta_0 = 0$, $Y_t = 0$ for all $t \geq 1$, regardless of the values of $\lambda_1$ and $\lambda_2$. So $Y = 0$ is the stationary point or long-run equilibrium value of {eq}`eq-9-21`.

If $\lambda_1$ and $\lambda_2$ are real, then $\lim_{t \to \infty} Y_t$ will equal zero if and only if both $|\lambda_1| < 1$ and $|\lambda_2| < 1$, regardless of the values of the parameters $\theta_0$ and $\eta_0$, so long as they are finite.

If $\lambda_1 > 1$, $|\lambda_2| < |\lambda_1|$, and $\theta_0 > 0$, then $\lim_{t \to \infty} Y_t = +\infty$. If $\lambda_1 > 1$, $|\lambda_2| < |\lambda_1|$, and $\theta_0 < 0$, then $\lim_{t \to \infty} Y_t = -\infty$. Thus, $Y$ will tend toward the stationary point zero as time passes provided that both $|\lambda_1| < 1$ and $|\lambda_2| < 1$. If one or both of the $\lambda$'s exceed one in absolute value, the behavior of $Y$ will eventually be "dominated" by the term in {eq}`eq-9-21` associated with the $\lambda$ that is larger in absolute value; i.e., eventually $Y$ will grow approximately as $\lambda_m^t$, where $\lambda_m$ is the $\lambda_j$ with the larger absolute value.

Now suppose that the roots are complex. If the roots are complex, they will occur as a complex conjugate pair, as the quadratic formula {eq}`eq-9-17` verifies. So assume that the roots are complex, and write them as

$$
\lambda_1 = re^{iw} = r(\cos w + i\sin w), \qquad \lambda_2 = re^{-iw} = r(\cos w - i\sin w)
$$

where the real part is $r\cos w$ and the imaginary part is $\pm r\sin w$. Notice that

$$
\lambda_1 - \lambda_2 = r(e^{iw} - e^{-iw}) = 2ri\sin w.
$$ (eq-9-22)

Equation {eq}`eq-9-21` becomes

$$
\begin{aligned}
Y_t &= \theta_0(re^{iw})^t + \eta_0(re^{-iw})^t \\
&= \theta_0(r^t e^{iwt}) + \eta_0(r^t e^{-iwt}) \\
&= \theta_0 r^t[\cos wt + i\sin wt] + \eta_0 r^t[\cos wt - i\sin wt] \\
&= (\theta_0 + \eta_0)r^t\cos wt + i(\theta_0 - \eta_0)r^t\sin wt.
\end{aligned}
$$ (eq-9-23)

Since $Y_t$ must be a real number for all $t$, it follows that $\theta_0 + \eta_0$ must be real and $\theta_0 - \eta_0$ must be imaginary. Therefore, $\theta_0$ and $\eta_0$ must be complex conjugates, say $\theta_0 = pe^{i\theta}$, $\eta_0 = pe^{-i\theta}$. Therefore we can write

$$
\begin{aligned}
Y_t &= pe^{i\theta}r^t e^{iwt} + pe^{-i\theta}r^t e^{-iwt} = pr^t[e^{i(wt+\theta)} + e^{-i(wt+\theta)}] \\
&= 2pr^t\cos(wt + \theta).
\end{aligned}
$$ (eq-9-24)

This is the solution of the "unforced" (i.e., $X_t = 0$ for all $t$) second-order difference equation with complex roots. The parameters $p$ and $\theta$ are chosen to satisfy two side conditions on the path of $Y_t$, say two initial conditions. The path of $Y_t$ oscillates with a frequency determined by $w$. The "damping factor" $r^t$ is determined by the amplitude $r$ of the complex roots. The value $Y_t = 0$ is the stationary point of the difference equation, which will be approached at $t \to \infty$ for arbitrary initial conditions if $r < 1$. If $r > 1$, the path of $Y_t$ displays explosive oscillations, unless the initial conditions are, say, $Y_0 = 0$, $Y_1 = 0$, so that $Y$ starts out at the stationary point for two successive values. If $r < 1$, the system displays damped oscillations provided $w \neq 0$, which is so when the roots are complex. If $r = 1$, $Y_t$ displays repeated oscillations of unchanging amplitude, and the solution is "periodic."

If $\lambda_1$ and $\lambda_2$ are complex, the distributed lag weights of {eq}`eq-9-19` are easily shown to oscillate. For we have

$$
\begin{aligned}
\frac{b}{\lambda_1 - \lambda_2}\sum_{j=0}^{\infty}(\lambda_1^{j+1} - \lambda_2^{j+1})X_{t-j} &= \frac{b}{re^{iw} - re^{-iw}}\sum_{j=0}^{\infty}[(re^{iw})^{j+1} - (re^{-iw})^{j+1}]X_{t-j} \\
&= \frac{b}{2ri\sin w}\sum_{j=0}^{\infty}2r^{j+1}i\{\sin w(j+1)\}X_{t-j} = b\sum_{j=0}^{\infty}r^j\frac{\sin w(j+1)}{\sin w}X_{t-j}.
\end{aligned}
$$

Since the roots are complex, $\sin w \neq 0$. Notice that the damping factor weighting the sine curve is $r^j$, so that the amplitude of the lag weights decreases as $j$ increases, provided that $r < 1$.

As noted above, the roots $\lambda_1$ and $\lambda_2$ are the reciprocals of the roots of the polynomial

$$
1 - t_1 z - t_2 z^2 = 0.
$$ (eq-9-25)

For we know that $1 - t_1 z - t_2 z^2 = (1 - \lambda_1 z)(1 - \lambda_2 z)$, with roots $1/\lambda_1$ and $1/\lambda_2$. Alternatively, multiply the above equation by $z^{-2}$ to obtain

$$
z^{-2} - z^{-1}t_1 - t_2 = 0 = (z^{-1} - \lambda_1)(z^{-1} - \lambda_2)
$$

or

$$
s^2 - t_1 s - t_2 = 0
$$ (eq-9-26)

where $s = z^{-1}$. Notice that the roots of {eq}`eq-9-26` are the reciprocals of the roots of {eq}`eq-9-25`. Thus, $\lambda_1$ and $\lambda_2$ are the roots of {eq}`eq-9-26`.

It is interesting to know what values of $t_1$ and $t_2$ yield complex roots. Using the quadratic formula we have that the roots of {eq}`eq-9-26` are

$$
\lambda_i = s = \frac{t_1 \pm \sqrt{t_1^2 + 4t_2}}{2}.
$$

For the roots to be complex, the term whose square root is taken must be negative, i.e.,

$$
t_1^2 + 4t_2 < 0,
$$ (eq-9-27)

which implies that $t_2 < 0$. In case {eq}`eq-9-27` is satisfied, the roots are

$$
\lambda_1 = \frac{t_1}{2} + \frac{i\sqrt{-(t_1^2 + 4t_2)}}{2} = a + bi, \qquad \lambda_2 = \frac{t_1}{2} - \frac{i\sqrt{-(t_1^2 + 4t_2)}}{2} = a - bi.
$$

To write $a + bi$ in polar form we recall that

$$
a + bi = r\cos w + ri\sin w = re^{iw}
$$

where $r = \sqrt{a^2 + b^2}$ and where $\cos w = a/r$. Thus we have that

$$
r = \sqrt{\left(\frac{t_1}{2}\right)^2 - \frac{(t_1^2 + 4t_2)}{4}} = \sqrt{-t_2}.
$$

We also have that

$$
\cos w = \frac{t_1}{2\sqrt{-t_2}} \qquad \text{or} \qquad w = \cos^{-1}\left(\frac{t_1}{2\sqrt{-t_2}}\right).
$$

For the oscillations to be damped we require that $r = \sqrt{-t_2} < 1$, which requires that $-t_2 < 1$.

The periodicity of the oscillations is $2\pi/\cos^{-1}(t_1/2\sqrt{-t_2})$; i.e., this is the number of periods from peak to peak in the oscillations.

If the roots are real, movements will be damped if both roots are less than one in absolute value. That requires

$$
-1 < \frac{t_1 + \sqrt{t_1^2 + 4t_2}}{2} < 1 \qquad \text{and} \qquad -1 < \frac{t_1 - \sqrt{t_1^2 + 4t_2}}{2} < 1.
$$

The condition $\frac{1}{2}(t_1 + \sqrt{t_1^2 + 4t_2}) < 1$ implies

$$
\sqrt{t_1^2 + 4t_2} < 2 - t_1, \qquad t_1^2 + 4t_2 < 4 + t_1^2 - 4t_1,
$$

$$
t_1 + t_2 < 1.
$$ (eq-9-28)

The condition $\frac{1}{2}(t_1 - \sqrt{t_1^2 + 4t_2}) > -1$ implies that

$$
-\sqrt{t_1^2 + 4t_2} > -2 - t_1, \qquad \sqrt{t_1^2 + 4t_2} < 2 + t_1,
$$

$$
t_1^2 + 4t_2 < t_1^2 + 4 + 4t_1,
$$

$$
t_2 < 1 + t_1.
$$ (eq-9-29)

Conditions {eq}`eq-9-28` and {eq}`eq-9-29` must be satisfied for the roots, if real, to be less than unity in absolute value.

Notice that both roots are negative and real if

$$
t_1^2 + 4t_2 > 0 \qquad \text{and} \qquad \frac{t_1 + \sqrt{t_1^2 + 4t_2}}{2} < 0,
$$

which implies

$$
t_1 < -\sqrt{t_1^2 + 4t_2}, \qquad t_1^2 > t_1^2 + 4t_2, \qquad 0 > t_2.
$$

Figure 1 depicts regions of the $t_1, t_2$ plane for which conditions {eq}`eq-9-27`, {eq}`eq-9-28`, or {eq}`eq-9-29` are or are not satisfied. The graph shows combinations of $t_1$ and $t_2$ that give rise to damped oscillations, explosive oscillations, etc.

```{figure} ../figures/fig-9-1.png
:name: fig-9-1
:width: 80%
:align: center

**Figure 1.** (Source: W. Baumol, *Economic Dynamics*, 2nd ed., p. 221, New York: Macmillan. Copyright 1959 by Macmillan Publishing Co., Inc.)
```

### A. An Example

Maybe the most famous second-order difference equation in economics is the one associated with Samuelson's multiplier accelerator model (see Samuelson, 1944). Samuelson posited the model

$$
\begin{aligned}
C_t &= cY_{t-1} + \alpha, &\qquad 1 > c > 0 \quad &\text{(consumption function)} \\
I_t &= \gamma(Y_{t-1} - Y_{t-2}), &\qquad \gamma > 0 \quad &\text{(accelerator)} \\
C_t + I_t &= Y_t,
\end{aligned}
$$

```{figure} ../figures/fig-9-2.png
:name: fig-9-2
:width: 80%
:align: center

**Figure 2.**
```

where $C_t$ is consumption and $I_t$ is investment. Substituting the first two equations into the third gives

$$
Y_t = (c + \gamma)Y_{t-1} - \gamma Y_{t-2} + \alpha
$$

or

$$
Y_t = t_1 Y_{t-1} + t_2 Y_{t-2} + \alpha
$$

where $t_1 = c + \gamma$, $t_2 = -\gamma$. Notice that $t_1 + t_2 = c$. So increases in the parameter $\gamma$ move the parameters $t_1$ and $t_2$ downward and to the right along the line

```{figure} ../figures/fig-9-3.png
:name: fig-9-3
:width: 80%
:align: center

**Figure 3.** Second-order difference equations.
```

$t_1 + t_2 = c$ in Figure 2. Using Figure 2, the values of $c$ and $\gamma$ compatible with damped oscillations, explosive oscillations, and so on, can easily be determined.

Figure 3 shows some realizations of second-order difference equations for various values of $t_1$ and $t_2$. In each case the values $Y_1$ and $Y_2$ were the two initial conditions which we specified arbitrarily. The reader can check that the behavior of these realizations matches that predicted by Figure 1.

(sec-9-3)=
## 3. Second-Order Difference Equations (Equal Roots)

The preceding treatment assumed that $\lambda_1 \neq \lambda_2$. (Notice that we divided by $\lambda_1 - \lambda_2$ to obtain {eq}`eq-9-19`.) If $\lambda_1 = \lambda_2$, then the polynomial we must study is

$$
\begin{aligned}
\frac{1}{(1 - \lambda L)(1 - \lambda L)} &= \frac{1}{1 - \lambda L}(1 + \lambda L + \lambda^2 L^2 + \cdots) \\
&= (1 + \lambda L + \lambda^2 L^2 + \cdots) + \lambda L(1 + \lambda L + \lambda^2 L^2 + \cdots) \\
&\quad + \lambda^2 L^2(1 + \lambda L + \lambda^2 L^2 + \cdots) + \cdots \\
&= 1 + 2\lambda L + 3\lambda^2 L^2 + \cdots,
\end{aligned}
$$

$$
\frac{1}{(1 - \lambda L)^2} = \sum_{i=0}^{\infty}(i + 1)\lambda^i L^i.
$$ (eq-9-30)

The lag distribution generated by the polynomial in {eq}`eq-9-30` is called a second-order Pascal lag distribution. It is the lag operator product or "convolution" of two geometric lag distributions with the same decay parameter[^fn-9-5] $\lambda$.

With the aid of {eq}`eq-9-30` we can study the solution to difference equations of the form

$$
(1 - \lambda L)^2 Y_t = a + bX_t.
$$

The general solution is

$$
Y_t = a\sum_{i=0}^{\infty}(i + 1)\lambda^i + b\sum_{i=0}^{\infty}(i + 1)\lambda^i X_{t-i} + c_1\lambda^t + c_2 t\lambda^t
$$ (eq-9-31)

where $c_1$ and $c_2$ are any real constants. To verify that {eq}`eq-9-31` is the solution, operate on both sides by $(1 - \lambda L)^2$ and verify that $(1 - \lambda L)^2\{c_1\lambda^t + c_2 t\lambda^t\} = 0$ for any choices of $c_1$ and $c_2$. To get a particular solution, we require two side conditions on the path of $Y_t$ to determine the two constants $c_1$ and $c_2$. For $a\sum_{i=0}^{\infty}(1 + i)\lambda^i$ to be finite either $|\lambda| < 1$ or $a = 0$ must be satisfied.

Assume that $a = 0$ and notice that {eq}`eq-9-31` can be rewritten for $t \geq 1$ as

$$
Y_t = b\sum_{i=0}^{t-1}(1 + i)\lambda^i X_{t-i} + b\sum_{i=t}^{\infty}(1 + i)\lambda^i X_{t-i} + c_1\lambda^t + c_2 t\lambda^t.
$$ (eq-9-32)

[^fn-9-5]: Let $w(L) = \sum_{j=-\infty}^{\infty}w_j L^j$, $v(L) = \sum_{j=-\infty}^{\infty}v_j L^j$. The *convolution* of the two sequences $\{w_j\}_{j=-\infty}^{\infty}$, $\{v_j\}_{j=-\infty}^{\infty}$ is the sequence whose $j$th element is $h_j = \sum_{k=-\infty}^{\infty}w_k v_{j-k}$. The reader should verify that $\sum_{j=-\infty}^{\infty}h_j L^j = h(L) = w(L)v(L) = v(L)w(L)$.


The second sum can be written as

$$
b \sum_{j=0}^{\infty} (j + 1 + t)\lambda^{t+j} X_{0-j} = b\lambda^t \sum_{j=0}^{\infty} (j + 1)\lambda^j X_{0-j} + bt\lambda^t \sum_{j=0}^{\infty} \lambda^j X_{0-j}.
$$

Therefore (32) becomes

$$
Y_t = b \sum_{i=0}^{t-1} (1 + i)\lambda^i X_{t-i} + \lambda^t \theta_0 + t\lambda^t \eta_0
$$ (eq-9-33)

where

$$
\theta_0 = \left\{ c_1 + b \sum_{j=0}^{\infty} (j + 1)\lambda^j X_{0-j} \right\} \quad \text{and} \quad \eta_0 = \left\{ c_2 + b \sum_{j=0}^{\infty} \lambda^j X_{0-j} \right\}.
$$

As for our earlier cases, (33) displays the solution for $t \geq 1$ as the sum of the distributed lag in $X_1, X_2, \ldots, X_t$ and two initial conditions. If $X_t = 0$ for $t \geq 1$, the sequence $Y_t$ will approach its stationary value of zero if $|\lambda| < 1$. If $|\lambda| > 1$, the sequence $Y_t$ will diverge from the stationary point of zero unless $\theta_0 = \eta_0 = 0$.

(sec-9-4)=
## 4. Nth-Order Difference Equations (Distinct Roots)

Consider a rational polynomial with $n$th-order denominator:

$$
A(L) = \frac{F(L)}{G(L)} = \frac{F(L)}{(1 - \lambda_1 L)(1 - \lambda_2 L) \cdots (1 - \lambda_n L)}.
$$

We assume that the degree of the numerator polynomial is less than the degree of the denominator polynomial. The zeros of $G(z) = 0$ are $z_1 = 1/\lambda_1$, $z_2 = 1/\lambda_2, \ldots, z_n = 1/\lambda_n$ since each of these values for $z$ satisfies the equation

$$
G(z) = (1 - \lambda_1 z)(1 - \lambda_2 z) \cdots (1 - \lambda_n z) = 0.
$$

Suppose the $n$ roots are distinct. Now the method of partial fractions enables us to express $A(L)$ as

$$
\frac{F(L)}{(1 - \lambda_1 L) \cdots (1 - \lambda_n L)} = \frac{A_1}{1 - \lambda_1 L} + \frac{A_2}{1 - \lambda_2 L} + \cdots + \frac{A_n}{1 - \lambda_n L}
$$ (eq-9-34)

where $A_1, A_2, \ldots, A_n$ are constants to be determined. To determine them, multiply both sides of the above equation by $(1 - \lambda_1 L) \cdots (1 - \lambda_n L)$ to get

$$
F(L) = A_1 (1 - \lambda_2 L) \cdots (1 - \lambda_n L) + A_2 (1 - \lambda_1 L)(1 - \lambda_3 L) \cdots (1 - \lambda_n L) + \cdots + A_n (1 - \lambda_1 L) \cdots (1 - \lambda_{n-1} L).
$$

Evaluating the above equation at $L = 1/\lambda_1$, the first root of $G(L) = 0$, gives

$$
A_1 = \frac{F(1/\lambda_1)}{(1 - (\lambda_2/\lambda_1)) \cdots (1 - (\lambda_n/\lambda_1))}.
$$

For the general term $A_i$, we obtain

$$
A_i = \frac{F(1/\lambda_i)}{(1 - (\lambda_1/\lambda_i)) \cdots (1 - (\lambda_{i-1}/\lambda_i))(1 - (\lambda_{i+1}/\lambda_i)) \cdots (1 - (\lambda_n/\lambda_i))}.
$$ (eq-9-35)

As an example, consider applying (34) to the second-order denominator polynomial

$$
A(L) = \frac{1}{(1 - \lambda_1 L)(1 - \lambda_2 L)} = \frac{A_1}{1 - \lambda_1 L} + \frac{A_2}{1 - \lambda_2 L}.
$$

We have

$$
A_1 = \frac{1}{1 - (\lambda_2/\lambda_1)} = \frac{\lambda_1}{\lambda_1 - \lambda_2} \quad \text{and} \quad A_2 = \frac{1}{1 - (\lambda_1/\lambda_2)} = \frac{\lambda_2}{\lambda_2 - \lambda_1}.
$$

Thus we obtain

$$
\frac{1}{(1 - \lambda_1 L)(1 - \lambda_2 L)} = \frac{1}{\lambda_1 - \lambda_2} \left( \frac{\lambda_1}{1 - \lambda_1 L} - \frac{\lambda_2}{1 - \lambda_2 L} \right),
$$

which we earlier used on the way to deriving (19).

Suppose we have an $n$th-order difference equation

$$
(1 - \lambda_1 L)(1 - \lambda_2 L) \cdots (1 - \lambda_n L) Y_t = b X_t.
$$ (eq-9-36)

The general solution to (36) is obtained by "dividing"[^fn-9-6] by $(1 - \lambda_1 L) \cdots (1 - \lambda_n L)$ to obtain

$$
Y_t = \frac{b}{(1 - \lambda_1 L) \cdots (1 - \lambda_n L)} X_t + c_1 \lambda_1^{\,t} + \cdots + c_n \lambda_n^{\,t},
$$

where $c_1, \ldots, c_n$ are any constants. That this is the solution can be verified by operating on both sides of the above equation with $(1 - \lambda_1 L) \cdots (1 - \lambda_n L)$. We now require $n$ side conditions on the path of $Y_t$ to determine the $n$ constants $c_1, \ldots, c_n$. We suppose that the $\lambda_j$ are all distinct. Then application of (34) to the above equation gives

$$
Y_t = b \sum_{r=1}^{n} \left( \frac{A_r}{1 - \lambda_r L} \right) X_t + \sum_{j=1}^{n} c_j \lambda_j^{\,t},
$$ (eq-9-37)

which shows that $Y_t$ can be expressed as the weighted sum of $n$ geometric distributed lags with decay coefficients $\lambda_1, \lambda_2, \ldots, \lambda_n$.

[^fn-9-6]: That is, "operating on both sides with the operator inverse of $(1 - \lambda_1 L)(1 - \lambda_2 L) \cdots (1 - \lambda_n L)$."

Given $n$ initial values of $Y$, and assuming $X_t = 0$ always, it is possible to start up difference equation (36) finitely far back in the past and to obtain a solution of the form

$$
Y_t = \lambda_1^{\,t} \eta_1 + \lambda_2^{\,t} \eta_2 + \cdots + \lambda_n^{\,t} \eta_n
$$

where $\eta_1, \ldots, \eta_n$ are constants chosen to satisfy the $n$ initial values. The above equation can be derived from (37) by applying calculations analogous to those applied above in the first- and second-order cases.

(sec-9-5)=
## 5. Nth-Order Difference Equations (N Equal Roots)

Consider the $n$th-order difference equation

$$
(1 - \lambda L)^n Y_t = b X_t,
$$ (eq-9-38)

which has the solution

$$
Y_t = \frac{b}{(1 - \lambda L)^n} X_t + c_1 \lambda^t + c_2 t \lambda^t + \cdots + c_n t^{n-1} \lambda^t.
$$

The polynomial $1/(1 - \lambda L)^n$ is the one associated with an $n$th-order Pascal lag distribution, which is formed by multiplying (convolving) $n$ geometric lag distributions with the same decay parameter $\lambda$. We have already studied the second-order Pascal distribution. The binomial expansion with negative exponent is

$$
(1 - x)^{-n} = 1 + nx + \frac{n(n + 1)}{2!} x^2 + \frac{n(n + 1)(n + 2)}{3!} x^3 + \cdots + \frac{n(n + 1) \cdots (n + i - 1)}{i!} x^i + \cdots, \qquad |x| < 1.
$$

Notice that $n(n + 1) \cdots (n + i - 1) = (n + i - 1)!/(n - 1)!$ The coefficient on $x$ in the above expansion is thus

$$
\binom{n + i - 1}{i} \equiv \frac{(n + i - 1)!}{i!(n - 1)!}.
$$

Therefore, the expansion can be written

$$
(1 - x)^{-n} = \sum_{i=0}^{\infty} \binom{n + i - 1}{i} x^i.
$$

Consequently, for a lag-generating function[^fn-9-7] with $n$ equal roots, we have

$$
\frac{1}{(1 - \lambda L)^n} = \sum_{i=0}^{\infty} \binom{n + i - 1}{i} \lambda^i L^i,
$$ (eq-9-39)

which agrees with our earlier formulas for the special cases $n = 1$ and $n = 2$.

To arrive at (39) in an alternative way, write

$$
f(\lambda L) = \frac{1}{1 - \lambda L} = \sum_{i=0}^{\infty} (\lambda L)^i.
$$

Differentiating with respect to $\lambda L$ gives

$$
f'(\lambda L) = \frac{1}{(1 - \lambda L)^2} = \sum_{i=0}^{\infty} i(\lambda L)^{i-1} = \sum_{i=0}^{\infty} (i + 1)(\lambda L)^i,
$$

which is (39) with $n = 2$. Differentiating again with respect to $\lambda L$ gives

$$
f''(\lambda L) = \frac{2}{(1 - \lambda L)^3} = \sum_{i=0}^{\infty} i(i - 1)(\lambda L)^{i-2}
$$

or

$$
\frac{1}{(1 - \lambda L)^3} = \sum_{i=0}^{\infty} \frac{(i + 1)(i + 2)}{2} (\lambda L)^i = \sum_{i=0}^{\infty} \binom{i + 2}{i} (\lambda L)^i,
$$

which is (39) with $n = 3$.

With the aid of (39), the solution to (38) can be written

$$
Y_t = b \sum_{i=0}^{\infty} \binom{n + i - 1}{i} \lambda^i X_{t-i} + c_1 \lambda^t + c_2 t \lambda^t + \cdots + c_n t^{n-1} \lambda^t.
$$

[^fn-9-7]: A term synonymous with "polynomial in the lag operator."

(sec-9-6)=
## 6. An Example of a First-Order System

Consider the following model studied by Cagan (1956). Let $m_t$ be the log of the money supply, $p_t$ the log of the price level and $p_{t+1}^{\,\mathrm{e}}$ the log of the price expected to prevail at time $t + 1$ given information available at time $t$. The model is

$$
m_t - p_t = \alpha(p_{t+1}^{\,\mathrm{e}} - p_t), \qquad \alpha < 0,
$$ (eq-9-40)

which is a portfolio equilibrium condition. The demand for real balances varies inversely with expected inflation $p_{t+1}^{\,\mathrm{e}} - p_t$. The variable $m_t$ is exogenous.

Suppose first that

$$
p_{t+1}^{\,\mathrm{e}} - p_t = \gamma(p_t - p_{t-1}),
$$ (eq-9-41)

so that the public expects inflation next period to be the current rate of inflation, $p_t - p_{t-1}$, multiplied by the constant $\gamma$. Then (40) becomes

$$
m_t - p_t = \alpha\gamma p_t - \alpha\gamma p_{t-1}.
$$

Using lag operators, this can be written as

$$
[(\alpha\gamma + 1) - \alpha\gamma L] p_t = m_t
$$

or

$$
\left[ 1 - \frac{\alpha\gamma}{1 + \alpha\gamma} L \right] p_t = \frac{1}{1 + \alpha\gamma} m_t.
$$

The solution can be written

$$
p_t = \frac{1}{1 + \alpha\gamma} \sum_{i=0}^{\infty} \left( \frac{\alpha\gamma}{1 + \alpha\gamma} \right)^i m_{t-i} + \left( \frac{\alpha\gamma}{1 + \alpha\gamma} \right)^t c
$$

where $c$ is a constant to be determined by, say, a given value of $p$ at some time. The solution for $p_t$ will be finite for the time path $m_t = \bar{m}$ for all $t$, provided that

$$
\left| \frac{\alpha\gamma}{1 + \alpha\gamma} \right| < 1.
$$

The above inequality is in the spirit of the "stability condition" developed by Cagan in his paper. It is a condition that delivers a bounded $p_t$ for all $t \geq 0$ for $\{m_t\}$ sequences that are bounded. Notice that

$$
\frac{1}{1 + \alpha\gamma} \sum_{i=0}^{\infty} \left( \frac{\alpha\gamma}{1 + \alpha\gamma} \right)^i = \frac{(1 + \alpha\gamma)^{-1}}{1 - \alpha\gamma(1 + \alpha\gamma)^{-1}} = 1.
$$

Thus, the long-run effect of a once-and-for-all jump in $m$ is to drive $p$ up by an equal amount (provided the above "stability condition" is met).

Returning to (40), let us abandon (41) and now assume perfect foresight:

$$
p_{t+1}^{\,\mathrm{e}} = p_{t+1}.
$$ (eq-9-42)

Substituting (42) into (40) gives

$$
m_t - p_t = \alpha p_{t+1} - \alpha p_t \quad \text{or} \quad \alpha p_{t+1} + (1 - \alpha) p_t = m_t.
$$

Write this as

$$
\left( L^{-1} + \frac{1 - \alpha}{\alpha} \right) p_t = \frac{1}{\alpha} m_t
$$

or

$$
\left( 1 - \frac{\alpha - 1}{\alpha} L \right) p_t = \frac{1}{\alpha} m_{t-1}.
$$ (eq-9-43)

Notice that since $\alpha < 0$, it follows that $(\alpha - 1)/\alpha > 1$. This fact is an invitation to solve (43) in the "forward" direction, i.e., to use (6). Dividing both sides of (43) by $(1 - ((\alpha - 1)/\alpha)L)$ gives

$$
p_t = \left( \frac{\alpha^{-1}}{1 - (\alpha - 1)\alpha^{-1} L} \right) m_{t-1} + c \left( \frac{\alpha - 1}{\alpha} \right)^t
$$

where $c$ is any constant. Using (6), this becomes

$$
\begin{aligned}
p_t &= \left( \frac{\alpha^{-1}(-\alpha(\alpha - 1)^{-1} L^{-1})}{1 - \alpha(\alpha - 1)^{-1} L^{-1}} \right) m_{t-1} + c \left( \frac{\alpha - 1}{\alpha} \right)^t \\
&= -\frac{1}{\alpha - 1} \left( \sum_{i=0}^{\infty} \left( \frac{\alpha}{\alpha - 1} \right)^i L^{-i} \right) m_t + c \left( \frac{\alpha - 1}{\alpha} \right)^t,
\end{aligned}
$$

$$
p_t = \frac{1}{1 - \alpha} \sum_{i=0}^{\infty} \left( \frac{\alpha}{\alpha - 1} \right)^i m_{t+i} + c \left( \frac{\alpha - 1}{\alpha} \right)^t.
$$ (eq-9-44)

Notice that since $\alpha < 0$, $0 < \alpha/(\alpha - 1) < 1$, so that the sum of the lag weights is finite. Equation (44) expresses the log of the current price as a moving sum of current and *future* values of the log of the money supply. Notice that

$$
\frac{1}{1 - \alpha} \sum_{i=0}^{\infty} \left( \frac{\alpha}{\alpha - 1} \right)^i = \frac{(1 - \alpha)^{-1}}{1 - \alpha(\alpha - 1)^{-1}} = 1,
$$

so that $p$ is a weighted *average* of current and future values of $m$. The solution for $p_t$ given by (44) will be finite for all finite $t$ if there is a constant $K > 0$ and an $x$, $1 \leq x < (\alpha - 1)/\alpha$ such that $|m_t| < Kx^t$ for all $t$. This is a condition that the money supply not grow too fast. For the solution given by (44) to be bounded for money supply paths that are bounded, we would require that $c = 0$. This amounts to an arbitrary terminal condition that rules out the occurrence of runaway inflations in the absence of runaway growth in the money supply.

(sec-9-7)=
## 7. An Example of a Second-Order System

Consider the following model studied by Muth (1961). Let $p_t$ be the price of a commodity at $t$, $C_t$ the demand for current consumption, $I_t$ the stock of inventories of the commodity, $Y_t$ the output of the commodity, and $p_t^{\,\mathrm{e}}$ the price previously expected to prevail at time $t$; $\{X_t\}$ is a bounded sequence of real numbers that represents the effects of the weather on supply. The model is

$$
\begin{aligned}
C_t &= -\beta p_t, & \beta &> 0 & \text{(demand curve)} \\
Y_t &= \gamma p_t^{\,\mathrm{e}} + X_t, & \gamma &> 0 & \text{(supply curve)} \\
I_t &= \alpha(p_{t+1}^{\,\mathrm{e}} - p_t), & \alpha &> 0 & \text{(inventory demand)} \\
Y_t &= C_t + (I_t - I_{t-1}) & & & \text{(market clearing)}.
\end{aligned}
$$

Let us suppose that there is perfect foresight so that $p_t^{\,\mathrm{e}} = p_t$ for all $t$. Making this assumption and substituting the first three equations into the fourth gives

$$
\gamma p_t + X_t = \alpha(p_{t+1} - p_t) - \alpha(p_t - p_{t-1}) - \beta p_t
$$

or

$$
\alpha p_{t+1} - (2\alpha + \beta + \gamma) p_t + \alpha p_{t-1} = X_t.
$$

Dividing by $\alpha$ gives

$$
p_{t+1} - \frac{2\alpha + \beta + \gamma}{\alpha} p_t + p_{t-1} = \alpha^{-1} X_t
$$

or

$$
(L^{-1} - \phi + L) p_t = \alpha^{-1} X_t
$$

where $\phi = ((\beta + \gamma)/\alpha) + 2 > 0$. Multiplying by $L$ gives

$$
(1 - \phi L + L^2) p_t = \alpha^{-1} X_{t-1}.
$$ (eq-9-45)

We need to factor the polynomial $1 - \phi L + L^2$ as

$$
1 - \phi L + L^2 = (1 - \lambda_1 L)(1 - \lambda_2 L) = 1 - (\lambda_1 + \lambda_2)L + \lambda_1 \lambda_2 L^2
$$

so that we require

$$
\lambda_1 + \lambda_2 = \phi, \qquad \lambda_1 \lambda_2 = 1.
$$

The second equality establishes that $\lambda_1 = 1/\lambda_2$, so that the two roots appear as a reciprocal pair. So we can write

$$
1 - \phi L + L^2 = (1 - \lambda L)(1 - \lambda^{-1} L)
$$

where $\lambda$ is chosen to satisfy $\lambda + \lambda^{-1} = \phi$. So (45) can be written

$$
(1 - \lambda L)(1 - \lambda^{-1} L) p_t = \alpha^{-1} X_{t-1}.
$$ (eq-9-46)

Since $(\beta + \gamma)/\alpha > 0$, it follows that $\phi = ((\beta + \gamma)/\alpha) + 2 > 2$. That implies that $\lambda$ does not equal 1 since $\lambda + \lambda^{-1} = \phi$. Notice that if $\lambda > 1$, $\lambda^{-1} < 1$. So one of our roots necessarily exceeds 1, the other necessarily is less than 1.

We divide both sides of (46) by $(1 - \lambda L)(1 - \lambda^{-1} L)$ to obtain

$$
p_t = \frac{1}{\alpha} \frac{1}{(1 - \lambda L)(1 - \lambda^{-1} L)} X_{t-1} + c_1 \lambda^t + c_2 \left( \frac{1}{\lambda} \right)^t
$$ (eq-9-47)

where $c_1$ and $c_2$ are any constants. Without loss of generality, suppose $\lambda < 1$ and let $\lambda_2 = 1/\lambda$. Use (6) and (19) to write the solution as

$$
\begin{aligned}
p_t &= \frac{\alpha^{-1} \lambda}{\lambda - \lambda_2} \left( \frac{1}{1 - \lambda L} \right) X_{t-1} - \frac{\alpha^{-1} \lambda_2}{(\lambda - \lambda_2)} \left( \frac{-(\lambda_2 L)^{-1}}{1 - (\lambda_2)^{-1} L^{-1}} \right) X_{t-1} + c_1 \lambda^t + c_2 \left( \frac{1}{\lambda} \right)^t \\
&= \frac{\alpha^{-1} \lambda}{\lambda - \lambda^{-1}} \left( \frac{1}{1 - \lambda L} \right) X_{t-1} + \frac{\alpha^{-1} \lambda^{-1}}{\lambda - \lambda^{-1}} \left( \frac{\lambda L^{-1}}{1 - \lambda L^{-1}} \right) X_{t-1} + c_1 \lambda^t + c_2 \lambda^{-t} \\
&= \frac{\alpha^{-1} \lambda}{\lambda - \lambda^{-1}} \left( \frac{1}{1 - \lambda L} \right) X_{t-1} + \frac{\alpha^{-1}}{\lambda - \lambda^{-1}} \left( \frac{1}{1 - \lambda L^{-1}} \right) X_t + c_1 \lambda^t + c_2 \lambda^{-t},
\end{aligned}
$$

$$
p_t = \frac{\alpha^{-1}}{\lambda - \lambda^{-1}} \sum_{i=1}^{\infty} \lambda^i X_{t-i} + \frac{\alpha^{-1}}{\lambda - \lambda^{-1}} \sum_{i=0}^{\infty} \lambda^i X_{t+i} + c_1 \lambda^t + c_2 \lambda^{-t},
$$

$$
p_t = \frac{\alpha^{-1}}{\lambda - \lambda^{-1}} \sum_{i=-\infty}^{\infty} \lambda^{|i|} X_{t-i} + c_1 \lambda^t + c_2 \lambda^{-t}.
$$ (eq-9-48)

The solution (48) expresses $p_t$ as a "two-sided" distributed lag of $X$, i.e., as a weighted sum of past, present, and future values of $X$. In this model the current price depends on the entire path of exogenous shock $X$ over the entire past *and* the entire future. As usual, the constants $c_1$ and $c_2$ are determined from two side conditions on the path of $p_t$. For example, if we were to impose the side conditions $\lim_{t \to -\infty} |p_t| < \infty$, $\lim_{t \to +\infty} |p_t| < \infty$, i.e., if we imposed boundedness on the $p_t$ path for all bounded $\{X_t\}$ sequences, then we would require $c_1 = c_2 = 0$.

(sec-9-8)=
## 8. An Optimization Example: Solving a System of Euler Equations

Thus far, our advice to solve stable roots backward and unstable roots forward has been to a certain extent arbitrary, being partly based on the desire to have solutions that are bounded. As it happens, when a linear difference equation emerges from a dynamic quadratic optimization problem, there are present additional necessary conditions for optimality which have precisely the effect of causing us to solve the stable root(s) backward and the unstable root(s) forward. We shall illustrate this for the case of a firm employing a single factor of production, labor, and facing a quadratic technology and quadratic costs of adjusting its labor force.

We consider a firm that at time $t$ chooses employment to maximize its present value

$$
v_t = \sum_{j=0}^{\infty} b^j \left\{ (f_0 + a_{t+j}) n_{t+j} - \frac{f_1}{2} n_{t+j}^2 - \frac{d}{2} (n_{t+j} - n_{t+j-1})^2 - w_{t+j} n_{t+j} \right\},
$$ (eq-9-49)


$n_{t-1}$ given. Here $f_0$, $f_1$, $d > 0$ and the discount factor $b$ obeys $0 < b < 1$. Here $w_{t+j}$ is the real wage faced by the firm at time $t + j$, while $n_{t+j}$ is the firm's employment at $t + j$. The firm faces a known sequence of $\{a_{t+j}\}_{j=0}^{\infty}$, $a_{t+j}$ being a shock to the marginal product of labor at time $t + j$. We assume for all $t$ that for some $K > 0$, $|w_t| < K(x)^t$, $|a_t| < K(x)^t$, where $1 \leq x < 1/\sqrt{b}$; sequences that satisfy these inequalities for some $K > 0$ and $1 \leq x < 1/\sqrt{b}$ will be termed of exponential order less than $1/\sqrt{b}$. The firm faces costs of adjusting its labor force rapidly, which we model by deducting the term $\frac{1}{2}d(n_{t+j} - n_{t+j-1})^2$ from real revenue each period. The firm faces known sequences $\{w_{t+j}\}_{j=0}^{\infty}$ and $\{a_{t+j}\}_{j=0}^{\infty}$ and is assumed to choose a sequence $\{n_{t+j}\}_{j=0}^{\infty}$ to maximize $v_t$.

In order to discover the appropriate solution to this problem it is convenient to consider the finite horizon problem: maximize

$$
v_t^T = \sum_{j=0}^{T} b^j \left\{ (f_0 + a_{t+j}) n_{t+j} - \frac{f_1}{2} n_{t+j}^2 - \frac{d}{2} (n_{t+j} - n_{t+j-1})^2 - w_{t+j} n_{t+j} \right\}
$$ (eq-9-50)

subject to $n_{t-1}$ being given. This problem matches our infinite horizon problem when we drive $T$ to infinity. If a sequence $\{n_{t+j}\}_{j=0}^{T}$ is to maximize {eq}`eq-9-50`, it must satisfy the following system of first-order necessary conditions, which we derive by differentiating {eq}`eq-9-50` with respect to $n_t, n_{t+1}, \ldots, n_{t+T}$:

$$
\begin{aligned}
f_0 + a_{t+j} - w_{t+j} - f_1 n_{t+j} &- d(n_{t+j} - n_{t+j-1}) \\
&+ db(n_{t+j+1} - n_{t+j}) = 0, \qquad j = 0, 1, \ldots, T-1,
\end{aligned}
$$ (eq-9-51)

$$
b^T [f_0 + a_{t+T} - w_{t+T} - f_1 n_{t+T} - d(n_{t+T} - n_{t+T-1})] = 0.
$$ (eq-9-52)

Equation {eq}`eq-9-51` is a system of second-order linear difference equations known as the "Euler equations," which we write as

$$
dbn_{t+j+1} - (f_1 + d(1+b)) n_{t+j} + dn_{t+j-1} = w_{t+j} - a_{t+j} - f_0
$$

or

$$
bn_{t+j+1} + \phi n_{t+j} + n_{t+j-1} = \frac{1}{d}(w_{t+j} - a_{t+j} - f_0), \qquad j = 0, 1, \ldots, T-1,
$$

$$
\phi = -\left( \frac{f_1}{d} + (1+b) \right).
$$ (eq-9-53)

To solve this second-order difference equation, we need two boundary conditions. One boundary condition is supplied by the historically given initial level of $n_{t-1}$, and the other by the terminal condition {eq}`eq-9-52`. That terminal condition is known as the transversality condition and is a necessary condition for optimality. To get the terminal condition in the infinite horizon problem, it is appropriate to take

$$
\lim_{T \to \infty} b^T [f_0 + a_{t+T} - w_{t+T} - f_1 n_{t+T} - d(n_{t+T} - n_{t+T-1})] n_{t+T} = 0
$$ (eq-9-54)

as the transversality condition. We obtain this condition by multiplying {eq}`eq-9-52` by the terminal stock of employment $n_{t+T}$, and taking the limit of the resulting equation as $T \to \infty$. Sufficient conditions for the transversality condition {eq}`eq-9-54` to hold are, first, that $\{a_{t+j}\}$ and $\{w_{t+j}\}$ be of exponential order less than $1/\sqrt{b}$, and, second, that the solution for $\{n_{t+j}\}$ be of exponential order less than $1/\sqrt{b}$. Thus, suppose that $\{w_t\}$, $\{a_t\}$, and $\{n_t\}$ are each of exponential order less than $1/\sqrt{b}$. Then notice that

$$
\begin{aligned}
|b^T [f_0 n_{t+T} &+ a_{t+T} n_{t+T} - w_{t+T} n_{t+T} - (f_1 + d) n_{t+T}^2 + d(n_{t+T} n_{t+T-1})]| \\
&\leq b^T f_0 |n_{t+T}| + b^T |a_{t+T} n_{t+T}| + b^T |w_{t+T} n_{t+T}| \\
&\quad + b^T (f_1 + d) |n_{t+T}^2| + b^T d |n_{t+T} n_{t+T-1}| \leq b^T f_0 K x^{t+T} \\
&\quad + b^T K x^{2(t+T)} + b^T K x^{2(t+T)} + b^T (f_1 + d) K x^{2(t+T)} \\
&\quad + b^T d K x^{2(t+T)-1} \\
&= f_0 K x^t b^T x^T + K(2 + f_1 + d) b^T x^{2T} x^{2t} + d K b^T x^{2T} x^{2t-1}.
\end{aligned}
$$

From the definition of a sequence of exponential order less than $1/\sqrt{b}$, we have that $0 < x\sqrt{b} < 1$. This condition implies that $xb < \sqrt{b} < 1$. Therefore the limit as $T \to \infty$ of the expression above equals zero. This proves that the conditions that $\{a_t\}$, $\{w_t\}$ and $\{n_t\}$ are of exponential order less than $1/\sqrt{b}$ are sufficient to imply that the transversality condition {eq}`eq-9-54` is satisfied. The necessary conditions for optimality for the infinite horizon problem are then satisfied if we can find a solution to the difference equation

$$
bn_{t+j+1} + \phi n_{t+j} + n_{t+j-1} = d^{-1}(w_{t+j} - a_{t+j} - f_0)
$$

subject to the transversality condition {eq}`eq-9-54` and the known initial value of $n_{t-1}$.

To solve the difference equation, write it as

$$
b \left( 1 + \frac{\phi}{b} L + \frac{1}{b} L^2 \right) n_{t+j+1} = \frac{1}{d}(w_{t+j} - a_{t+j} - f_0).
$$

We seek a factorization

$$
\left( 1 + \frac{\phi}{b} L + \frac{1}{b} L^2 \right) = (1 - \lambda_1 L)(1 - \lambda_2 L) = 1 - (\lambda_1 + \lambda_2) L + \lambda_1 \lambda_2 L^2.
$$

Equating powers of $L$ gives

$$
-\phi/b = (\lambda_1 + \lambda_2), \qquad 1/b = \lambda_1 \lambda_2, \qquad \text{or} \qquad 1/\lambda_1 b = \lambda_2.
$$

Thus we have that $\lambda_1$ must satisfy

$$
-\frac{\phi}{b} = \left( \lambda_1 + \frac{1}{\lambda_1 b} \right) \qquad \text{or} \qquad \frac{f_1}{d} + (1+b) = -\phi = \lambda_1 b + \frac{1}{\lambda_1}.
$$

```{figure} ../figures/fig-9-4.png
:name: fig-9-4
:width: 80%
:align: center

**Figure 4.**
```

Figure 4 depicts the determination of $\lambda_1$ and $\lambda_2$. The function $\lambda b + \lambda^{-1}$ attains a minimum at $\lambda = \sqrt{1/b}$, attaining a value of $2\sqrt{b}$ there. For $0 < b < 1$, we have that $1 + b \geq 2\sqrt{b}$ with equality at $b = 1$, so that even with $f_1 = 0$, in which case $\phi = -(1+b)$, we have $-\phi \geq 2\sqrt{b}$. This assures that with $f_1 > 0$ the solutions for $\lambda_1$ and $\lambda_2$ as depicted in Figure 4 are real and distinct. Without loss of generality, let $\lambda_1$ be the smaller root. Then notice that the preceding implies that $\lambda_1 < 1/\sqrt{b} < \lambda_2$. By using Figure 4 we can establish that $\lambda_1 < 1 < 1/b < \lambda_2$. From $-\phi = ((f_1/d) + 1 + b)$, we have that $-\phi > 1 + b$, since $f_1/d > 0$. Now if $-\phi$ were equal $1 + b$, we would have $\lambda_1$ and $\lambda_2$ being determined from

$$
1 + b = \lambda b + \lambda^{-1}
$$

or

$$
b\lambda^2 - (1+b)\lambda + 1 = 0, \qquad (\lambda b - 1)(\lambda - 1) = 0,
$$

so that $\lambda_1 = 1$, $\lambda_2 = 1/b$. Since $-\phi > 1 + b$, inspection of Figure 4 shows that unity must be an upper bound for $\lambda_1$ and $1/b$ a lower bound for $\lambda_2$.[^fn-9-8]

Having achieved this factorization, we write the difference equation as

$$
b(1 - \lambda_1 L)(1 - \lambda_2 L) n_{t+j+1} = d^{-1}(w_{t+j} - a_{t+j} - f_0).
$$

[^fn-9-8]: A version of Figure 4 with $b = 1$ determines the intriguing golden ratio or golden section. The golden ratio is the unique positive number $\lambda$ whose reciprocal equals itself plus unity: $\lambda^{-1} = 1 + \lambda$. This equation can be rearranged to read $\lambda^{-1} + \lambda = 1 + 2\lambda$. From the quadratic formula, the golden ratio equals $(\sqrt{5} - 1)/2 = 0.618034$ and is found as the intersection in the positive quadrant of the line $1 + 2\lambda$ with the curve $\lambda + \lambda^{-1}$. The golden ratio, which appears repeatedly in nature and mathematics, fascinated ancient people, and is reflected in the design of the Parthenon. One place that the number occurs in mathematics is as the limit as $t$ goes to infinity of the ratio of successive terms in a Fibonacci sequence. A Fibonacci sequence is generated by the difference equation $x_{t+1} = x_t + x_{t-1}$ with initial conditions $x_0 = 1$, $x_{-1} = 0$. The characteristic polynomial $(1 - L - L^2)$ associated with this equation can be factored as $(\lambda - L)(\lambda^{-1} + L)$ where $\lambda$ is the golden ratio. For more on the golden ratio, see "Math and Music: The Deeper Links," *The New York Times*, Sunday, August 29, 1982.

To satisfy the transversality condition, operate on both sides of this equation with the "forward" inverse of $1 - \lambda_2 L$ to get

$$
(1 - \lambda_1 L) n_{t+j+1} = \frac{-(b\, d\lambda_2)^{-1} L^{-1}}{1 - \lambda_2^{-1} L^{-1}} (w_{t+j} - a_{t+j} - f_0) + c\lambda_2^t
$$

where $c$ is a constant. However, since $\lambda_2 > 1/b > 1/\sqrt{b}$, we must set $c = 0$ in order to satisfy the transversality condition, for otherwise $\{n_{t+j}\}$ would fail to be of exponential order less than $1/\sqrt{b}$. Setting $c = 0$ and observing that $1/\lambda_2 = b\lambda_1$, we have

$$
n_{t+j+1} = \lambda_1 n_{t+j} - \frac{\lambda_1}{d} \sum_{i=0}^{\infty} \left( \frac{1}{\lambda_2} \right)^i (w_{t+j+1+i} - a_{t+j+1+i} - f_0), \qquad j = 0, 1, \ldots.
$$ (eq-9-55)

This is a demand schedule for employment that expresses current employment as a function of once-lagged employment and current and all future values of the real wage and the shock to productivity. Since $\lambda_2 > 1/b$, we have that $1/\lambda_2 < b$, so that the infinite weighted sum on the right-hand side of {eq}`eq-9-55` converges for $\{a_{t+j}\}$ and $\{w_{t+j}\}$ sequences that are of exponential order of less than $1/b$. We have assumed the stronger condition that $\{a_t\}$, $\{w_t\}$ are sequences of exponential order less than $1/\sqrt{b}$, so that the infinite weighted sum on the right side of {eq}`eq-9-55` converges.

It remains to show that {eq}`eq-9-55` holds for $j = -1$ as well. To show this, take the solution {eq}`eq-9-55` for $j = 0$ and use it to eliminate $n_{t+1}$ from the Euler equation {eq}`eq-9-53` for $j = 0$, thereby obtaining

$$
\{b\lambda_1 + \phi\} n_t + n_{t-1} = z_t + \frac{b\lambda_1}{1 - \lambda_2^{-1} L^{-1}} z_{t+1}
$$

where $z_t \equiv d^{-1}(w_t - a_t - f_0)$. Since $\lambda_1 + \lambda_2 = -\phi/b$ and $b\lambda_1 = \lambda_2^{-1}$, we can rearrange the above equation to read

$$
n_t = \lambda_1 n_{t-1} - \frac{\lambda_1}{1 - \lambda_2^{-1} L^{-1}} z_t,
$$

which asserts that {eq}`eq-9-55` holds for $j = -1$. Thus, {eq}`eq-9-55` gives the firm's optimal plan for setting $n_{t+j+1}$ for $j = -1, 0, 1, \ldots$.

We have in effect imposed the transversality condition first by using the "unstable" root $\lambda_2$ to solve in the forward direction and second by setting $c\lambda_2^t = 0$ in {eq}`eq-9-55`. To see this, write {eq}`eq-9-55` as

$$
(1 - \lambda_1 L) n_{t+j+1} = \varepsilon_{t+j+1}, \qquad j = -1, 0, 1, 2, \ldots,
$$ (eq-9-56)

where

$$
\varepsilon_{t+j+1} \equiv -\frac{\lambda_1}{d} \frac{1}{1 - \lambda_2^{-1} L^{-1}} (w_{t+j+1} - a_{t+j+1} - f_0).
$$

Since $\{w_{t+j}\}$ and $\{a_{t+j}\}$ are of exponential order less than $1/\sqrt{b}$, it follows[^fn-9-9] that the sequence $\{\varepsilon_{t+j+1}\}$ is of exponential order less than $1/\sqrt{b}$ and is finite for finite $t + j$. The solution of the difference equation {eq}`eq-9-55` with $n_{t-1}$ given is

$$
n_{t+j+1} = \lambda_1^{j+2} n_{t-1} + \sum_{i=0}^{j+1} \lambda_1^i \varepsilon_{t+j+1-i}, \qquad j = -2, -1, 0, 1, \ldots.
$$ (eq-9-57)

This is the unique solution to the Euler equation that satisfies the initial condition and the transversality condition {eq}`eq-9-54`. It is straightforward to verify that both the initial condition and the transversality condition are satisfied by this solution, the latter condition being satisfied by virtue of the exponential order being less than $1/\sqrt{b}$ for the $\{\varepsilon_{t+j}\}$ sequence and the fact that $|\lambda_1| < 1$. For notice that

$$
\left| \sum_{i=0}^{j+1} \lambda_1^i \varepsilon_{t+j+1-i} \right| < \sum_{i=0}^{j+1} \lambda_1^i K (x)^{t+j+1-i} \leq K(x)^{t+j+1} \left( \frac{1}{1 - (\lambda_1/x)} \right),
$$

which proves that $\sum_{i=0}^{j+1} \lambda_1^i \varepsilon_{t+j+1-i}$ is of exponential order less than $1/\sqrt{b}$. It follows that {eq}`eq-9-57` implies that $n_{t+j}$ is of exponential order less than $1/\sqrt{b}$, and therefore that the transversality condition {eq}`eq-9-54` is satisfied.

It should be emphasized that the transversality condition compels us to solve the unstable root forward. If we had solved the unstable root $\lambda_2$ "backward," then eventually the solution for $n_{t+j+1}$ would come to be dominated by terms of exponential order $\lambda_2$. Since $\lambda_2 > 1/\sqrt{b}$, that solution would violate the transversality condition {eq}`eq-9-54`.

In conclusion, the solution of our infinite horizon maximum problem is for the firm to set its employment according to the demand schedule or decision rule

$$
n_{t+j+1} = \lambda_1 n_{t+j} - \frac{\lambda_1}{d} \frac{1}{1 - \lambda_2^{-1} L^{-1}} (w_{t+j+1} - a_{t+j+1} - f_0)
$$

for $j = -1, 0, 1, 2, \ldots$.

[^fn-9-9]: For example, since $|w_{t+j+1}| < K(x)^{t+j+1}$ for some $K > 0$ and for $1 \leq x < 1/\sqrt{b}$, we have

    $$
    F_{t+j+1} \equiv \left| \sum_{i=0}^{\infty} \left( \frac{1}{\lambda_2} \right)^i w_{t+j+i+1} \right| < K \sum_{i=0}^{\infty} \left( \frac{1}{\lambda_2} \right)^i (x)^{t+j+i+1} = K \left( \frac{1}{1 - \lambda_2^{-1} x} \right) x^{t+j+1}
    $$

    where $1/\lambda_2 < \sqrt{b}$ and $x < 1/\sqrt{b}$. Therefore, $F$ is of exponential order less than $1/\sqrt{b}$.

(sec-9-9)=
## 9. A More General Univariate Optimization Problem

The problem of the preceding section is a special case of a more general quadratic optimization problem, one version of which follows. We define the polynomial in the lag operator,

$$
d(L) = d_0 + d_1 L + \cdots + d_m L^m,
$$ (eq-9-58)

where $d_0 \neq 0$, $d_m \neq 0$. We assume that $g_t$ is a sequence of exponential order less than $1/\sqrt{b}$, where $0 < b < 1$. Then the problem is to choose a sequence $\{y_t, t \geq 0\}$ to maximize

$$
\sum_{t=0}^{\infty} b^t \{ g_t y_t - \tfrac{1}{2} h y_t^2 - \tfrac{1}{2} [d(L) y_t]^2 \}
$$ (eq-9-59)

where $h > 0$, subject to $y_{-1}, y_{-2}, \ldots, y_{-m}$ given. The problem of Section 8 is a version of problem {eq}`eq-9-59` with $g_t = (f_0 + a_t - w_t)$, $m = 1$, $y_t = n_t$, and $d(L) = \sqrt{d}(1 - L)$, $h = f_1$.

The Euler equation for this problem is

$$
[d(bL^{-1}) d(L) + h] y_t = g_t.
$$ (eq-9-60)

We invite the reader to verify that this is the Euler equation by differentiating {eq}`eq-9-59` with respect to $y_t$ and rearranging. Note that {eq}`eq-9-53` is a special case of {eq}`eq-9-60` with $d(L) = (1 - L)\sqrt{d}$, $h = f_1$, $g_t = (f_0 + a_t - w_t)$.

In addition to the Euler equation {eq}`eq-9-60`, the necessary and sufficient conditions for maximization of {eq}`eq-9-59` are completed by the condition

$$
\frac{h}{2} \sum_{t=0}^{\infty} b^t y_t^2 < +\infty.
$$ (eq-9-61)

Inspection of {eq}`eq-9-59` reveals that any $\{y_t\}$ path that violates {eq}`eq-9-61`, even if it satisfies the Euler equation, gives a very bad outcome for the criterion function. Imposing conditions {eq}`eq-9-61` is equivalent to imposing transversality conditions. Transversality conditions could be derived as in Section 8, by differentiating a finite $T$ version of {eq}`eq-9-59` with respect to $y_T, \ldots, y_{T-m+1}$, multiplying the resulting $m$ expressions by $y_T, \ldots, y_{T-m+1}$, respectively, and setting the limits as $T \to \infty$ to zero. As do the transversality conditions, {eq}`eq-9-61` forces the solution $\{y_t\}$ to be of exponential order less than $1/\sqrt{b}$.

We briefly describe how {eq}`eq-9-60` is solved subject to {eq}`eq-9-61`. We first note a simple but important feature of the characteristic polynomial $h + d(bL^{-1}) d(L)$ that appears in {eq}`eq-9-60`. This polynomial evaluated at any value $z_0$ equals the polynomial evaluated at $b z_0^{-1}$. In particular, if $z_0$ is a zero of this polynomial, then so is $b z_0^{-1}$. To prove this, suppose that $z_0$ is a zero, which means setting $L = z_0$ gives $h + d(b z_0^{-1}) d(z_0) = 0$. The claim is that this implies that setting $L = b z_0^{-1}$ will also set the polynomial to zero. Making this substitution gives $h + d(z_0) d(b z_0^{-1})$, which equals the characteristic polynomial evaluated at $z_0$, which equals zero by assumption. Thus, the zeros of the characteristic polynomial of the Euler equation come in pairs of the form $z_k, b z_k^{-1}$, $k = 1, \ldots, m$.

Let us assume that the zeros of the characteristic polynomial are distinct. (The reader can apply the methods of Section 5 to analyze the case of repeated zeros, if this becomes necessary.) There are $2m$ zeros of the polynomial $h + d(b z^{-1}) d(z)$. Denote these zeros $z_1, z_2, \ldots, z_{2m}$. Without loss of generality, let the zeros be arranged in descending order according to their absolute value. Then we have that the zeros are ordered as $|z_1| \geq |z_2| \geq \cdots \geq |z_m| > |b z_m^{-1}| \geq |b z_{m-1}^{-1}| \geq \cdots \geq |b z_1^{-1}|$. (The weak inequalities leave open the possibility that there are complex roots.) Let us define $\lambda_k = 1/z_k$, $k = 1, \ldots, m$. Notice that the above ordering, and in particular $|z_m| > |b z_m^{-1}|$, implies that $|\lambda_k| < 1/\sqrt{b}$ for $k = 1, \ldots, m$. Then it is possible to factor the characteristic polynomial of the Euler equation as

$$
h + d(bL^{-1}) d(L) = c(bL^{-1}) c(L)
$$ (eq-9-62)

where

$$
c(L) = c_0 (1 - \lambda_1 L)(1 - \lambda_2 L) \ldots (1 - \lambda_m L)
$$ (eq-9-63)

$$
c_0 = [(-1)^m \lambda_0 / (\lambda_1 \lambda_2 \ldots \lambda_m)]
$$ (eq-9-64)

and where $\lambda_0$ is a constant that is uniquely determined by $h$ and $d(L)$.

Since $|\lambda_k| < 1/\sqrt{b}$ for $k = 1, \ldots, m$, {eq}`eq-9-62`–{eq}`eq-9-63` imply that $h + d(bL^{-1}) d(L)$ has been factored into $c(bL^{-1}) c(L)$, where the zeros of $c(z)$ exceed $\sqrt{b}$ in modulus, (from {eq}`eq-9-63` the zeros of $c(z)$ are the $1/\lambda_k$'s) and where the zeros of $c(bz^{-1})$ are *less* than $\sqrt{b}$ in modulus.

Using factorization {eq}`eq-9-62`, the Euler equation {eq}`eq-9-60` can be represented as

$$
c(bL^{-1}) c(L) y_t = g_t.
$$ (eq-9-65)

The unique solution of the Euler equation that satisfies {eq}`eq-9-61` is given by

$$
c(L) y_t = c(bL^{-1})^{-1} g_t.
$$ (eq-9-66)

or

$$
(1 - \lambda_1 L) \ldots (1 - \lambda_m L) y_t = c_0^{-2} [(1 - \lambda_1 b L^{-1}) \ldots (1 - \lambda_m b L^{-1})]^{-1} g_t.
$$

A partial fraction representation of $c_0 c(bL^{-1})^{-1}$ is given by

$$
[(1 - \lambda_1 b L^{-1}) \ldots (1 - \lambda_m b L^{-1})]^{-1} = \sum_{k=1}^{m} \frac{A_k}{(1 - \lambda_k b L^{-1})}
$$ (eq-9-67)

where

$$
A_k = \lim_{z \to \lambda_k b} (1 - \lambda_k b z^{-1}) c_0 / c(bz^{-1}).
$$ (eq-9-68)

Using {eq}`eq-9-67`, solution {eq}`eq-9-66` can be represented as

$$
(1 - \lambda_1 L) \ldots (1 - \lambda_m L) y_t = c_0^{-2} \sum_{k=1}^{m} \left( \frac{A_k}{1 - \lambda_k b L^{-1}} \right) g_t.
$$ (eq-9-69)

or

$$
(1 - \lambda_1 L) \ldots (1 - \lambda_m L) y_t = c_0^{-2} \sum_{k=1}^{m} A_k \sum_{j=0}^{\infty} (b\lambda_k)^j g_{t+j}.
$$ (eq-9-70)

Since $|\lambda_k| < 1/\sqrt{b}$, $|b\lambda_k| < \sqrt{b}$. This condition, together with our having assumed that $\{g_t\}$ is of exponential order less than $1/\sqrt{b}$ guarantees that the right hand side converges, and that as a function of $t$ it is itself of exponential order less than $1/\sqrt{b}$. It follows that the solution {eq}`eq-9-70` starting from the given starting values $y_{-1}, \ldots, y_{-m}$ satisfies condition {eq}`eq-9-61` and is therefore optimal.

We close this section with a remark on terminology. Consider solution {eq}`eq-9-66`, $c(L) y_t = c(bL^{-1})^{-1} g_t$. In the control literature, $c(L) y_t$ is called the "feedback part" of the solution for $y$, while $c(bL^{-1})^{-1} g_t$ is called the "feedforward part."

(sec-9-10)=
## 10. Introduction to Multivariate Dynamic Optimization

We now briefly describe a multivariate generalization of the problem of Section 9, and its solution. This section is a brief introduction to ideas described more fully in Hansen and Sargent [1981]. We define the matrix polynomial in the lag operator

$$
D(L) = D_0 + D_1 L + \cdots + D_m L^m
$$ (eq-9-71)

where $D_j$ is an $(n \times n)$ matrix. We let $\{G_t\}$ be a sequence of $(n \times 1)$ vectors, each component of which is a sequence of exponential order less than $1/\sqrt{b}$. We let $H$ be an $(n \times n)$ positive definite matrix. Finally, we let $\{Y_t\}$ be a sequence of $(n \times 1)$ vectors of variables that are to be chosen for $t \geq 0$, with given initial values $Y_{-1}, Y_{-2}, \ldots, Y_{-m}$. The problem we are interested in is to choose $\{Y_t, t \geq 0\}$ to maximize

$$
\sum_{t=0}^{\infty} b^t \{ G_t' Y_t - \tfrac{1}{2} Y_t' H Y_t - \tfrac{1}{2} [D(L) Y_t]' [D(L) Y_t] \}
$$ (eq-9-72)

given $\{G_t, t \geq 0\}$ and $Y_{-1}, \ldots, Y_{-m}$. In {eq}`eq-9-72`, the prime $(')$ denotes matrix transposition.

Necessary and sufficient conditions for a maximum of {eq}`eq-9-72` are

$$
\sum_{t=0}^{\infty} b^t Y_t' H Y_t < +\infty
$$ (eq-9-73)

and the Euler equations

$$
\{ H + D(bL^{-1})' D(L) \} Y_t = G_t.
$$ (eq-9-74)

That the Euler equations assume the form {eq}`eq-9-74` can be proved using the method of Sections 8 and 9, namely by differentiating {eq}`eq-9-72` with respect to $Y_t$, equating the result to a vector of zeros, and rearranging. Condition {eq}`eq-9-73` is the correct generalization of {eq}`eq-9-61`, and is justified by similar reasoning.

The related-pairs property of the zeros of the characteristic polynomial of the Euler equation that held in the univariate case generalizes as follows.


If $z_0$ is a zero of $\det\{H + D(bz^{-1})'D(z)\}$ then so is $bz_0^{-1}$. Here "det" denotes the determinant of a matrix. The appropriate matrix analogue of the scalar polynomial factorization (62) is a matrix polynomial factorization $C(bL^{-1})'C(L)$, where the zeros of $\det C(z)$ exceed $\sqrt{b}$ in absolute value, while those of $\det C(bz^{-1})$ are less than $\sqrt{b}$ in absolute value. By a theorem on matrix polynomials of the form that appears in (74), there always exists a matrix factorization

$$
H + D(bL^{-1})'D(L) = C(bL^{-1})'C(L)
$$ (eq-9-75)

where

$$
C(L) = C_0 + C_1 L + \cdots + C_m L^m
$$ (eq-9-76)

and where the zeros of $\det C(z)$ exceed $\sqrt{b}$ in absolute value.

The solution of (74) satisfying (73) is then given by

$$
C(L)Y_t = C(bL^{-1})'^{-1}G_t.
$$ (eq-9-77)

By the above mentioned property of the location of the zeros of $\det C(bz^{-1})$, we can represent $\det C(bL^{-1})'$, as

$$
\det C(bL^{-1})' = \lambda_0(1 - \lambda_1 bL^{-1})\ldots(1 - \lambda_k bL^{-1})
$$

where $|\lambda_j| < 1/\sqrt{b}$, $j = 1, \ldots, k$, where $k = mn$, and where we have assumed distinct zeros of $\det C(bz^{-1})$. Using a matrix version of partial fractions, $C(bL^{-1})'^{-1}$ can be represented as

$$
C(bL^{-1})'^{-1} = \sum_{h=1}^{k} \frac{A_h}{(1 - \lambda_h bL^{-1})}
$$ (eq-9-78)

where

$$
A_h = \lim_{z \to \lambda_h b} (1 - \lambda_h bz^{-1})C(bz^{-1})'^{-1}.
$$ (eq-9-79)

Note that each $A_h$ is an $(n \times n)$ matrix. Using (78), (77) can be represented

$$
C(L)Y_t = \sum_{h=1}^{k} A_h \sum_{j=0}^{\infty} (\lambda_h b)^j G_{t+j}.
$$ (eq-9-80)

Representation (80) is the vector generalization of (70).

Note in (77) a sort of symmetry in form between the "feedback part" $C(L)Y_t$ and the "feedforward part" $C(bL^{-1})'^{-1}G_t$ that generalizes a similar symmetry that we observed in the univariate problems.

A key step in solving problems of this sort is achieving the matrix factorization (75). Hansen and Sargent (1981) describe several methods for accomplishing this.

### Example (i): Interrelated Factor Demand

Consider the problem of a firm that maximizes

$$
\sum_{t=0}^{\infty} b^t\left\{q_t - w_t n_t - J_t k_t - \frac{d_1}{2}\left[(1 - L)n_t\right]^2 - \frac{d_2}{2}\left[(1 - L)k_t\right]^2\right\}, \quad 0 < b < 1,
$$ (eq-9-81)

subject to

$$
q_t = f_1'\binom{n_t}{k_t} - \frac{1}{2}\binom{n_t}{k_t}' F \binom{n_t}{k_t},
$$ (eq-9-82)

where $k_t$ is the stock of capital, $n_t$ is the stock of labor, $w_t$ is the real rental on labor, and $J_t$ is the real rental on capital; $f_1$ is a $(2 \times 1)$ vector of positive constants, $F$ is a positive definite matrix, and $d_1$ and $d_2$ are positive constants measuring adjustment costs. We assume that $\{J_t\}$ and $\{w_t\}$ for $t \geq 0$ are known sequences of exponential order less than $1/\sqrt{b}$. The problem is to choose sequences $\{k_t, n_t, t \geq 0\}$ to maximize (81) subject to (82), given initial values $n_{-1}$, $k_{-1}$ and given sequences for $w_t$ and $J_t$.

Problem (81)–(82) is a special case of (72) with

$$
G_t = f_1 - \begin{bmatrix} w_t \\ J_t \end{bmatrix},
$$

$$
D(L) = \begin{bmatrix} \sqrt{d_1}\,(1 - L) & 0 \\ 0 & \sqrt{d_2}\,(1 - L) \end{bmatrix},
$$

and $H = F$. The solution (80) is an interrelated pair of decision rules for $(n_t, k_t)$ of the form

$$
C(L)\begin{bmatrix} n_t \\ k_t \end{bmatrix} = \sum_{h=1}^{2} A_h \sum_{j=0}^{\infty} (\lambda_h b)^j\left\{f_1 - \begin{bmatrix} w_{t+j} \\ J_{t+j} \end{bmatrix}\right\}
$$ (eq-9-83)

where $C(L) = C_0 + C_1 L$. In (83) the decision rules for capital and labor interact in the sense that each of $(n_t, k_t)$ depends on lagged values of the other, and that each depends on future rental rates for the other. This interdependence occurs so long as either $F$ (or $H$) or $D(L)$ is not diagonal.

Versions of this problem were studied and utilized by Nadiri and Rosen (1973), Hansen and Sargent (1981), and Eichenbaum (1983).

### Example (ii): A Dynamic Nash Equilibrium[^fn-9-10]

We consider a duopoly in an industry producing a single good. Demand is governed by a linear demand schedule,

$$
p_t = A_0 - \frac{A_1}{2}(q_{1t} + q_{2t}) + u_t, \quad A_0, A_1 > 0,
$$ (eq-9-84)

[^fn-9-10]: A variety of linear quadratic dynamic games is analyzed by Hansen, Epple, and Roberds (1985). These authors use lag operator methods to obtain explicit solutions of their games, and to extract the econometric restrictions implied by their games.

where $q_{it}$ is output of firm $i$ at $t$, $p_t$ is the price of output at $t$, and $u_t$ is a given sequence of disturbances to demand of exponential order less than $1/\sqrt{b}$. Firm $i$ maximizes

$$
\sum_{t=0}^{\infty} b^t\{p_t q_{it} - q_{it}s_{it} - \tfrac{1}{2}[d_i(L)q_{it}]^2\}, \quad 0 < b < 1,
$$ (eq-9-85)

where $d_i(L) = d_{i0} + \cdots + d_{im}L^m$, and where $s_{it}$ is a sequence of shocks to costs of production of firm $i$, assumed to be a known sequence of exponential order less than $1/\sqrt{b}$. Here $[d_i(L)q_{it}]^2$ stands for costs of adjusting production rapidly. The maximization of (85) by $i$ takes $\{s_{it}\}$, $\{u_t\}$, and $\{q_{jt}, j \neq i\}$ given for $t \geq 0$, and $q_{it}$ given for $\{t = -m, \ldots, -1\}$, $i = 1, 2$. In particular, firm $i$ is supposed to regard firm $j$'s output sequence as given and beyond its influence.

Substituting (84) into (85) gives

$$
\sum_{t=0}^{\infty} b^t\left\{\left[A_0 - \frac{A_1}{2}(q_{1t} + q_{2t}) + u_t\right]q_{it} - q_{it}s_{it} - \tfrac{1}{2}[d_i(L)q_{it}]^2\right\}.
$$ (eq-9-86)

The Euler equations for this problem for firms 1 and 2 are

$$
\begin{bmatrix} (A_1 + d_1(bL^{-1})d_1(L)) & \dfrac{A_1}{2} \\[2ex] \dfrac{A_1}{2} & (A_1 + d_2(bL^{-1})d_2(L)) \end{bmatrix} \begin{bmatrix} q_{1t} \\ q_{2t} \end{bmatrix}
$$

$$
= \begin{bmatrix} A_0 + u_t - s_{1t} \\ A_0 + u_t - s_{2t} \end{bmatrix}.
$$ (eq-9-87)

We define a Nash equilibrium in the space of sequences of quantities $q_{1t}$, $q_{2t}$ as a pair of sequences for $q_{1t}$, $q_{2t}$[^fn-9-11] that solve the interrelated Euler equations (87) and satisfy the boundary conditions

$$
\sum_{t=0}^{\infty} b^t q_{it}^2 < +\infty \quad \text{for } i = 1, 2.
$$ (eq-9-88)

Equivalently, the Nash equilibrium is the pair of $(q_{1t}, q_{2t})$ sequences that satisfies the following conditions:

(i) Firm $i$'s quantity sequence maximizes its present value (85), given firm $j$'s sequence, for $(i, j) = (1, 2)$ and $(2, 1)$.

(ii) The output market clears, in the sense that (84) holds for all $t$.

Equation (87) is evidently in the form of a vector Euler equation in $(q_{1t}, q_{2t})$. The matrix polynomial on the left side of (87) can be factored into the form

[^fn-9-11]: It is necessary to add the qualifier "in the space of sequences of quantities $q_{1t}$, $q_{2t}$" because different definitions of strategy spaces in general give rise to distinct Nash equilibria.

$C(bL^{-1})'C(L)$, where $C(L)$ is a $(2 \times 2)$ matrix polynomial with the zeros of $\det C(z)$ exceeding $\sqrt{b}$ in absolute value. Then the Nash equilibrium can be represented

$$
C(L)\begin{bmatrix} q_{1t} \\ q_{2t} \end{bmatrix} = C(bL^{-1})'^{-1}\begin{bmatrix} A_0 + u_t - s_{1t} \\ A_0 + u_t - s_{2t} \end{bmatrix}
$$

where

$$
C(L) = C_0 + C_1 L + \cdots + C_m L^m.
$$

Note that the vector Euler equation (87) was discovered by our having solved an interrelated pair of univariate dynamic optimization problems. The resulting system of Euler equations (87) can itself be viewed as solving some vector optimization problem, since the matrix characteristic polynomial is expressible in the form $H + D(bL^{-1})'D(L)$. To seek the vector dynamic optimization problem that is implicitly solved by (87) is to pose a version of an "inverse optimal control problem:" given a system of difference equations, attempt to synthesize an optimum problem for which they are necessary conditions.[^fn-9-12][^fn-9-13]

[^fn-9-12]: We will encounter another "inverse optimal control problem" in Chapter XIV. There we will study the implicit social welfare criteria that are maximized by various versions of an equilibrium model of investment under uncertainty. In Chapter XI, we will study a technically related problem called the "inverse optimal predictor problem." Such a problem was posed by John F. Muth (1960). Given the expectations formation scheme utilized by Cagan (1956), Muth sought to discover a random environment which would render such a forecasting scheme optimum.

[^fn-9-13]: In Chapter XV, we shall study a dominant player dynamic game in some detail.

(sec-9-11)=
## 11. Another Example: Learning by Doing

The following example is interesting from both a substantive and a technical point of view. From a technical point of view, the example will introduce a control problem in which a rational characteristic polynomial must be factored in order to solve the Euler equation.

We consider a monopolist who faces the linear law of demand

$$
p_t = A_0 - (A_1/2)Q_t + u_t, \quad A_0, A_1 > 0,
$$ (eq-9-89)

where $p_t$ is price, $Q_t$ is output, and $u_t$ is a known sequence of shocks to demand, assumed to be of exponential order less than $1/\sqrt{b}$. The monopolist's total costs at $t$ are given by

$$
C(Q_t) = c_0 + c_1 Q_t + \frac{c_2}{2}Q_t^2 + c_3 Q_t s_t - c_4[h(L)Q_t]Q_t, \quad c_0, c_1, c_2, c_3, c_4 > 0,
$$ (eq-9-90)

where $s_t$ is a known sequence of shocks to marginal cost, assumed to be of exponential order less than $1/\sqrt{b}$, where $h(L) = \sum_{j=0}^{\infty} h_j L^j$. As an example of the sort of $h(L)$ that we have in mind, we shall later consider the special case $h(L) = 1/(1 - \rho L)$ where $\rho < 1$ but where $\rho$ is close to one. Then (90) captures the notion that marginal costs of current output fall with cumulated past output. This is one version of a "learning-by-doing" cost structure.

The firm maximizes

$$
\sum_{t=0}^{\infty} b^t\{p_t Q_t - C(Q_t)\}, 0 < b < 1,
$$

or using (89) and (90)

$$
\sum_{t=0}^{\infty} b^t\left\{\left[A_0 - \frac{A_1}{2}Q_t + u_t\right]Q_t - \left[c_0 + c_1 Q_t + \frac{c_2}{2}Q_t^2 + c_3 Q_t s_t - c_4[h(L)Q_t]Q_t\right]\right\}.
$$ (eq-9-91)

The firm chooses a sequence $\{Q_t\}_{t=0}^{\infty}$ to maximize (91), taking as given $Q_t$ for $t < 0$ and the sequences $\{u_t\}_{t=0}^{\infty}$, $\{s_t\}_{t=0}^{\infty}$. By using the preceding method, the Euler equation for the firm's problem is found to be, after some rearrangement,

$$
[(A_1 + c_2) - c_4 h(L) - c_4 h(bL^{-1})]Q_t = A_0 + u_t - c_1 - c_3 s_t.
$$ (eq-9-92)

For convenience, set $c_1 = 0$. The boundary condition is

$$
A_1 \sum_{t=0}^{\infty} b^t Q_t^2 < +\infty.
$$ (eq-9-93)

We shall now consider the special case of the model that emerges when we set $h(L) = 1/(1 - \rho L)$, $0 < \rho < 1$. In this case, the characteristic polynomial on the left side of the Euler equation (92) becomes the rational polynomial

$$
k - \frac{c_4}{1 - \rho L} - \frac{c_4}{1 - \rho b L^{-1}}
$$

where $k \equiv (A_1 + c_2)$. To solve the Euler equation (92) subject to boundary condition (93), the first step is to express the characteristic polynomial in terms of a common denominator, which gives

$$
\frac{[k(1 + \rho^2 b) - 2c_4] - (k - c_4)\rho L - (k - c_4)b\rho L^{-1}}{(1 - \rho b L^{-1})(1 - \rho L)}.
$$ (eq-9-94)

Notice that the denominator is already factored, and that the zeros of the numerator come in the familiar pair $(z_0, bz_0^{-1})$. Our next step on the way to solving our problem is to factor the numerator. Note that the numerator is proportional to

$$
\left[-L^{-1} + \left[\frac{k(1 + \rho^2 b) - 2c_4}{(k - c_4)\rho b}\right] - 1/bL\right] \equiv \frac{1}{\alpha b}(1 - \alpha b L^{-1})(1 - \alpha L)
$$ (eq-9-95)

where $1/\alpha$ is the zero of the characteristic polynomial that exceeds $1/\sqrt{b}$. Equating powers of $L$, as in section 8, shows that $\alpha$ must solve

$$
\frac{1}{\alpha b} + \alpha = \left[\frac{k(1 + \rho^2 b) - 2c_4}{(k - c_4)\rho b}\right].
$$ (eq-9-96)

We assume that the parameter $k \equiv A_1 + c_2$, $c_4$, $\rho$, $b$ are such that the right side exceeds $2\sqrt{b}$ in absolute value. This guarantees the existence of a real number $\alpha$ that solves (96). Note that for $\rho = b = 1$, the above equation has the solution $\alpha = 1/\alpha = 1$. By continuity of the solution in the arguments on the right hand side, for values of $\rho$ and $b$ sufficiently close to 1, $\alpha$ will be close to one.

Using (94) and (95), the Euler equation can be expressed

$$
\frac{(1 - \alpha b L^{-1})(1 - \alpha L)}{(1 - \rho b L^{-1})(1 - \rho L)}Q_t = \frac{\alpha b}{(k - c_4)\rho b}[A_0 + u_t - c_3 s_t].
$$ (eq-9-97)

The solution of the Euler equation (97) that satisfies boundary condition (93) can be expressed in "feedback-feedforward" form

$$
\frac{(1 - \alpha L)}{(1 - \rho L)}Q_t = \frac{\alpha b}{(A_1 + c_2 - c_4)\rho b} \cdot \frac{(1 - \rho b L^{-1})}{(1 - \alpha b L^{-1})}[A_0 + u_t - c_3 s_t].
$$ (eq-9-98)

Unless $\alpha = \rho$, $Q_t$ feeds back on an infinite number of its own past values, reflecting the dynamics of the firm's optimally coping with the learning-by-doing cost structure. In general, $\rho \neq \alpha$. However, in the special limiting case $\rho = b = 1$, it can be verified that $\rho = \alpha$. In this special case, (98) collapses to the static decision rule

$$
Q_t = \frac{1}{A_1 + c_2 - c_4}[A_0 + u_t - c_3 s_t],
$$

despite the presence of the learning-by-doing cost structure.

(sec-9-12)=
## 12. A Consumption Example

We begin with a technical remark. Suppose that we wish to differentiate with respect to $x_t$ the expression

$$
V = \sum_{t=0}^{\infty} b^t[a(L)z_t][d(L)x_t]
$$

where $a(L) = \sum_{j=-\infty}^{\infty} a_j L^j$, $d(L) = \sum_{j=-\infty}^{\infty} d_j L^j$. We invite the reader to verify that

$$
\frac{\partial V}{\partial x_t} = b^t d(bL^{-1})a(L)z_t.
$$ (eq-9-99)

We shall apply this formula in studying an optimal consumption problem.

A consumer is assumed to maximize

$$
\sum_{t=0}^{\infty} b^t\left\{u_0 + u_1 c_t - \frac{u_2}{2}c_t^2\right\}, \quad 0 < b < 1, \quad u_0, u_1, u_2 > 0,
$$ (eq-9-100)

subject to

$$
A_{t+1} = R[A_t + y_t - c_t], \quad R \geq 1,
$$ (eq-9-101)

and

$$
A_t > B \quad \text{for all } t \geq 0,
$$ (eq-9-102)

$$
A_0 \quad \text{given.}
$$ (eq-9-103)

Here $c_t$ is consumption, $y_t$ is labor income, $A_t$ is nonhuman assets available at the beginning of period $t$, and $R$ is the gross rate of return on nonhuman assets. We assume that $0 > B > -\infty$, so that (102) requires that assets be uniformly bounded from below. This condition rules out a strategy of consuming large amounts and financing it by borrowing along an increasing and unbounded path.

Use (101) to express $c_t$ as

$$
c_t = y_t - (R^{-1}L^{-1} - 1)A_t.
$$ (eq-9-104)

Substituting (104) into (100) gives the criterion

$$
J = \sum_{t=0}^{\infty} b^t\left\{u_0 + u_1(y_t - (R^{-1}L^{-1} - 1)A_t) - \frac{u_2}{2}\left(y_t^2 - 2y_t(R^{-1}L^{-1} - 1)A_t + [(R^{-1}L^{-1} - 1)A_t]^2\right)\right\},
$$ (eq-9-105)

which is to be maximized by selecting a sequence for $A_t$, subject to the boundary conditions (102) and (103). Repeatedly applying formula (99), we find the Euler equation

$$
-\frac{u_1}{u_2}(R^{-1}b^{-1} - 1) + (R^{-1}b^{-1}L - 1)y_t - (R^{-1}b^{-1}L - 1)(R^{-1}L^{-1} - 1)A_t = 0
$$ (eq-9-106)

$$
\text{for } t = 1, 2, \ldots.
$$

The zeros of the characteristic polynomial in $A_t$ are at $L = R^{-1} < 1$ and $L = Rb$. In order to satisfy the boundary conditions (102) and (103), the polynomial $(R^{-1}L^{-1} - 1)$ must be solved forwards, and the polynomial $(R^{-1}b^{-1}L - 1)$ solved backwards. Dividing (106) by $(R^{-1}L^{-1} - 1)$ gives

$$
(R^{-1}b^{-1}L - 1)A_t = \frac{-u_1}{u_2}\frac{(R^{-1}b^{-1} - 1)}{(R^{-1} - 1)} + \left(\frac{R^{-1}b^{-1}L - 1}{R^{-1}L^{-1} - 1}\right)y_t
$$

or

$$
A_t = R^{-1}b^{-1}A_{t-1} + \frac{u_1}{u_2}\frac{(R^{-1}b^{-1} - 1)}{(R^{-1} - 1)} + \left(\frac{1 - R^{-1}b^{-1}L}{R^{-1}L^{-1} - 1}\right)y_t.
$$ (eq-9-107)

Subtracting (107) from (101) lagged one period, namely,

$$
A_t = R[A_{t-1} + y_{t-1} - c_{t-1}],
$$

gives

$$
0 = (R - R^{-1}b^{-1})A_{t-1} - \frac{u_1}{u_2}\frac{(R^{-1}b^{-1} - 1)}{R^{-1} - 1} + \frac{R(R^{-1}L^{-1} - 1)}{(R^{-1}L^{-1} - 1)}Ly_t - \frac{(1 - R^{-1}b^{-1}L)}{(R^{-1}L^{-1} - 1)}y_t - Rc_{t-1}
$$

or

$$
c_{t-1} = (1 - R^{-2}b^{-1})A_{t-1} - \frac{u_1}{u_2}\frac{(R^{-1}b^{-1} - 1)}{1 - R} + \frac{(1 - R^{-2}b^{-1})}{1 - L^{-1}R^{-1}}y_{t-1}.
$$ (eq-9-108)

Writing out the polynomial in $L$ explicitly and shifting forward one period gives

$$
c_t = (1 - R^{-2}b^{-1})\left[A_t + \sum_{j=0}^{\infty} R^{-j}y_{t+j}\right] - \frac{u_1}{u_2}\left(\frac{R^{-1}b^{-1} - 1}{1 - R}\right).
$$ (eq-9-109)

This is a certainty version of the permanent income theory of consumption. According to (109), consumption is a linear function of wealth, defined as the sum of nonhuman wealth $A_t$ and the discounted present value of labor income, $\sum_{j=0}^{\infty} R^{-j}y_{t+j}$. In the special case that $bR = 1$, (109) collapses to

$$
c_t = (1 - R^{-1})\left[A_t + \sum_{j=0}^{\infty} R^{-j}y_{t+j}\right].
$$

It is also true that when $Rb = 1$, (109) and (101) imply that $c_t = c_0$ for all $t$, so that the model implies complete consumption smoothing. To convince oneself of this, the reader is invited to verify that the Euler equation of this model can be represented as

$$
c_t = (Rb)^{-1}c_{t-1}.
$$

Robert Hall [1978] pursued the empirical implications of a stochastic version of such an Euler equation. In Chapter XII below, we shall study this consumption example further.


(sec-9-13)=
## 13. Nonanticipative Representations and Lucas's Critique

We return to a version of the difference equation associated with Cagan's model of hyperinflation, which we studied in Section 6 above. In particular, we have

$$
y_t = \lambda y_{t+1} + x_t, \qquad |\lambda| < 1,
$$

where $\{x_t\}$ is a sequence of exponential order less than $1/\lambda$. The unique bounded solution of this equation is given by

$$
y_t = \sum_{j=0}^{\infty} \lambda^j x_{t+j}.
$$ (eq-9-110)

This equation represents the solution for $y_t$ in terms of an infinite number of future values of the "forcing variable" $x_t$. In practice, one often prefers a "nonanticipative" representation of a solution, in which $y_t$ is expressed as a solution only of $x$'s dated $t$ and earlier. We now briefly show how representation (110) can be converted to an alternative representation expressing $y_t$ as such a function of current and past $x$'s only. Given (110), such a representation in terms of current and past $x$'s can be attained by positing a specific difference equation for $\{x_t\}$ itself, using this difference equation to express $x_{t+j}$ for $j \geq 1$ as functions of $\{x_t, x_{t-1}, \ldots\}$, and then using these functions to eliminate $x_{t+j}$ for $j \geq 1$ from (110). The parameters of the representation for $y_t$ as a function of current and past $x$'s will depend on the parameters of the difference equation for $x_t$. This is so because the representation expressing $y_t$ as a function of $\{x_t, x_{t-1}, \ldots\}$ incorporates the solution of the difference equation expressing $x_{t+j}$, $j \geq 1$, in terms of $\{x_t, x_{t-1}, x_{t-2}, \ldots\}$. For example, suppose that $x_t$ is governed by the difference equation

$$
x_t = \rho x_{t-1}, \qquad |\rho \lambda| < 1.
$$

Then we have that

$$
x_{t+j} = \rho^j x_t, \qquad j \geq 1.
$$

Substituting this into (110) gives

$$
y_t = \left( \frac{1}{1 - \lambda \rho} \right) x_t,
$$ (eq-9-111)

which is an expression for $y_t$ in terms of $x_t$ in which the parameter $\rho$ of the difference equation for $x_t$ appears.

More generally, assume that $x_t$ is governed by the difference equation

$$
a(L) x_t = 0
$$ (eq-9-112)

where $a(L) = 1 - a_1 L - a_2 L^2 - \cdots - a_r L^r$ and where the zeros of $a(z)$ exceed $1/\lambda$ in absolute value. Then Hansen and Sargent (1980) show that the solution of (110) can be represented as

$$
y_t = a(\lambda)^{-1} \left[ 1 + \sum_{j=1}^{r-1} \left( \sum_{k=j+1}^{r} \lambda^{k-j} a_k \right) L^j \right] x_t.
$$ (eq-9-113)

(We shall show how to derive this formula in Chapter XI.) Equation (113) is evidently a generalization of (111). That is, (111) is (113) for the case $r = 1$, $a_0 = 1$, $a_1 = \rho$.

Equations (112) and (113) display a hallmark of economic models incorporating foresight as in (110), namely the presence of cross-equation restrictions between, on the one hand, the parameters of the process $\{x_t\}$ of the forcing variables, and, on the other hand, the solution for $y_t$ in terms of $\{x_t, x_{t-1}, \ldots\}$. This characteristic of models with foresight means that it will not be possible to find a representation expressing $y_t$ as a function of $\{x_t, x_{t-1}, \ldots\}$ that is independent of the law of motion (difference equation) governing the forcing process $\{x_t\}$. Alterations in the law of motion for $x_t$ will alter the *function* describing the dependence of $y_t$ on current and past $x$'s.

The presence of such cross-equation restrictions in models with foresight is the source of Robert E. Lucas's (1976) criticism of some econometric policy evaluation procedures of the past. Lucas criticized methods that took representations of a variable like $y_t$ as a function of $\{x_t, x_{t-1}, \ldots\}$ and treated them as being invariant with respect to alterations in the difference equation generating the $\{x_t\}$ sequence.

(sec-9-exercises)=
## Exercises

**1.** Verify that the presence of both *lagged* employment and *future* values of $w_t$ and $a_t$ on the right-hand side of the employment decision rule (55) (demand schedule) depends on having the adjustment cost parameter $d$ strictly positive. (Set $d = 0$ and rework the firm's optimum problem.)

**2.** Determine the effect of an increase in $d$ on the speed of adjustment parameters $\lambda_1$ and $\lambda_2$. Does a firm facing a small $d$ adjust its labor force more or less quickly in response to current conditions than a firm facing a larger value of $d$? (*Hint:* use Figure 4.)

**3.** (*A Keynesian investment schedule*) A firm chooses a sequence of capital $\{k_{t+j}\}_{j=0}^{\infty}$ to maximize

$$
v_t = \sum_{j=0}^{\infty} b^j \left\{ a_0 k_{t+j} - \frac{a_1}{2} k_{t+j}^2 - J_{t+j}(k_{t+j} - k_{t+j-1}) - \frac{d}{2}(k_{t+j} - k_{t+j-1})^2 \right\},
$$

given $k_{t-1} > 0$, where $a_0, a_1, d > 0$ and where $\{J_{t+j}\}_{j=0}^{\infty}$ is a known sequence of the price of capital relative to the price of the firm's output. The sequence $\{J_{t+j}\}$ is of exponential order less than $1/b$, where $0 < b < 1$ is the discount factor (the reciprocal of one plus the real rate of interest).

&nbsp;&nbsp;&nbsp;&nbsp;**A.** Derive the Euler equations and the transversality condition.

&nbsp;&nbsp;&nbsp;&nbsp;**B.** Show that the optimum decision rule of the firm is of the form

$$
k_{t+j+1} = \lambda_1 k_{t+j} + c_0 + \frac{c_1}{1 - \lambda_2^{-1} L^{-1}} (b J_{t+j+2} - J_{t+j+1}),
$$

where $0 < \lambda_1 < 1 < \lambda_2$, and $c_0$ and $c_1$ are constants. Find $c_0$ and $c_1$; show how to find $\lambda_1$ and $\lambda_2$; and prove that the $\lambda$'s obey the inequalities just stated.

**4.** (*Keynesian stabilization policy*) The reduced form for GNP ($Y$) is

$$
Y_t = \alpha + B S_t + c g_t, \qquad B > 0, \quad c > 0,
$$

where $g$ is government purchases and $S_t$ is exports, an exogenous variable outside the government's control. The sequence of exports $\{S_t\}_{t=0}^{\infty}$ is of exponential order less than $1/b$ where $0 < b < 1$. Suppose that the government sets $g_t$ to minimize the loss function

$$
T = \sum_{t=0}^{\infty} b^t \{ (Y_t - Y_t^*)^2 + d(g_t - g_{t-1})^2 \}, \qquad d \geq 0,
$$

where $\{Y_t^*\}_{t=0}^{\infty}$ is a target sequence of GNPs that is of exponential order less than $1/b$, and $d$ is non-negative and measures the cost of changing the setting of $g$ from its previous value.

&nbsp;&nbsp;&nbsp;&nbsp;**A.** Derive an optimal rule for setting $g_t$ under the assumption that $d = 0$.

&nbsp;&nbsp;&nbsp;&nbsp;**B.** Derive the optimal rule for setting $g_t$ under the assumption that $d > 0$.

**5.** (*Cass-Koopmans optimum growth problem*) A planner wants to choose (consumption, capital) sequences that maximize

$$
\sum_{t=0}^{\infty} \beta^t \left( u_0 c_t - \frac{u_1}{2} c_t^2 \right), \qquad u_0, u_1 > 0
$$

subject to $c_t + k_{t+1} \leq f_0 k_t$, where $k_0$ is given and satisfies $0 < k_0 < (u_0/u_1)(1/(f_0 - 1))$. Here $c_t$ is per capita consumption and $k_t$ is per capita capital. Assume that $f_0 > 1/\beta > 1$, where $\beta$ is the discount factor. The planner imposes the side condition that the $\{k_t\}$ sequence be bounded. (Actually, all we need is that $k_t$ be required to be nonnegative, but boundedness does the job and is easier to handle technically.)

&nbsp;&nbsp;&nbsp;&nbsp;**A.** Find the Euler equation and the transversality condition. Interpret the transversality condition. (*Hint:* use the constraint to eliminate $c_t$ from the objective function and differentiate with respect to successive $k$'s.)

&nbsp;&nbsp;&nbsp;&nbsp;**B.** Solve the Euler equation for the steady state value of $k$, say $\bar{k}$. Find the steady-state value of $\bar{c}$. (*Hint:* you should get

$$
\bar{k} = \left( \frac{1}{f_0 - 1} \right) \frac{u_0}{u_1}, \qquad \bar{c} = \frac{u_0}{u_1}.)
$$

Interpret the steady-state value of $\bar{c}$.

&nbsp;&nbsp;&nbsp;&nbsp;**C.** Show that the optimal feedback rule for setting $k$ is

$$
k_{t+1} = \frac{1}{f_0 \beta} k_t - \left( \frac{1}{f_0 - 1} \right) \left( \frac{u_0}{u_1} \left( \frac{1}{\beta f_0} - 1 \right) \right).
$$ (*)

(*Hint:* Solve the Euler equation.) Does this solution satisfy the transversality condition? If you had solved the "other" root backward, rather than forward, would the transversality condition be satisfied?

&nbsp;&nbsp;&nbsp;&nbsp;**D.** Prove that as $t \to \infty$ the solution for $k_{t+1}$ in (*) converges to $\bar{k}$.

&nbsp;&nbsp;&nbsp;&nbsp;**E.** Prove that consumption increases as capital increases toward its steady-state value.

**6.** A consumer is assumed to face the sequence of budget constraints

$$
A_{t+1} = (1 + r) A_t + (1 + r)(y_t - c_t), \qquad t = 0, 1, 2, \ldots,
$$ (†)

where $A_t$ is his asset holdings at the beginning of period $t$, $y_t$ is exogenous income, and $c_t$ is consumption at $t$. Here $r > 0$ is the real rate of return on assets, assumed constant over time. The consumer earns (pays) each period $1 + r$ times his initial assets plus $1 + r$ times the addition (subtraction) to his assets made by consuming less (more) than his income. Assume that both $c_t$ and $y_t$ are of exponential order less than $1 + r$.

Suppose we impose upon the consumer the boundary condition

$$
\lim_{t \to \infty} (1 + r)^{-t} A_t = 0.
$$

Show that then the sequence of budget constraints (†) implies

$$
A_t + \sum_{t=0}^{\infty} \frac{y_{t+i}}{(1+r)^i} = \sum_{t=0}^{\infty} \frac{c_{t+i}}{(1+r)^i}, \qquad t = 0, 1, 2, \ldots.
$$

Interpret this result.

Assume that the consumer starts out with initial assets of $A_0$ at time 0. Show that the sequence of budget constraints (†) implies

$$
A_t = \sum_{i=0}^{t-1} (1 + r)^{i+1}(y_{t-1-i} - c_{t-1-i}) + (1 + r)^t A_0.
$$

Interpret this result.

**7.** Suppose that portfolio equilibrium is described by Cagan's equation:

$$
m_t - p_t = -1(p_{t+1}^e - p_t), \qquad t = 0, 1, 2, \ldots,
$$

where $m_t$ is the log of the money supply, $p_t$ the log of the price level at $t$, and $p_{t+1}^e$ the public's expectation of $p_{t+1}$, formed at time $t$. Suppose that expectations are "rational," so that

$$
p_{t+1}^e = p_{t+1}.
$$

&nbsp;&nbsp;&nbsp;&nbsp;**A.** Suppose that $\{m_t\}$, $t = 0, 1, \ldots$, is given by

$$
m_t = 10 \lambda^t, \qquad t = 0, 1, 2, \ldots.
$$

Compute the equilibrium value of $p_t$ for $t = 0, 1, 2$, and $t = 5, 6$ for the following values of $\lambda$:

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;(i) $\lambda = 1$;

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;(ii) $\lambda = 1.5$;

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;(iii) $\lambda = 2.0$.

&nbsp;&nbsp;&nbsp;&nbsp;**B.** Suppose that $m_t$ follows the path

$$
m_t = 10 \lambda^t, \qquad t = 0, 1, 2, 3, 4
$$

$$
m_t = 11 \lambda^t, \qquad t = 5, 6, 7, 8, \ldots.
$$

Compute the equilibrium values of $p_t$ for $t = 0, 1, 2, 3, 4, 5, 6$ assuming that $\lambda = 1.5$. How does this time path for $\{p_t\}$ compare with that computed in A(ii)? Graph the two paths.

**8.** (Advertising) A monopolist faces the following demand curve for his product,

$$
p_t = A_0 - A_1 Q_t + g(L) a_t + u_t, \qquad A_0, A_1 > 0,
$$

where $p_t$ is price, $Q_t$ is output, $a_t$ is advertising, $u_t$ is a sequence of shocks to demand, and $g(L) = g_0 + g_1 L + \cdots + g_m L^m$, where $g_j > 0$ for $j = 0, \ldots, m$. The firm maximizes

$$
\sum_{t=0}^{\infty} \beta^t \left\{ p_t Q_t - Q_t s_t - (1/2)[d(L) Q_t]^2 - \frac{\gamma}{2} a_t^2 - a_t w_t \right\}, \qquad 0 < \beta < 1,
$$ (1)

where $d(L) = \sum_{j=0}^{n} d_j L^j$. In (1), $s_t$ is a shock to costs, $\frac{1}{2}[d(L) Q_t]^2$ represents costs of rapid adjustment, and the marginal costs of advertising at $t$ are $(w_t + \gamma a_t)$, where $w_t$ is a known sequence. We assume that $(u_t, s_t, w_t)$ are known sequences of exponential order less than $1/\sqrt{\beta}$. The criterion (1) is to be maximized over sequences for $\{Q_s, a_s, s \geq 0\}$ taking as given $\{Q_s, a_s, s < 0\}$.

&nbsp;&nbsp;&nbsp;&nbsp;**A.** Find the Euler equations for this problem.

&nbsp;&nbsp;&nbsp;&nbsp;**B.** Argue that the solution will be linear laws of motion for $(Q_t, a_t)$ in which each of $(Q_t, a_t)$ depends on lagged values of both $Q$ and $a$, and current and future values of all of $(u, s, w)$.

**9.** (Time to build with two processes)

Consider a monopolist whose output satisfies

$$
Q_t = f(L) n_{1t} + g(L) n_{2t}
$$ (1)

where

$$
f(L) = \sum_{j=0}^{m} f_j L^j, \quad g(L) = \sum_{j=0}^{r} g_j L^j, \qquad f_j \geq 0, \quad g_j \geq 0,
$$

for all $j$. In (1), $n_{1t}$ is the amount of labor at $t$ that is assigned to process 1, while $n_{2t}$ is the amount that is assigned to process 2. The idea is that output can be produced via two processes, with different timing characteristics, e.g., to represent the notion that the first process is fast but wasteful, while the other is efficient but time consuming, we might set $f(L) = L$, $g(L) = (\frac{1}{4})[L + L^2 + L^3 + L^4]$. The firm faces the demand curve

$$
p_t = A_0 - A_1 Q_t + u_t, \qquad A_0, A_1 > 0,
$$ (2)

where $u_t$ is a known sequence of exponential order less than $1/\sqrt{\beta}$. The firm hires labor at the wage rate $w_t$, where $w_t$ is a known sequence of exponential order less than $1/\sqrt{\beta}$. The firm's problem is to maximize

$$
\sum_{t=0}^{\infty} \beta^t \{ p_t Q_t - w_t(n_{1t} + n_{2t}) \}, \qquad 0 < \beta < 1,
$$ (3)

subject to (1) and (2), with $\{n_{1s}, s = -1, \ldots, -m\}$, $\{n_{2s}, s = -1, \ldots, -r\}$ given.

&nbsp;&nbsp;&nbsp;&nbsp;**A.** Find the Euler equations for this problem.

&nbsp;&nbsp;&nbsp;&nbsp;**B.** Indicate the form of the optimum decision rules for $(n_{1t}, n_{2t})$.

**10.** At time $t$, a farmer plants $A_t$ units of land from which he produces $y_{t+1} = f_1 A_t$ units of output available to be sold at time $(t + 1)$ at price $p_{t+1}$ where $f_1 > 0$. If $A_t$ units of land are planted, the farmer's costs at time $t$ are $c_0 A_t + (c_1/2) A_t^2 + c_2 A_t A_{t-1}$, where $c_0, c_1, c_2 > 0$. The term $c_2 A_t A_{t-1}$ reflects the wearing out of the land from two successive heavy plantings. The price of output $p_t$ obeys the first-order difference equation

$$
p_{t+1} = \alpha p_t
$$

with $p_0$ given and $|\alpha| < 1/\sqrt{\beta}$. The farmer faces the $p_t$ sequence as a price-taker. The farmer chooses $A_0, A_1, \ldots$ to maximize

$$
\sum_{t=0}^{\infty} \beta^{t+1} p_{t+1}(f_1 A_t) - \sum_{t=0}^{\infty} \beta^t \left[ c_0 A_t + \frac{c_1}{2} A_t^2 + c_2 A_t A_{t-1} \right], \qquad 0 < \beta < 1,
$$

subject to $A_{-1}$ given.

&nbsp;&nbsp;&nbsp;&nbsp;**A.** Obtain the first-order necessary conditions for optimization.

&nbsp;&nbsp;&nbsp;&nbsp;**B.** Assume that $c_1 = 5/2$, $c_2 = 1$, $\beta = 1$, $f_1 = 10$, and $c_0 = 1$. Find an optimal acreage plan for the farmer of the form

$$
A_t = L(A_{t-1}, p_t, \text{constant})
$$

where $L$ is a linear function. Compute actual numerical values for the linear coefficients under the alternative assumptions $\alpha = 0.5$, $\alpha = 0$, and $\alpha = -0.5$.

&nbsp;&nbsp;&nbsp;&nbsp;**C.** Using your answer to B, describe how you would expect the farmer's supply of output next period, $y_{t+1} = f_1 A_t$, to vary with this period's price, $p_t$. When is supply $y_{t+1}$ to vary inversely with $p_t$?

&nbsp;&nbsp;&nbsp;&nbsp;**D.** Use this example to illustrate the general observation that decision rules change whenever the economic environment changes.

**11.** A social planner wants to maximize the criterion

$$
\sum_{t=0}^{\infty} b^t \{ u_1(c_t - a) - (\tfrac{1}{2}) u_2(c_t - a)^2 - u_3 n_t - (\tfrac{1}{2}) u_4 n_t^2 \},
$$ (1)

$$
0 < b < 1, \qquad u_1, u_2, u_3, u_4 > 0, \quad a > 0,
$$

subject to

$$
\text{(a)} \quad c_t + n_{t+1} + g = f n_t, \qquad g > 0, \quad f > 1, \quad fb > 1,
$$

$$
\text{(b)} \quad n_0 > 0, \quad \text{given}.
$$ (2)

Here $c_t$ is consumption, $n_t$ is labor supplied, $g$ is government purchases, which are fixed and outside the control of the planner, and $a$ is a "bliss level" of consumption. (In the background in (1)–(2), it is understood that output at $t$ equals $f \min(k_t, n_t)$ where $k_t$ is capital at $t$. This fixed proportions technology makes it optimal to set $n_t = k_t$. We have substituted this into an original set of constraints to get (2).)

&nbsp;&nbsp;&nbsp;&nbsp;**A.** Show that the optimal path for $\{n_t\}$ converges as $t \to \infty$. Call the limit point the "stationary optimal value of $n$."

&nbsp;&nbsp;&nbsp;&nbsp;**B.** Give an explicit formula for the stationary optimal value of $n$, call it $\bar{n}$.

&nbsp;&nbsp;&nbsp;&nbsp;**C.** Stationary values of $n_t$ must satisfy (2a) with $n_{t+1} = n_t = n$:

$$
c + n + g = f n.
$$ (3)

Consider the following "Golden rule" problem: to maximize

$$
u_1(c - a) - (\tfrac{1}{2}) u_2(c - a)^2 - u_3 n - (\tfrac{1}{2}) u_4 n^2
$$

subject to (3), by choice of $c$ and $n$. (In words, find the utility maximizing sustainable levels of consumption and labor.) Solve this problem, finding an explicit formula for the "Golden rule" employment level, $n_g$.

&nbsp;&nbsp;&nbsp;&nbsp;**D.** Does $n_g = \bar{n}$, in general? If not, describe special settings of the parameter values for which $n_g = \bar{n}$. Interpret these special settings.

## References

- Allen, R. G. D. (1960). *Mathematical Economics*, 2nd ed., London: Macmillan.
- Baumol, W. J. (1959). *Economic Dynamics*, New York: Macmillan.
- Cagan, P. (1956). "The monetary dynamics of hyperinflation." *Studies in the Quantity Theory of Money* (M. Friedman, ed.), Chicago: University of Chicago Press.
- Eichenbaum, M. S. (1983). "A rational expectations equilibrium model of inventories of finished goods and employment." *Journal of Monetary Economics*, Vol. 29, No. 3, pp. 259–277.
- Gabel, R. A. and Roberts, R. A. (1973). *Signals and Linear Systems*, New York: Wiley.
- Hall, R. E. (1978). "Stochastic implications of the life cycle-permanent income hypothesis: theory and evidence." *Journal of Political Economy*, Vol. 86, No. 6, pp. 971–987, reprinted in *Rational Expectations and Econometric Practice* (R. E. Lucas, Jr. and T. J. Sargent, eds.), Minneapolis: University of Minnesota Press, 1981.
- Hansen, L. P. and Sargent, T. J. (1980). "Formulating and estimating dynamic linear rational expectations models." *Journal of Economic Dynamics and Control*, Vol. 2, No. 1, pp. 7–46, reprinted in *Rational Expectations and Econometric Practice* (R. E. Lucas, Jr. and T. J. Sargent, eds.), Minneapolis: University of Minnesota Press, 1981.
- Hansen, L. P. and Sargent, T. J. (1981). "Linear rational expectations models for dynamically interrelated variables." *Rational Expectations and Econometric Practice*, (R. E. Lucas, Jr. and T. J. Sargent, eds.), Minneapolis: University of Minnesota Press.
- Hansen, L. P., Epple, D., and Roberds, W. (1985). "Linear quadratic models of resource depletion." *Energy, Foresight, and Strategy* (T. J. Sargent, ed.), Washington, DC: Resources for the Future.
- Lucas, R. E. Jr. (1976). "Econometric policy evaluation: a critique." *The Phillips Curve and Labor Markets* (K. Brunner and A. Meltzer, eds.), Vol. 1, Carnegie-Rochester Conference Series on Public Policy, Amsterdam: North-Holland.
- "Math and music: the deeper links." *New York Times*, Sunday, August 29, 1982.
- Muth, J. F. (1960). "Optimal properties of exponentially weighted forecasts." *Journal of the American Statistical Association*, Vol. 55, No. 290, pp. 299–321.
- Muth, J. F. (1961). "Rational expectations and the theory of price movements." *Econometrica*, Vol. 29, No. 3, pp. 315–335.
- Nadiri, M. and Rosen, S. (1973). *A Disequilibrium Model of Demand for Factors of Production*, New York: Columbia University Press, for the NBER.
- Samuelson, P. (1944). "Interaction between the multiplier analysis and the principle of acceleration." *Readings in Business Cycle Theory*, American Economic Association, New York: McGraw-Hill.


