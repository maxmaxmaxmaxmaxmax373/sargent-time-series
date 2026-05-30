# The Effects of Filtering on One-Sided and Two-Sided Projections

Consider a jointly covariance stationary, linearly indeterministic process $(y_t,x_t)$. Consider the projections

```{math}
:label: eq-141
y_t = \sum_{j=-\infty}^{\infty} b_j x_{t-j} + \epsilon_t, \qquad E\epsilon_t x_{t-j} = 0 \quad \text{for all integer } j
```

```{math}
:label: eq-142
y_t = \sum_{j=0}^{\infty} h_j x_{t-j} + v_t, \qquad Ev_t x_{t-j} = 0 \quad \text{for all integer } j \geq 0
```

Assume that $x_t$ has Wold representation $x_t=d(L)\eta_t$, where $E\eta_t^2 = 1$. We have seen that the generating functions $b(L)$ and $h(L)$ are given by the formulas

```{math}
:label: eq-143
b(z) = \frac{g_{yx}(z)}{d(z^{-1})d(z)}
```

```{math}
:label: eq-144
h(z) = \left[\frac{g_{yx}(z)}{d(z^{-1})}\right]_+\frac{1}{d(z)}
```

where $x_t=d(L)\eta_t$ is a Wold representation for $x_t$. Now consider the filtered series

$$
y_t^f=f(L)y_t, \qquad x_t^f=f(L)x_t
$$

where

$$
f(L) = \sum_{j=-\infty}^\infty f_jL^j, \qquad \sum_{j=-\infty}^\infty f_j^2<+\infty.
$$

We have that the covariance generating functions $g_{y^fx^f}(z),g_{x^f}(z),g_{y^f}(z)$ are given by

```{math}
:label: eq-145
\begin{aligned}
g_{y^f}(z) &= f(z)f(z^{-1})g_{y}(z) \\
g_{x^f}(z) &= f(z)f(z^{-1})g_{x}(z) \\
g_{y^fx^f}(z) &= f(z)f(z^{-1})g_{yx}(z)
\end{aligned}
```

Using these formulas, we consider the projections analogous to {eq}`eq-141` and {eq}`eq-142` with the filtered series, $y_t^f,x_t^f$:

```{math}
:label: eq-146
y_t^f = \sum_{j=-\infty}^{\infty} b^f_j x^f_{t-j} + \tilde{\epsilon}_t, \qquad E\tilde{\epsilon}_t x^f_{t-j} = 0 \quad \text{for all integer } j
```

```{math}
:label: eq-147
y_t^f = \sum_{j=0}^{\infty} h^f_j x^f_{t-j} + \tilde{v}_t, \qquad E\tilde{v}_t x^f_{t-j} = 0 \quad \text{for all integer } j \geq 0
```

Let $b^f(L) =\sum_{j=-\infty}^{\infty}b^f_jL^j$, $h^f(L) =\sum_{j=0}^{\infty}h^f_jL^j$. To implement the filtering formula {eq}`eq-144`, we require a Wold representation for the filtered series $x_t^f$. We can obtain such a representation

$$
x_t^f=c(L)a_t
$$

where $c(z)$ has no zeros inside the unit circle and $Ea_t^2=1$; $c(z)$ satisfies

```{math}
:label: eq-148
c(z)c(z^{-1})=g_{x^f}(z) = f(z)f(z^{-1})d(z)d(z^{-1}).
```

However, in general, $f(z)d(z) \neq c(z)$. Now applying the two-sided filtering formula {eq}`eq-143` gives

```{math}
:label: eq-149
b^f(z)=\frac{f(z)f(z^{-1})g_{yx}(z)}{f(z)f(z^{-1})g_{x}(z)}
```

or

$$
b^f(z)=\frac{g_{yx}(z)}{d(z)d(z^{-1})} = b(z)
$$

Applying the one-sided filtering formula {eq}`eq-144` gives

```{math}
:label: eq-150
h^f(z)=\left[\frac{f(z)f(z^{-1})g_{yx}(z)}{c(z^{-1})}\right]_+\frac{1}{c(z)}.
```

Comparison of {eq}`eq-143` and {eq}`eq-149` shows that filtering both series by a common filter leaves unaltered the two-sided infinite projection of $y$ on $x$. Comparison of {eq}`eq-144` and {eq}`eq-150` shows that in general, filtering $y$ and $x$ with a common filter alters the projection of $y$ on the one-sided infinite past of $x$. However, there is an important special case in which the one-sided projection is unaltered by filtering both $y$ and $x$. This special case occurs when the one-sided projection equals the two-sided projection, which we know to be unaltered by filtering. Sims's theorem establishes that this case obtains when the $y$ process fails to Granger cause $x$. We thus have established the following proposition:

If $y$ fails to Granger cause $x$, filtering fails to alter the projection of $y$ on $x$.

A converse of this proposition also holds:

Consider the one-sided projection of $y_t$ on $[x_t,x_{t-1},\ldots]$,

```{math}
:label: eq-151
y_t=\sum_{j=0}^\infty h_jx_{t-j} + v_t = h(L)x_t + v_t
```

where $Ev_tx_{t-j}=0$ for $j\geq0$. Let $F$ be the space of one-sided, square summable sequences with elements $f=\{f_j\}_{j=0}^\infty$, $\sum_{j=0}^\infty f_j^2 < +\infty$. Assume that for every $f \in F$, the one-sided projection of $f(L)y_t$ on $[f(L)x_t,f(L)x_{t-1},\ldots]$ can be represented as

```{math}
:label: eq-152
f(L)y_t=h(L)[f(L)x_t] + f(L)v_t
```

where $E[f(L)v_t]\cdot[f(L)x_{t-k}] = 0$ for $k\geq0$, and where $h(L)$ is the same operator $h(L)$ that appears in {eq}`eq-151`. Then $\{y_t\}$ fails to be Granger-caused by $\{x_t\}$.

To prove the proposition, we begin by evaluating $E[f(L)v_t]\cdot[f(L)x_{t-k}]$, which is assumed to be zero for all $j\geq0$ and all $f\in F$. For $j=0$, we obtain

$$
E(f_0v_t + f_1v_{t-1} + f_2v_{t-2} + \cdots)(f_0x_t + f_1x_{t-1}+ \cdots) = \sum_{h=0}^\infty f_h \sum_{k=h+1}^\infty f_kEx_{t-h}v_{t-k}.
$$

Evidently, unless $Ex_t v_{t-k}=0$ for all $k$ (that is, $k>0$ as well as $k\leq 0$), $f$ can be chosen to make the above expression unequal to zero. Therefore, the condition that the one-sided projection of $y_t^f \equiv f(L)y_t$ on $x_t^f \equiv f(L)x_t$ have the same generating function $h(L)$ for all $f \in F$ implies that $v_t$ is orthogonal to all *future* $x_t$'s $\{x_{t+k}, k\geq 0\}$, as well as current and past $x_t$'s. Therefore, the one-sided projection of $\{y_t\}$ on current and past $x_t$ equals the two-sided projection of $\{y_t\}$ on past, present, and future $\{x_t\}$. It follows from Sims's theorem 2 that $\{y_t\}$ fails to Granger-cause $\{x_t\}$.

These propositions provide the foundation of the result in econometrics that generalized least squares estimates of projection equations like {eq}`eq-142` are consistent if and only if $x_t$ fails to be Granger-caused by $y_t$. Under standard regularity conditions, generalized least squares estimates of $h(z)$ converge to $h^f(z)$ given by {eq}`eq-150`, where $f(z)$ is a filter that is one-sided in nonnegative powers of $z$ and that whitens $v_t$.
