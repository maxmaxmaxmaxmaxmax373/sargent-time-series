# Evaluating the Inverse $z$ Transform

Given a square-summable sequence $\{c_j\}$, the $z$-transform is

$$
g(z) = \sum_{j=-\infty}^{\infty} c_j\, z^j,
$$

well-defined at least on the unit circle $z = e^{i\omega}$, $\omega \in [0, 2\pi]$.
The sequence $\{c_k\}$ can be recovered by the inversion formula

```{math}
:label: eq-21
c_k = \frac{1}{2\pi i} \int_\Gamma g(z)\, z^{-k-1}\, dz,
\qquad k = 0, \pm 1, \pm 2, \ldots
```

where $\Gamma$ denotes the unit circle and the integral is a contour integral.[^fn-izl-1]
In this section we give simple formulas for evaluating {eq}`eq-21` using **residues**.
Virtually no knowledge of complex analysis is required to use these formulas.

[^fn-izl-1]: Churchill, Brown, and Verhey (1974) is a good reference on the complex
analysis used here. Gabel and Roberts (1973) is a good reference on operational methods
for linear dynamic systems.

## Poles and Residues

A **pole** of $g(z)$ is a point $z_0$ in the complex plane where $g(z) \to \infty$
as $z \to z_0$. In this book, $g(z)$ is almost always a *rational function*—a ratio
of finite-order polynomials in $z$—so the poles are simply the zeros of the denominator.

*Test for poles of order $m$:* If for some positive integer $m$ the function

```{math}
:label: eq-22
\phi(z) = (z - z_0)^m\, g(z)
```

can be defined[^fn-izl-2] with $\phi(z_0) \neq 0$, then $g(z)$ has a pole of order
$m$ at $z_0$.

[^fn-izl-2]: The function $\phi(z)$ must also be "analytic" at $z_0$, meaning its
derivative exists at $z_0$ and nearby. In our examples this requirement is routinely
satisfied.

For example, if $g(z) = 1/(1-\lambda z)^r$, then $g(z)$ has a pole of order $r$ at
$z = \lambda^{-1}$.

*Definition of residue:* Suppose $g(z)$ has a pole of order $m$ at $z = z_0$, and
define $\phi(z) = (z-z_0)^m g(z)$. The **residue** at $z_0$ is

```{math}
:label: eq-23
\operatorname{res}(z_0) = \frac{\phi^{(m-1)}(z_0)}{(m-1)!},
```

where $\phi^{(m-1)}$ denotes the $(m-1)$th derivative. When $m = 1$ (a simple pole),
this reduces to

```{math}
:label: eq-24
\operatorname{res}(z_0) = \lim_{z \to z_0} (z - z_0)\, g(z).
```

## The Residue Formulas

The inversion integral {eq}`eq-21` can be evaluated by either of two equivalent
formulas:[^fn-izl-3]

```{math}
:label: eq-25
c_j = \frac{1}{2\pi i}\int_\Gamma g(z)\, z^{-j-1}\, dz
    = \begin{cases}
        \text{sum of residues of } g(z^{-1})z^{-j-1}
        \text{ at poles inside unit circle} \\[4pt]
        \text{sum of residues of } g(z)z^{-j-1}
        \text{ at poles inside unit circle.}
      \end{cases}
```

Choose whichever branch avoids poles of order greater than one at $z = 0$.

[^fn-izl-3]: The second representation is standard (e.g., Churchill, Brown, and Verhey
1974). The first follows by noting that $g(z^{-1}) = \sum_j c_{-j} z^j \equiv \sum_j d_j z^j$
with $d_j = c_{-j}$; applying the second representation to compute $d_j$ then gives
$c_j = \text{sum of residues of } g(z^{-1})z^{-j-1}$ inside the unit circle.

## Example 1: $g(z) = 1/(1-\lambda z)$, $|\lambda| < 1$

Using the first branch of {eq}`eq-25`, the function
$g(z^{-1})z^{-j-1} = z^{j-1}/(1-\lambda z^{-1})$ has:

- For $j > 0$: a single simple pole at $z = \lambda$ (inside the unit circle), with
  residue $\lambda^j$. So $c_j = \lambda^j$.
- For $j < 0$: using the second branch, $z^{-j-1}/(1-\lambda z)$ has no poles inside
  the unit circle, so $c_j = 0$.
- For $j = 0$: using the second branch, $z^{-1}/(1-\lambda z)$ has a simple pole at
  $z = 0$ with residue $1$, so $c_0 = 1$.

Thus,

$$
c_j = \begin{cases} \lambda^j & j \geq 0, \\ 0 & j < 0. \end{cases}
$$

(These results are more easily obtained by expanding $1/(1-\lambda z)$ as a geometric
series, but residue calculations are faster for more complex examples.)

## Example 2: The $n$th-Order AR Covariance Generating Function

Consider the covariance generating function

$$
g_y(z) = \frac{1}{\prod_{j=1}^n(1-\lambda_j z)\,\prod_{j=1}^n(1-\lambda_j z^{-1})}
= \frac{z^n}{\prod_{j=1}^n(1-\lambda_j z)\,\prod_{j=1}^n(z-\lambda_j)},
$$

where $|\lambda_j| < 1$ for $j = 1, \ldots, n$. Applying {eq}`eq-25`,

```{math}
:label: eq-26
c_y(\tau) =
\begin{cases}
  \text{sum of residues of }
  \dfrac{z^{n-\tau-1}}{\prod_{j=1}^n(1-\lambda_j z)\,\prod_{j=1}^n(z-\lambda_j)}
  \text{ inside } \Gamma, & \tau \leq 0, \\[10pt]
  \text{sum of residues of }
  \dfrac{z^{n+\tau-1}}{\prod_{j=1}^n(1-\lambda_j z)\,\prod_{j=1}^n(z-\lambda_j)}
  \text{ inside } \Gamma, & \tau \geq 0.
\end{cases}
```

Use the first line for $\tau \leq 0$ and the second for $\tau \geq 0$ to avoid
higher-order poles at $z = 0$. In each case, the only poles inside the unit circle are
the simple poles at $\lambda_1, \ldots, \lambda_n$. For $\tau \geq 0$ the residue at
$\lambda_j$ is

$$
\operatorname{res}(\lambda_j) = \frac{\lambda_j^{n+\tau-1}}
  {\prod_{k=1}^n(1-\lambda_k\lambda_j)\,\prod_{k=1,\, k\neq j}^n(\lambda_j-\lambda_k)},
\qquad \tau \geq 0.
$$

Summing over $j$ and using the symmetry $c_y(\tau) = c_y(-\tau)$:

```{math}
:label: eq-27
c_y(\tau) = \sum_{j=1}^{n}
  \frac{\lambda_j^{n+|\tau|-1}}
  {\prod_{k=1}^{n}(1-\lambda_k\lambda_j)\,\prod_{k=1,\, k\neq j}^n(\lambda_j-\lambda_k)}.
```
