# A Mathematical Digression: Fourier and $z$ Transforms

```{note}
This section provides the mathematical foundations for the $z$-transform and Fourier
methods used throughout the chapter. It **may be omitted on a first reading** without
loss of continuity.
```

The following theorem provides the foundation for the $z$-transform, Fourier transform,
and lag operator methods used repeatedly in these pages.

---

**Theorem (Riesz-Fischer).**[^fn-rf-1]
Let $\{c_n\}_{n=-\infty}^{\infty}$ be a sequence of complex numbers for which
$\sum_{n=-\infty}^{\infty} |c_n|^2 < \infty$.
Then there exists a complex-valued function $f(\omega)$ defined for
$\omega \in [-\pi, \pi]$ such that

```{math}
:label: eq-19
f(\omega) = \sum_{j=-\infty}^{\infty} c_j\, e^{-i\omega j},
```

where the series converges in mean square:

$$
\lim_{n\to\infty} \int_{-\pi}^{\pi}
  \left|\sum_{j=-n}^{n} c_j e^{-i\omega j} - f(\omega)\right|^2 d\omega = 0.
$$

The function $f(\omega)$ is called the **Fourier transform of the sequence $\{c_k\}$**
and satisfies $\int_{-\pi}^{\pi} |f(\omega)|^2 d\omega < \infty$
(i.e., $f$ belongs to $L_2[-\pi,\pi]$, a Lebesgue integral).

Given $f(\omega)$, the sequence $\{c_k\}$ can be recovered via the **inversion formula**:

```{math}
:label: eq-20
c_k = \frac{1}{2\pi}\int_{-\pi}^{\pi} f(\omega)\, e^{+i\omega k}\, d\omega.
```

Finally, $f(\omega)$ and $\{c_k\}$ satisfy **Parseval's relation**:

$$
\frac{1}{2\pi}\int_{-\pi}^{\pi} |f(\omega)|^2\, d\omega
= \sum_{j=-\infty}^{\infty} |c_j|^2.
$$

---

[^fn-rf-1]: For a proof, see Apostol (1974, Chapter 11).

## The $l_2$ and $L_2$ Spaces

Consider the space of all doubly infinite square-summable sequences
$\{x_k\}_{k=-\infty}^{\infty}$ with $\sum_{k=-\infty}^{\infty} |x_k|^2 < \infty$.
We denote this space $l_2(-\infty,\infty)$.[^fn-rf-2]
It is a *linear space*:

1. If $\{x_t\} \in l_2$ and $\alpha$ is a scalar, then
   $\{\alpha x_k\} \in l_2$ (i.e., $\sum|\alpha x_k|^2 < \infty$).
2. If $\{x_t\}, \{y_t\} \in l_2$, then $\{x_k + y_k\} \in l_2$.

[^fn-rf-2]: For properties of linear spaces, see Naylor and Sell (1971).

Similarly, $L_2[-\pi,\pi]$ is the space of functions $f(\omega)$ with
$\int_{-\pi}^{\pi} |f(\omega)|^2 d\omega < \infty$, which is also a linear space.

Both spaces are **metric spaces**. On $l_2$ the metric is

$$
d_2(x, y) = \left[\sum_{k=-\infty}^{\infty} |x_k - y_k|^2\right]^{1/2},
$$

and on $L_2[-\pi,\pi]$ it is[^fn-rf-3]

$$
D_2(f, g) = \left\{\frac{1}{2\pi}\int_{-\pi}^{\pi} |f(\omega) - g(\omega)|^2\, d\omega\right\}^{1/2}.
$$

[^fn-rf-3]: We adopt the convention that $f = g$ whenever they differ only on a set of
Lebesgue measure zero. See Naylor and Sell (1971).

## The Fourier Transform as an Isometric Isomorphism

The Fourier transform {eq}`eq-19` defines a mapping from $l_2(-\infty,\infty)$ to
$L_2[-\pi,\pi]$. Let $\{x_k\}, \{y_k\} \in l_2$ with Fourier transforms
$x(\omega) = \sum_k x_k e^{-i\omega k}$ and $y(\omega) = \sum_k y_k e^{-i\omega k}$.
Then:

$$
x(\omega) + y(\omega) = \sum_k (x_k + y_k)\, e^{-i\omega k}, \qquad
\alpha\, x(\omega) = \sum_k \alpha x_k\, e^{-i\omega k},
$$

so the Fourier transform is an *isomorphism* (it preserves linear structure). Moreover,
by Parseval's relation,

$$
D_2(x(\omega), y(\omega)) = d_2(x, y),
$$

so it is also *isometric* (it preserves distance). The converse of the Riesz-Fischer
theorem assures the mapping is *onto* and one-to-one, giving us an
**isometric isomorphism** from $l_2(-\infty,\infty)$ to $L_2[-\pi,\pi]$.

## Convolution Becomes Multiplication

The Fourier transformation {eq}`eq-19` puts square-summable sequences $\{x_k\}$ into
one-to-one correspondence with square-integrable functions $f(\omega)$ on $[-\pi,\pi]$.
The transformation preserves linear structures and a measure of distance, as we have seen.
The benefit from using the transformation is that operations that are complicated in one
space are sometimes the counterparts of simple operations in another space. In particular,
consider the *convolution* of two sequences $\{x_k\}$ and $\{y_k\}$, defined to be
the new sequence

$$
\{y * x_k\}_{k=-\infty}^{\infty} \equiv \left\{\sum_{s=-\infty}^{\infty} y_s x_{k-s}\right\}_{k=-\infty}^{\infty}.
$$

The Fourier transform of $(y * x_k)_k$ is given by

$$
\sum_{k=-\infty}^{\infty} \sum_{s=-\infty}^{\infty} y_s x_{k-s} e^{-i\omega k}
= \sum_{s=-\infty}^{\infty} y_s e^{-i\omega s}
  \sum_{k=-\infty}^{\infty} x_{k-s} e^{-i\omega(k-s)}
= y(\omega)\, x(\omega),
$$

where $y(\omega) = \sum_k y_k e^{-i\omega k}$ and $x(\omega) = \sum_k x_k e^{-i\omega k}$.

Thus: **the Fourier transform of the convolution of $\{x_k\}$ and $\{y_k\}$ is the
product of their Fourier transforms.** The complicated convolution operation corresponds
simply to multiplication of Fourier transforms.

All transform techniques exploit properties like the preceding one. The aim is to
transform a problem from one space where it appears complicated to another isometrically
isomorphic space where the operations are simpler, then to transform back to the original
space using the inversion mapping such as {eq}`eq-20` after the calculations have been
performed.

## The $z$-Transform Corollary

Making the change of variable $z = e^{i\omega}$ in the Riesz-Fischer theorem yields:

---

**Corollary.**
Let $\{c_n\}_{n=-\infty}^{\infty}$ be a complex sequence with
$\sum_{n=-\infty}^{\infty} |c_n|^2 < \infty$.
Then there exists a complex-valued function $g(z)$ such that

$$
g(z) = \sum_{j=-\infty}^{\infty} c_j\, z^j,
$$

converging in mean square over the unit circle $\Gamma$:

$$
\lim_{n\to\infty} \int_\Gamma
  \left|\sum_{j=-n}^n c_j z^j - g(z)\right|^2 \frac{dz}{z} = 0.
$$

The function $g(z)$—called the **$z$-transform** of $\{c_k\}$—is defined at least on
the unit circle and satisfies
$\left|(2\pi i)^{-1}\int_\Gamma |g(z)|^2 dz/z\right| < \infty$.
The coefficients can be recovered by

$$
c_k = \frac{1}{2\pi i} \int_\Gamma g(z)\, z^{-k-1}\, dz.
$$

---

## Orthogonality of the Exponential Basis

The functions $e^{i\omega j}$, $j = 0, \pm 1, \pm 2, \ldots$ are orthogonal on
$[-\pi, \pi]$. For $n \neq m$:

$$
\frac{1}{2\pi}\int_{-\pi}^{\pi} e^{i\omega n}\, e^{-i\omega m}\, d\omega
= \frac{1}{\pi(n-m)}\sin\pi(n-m) = 0,
$$

since $n - m$ is a nonzero integer. This orthogonality underlies the inversion formula
{eq}`eq-20` and all of the spectral analysis that follows.

```{note}
For the most part, the Riesz-Fischer theorem and its corollary are sufficient for our
needs. We will briefly encounter a deterministic process for which
$\sum|c_k|^2 < \infty$ is violated (where the $c_k$ depict the covariogram), so
that the theorem does not suffice to define the Fourier transform. In that case the
Fourier transform is still defined in the sense of a *generalized function* (distribution).
```
