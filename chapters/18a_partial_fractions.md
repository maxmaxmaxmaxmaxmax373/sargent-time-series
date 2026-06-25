# The Residue Theorem Behind Partial Fractions

Several times now we have broken a rational function of the lag operator into a sum of
first-order pieces of the form $\dfrac{A_j}{1-\lambda_j L}$. We did it to invert a moving
average operator, writing $d(L)^{-1} = \sum_j A_j/(1-\lambda_j L)$ in {eq}`eq-70`; we will do
it again, implicitly, every time we apply the annihilation operator $(\cdot)_+$ of the
Wiener–Kolmogorov prediction formula {eq}`eq-62`. This "partial-fractions" maneuver is not a
bag of tricks. Its coefficients are **residues**, and the legitimacy of operating on the
pieces one at a time rests on the residue theorem and the Laurent expansions of
{doc}`04_fourier_z_transforms` and {doc}`05_inverse_z_transform`. This short section makes
the connection explicit and shows why first-order building blocks keep reappearing.

## Partial-fraction coefficients are residues

Let $g(z) = N(z)/D(z)$ be a ratio of polynomials with $\deg N < \deg D$ (if not, divide out
the polynomial part first), and suppose $D$ has distinct zeros $p_1,\dots,p_n$. Then $g$ has a
partial-fraction expansion

```{math}
:label: eq-pf-1
g(z) = \sum_{j=1}^{n} \frac{r_j}{z - p_j},
\qquad
r_j = \operatorname*{Res}_{z=p_j} g(z) = \lim_{z\to p_j}(z-p_j)\,g(z) = \frac{N(p_j)}{D'(p_j)}.
```

The coefficient $r_j$ is exactly the residue of $g$ at the simple pole $p_j$ defined in
{eq}`eq-24`. The "cover-up rule" of elementary calculus — multiply by $(z-p_j)$, set
$z=p_j$ — *is* the residue evaluation. The residue theorem supplies the reason the two
agree: integrating {eq}`eq-pf-1` around a small circle $\gamma_j$ enclosing only $p_j$,

```{math}
:label: eq-pf-2
\frac{1}{2\pi i}\oint_{\gamma_j} g(z)\,dz = r_j,
```

since each term $\frac{r_k}{z-p_k}$ integrates to $r_k$ if $p_k$ is inside $\gamma_j$ and to
$0$ otherwise. The residue is at once the coefficient in the expansion and the contour
integral; that is what lets us pass between the two views.

For time-series work it is more convenient to factor the denominator into the
lag-operator-friendly form $D(z) = \prod_{j}(1-\lambda_j z)$, whose poles sit at
$z = 1/\lambda_j$. The expansion becomes

```{math}
:label: eq-pf-3
g(z) = \frac{N(z)}{\prod_{j=1}^{n}(1-\lambda_j z)} = \sum_{j=1}^{n} \frac{A_j}{1-\lambda_j z},
\qquad
A_j = \Big[(1-\lambda_j z)\,g(z)\Big]_{z=1/\lambda_j}
    = \frac{N(1/\lambda_j)}{\prod_{k\neq j}\big(1-\lambda_k/\lambda_j\big)} .
```

This $A_j$ is the same residue in different clothing: because $1-\lambda_j z = -\lambda_j(z -
1/\lambda_j)$, the two coefficients are related by $A_j = -\lambda_j\,
\operatorname*{Res}_{z=1/\lambda_j} g(z)$. The cover-up rule {eq}`eq-pf-3` is the form we use
in practice.

**Repeated poles.** When $D$ has a zero of order $m$, say a factor $(1-\lambda z)^m$, the
expansion carries terms $\dfrac{A^{(1)}}{1-\lambda z} + \cdots + \dfrac{A^{(m)}}{(1-\lambda
z)^m}$, and the coefficients are the higher-order residues of {eq}`eq-23`, obtained by
differentiating $(1-\lambda z)^m g(z)$. The order-$m$ residue formula of
{doc}`05_inverse_z_transform` is precisely what a repeated pole demands.

## Why this is the natural tool: inversion pole by pole

The payoff is that each simple term is the $z$-transform of a *geometric* sequence. For
$|\lambda_j|<1$,

```{math}
:label: eq-pf-4
\frac{A_j}{1-\lambda_j z} = A_j\sum_{k=0}^{\infty}\lambda_j^{k} z^{k}
\quad\Longleftrightarrow\quad
\{A_j\lambda_j^{k}\}_{k\ge 0},
```

one-sided in nonnegative powers — exactly Example 1 of {doc}`05_inverse_z_transform`,
computed there by residues. A partial-fraction expansion therefore turns the inversion of a
complicated rational function into the inversion of a handful of geometrics, one per pole,
which are then added. This is the same calculation as the inverse-$z$-transform-by-residues
of {eq}`eq-25`, simply organized by pole rather than performed as a single contour integral.

Two facts we have already used are instances of it.

- Inverting an $m$th-order moving average. The expansion $d(L)^{-1}=\sum_j A_j/(1-\lambda_j
  L)$ of {eq}`eq-70` gives the autoregressive weights as a sum of geometrics, $a_k = \sum_j
  A_j\lambda_j^{k}$. Convergence of {eq}`eq-71` is just the requirement $|\lambda_j|\le 1$
  that makes each geometric {eq}`eq-pf-4` summable — the fundamental-versus-nonfundamental
  distinction read off the poles.

- The autocovariances of an autoregression. The closed form {eq}`eq-27` for $c_y(\tau)$ —
  a sum over the poles $\lambda_j$ of residue terms — is nothing but the partial-fraction
  expansion of the covariance generating function $g_y(z)$ inverted term by term. The "sum of
  residues inside $\Gamma$" of {eq}`eq-26` and the partial-fraction sum are the same object.

## Example 1: inverting a moving average operator

Let $d(L) = (1-\tfrac12 L)(1-\tfrac13 L)$, so that $\lambda_1=\tfrac12,\ \lambda_2=\tfrac13$,
and we want $d(L)^{-1} = \sum_j A_j/(1-\lambda_j L)$. The cover-up rule {eq}`eq-pf-3` with
$N\equiv 1$ gives

$$
A_1 = \frac{1}{1-\lambda_2/\lambda_1} = \frac{1}{1-\tfrac{2}{3}} = 3,
\qquad
A_2 = \frac{1}{1-\lambda_1/\lambda_2} = \frac{1}{1-\tfrac{3}{2}} = -2 ,
$$

so that

$$
d(L)^{-1} = \frac{3}{1-\tfrac12 L} - \frac{2}{1-\tfrac13 L},
\qquad
a_k = 3\left(\tfrac12\right)^{k} - 2\left(\tfrac13\right)^{k}.
$$

As a check, $a_0 = 3-2 = 1$ (consistent with $d_0=1$), and $a_1 = \tfrac32 - \tfrac23 =
\tfrac56$, which matches the coefficient of $L$ in the direct inverse of $d(L) = 1 - \tfrac56
L + \tfrac16 L^2$. The weights decay at the rate of the *largest* pole reciprocal,
$\lambda_1=\tfrac12$, because that geometric dominates the tail — a fact that the
partial-fraction form displays at a glance but the recursive inversion hides.

If instead the pole is repeated, $d(L) = (1-\tfrac12 L)^2$, the order-two residue
{eq}`eq-23` gives $d(L)^{-1} = \sum_{k\ge 0}(k+1)\left(\tfrac12\right)^{k} L^{k}$: a geometric
multiplied by the polynomial $k+1$, the signature of a double pole.

## Example 2: the annihilation operator is diagonal in this basis

The prediction formula {eq}`eq-62` hinges on the annihilation operator $(\cdot)_+$, which
discards negative powers of $L$. Partial fractions is what makes it computable, because
$(\cdot)_+$ acts on each first-order piece independently. Consider the mixed process of
{doc}`15_prediction_examples`,

$$
x_t = \frac{1+aL}{1-\beta L}\,\epsilon_t, \qquad |a|<1,\ |\beta|<1,
$$

with $\epsilon_t$ fundamental. To form $P_{t-1}x_t = \big(\tfrac{d(L)}{L}\big)_+ d(L)^{-1}
x_{t-1}$ we expand the inside operator in partial fractions before annihilating:

$$
\frac{L^{-1}(1+aL)}{1-\beta L}
= \underbrace{\frac{L^{-1}}{1-\beta L}}_{\text{has an }L^{-1}\text{ term}}
  + \frac{a}{1-\beta L}.
$$

The second piece is already a sum of nonnegative powers and survives untouched; the first
contributes its negative power $L^{-1}$ (annihilated) and leaves the geometric tail
$\beta/(1-\beta L)$. Hence

$$
\left(\frac{L^{-1}(1+aL)}{1-\beta L}\right)_+ = \frac{\beta + a}{1-\beta L},
\qquad
P_{t-1}x_t = \frac{\beta+a}{1-\beta L}\cdot\frac{1-\beta L}{1+aL}\,x_{t-1}
           = \frac{a+\beta}{1+aL}\,x_{t-1},
$$

reproducing the result of {doc}`15_prediction_examples`. The lesson is structural:
**$(\cdot)_+$, like inversion and like filtering, is diagonal in the partial-fraction
basis** — it can be applied to one pole at a time. That is the real reason the first-order
operator $1/(1-\lambda L)$ is the recurring atom of this theory.

## Echoes and influences

This circle of ideas — *factor, expand into one term per pole, operate pole by pole, add* —
runs through the neighboring sections:

- **The inverse $z$-transform ({doc}`05_inverse_z_transform`).** Partial fractions is the
  bookkeeping that organizes the residue formula {eq}`eq-25`: each pole contributes one
  geometric, and the inversion is their sum.

- **Finding Wold representations ({doc}`17_wold_ma`, {doc}`18_wold_arma`).** Factoring the
  covariance generating function locates the poles and zeros; partial fractions then inverts
  the chosen $d(L)$, and the position of each pole relative to the unit circle decides
  whether the corresponding geometric is one-sided in the past (fundamental) or the future.

- **Prediction ({doc}`14_linear_prediction`, {doc}`15_prediction_examples`).** The
  Wiener–Kolmogorov operator $(\cdot)_+$ keeps the nonnegative-power part of each
  partial-fraction term, as in Example 2; every worked forecast there is a pole-by-pole
  calculation.

- **Signal extraction ({doc}`19_signal_extraction`).** The whole apparatus — common
  denominator, factor the numerator, invert — is partial fractions applied to a rational
  spectral density; the Muth example is the first-order case.

A unifying way to see it: partial fractions **diagonalizes** a rational lag operator into a
sum of one-pole operators, in the same way that an eigen-decomposition diagonalizes a matrix.
Operations that look formidable on the original operator — inverting it, predicting with it,
filtering through it, annihilating its negative powers — become trivial on each $1/(1-\lambda
L)$ piece and are simply summed back. The residues are the weights in that decomposition, and
the residue theorem is the guarantee that the pieces fit back together.

## Exercises

```{admonition} Exercise 1: residues, inversion, and a repeated pole
:class: tip

1. Let $d(L) = (1-0.6L)(1-0.2L)$. Use the cover-up rule {eq}`eq-pf-3` to find $A_1,A_2$ in
   $d(L)^{-1} = A_1/(1-0.6L) + A_2/(1-0.2L)$, and write the autoregressive weights $a_k$.
   Verify that $a_0 = 1$ and that $a_1$ equals the coefficient of $L$ in the direct
   expansion of $d(L)^{-1}$. Which pole governs the rate at which $a_k\to 0$?
2. Confirm the residue/cover-up coefficient $A_j$ equals $-\lambda_j$ times the ordinary
   residue $\operatorname*{Res}_{z=1/\lambda_j} d(z)^{-1}$ of {eq}`eq-24`, for this $d$.
3. Now take the repeated pole $d(L) = (1-0.6L)^2$. Using the order-two residue {eq}`eq-23`,
   show $d(L)^{-1} = \sum_{k\ge0}(k+1)(0.6)^k L^k$, and explain in one line why a double pole
   produces a geometric multiplied by a polynomial in $k$.
```

```{admonition} Exercise 2: annihilation pole by pole
:class: tip

Let $x_t = \dfrac{1}{(1-\tfrac12 L)(1-\tfrac13 L)}\,\epsilon_t$ with $\epsilon_t$ fundamental
white noise.

1. Expand $\dfrac{1}{(1-\tfrac12 L)(1-\tfrac13 L)}$ in partial fractions (you computed the
   coefficients in Example 1) and write the moving-average weights $d_k$.
2. Use the Wiener–Kolmogorov formula {eq}`eq-62` with $k=1$,
   $P_{t-1}x_t = \big(\tfrac{d(L)}{L}\big)_+\, d(L)^{-1}\,x_{t-1}$, applying $(\cdot)_+$ to
   each partial-fraction term separately, to show that the one-step-ahead forecast is a
   geometric distributed lag of past $x$'s. Identify its decay rate.
3. *Optional.* Compute the autocovariances $c_x(\tau)$ from the residue formula {eq}`eq-27`
   and confirm that they are a sum of two geometrics in $|\tau|$, one per pole — the
   covariance-generating-function counterpart of the partial-fraction expansion in part 1.
```

## References

- Churchill, R. V., J. W. Brown, and R. F. Verhey (1974). *Complex Variables and
  Applications*, 3rd ed. McGraw-Hill.
- Hansen, L. P., and T. J. Sargent (1980). Formulating and estimating dynamic linear rational
  expectations models. *Journal of Economic Dynamics and Control* 2, 7–46.
- Whittle, P. (1983). *Prediction and Regulation by Linear Least-Square Methods*, 2nd ed.
  University of Minnesota Press.
