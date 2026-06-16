# Nonlinear (Volterra) Moving-Average Representations

```{note}
This chapter is a modern postscript to the linear theory developed in
Sections 13–18. Everything in Chapter XI rests on the
**Wold representation** $x_t=d(L)\epsilon_t$: the statement that a covariance
stationary, linearly indeterministic process is a *linear* moving average of its
own one-step prediction errors. We now ask what happens when the map from the
innovations $\{\epsilon_t\}$ to $x_t$ is allowed to be **nonlinear**, and we
develop the corresponding moving-average theory — the Volterra and Wiener–Itô
expansions, their frequency-domain counterparts (the polyspectra), the canonical
model classes, and the conditions under which the representation exists. The
companion chapter {doc}`40_nonlinear_lab` carries out the constructions
numerically.
```

(nl-intro)=
## 1. From the Linear Wold Theorem to Nonlinear Functionals

The spine of this chapter has been the linear moving-average representation
{eq}`eq-59`,

$$
x_t = d(L)\epsilon_t = \sum_{j=0}^\infty d_j\,\epsilon_{t-j},
\qquad \epsilon_t = x_t - P[x_t\mid x_{t-1},x_{t-2},\ldots],
$$

in which $x_t$ is a *linear* functional of present and past innovations. The Wold
theorem guarantees that such a representation always exists for the class of
processes we have studied, and the whole apparatus — spectra, prediction,
filtering, Granger causality — is built on it.

But the Wold innovations $\epsilon_t$ are only **linearly** unpredictable: by
construction $\epsilon_t\perp\{x_{t-1},x_{t-2},\ldots\}$ in the sense of *linear*
projection. They need not be independent, and $x_t$ need not be a linear function
of them. The simplest illustration is a process that is white (zero
autocovariance at every nonzero lag) yet strongly dependent: $x_t=\epsilon_t
\epsilon_{t-1}$ with $\epsilon$ i.i.d. has $c_x(\tau)=0$ for $\tau\neq0$, so its
spectrum is flat — second-order analysis sees pure noise — while $x_t$ is plainly
forecastable in mean square from its past once we allow *nonlinear* functions.
This gap is exactly what the present chapter addresses.

The governing object is now a general **causal, time-invariant functional**

```{math}
:label: eq-nl-functional
x_t = F(\ldots,\epsilon_{t-1},\epsilon_t),
```

with $\{\epsilon_t\}$ i.i.d., mean zero, variance $\sigma^2$. The aim is a
constructive moving-average expansion of $F$ analogous to $d(L)$, a frequency
domain analogous to the spectrum, and existence conditions analogous to linear
indeterminism and invertibility.

A theme worth stating at the outset distinguishes the linear and nonlinear worlds:

> A linear filter **never creates new frequencies** — it only reweights and
> phase-shifts the frequencies already present in its input. A genuinely
> nonlinear functional **does** create new frequencies, mixing a pair
> $(\omega_1,\omega_2)$ of input frequencies into sum and difference components
> $\omega_1\pm\omega_2$ at the output.

Frequency mixing is the operational signature of nonlinearity, and it is precisely
what the second-order spectrum cannot detect but the higher-order spectra can.

(nl-volterra)=
## 2. The Volterra Moving-Average Series

The natural nonlinear generalization of $x_t=\sum_j d_j\epsilon_{t-j}$ adds
quadratic, cubic, and higher products of the innovations. The **Volterra
expansion** Volterra 1930; Schetzen 1980 is

```{math}
:label: eq-nl-volterra
x_t = h_0
  + \sum_{j\ge0} h_1(j)\,\epsilon_{t-j}
  + \sum_{j_1,j_2\ge0} h_2(j_1,j_2)\,\epsilon_{t-j_1}\epsilon_{t-j_2}
  + \sum_{j_1,j_2,j_3\ge0} h_3(j_1,j_2,j_3)\,\epsilon_{t-j_1}\epsilon_{t-j_2}\epsilon_{t-j_3}
  + \cdots
```

The deterministic functions $h_k:\mathbb{Z}_+^k\to\mathbb{R}$ are the **Volterra
kernels**. Without loss of generality each $h_k$ may be taken symmetric in its $k$
arguments, since only the symmetric part survives the symmetric sum. The first
kernel $h_1(\cdot)$ is exactly the ordinary linear MA filter — $h_1(j)=d_j$ — so
{eq}`eq-nl-volterra` contains the Wold representation as its $k=1$ term; the
higher kernels $h_2,h_3,\ldots$ encode genuine nonlinearity.

Convergence of {eq}`eq-nl-volterra` is required in $L^2$ (mean square), or almost
surely depending on the formulation, and this imposes summability conditions on
the kernels. For the quadratic truncation, for instance, $L^2$ convergence
requires $\sum_{j_1,j_2}h_2(j_1,j_2)^2<\infty$ together with the linear condition
$\sum_j h_1(j)^2<\infty$; the general statement is deferred to {ref}`Section 7
<nl-existence>`, where combinatorial moment growth makes the polynomial
conditions delicate.

(nl-wienerito)=
## 3. Orthogonalization: the Wiener–Itô (Hermite) Representation

The Volterra series {eq}`eq-nl-volterra` is the nonlinear analogue of a *raw*
moving average, but it is **not** the analogue of the *Wold* moving average,
because its building blocks are not orthogonal. This is the conceptual heart of
the chapter.

### The non-orthogonality problem

In the linear theory the terms $d_j\epsilon_{t-j}$ are mutually orthogonal because
$E[\epsilon_{t-i}\epsilon_{t-j}]=0$ for $i\neq j$, and this orthogonality is what
makes $d_j$ the coefficients of an *orthogonal projection* (and gives the clean
variance decomposition $Ex_t^2=\sigma^2\sum_j d_j^2$). For the nonlinear series
the monomials

$$
\epsilon_{t-j_1}\cdots\epsilon_{t-j_k}
$$

of different orders $k$ are generally **not** mutually orthogonal in $L^2$ when
the input is i.i.d.: even-order monomials have nonzero means
($E[\epsilon_{t-j}^2]=\sigma^2\neq0$) and nonzero cross-moments with lower-order
terms. Consequently the raw Volterra kernels $h_k$ are *not* the orthogonal
projections of $x_t$ onto the "$k$-th order nonlinearity," and they do not deliver
a clean variance decomposition.

### Wiener's fix: Hermite chaos

Wiener's solution (Wiener 1958), made rigorous by Itô (1951),
replaces the monomials by **Hermite polynomials** of the i.i.d. Gaussian input.
Recall that the Hermite polynomials $H_k$ are orthogonal with respect to the
standard Gaussian measure,

$$
E\!\left[H_m(\xi)\,H_n(\xi)\right]=m!\,\delta_{mn},\qquad \xi\sim\mathcal N(0,1),
$$

with $H_0=1$, $H_1(\xi)=\xi$, $H_2(\xi)=\xi^2-1$, $H_3(\xi)=\xi^3-3\xi$, and so
on. The orthogonalized expansion is

```{math}
:label: eq-nl-hermite
x_t = \sum_{k\ge0}\ \sum_{j_1\le\cdots\le j_k}
        g_k(j_1,\ldots,j_k)\,
        H_k\!\left(\epsilon_{t-j_1},\ldots,\epsilon_{t-j_k}\right),
```

where the terms — the multiple **Wiener–Itô integrals** — are pairwise orthogonal
*across* orders $k$ as well as within. For Gaussian i.i.d. input,
{eq}`eq-nl-hermite` is the rigorous nonlinear analogue of the Wold decomposition
(Wold 1938):

```{admonition} The nonlinear Wold theorem
:class: tip
For an $L^2$ functional of a Gaussian i.i.d. input, {eq}`eq-nl-hermite` is the
**unique orthogonal decomposition** of $x_t$ into homogeneous *chaoses* (the
order-$k$ subspaces). The space $L^2$ of the input splits as an orthogonal direct
sum $\bigoplus_{k\ge0}\mathcal{C}_k$, with $\mathcal{C}_0$ the constants,
$\mathcal{C}_1$ the linear (Wold) part, $\mathcal{C}_2$ the quadratic chaos, and so
on.
```

The **Hermite kernels** $g_k$ are recovered, exactly as in the linear case, by
*projection*: $g_k$ is the orthogonal projection of $x_t$ onto the $k$-th chaos,
equal to the conditional expectation of $x_t$ onto that subspace divided by $k!$.
The order-one kernel $g_1$ coincides with the Wold filter $d_j$; the constant
$g_0=Ex_t$. Because the chaoses are orthogonal, one again obtains a clean variance
decomposition,

$$
\operatorname{Var}(x_t)=\sum_{k\ge1} k!\sum_{j_1\le\cdots\le j_k} g_k(j_1,\ldots,j_k)^2,
$$

the nonlinear counterpart of $\sigma^2\sum_j d_j^2$.

### Relation between the two sets of kernels

The polynomial Volterra kernels $h_k$ and the Hermite kernels $g_k$ are related by
a **triangular linear transformation** whose entries involve the moments of
$\epsilon$. Triangularity reflects the fact that $H_k$ is a degree-$k$ polynomial:
expressing a monomial of degree $k$ in the Hermite basis involves only Hermite
polynomials of degree $\le k$, and conversely. For Gaussian input only even
corrections appear (odd central moments vanish); for example $\epsilon^2=H_2+ 1$
moves variance from the quadratic monomial into the order-zero (mean) term — which
is exactly the source of the non-orthogonality the Hermite basis repairs.

(nl-polyspectra)=
## 4. The Frequency Domain: Generalized Transfer Functions and Polyspectra

The linear filter $d(L)$ has a transfer function $d(e^{-i\omega})$, and the
spectrum of the output is $|d(e^{-i\omega})|^2\sigma^2$ (Section 6). Each Volterra kernel has its own transfer function.

### Generalized transfer functions

Define the order-$k$ **generalized transfer function** as the multidimensional
Fourier transform of the $k$-th kernel,

```{math}
:label: eq-nl-Hk
H_k(\omega_1,\ldots,\omega_k)
= \sum_{j_1,\ldots,j_k\ge0} h_k(j_1,\ldots,j_k)\,
   e^{-i(\omega_1 j_1+\cdots+\omega_k j_k)}.
```

The linear theory has only $H_1(\omega)=d(e^{-i\omega})$. The first genuinely
nonlinear object is the **quadratic transfer function** $H_2(\omega_1,\omega_2)$,
which describes how a pair of input frequencies $(\omega_1,\omega_2)$ is converted
into output power at the sum and difference frequencies $\omega_1\pm\omega_2$. This
is the precise content of the "frequency-mixing" statement of
{ref}`Section 1 <nl-intro>`: a sinusoidal pair fed through a quadratic term
produces, via $\cos\omega_1 t\,\cos\omega_2 t=\tfrac12[\cos(\omega_1+\omega_2)t
+\cos(\omega_1-\omega_2)t]$, energy at frequencies that were not present at the
input. Linear filters cannot do this.

### Polyspectra: the bispectrum and bicoherence

The second-order spectrum $g_x(\omega)$ is the Fourier transform of the
autocovariance $c_x(\tau)=E[x_t x_{t+\tau}]$. The natural higher-order analogues
are the Fourier transforms of the higher **cumulants**. The third-order cumulant
of a zero-mean process,

$$
c_3(\tau_1,\tau_2)=E[x_t\,x_{t+\tau_1}\,x_{t+\tau_2}],
$$

has as its double Fourier transform the **bispectrum**

```{math}
:label: eq-nl-bispectrum
B_x(\omega_1,\omega_2)=\sum_{\tau_1,\tau_2} c_3(\tau_1,\tau_2)\,
   e^{-i(\omega_1\tau_1+\omega_2\tau_2)}.
```

The bispectrum measures the **phase coupling** between the frequency components at
$\omega_1$, $\omega_2$, and $\omega_1+\omega_2$. For a process that is linear in
Gaussian innovations the third cumulant vanishes identically, so $B_x\equiv0$: a
nonzero bispectrum is a direct fingerprint of nonlinearity (or of non-Gaussian
shocks). The normalized, scale-free version is the **bicoherence**

```{math}
:label: eq-nl-bicoherence
b_x(\omega_1,\omega_2)=
\frac{|B_x(\omega_1,\omega_2)|^2}
     {g_x(\omega_1)\,g_x(\omega_2)\,g_x(\omega_1+\omega_2)},
```

a quantity in $[0,1]$ that plays, for quadratic coupling, the role the coherence
of Section 7 played for linear cross-dependence.
Both objects are computed and visualized in {doc}`40_nonlinear_lab`.

(nl-classes)=
## 5. Three Canonical Model Classes

Three structured subclasses of the Volterra series recur throughout the applied
literature; each has a distinctive kernel and therefore a distinctive bispectral
signature. The laboratory chapter shows that these signatures are visible in
simulated data.

### Hammerstein: static nonlinearity, then linear filter

$$
x_t=\sum_j \psi_j\, g(\epsilon_{t-j}).
$$

Expanding the static nonlinearity $g$ in a power (or Hermite) series gives
**diagonal** Volterra kernels: $h_k(j_1,\ldots,j_k)$ is nonzero only when
$j_1=\cdots=j_k$, because each $\epsilon_{t-j}$ enters $g(\epsilon_{t-j})$ at the
*same* lag. In the frequency domain this collapses to

$$
H_k(\omega_1,\ldots,\omega_k)=c_k\,\Psi(\omega_1+\cdots+\omega_k),
$$

so $H_2$ depends only on the **sum frequency** $\omega_1+\omega_2$. The
consequence is that bispectral mass concentrates along the anti-diagonal *ridges*
$\omega_1+\omega_2=\text{const}$. Notably the Hammerstein model leaves the
*shape* of the second-order power spectrum essentially unchanged (it merely
rescales it), so it is **invisible** to ordinary spectral analysis — the
bispectrum is the only linear-invariant that detects it.

### Wiener cascade: linear filter, then static nonlinearity

$$
u_t=\sum_j \psi_j\,\epsilon_{t-j},\qquad x_t=g(u_t).
$$

Expanding $g$ in a Taylor series gives **separable** kernels,

$$
h_k(j_1,\ldots,j_k)=\frac{g^{(k)}(0)}{k!}\,\psi_{j_1}\cdots\psi_{j_k},
\qquad
H_k(\omega_1,\ldots,\omega_k)=\frac{g^{(k)}(0)}{k!}\,\Psi(\omega_1)\cdots\Psi(\omega_k),
$$

a product of one-dimensional transfer functions. The quadratic kernel
$H_2(\omega_1,\omega_2)=\tfrac12 g''(0)\Psi(\omega_1)\Psi(\omega_2)$ is a
**rank-one** (separable) array, so bispectral mass concentrates at *point
clusters* $(\pm\omega_0,\pm\omega_0)$ where $\Psi$ peaks, rather than along
ridges. The Wiener cascade also generates a *harmonic* at $2\omega_0$ in the
ordinary power spectrum (through the self-convolution of $u_t$), so unlike the
Hammerstein model it is partly visible at second order.

### Bilinear processes

The bilinear model (Granger and Andersen 1978) couples the state to the noise,

$$
x_t=\sum_i a_i\,x_{t-i}+\epsilon_t+\sum_{i,j} c_{ij}\,x_{t-i}\,\epsilon_{t-j}.
$$

Solving recursively for $x_t$ in terms of past $\epsilon$'s generates a Volterra
series of **infinite order**: the cross-product term propagates into every kernel.
The bispectrum is correspondingly diffuse and asymmetric, dominated by neither
ridges nor point masses. Bilinear models are valued in practice because they
admit finite-dimensional state-space representations even though their Volterra
expansion is infinite.

```{admonition} GARCH and stochastic volatility
:class: note
GARCH and stochastic-volatility processes admit Volterra-type expansions of the
*level* $x_t$, but the more natural Hermite expansion is for the *squared* process
$x_t^2$. Diagnostically, the bicoherence is nearly flat for GARCH **levels** but
distinctly non-flat for the **squares** — a reminder that the right object to
spectrum-analyse depends on the question being asked.
```

(nl-identification)=
## 6. Constructing the Kernels from Cumulants and Polyspectra

The representation is not merely descriptive: the kernels can be recovered from
data through higher-order spectral analysis (Brillinger 1975; Priestley 1988; Nikias and Petropulu 1993).

Start from the cumulants of $x$. The autocovariance receives a linear and
nonlinear contribution,

```{math}
:label: eq-nl-autocov
c_x(\tau)=\operatorname{cov}(x_t,x_{t+\tau})
=\sigma^2\sum_j h_1(j)\,h_1(j+\tau)
 \;+\;\text{contributions from } h_2,h_3,\ldots,
```

so the power spectrum alone cannot separate the linear filter from the nonlinear
kernels. The third cumulant breaks the degeneracy. For **Gaussian** $\epsilon$
(so the third moment $\mu_3=0$ and all odd cumulants of the input vanish), the
third cumulant of the *output* comes purely from $h_2$ interacting with $h_1$
through combinatorial Isserlis/Wick sums. Fourier transforming yields a relation
of the form

```{math}
:label: eq-nl-bisp-factor
B_x(\omega_1,\omega_2)\ \propto\
H_1(\omega_1)\,H_1(\omega_2)\,H_2^{*}(\omega_1,\omega_2)\ +\ \text{permutations}.
```

This is the key identification equation. Given Gaussian i.i.d. innovations and the
linear filter $H_1$ already identified from the spectrum — up to the usual
Blaschke (phase) flips, as in Lii and Rosenblatt (1982) — the quadratic kernel $H_2$ is
recovered from the bispectrum by dividing out the $H_1$ factors,

$$
\hat H_2(\omega_1,\omega_2)\ \propto\
\frac{B_x(\omega_1,\omega_2)}{H_1(\omega_1)H_1(\omega_2)H_1^{*}(\omega_1+\omega_2)}.
$$

Higher kernels $H_k$ are identified, in the same way, from the order-$(k{+}1)$
polyspectrum. The procedure is constructive: estimate the polyspectra from data,
then solve a sequence of (in principle linear) inversion problems for the $H_k$.
Whether to attempt this at all is decided by the **linearity and Gaussianity
tests** of Subba Rao and Gabr (1980) and Hinich (1982), which test precisely
whether the bispectrum is flat (linear/Gaussian) or structured (nonlinear).

```{admonition} Classifying a model from its quadratic kernel
:class: tip
The three classes of {ref}`Section 5 <nl-classes>` are separated by the algebraic
*structure* of $\hat H_2$:

- **separable / rank-one** $\Rightarrow$ Wiener cascade;
- **diagonal** (mass on the anti-diagonals) $\Rightarrow$ Hammerstein;
- **neither** $\Rightarrow$ bilinear or higher-order interaction.

A singular-value decomposition of $|\hat H_2|$ makes this operational: the
fraction of singular-value mass in the leading singular value is a
**separability score** close to $1$ for the Wiener cascade and small for the
Hammerstein model. This SVD test is implemented in {doc}`40_nonlinear_lab`.
```

(nl-existence)=
## 7. Existence, Uniqueness, and Caveats

When does a process actually admit a Volterra/Wiener moving-average representation
in i.i.d. shocks? Three nontrivial conditions are needed; they are the nonlinear
counterparts of *linear indeterminism* and *invertibility* from
Sections 13–14.

1. **Bernoulli-shift property.** The process must be measurably generated by an
   i.i.d. sequence, $x_t=F(\ldots,\epsilon_{t-1},\epsilon_t)$ for some measurable
   $F$ — the exact analogue of *linear indeterminism*. Not every stationary
   ergodic process qualifies: Ornstein's theory (Ornstein 1974) distinguishes
   Bernoulli shifts from the broader class of $K$-automorphisms. Most parametric
   models used in practice are Bernoulli shifts.

2. **Causal nonlinear invertibility.** The $\sigma$-algebra generated by past $x$
   should equal that generated by past $\epsilon$ — the analogue of invertibility
   of $d(z)$, i.e. of $\epsilon_t$ lying in the space spanned by current and past
   $x$. The nonlinear version is harder and usually verifiable only model by model
   (bilinear and GARCH models have explicit parameter conditions).

3. **Summability/convergence of the kernels in $L^2$.** For the Hermite expansion
   the clean sufficient condition is $\sum_k k!\,\|g_k\|_2^2<\infty$; for the raw
   polynomial form the conditions are stricter, because the combinatorial growth
   of the moments of $\epsilon$ inflates the cross-terms.

When (1)–(3) hold, the nonlinear moving average **exists and is unique**, the
kernels are identified from the polyspectra given a non-degenerate cumulant
structure of $\epsilon$ (Rosenblatt 1985), and one has a complete nonlinear
analogue of the Wold representation.

When they fail — most commonly when only the *linear* projection structure is
invertible — one still has a Wold MA $x_t=d(L)\epsilon_t$ with **white-but-
dependent** innovations, but **no** genuine i.i.d.-driven Volterra MA. This is
the intermediate regime in which $\epsilon_t$ is uncorrelated yet not independent
(the $x_t=\epsilon_t\epsilon_{t-1}$ example of {ref}`Section 1 <nl-intro>`), and
the **bicoherence is the natural diagnostic** that one is in it: a flat spectrum
combined with a structured bicoherence is the signature of nonlinear dependence
hiding behind a white second-order structure.

## 8. Summary

```{list-table} What each order of analysis can detect
:header-rows: 1
:name: tab-nl-summary

* - Object
  - Linear / Gaussian
  - Hammerstein
  - Wiener cascade
  - Bilinear
* - Power spectrum $g_x(\omega)$
  - Full description
  - Same shape (rescaled) — **invisible**
  - Extra harmonic at $2\omega_0$
  - Altered shape
* - Bispectrum $B_x(\omega_1,\omega_2)$
  - Zero
  - Ridge on $\omega_1+\omega_2=\text{const}$
  - Point peaks at $(\pm\omega_0,\pm\omega_0)$
  - Diffuse, asymmetric
* - $|H_2|$ structure (SVD)
  - —
  - Diagonal (low rank-one score)
  - Separable (high rank-one score)
  - Neither
```

The linear Wold representation of Chapter XI is one term — the order-one chaos —
of a richer nonlinear decomposition. The second-order spectrum that has occupied
us throughout sees only that one term; the polyspectra see the rest. The
constructive route is: estimate the spectrum to get $H_1$, estimate the bispectrum
to get $H_2$ (and classify the model by the structure of $\hat H_2$), and proceed
up the hierarchy of polyspectra for the higher kernels — always subject to the
existence conditions of {ref}`Section 7 <nl-existence>`. The next chapter turns
each of these statements into running code.

## References

- Brillinger, D. R. (1975). *Time Series: Data Analysis and Theory.* Holt, Rinehart and Winston.
- Granger, C. W. J., and A. P. Andersen (1978). *An Introduction to Bilinear Time Series Models.* Vandenhoeck and Ruprecht.
- Hinich, M. J. (1982). Testing for Gaussianity and linearity of a stationary time series. *Journal of Time Series Analysis* 3(3), 169–176.
- Itô, K. (1951). Multiple Wiener integral. *Journal of the Mathematical Society of Japan* 3(1), 157–169.
- Lii, K.-S., and M. Rosenblatt (1982). Deconvolution and estimation of transfer function phase and coefficients for non-Gaussian linear processes. *Annals of Statistics* 10(4), 1195–1208.
- Nikias, C. L., and A. P. Petropulu (1993). *Higher-Order Spectra Analysis: A Nonlinear Signal Processing Framework.* Prentice Hall.
- Ornstein, D. S. (1974). *Ergodic Theory, Randomness, and Dynamical Systems.* Yale University Press.
- Priestley, M. B. (1988). *Non-linear and Non-stationary Time Series Analysis.* Academic Press.
- Rosenblatt, M. (1985). *Stationary Sequences and Random Fields.* Birkhäuser.
- Schetzen, M. (1980). *The Volterra and Wiener Theories of Nonlinear Systems.* John Wiley & Sons.
- Subba Rao, T., and M. M. Gabr (1980). A test for linearity of stationary time series. *Journal of Time Series Analysis* 1(2), 145–158.
- Volterra, V. (1930). *Theory of Functionals and of Integral and Integro-Differential Equations.* Blackie and Son.
- Wiener, N. (1958). *Nonlinear Problems in Random Theory.* MIT Press.
- Wold, H. O. A. (1938). *A Study in the Analysis of Stationary Time Series.* Almqvist and Wiksell.
