# Index Models

Let $y_t$ be an $(n \times 1)$ covariance stationary, linearly indeterministic stochastic process. The process $y_t$ is said to satisfy an unobservable *index model*, in the sense of Sargent and Sims (1977), if it possesses a representation:

```{math}
:label: eq-55
y_t = \Lambda(L) f_t + D(L) \epsilon_t
```

where $f_t$ is a $(k \times 1)$ vector white noise, $\epsilon_t$ is an $(n \times 1)$ vector white noise, $\Lambda(L)$ is an $(n \times k)$ matrix of square summable polynomials in the lag operator, $D(L)$ is an $(n \times n)$ diagonal matrix of square summable polynomials in the lag operator, and the white noises $f_t$, $\epsilon_t$, satisfy the orthogonality conditions

```{math}
:label: eq-56
\begin{aligned}
E \begin{bmatrix} f_t \\ \epsilon_t \end{bmatrix} \begin{bmatrix} f_{t-s} \\ \epsilon_{t-s}\end{bmatrix}^T &= \begin{bmatrix} \Sigma_f & 0 \\ 0 & \Sigma_{\epsilon}\end{bmatrix}, \qquad s = 0 \\
&= [0] \qquad s \neq 0
\end{aligned}
```

where $\Sigma_f$ and $\Sigma_{\epsilon}$ are each diagonal matrices. Thus, each component of $f_t$ and $\epsilon_t$ is orthogonal to every other component at all times, and to itself at all nonzero leads and lags. The model is restrictive when $k$ is sufficiently smaller than $n$. Usually, it has been assumed that $k=1$ or $2$ in applied work.

Equations ({eq}`eq-55` and {eq}`eq-56`) imply that the covariance generating function of $y_t$ satisfies

```{math}
:label: eq-57
S_y(z) = \Lambda(z)\Sigma_f \Lambda(z^{-1})^T + D(z)\Sigma_{\epsilon} D(z^{-1})^T
```

where $T$ denotes transposition. Equation {eq}`eq-57` states that the covariance generating function of $y$ (or the spectral density matrix of $y$) is the sum of a matrix of rank $k$, namely $\Lambda(z)\Sigma_f \Lambda(z^{-1})^T$, and a diagonal matrix of rank $n$, namely $D(z)\Sigma_{\epsilon} D(z^{-1})^T$. In other words, all of the covariance among distinct components of y is mediated through their common dependence on the $k$ indexes $f_t$. Equation {eq}`eq-57` is a frequency domain version of a "factor analysis" model, being a factor analytic model for each value of $z$. Geweke and Singleton (l981a,b) call {eq}`eq-57` the "dynamic factor" model.

Following T. C. Koopmans (1947), Sargent and Sims (1977) proposed {eq}`eq-55` with $k=1$ as a statistical model that seemed consistent with the conception of the business cycle used by Arthur Burns and Wesley Claire Mitchell (1946) at the National Bureau of Economic Research. With $k=1$, the covariances across distinct series entirely reflect their common dependence on the one-dimensional shock $f_t$

Sargent and Sims (1977) described some macroeconomic models that would lead an index model like {eq}`eq-55` to fit well. More recently, several researchers have applied stochastic optimal growth models which have assumed the form of a one unobservable index model. For example, see Kydland and Prescott (1982) and Altug (1985). The one index in these models is the innovation to a technology shock, while the components of $\epsilon_t$, are interpreted simply as mutually orthogonal measurement errors. Brock (1982) described a specification of preferences and technology in a stochastic growth model which could lead stock prices to exhibit an unobservable index structure.

A simple example of such a model can be obtained by adding mutually orthogonal measurement errors to the true series for $y_{nt}, c_t, $ and $k_{t+1} - k_t$, in the linear stochastic growth model of Chapter XII, Section 9. The one index in this model is $\theta_t$ the innovation to the technology shock.
