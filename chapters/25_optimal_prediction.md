# Optimal Prediction: Compact Notation

Using the fact that $\epsilon_t$ in {eq}`eq-104` is a serially uncorrelated vector process, it is straightforward to deduce from {eq}`eq-105` that the projection of $x_{t+\tau}$ against $x_t$ is given by

```{math}
:label: eq-108
P[x_{t+\tau}|x_t] = A^{\tau}x_t.
```

This is a compact formula for linear least squares predictors of a vector governed by a finite order stochastic difference equation. As an example illustrating the use of this formula, return to the portfolio balance example {eq}`eq-93`, which leads to a solution for the log of the price level of the form

```{math}
:label: eq-109
p_t = (1 - \lambda) \sum_{j=0}^\infty \lambda^j P_t m_{t+j}
```

where $\lambda = -\alpha/(1-\alpha)$, and where $m_t$ is the log of the money supply. Suppose that $m_t$ is the first element of a vector $x_t$ that evolves according to $x_t = A x_{t-1} + \epsilon_t$ where $\epsilon_t$ is a vector white noise. Let $e$ be the unit vector that validates our writing $m_t = e x_t$. Then substituting {eq}`eq-108` into the above solution for $p_t$ gives

$$
p_t = (1 - \lambda) e\left(\sum_{j=0}^\infty \lambda^j A^j\right)x_t.
$$

If the eigenvalues of $A$ are bounded by $1/\lambda$ in modulus,[^fn-opt-1] then we have that $\sum_{j=0}^\infty \lambda^j A^j = (I -\lambda A)^{-1}$. Therefore our solution can be represented

```{math}
:label: eq-110
\begin{aligned}
p_t &= (1 - \lambda)e(I - \lambda A)^{-1}x_t \\
x_t &= A x_{t-1} + \epsilon_t
\end{aligned}
```

Two comments about this derivation are in order. First, in the special case in which only lagged $m$'s appear in the $x_{t-1}$ vector, the above formula is equivalent with formula {eq}`eq-90`. In fact formula {eq}`eq-90` could be derived from the above one simply by explicitly inverting $(I - \lambda A)$.

Second, we notice from {eq}`eq-110` that not only lagged $m$'s but also any other variables that appear in the vector $x_t$ also enter the equation {eq}`eq-110` for $p_t$. Thus, any variables that help predict future $m$'s end up in the equation {eq}`eq-110` expressing $p_t$ as a function of current and lagged variables.

[^fn-opt-1]: Nonstationarity of $\{m_t\}$ in the form of $\max_i |\lambda_i|>1$ still leaves the proposed solution valid as long as $\max_i[(-\alpha/(1-\alpha))\lambda_i] < 1$.
