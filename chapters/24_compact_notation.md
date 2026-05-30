# A Compact Notation

It is always possible to write an $m$th-order difference equation in terms of a vector first-order system. For example, consider the bivariate system

$$
x_{1,t+1} = \alpha_1 x_{1,t} + \cdots + \alpha_m x_{1,t-m+1} + \alpha_{m+1}x_{2,t} + \cdots + \alpha_{2 m}x_{2,t-m+1} + \epsilon_{1,t+1}
$$

```{math}
:label: eq-103
x_{2,t+1} = \beta_1 x_{1,t} + \cdots + \beta_m x_{1,t-m+1} + \beta_{m+1}x_{2,t} + \cdots + \beta_{2m}x_{2,t-m+1} + \epsilon_{2,t+1}
```

where $(\epsilon_{1,t+1},\epsilon_{2,t+1})$ are two serially uncorrelated white-noise processes. Equations {eq}`eq-103` can be written as

```{math}
:label: eq-104
x_{t+1} = A x_t + \epsilon_{t+1}
```

$$
A =
\begin{bmatrix}
    \alpha_1 & \alpha_2 & \cdots & \alpha_m & \alpha_{m+1} & & & \cdots & \alpha_{2 m} \\
    1 & 0 & \cdots & 0 & 0 & & & \cdots & 0 \\
    \vdots & & & & & & & & \\
    0 & 0 & 1 & 0 & 0 & & & & \\
    \beta_1 & \beta_2 & \cdots & \beta_m & \beta_{m+1} & & & \cdots & \beta_{2 m} \\
    0 & & & 0 & 1 & & 0 & \cdots & 0 \\
    \vdots & & & & & & & & \\
    0 & 0 & & 0 & 0 & \cdots & 0 & 1 & 0
\end{bmatrix}
\leftarrow(m+1)\text{th row}\\
\uparrow(m+1)\text{th column}
$$

$$
x_{t+1} =
\begin{bmatrix}
    x_{1, t+1} \\
    x_{1, t} \\
    \vdots \\
    x_{1, t-m+2} \\
    x_{2, t+1} \\
    x_{2, t} \\
    \vdots \\
    x_{2, t-m+2}
\end{bmatrix}
\qquad
\epsilon_{t+1} =
\begin{bmatrix}
\epsilon_{1, t+1} \\
0 \\
0 \\
\vdots \\
0 \\
\epsilon_{2, t+1} \\
0 \\
0 \\
\vdots \\
0
\end{bmatrix}
\leftarrow(m+1)\text{th row}
$$

The solution of the vector difference equation {eq}`eq-104` can be written

```{math}
:label: eq-105
x_{t+\tau} = A^\tau x(t) + \epsilon(t + \tau) + A \epsilon(t + \tau - 1) + \cdots + A^{\tau-1}\epsilon(t + 1).
```

Since $E \epsilon(t + \tau)x(t)' = 0_{2m \times 2m}$ for all $\tau \geq 1$, multiplying the solution {eq}`eq-105` through by $x(t)$'s and taking expected values gives the matrix Yule-Walker equation

```{math}
:label: eq-106
E x_{t+\tau}x_t' = A^{\tau} E x_t x_t', \quad \tau \geq 1, \qquad \text{or} \qquad C_x({\tau}) = A^{\tau} C_x(0), \quad \tau \geq 1,
```

where $C_x(\tau) = E x(t + \tau)x(t)'$. As before, we have the result that the covariogram (this time the matrix covariogram) obeys the deterministic part of the difference equation with initial conditions given by the lagged covariances that are in $C_x(0)$.

Using the compact notation {eq}`eq-104`, it is straightforward to show that the cross-spectral density matrix $\Omega(e^{i\omega})$ of the vector $x$ process is given by

```{math}
:label: eq-107
\Omega(e^{-i \omega}) = (e^{i \omega} I - A)^{-1} V(I e^{-i\omega} - A')^{-1}
```

where $V = E \epsilon_t \epsilon_t'$ and where it is assumed that the process is stationary, which requires that the eigenvalues of $A$ have absolute values less than unity.

Assuming that the eigenvalues of $A$ are distinct, it is possible to represent $A$ in the form $A = P \Lambda P^{-1}$ where the columns of $P$ are the eigenvectors of $A$ while $\Lambda$ is the diagonal matrix whose diagonal entries are the eigenvalues of $A$. Then we have

$$
A^t = P \Lambda^t P^{-1}
$$

so that the solution {eq}`eq-106` can be written

$$
C_x(\tau) = P \Lambda^{\tau} P^{-1} C_x(0).
$$

This expression shows how the eigenvalues of $A$ govern the behavior of the solution. It also illustrates how increasing the number of variables in the system or increasing the number of lags in any particular equation increases the order of the $A$ matrix, and thereby contributes to the potential for generating complicated covariograms. Reference to this point can be used to show, for example, that while a one variable, first order difference equation cannot deliver a covariogram with damped oscillations of period greater than two periods (the periodicity if the single root is negative), a multivariate, first-order (i.e., single-lag) system can have complex roots and may therefore generate oscillatory covariograms.

Formula {eq}`eq-102` or {eq}`eq-107` has been used to summarize and analyze the stochastic properties of linear macroeconometric models. For interesting examples of such work, the reader is referred to articles by Chow and Levitan (1969) and by Howrey (1971).
