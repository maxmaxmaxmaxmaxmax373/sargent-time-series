# Theories Restricting One-Sided Projections

There is a class of linear rational expectation models that give rise to a theoretical restriction on the one-sided projection of a process $y$ on a process $x$, but which leave the serial correlation of the residual process unrestricted. These specifications then fit naturally into the setting analyzed in the preceding section.

These specifications are motivated by an idea of Shiller (1972), who used an iterated projections argument to show that when private agents are assumed to observe more information relevant for forecasting than does the econometrician, the error thereby induced in the equation of interest of the econometrician remains orthogonal to part of the econometrician's information set. These orthogonality conditions can restrict the variables observable to the econometrician, and provide a sufficient basis for estimating and testing the model. We illustrate Shiller's idea in the context of a rational expectations version of Cagan's model of hyperinflation.

The log of the price level $p_t$ is assumed to be determined as the average of expected future values of the log of the money supply:

```{math}
:label: eq-160
p_t = (1-\lambda)\sum_{j=0}^\infty \lambda^jP[m_{t+j}|\Omega_t],\quad |\lambda|<1
```

where $\Omega_t \supset \Omega_{t-1} \supset \Omega_{t-2}, \ldots$, and where $\Omega_t$ is the information set available to private agents. We assume that there is an econometrician who observes a proper subset $\phi_t \subset \Omega_t$ at each date, where $\phi_t \supset \phi_{t-1} \supset \phi_{t-2}, \ldots$. Thus, private agents are assumed to base their forecasts of $m_{t+j}$ on more information than is available to the econometrician.

Shiller's idea was to define the error term

```{math}
:label: eq-161
s_t = (1-\lambda)\sum_{j=0}^\infty \lambda^j[P[m_{t+j}|\Omega_t]-P[m_{t+j}|\phi_t]]
```

From the law of iterated projections, it follows that

$$
P[s_t|\phi_t] = 0
$$

From the orthogonality principle that states that a projection error is orthogonal to every component in the information set used to form the projection, it follows that

$$
Es_tx_{t-k}=0,\quad k\geq0
$$

for any $x_t \in \phi_t$. Therefore, we can represent {eq}`eq-160` as

```{math}
:label: eq-162
p_t = (1-\lambda)\sum_{j=0}^\infty \lambda^jP[m_{t+j}|\phi_t] + s_t,
```

where $s_t$ is orthogonal to $\phi_t$. To be more concrete, let us suppose that the econometrician has data on $(m_t,m_{t-1},\ldots)$ as his $\phi_t$. Suppose that $m_t$ has the moving average representation

$$
m_t=c(L)\epsilon_t
$$

where $c(L)=\sum_{j=0}^\infty c_jL^j$, $\sum_{j=0}^\infty c_j^2<+\infty$, $\epsilon_t=m_t-P[m_t|m_{t-1},m_{t-2},\ldots]$. Then from formula {eq}`eq-89` we have that

```{math}
:label: eq-163
p_t = (1-\lambda)\left[\frac{1-\lambda c(\lambda)c(L)^{-1}L^{-1}}{1-\lambda L^{-1}} \right]m_t + s_t,
```

```{math}
:label: eq-164
c(L)^{-1}m_t = \epsilon_t
```

Equation {eq}`eq-163` is a projection equation, i.e., $Es_tm_{t-j}=0$ for $j\geq0$. However, in general $Es_tm_{t-j} \neq 0$ for $j<0$. Intuitively, the reason can be understood from equation {eq}`eq-161` that defines $s_t$. The random variable $s_t$ reflects information in $\Omega_t$ that is useful for predicting future variables of $m$ and that is not reflected in $\phi_t$. For this reason, $s_t$ will in general be correlated with $m_{t+k}$, $k>0$.

In addition, $s_t$ is in general serially correlated. Equation {eq}`eq-161` defines $s_t$ as a linear least squares error in projecting $z_t \equiv (1-\lambda)\sum_{j=0}^\infty\lambda^jP[m_{t+j}|\Omega_t]$ on $\phi_t$. Because lagged values of the random variable $z_t$ being projected are *not* included in the information set $\phi_t$ being projected onto, there is no reason to expect $s_t$ to be serially uncorrelated. Compare this situation with the argument in Chapter XVI, p. 445 that implies that a forecast error is serially uncorrelated in the case in which lagged values of the variable being projected are in the information set.

Thus, {eq}`eq-163` is an example that fits the situation analyzed in the preceding section: it is a projection equation, but the disturbance $s_t$ is both serially correlated and correlated with future $m_{t+k}$'s.

Using $a(L) = c(L)^{-1} = (1- a_1L - a_2L^2 - \cdots)$, we can rewrite {eq}`eq-163` and {eq}`eq-164` as

```{math}
:label: eq-165
p_t = \sum_{j=0}^\infty h_jm_{t-j} + s_t
```

```{math}
:label: eq-166
m_t = \sum_{j=1}^\infty a_jm_{t-j} + \epsilon_t
```

where

```{math}
:label: eq-167
h(L)=(1-\lambda)\frac{[1-\lambda a(\lambda)^{-1}a(L)L^{-1}]}{1-\lambda L^{-1}}
```

Using the orthogonality conditions $Es_tm_{t-j} = E\epsilon_t m_{t-j} = 0$ for $j\geq1$, {eq}`eq-165` and {eq}`eq-166` imply

```{math}
:label: eq-168
Ep_tm_{t-k} = \sum_{j=0}^\infty h_jEm_{t-j}m_{t-k},\quad k\geq0
```

```{math}
:label: eq-169
Em_tm_{t-k} = \sum_{j=1}^\infty a_jEm_{t-j}m_{t-k},\quad k\geq0
```

Equations {eq}`eq-167`, {eq}`eq-168`, and {eq}`eq-169` supply an extensive set of restrictions on the moments of the data that are observable to the econometrician. The restrictions {eq}`eq-169` are sufficient to uniquely determine the $a_j$'s (equations {eq}`eq-169` are the univariate Yule-Walker equations for the $m$ process). From {eq}`eq-167`, the parameters $h_j$ in {eq}`eq-168` are not all free to vary independently. Rather, given the $a_j$'s, {eq}`eq-167` restricts the $h_j$'s to vary together only as a function of the single parameter $\lambda$.

A way of looking at the idea of Hayashi and Sims (1983) is to note that equations {eq}`eq-168` and {eq}`eq-169` fail to exhaust all of the implications of {eq}`eq-165`, {eq}`eq-166`, and {eq}`eq-167`. Equation {eq}`eq-168` is the condition that $s_t$ is orthogonal to $\{m_t,m_{t-1},m_{t-2},\ldots\}$. "Backward filtering" destroys the orthogonality conditions because in general $s_{t-1}$ is *not* orthogonal to $\{m_t,m_{t-1},m_{t-2},\ldots\}$. However, for $k\geq0$, $s_{t+k}$ *is* orthogonal to $\{m_t,m_{t-1},m_{t-2},\ldots\}$. Therefore, "forward filtering" leaves the orthogonality conditions intact. Forward filtering can be interpreted as a measure designed to impose upon the covariance structure more of the restrictions implied by the theory. By imposing more of these restrictions, more efficient estimation procedures can be constructed.[^fn-osp-1]

[^fn-osp-1]: See Hansen and Sargent (1982a) for a formal development of this argument.
