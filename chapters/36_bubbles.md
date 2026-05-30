# Bubbles

Following Blanchard and Watson (1982), consider the stochastic expectational difference equation

```{math}
:label: eq-203
y_t = \lambda P_t y_{t+1} + x_t,\quad |\lambda|<1
```

where $x_t$ is a stationary autoregressive process

```{math}
:label: eq-204
x_t = \rho x_{t-1} + \epsilon_t,\quad |\rho|<1
```

where $\epsilon_t$ is a white noise that is fundamental for $x_t$. In {eq}`eq-203`, $P_t$ is the linear least squares projection operator, conditioned on information known at $t$. An application of formula {eq}`eq-90` shows that the *stationary* solution of {eq}`eq-203` given {eq}`eq-204` is

```{math}
:label: eq-205
y_t = \frac{1}{1-\lambda\rho}x_t
```

Blanchard and Watson noted that in addition to the stationary solution {eq}`eq-205` there are many nonstationary solutions. These nonstationary solutions can be characterized as follows. Let $c_t$ be any martingale, that is, let $c_t$ be any stochastic process that satisfies $P_t c_{t+1} = c_t$. Then a solution of {eq}`eq-203` is

```{math}
:label: eq-206
y_t = \frac{1}{1-\lambda\rho}x_t + \left(\frac{1}{\lambda}\right)^t c_t
```

That {eq}`eq-206` is a solution of {eq}`eq-203` can be verified directly.

Three examples of martingales $c_t$ can usefully be given. A first is the constant $c_t = c$ for all $t$. This is the sort of deterministic bubble encountered in Chapter IX. A second is the one proposed by Blanchard and Watson, namely, the process

$$
c_{t+1} = \begin{cases}
  c_t/\pi & \text{with prob } \pi, \quad 0<\pi<1 \\
  0 & \text{with prob } 1-\pi
\end{cases}
$$

The process $c_{t+1}$ is readily verified to be a martingale. A third example is generated from the $x_t$ process itself. Simply set $c_t = \rho^{-1}x_t$, which is a martingale in light of {eq}`eq-204`.

In Chapter XIV, we study a model in which a transversality condition serves to make setting $c_t= 0$ the only admissible solution for a version of equation {eq}`eq-203`. There has recently been work designed to estimate and test models in which there are insufficient boundary conditions to justify setting $c_t = 0$. (See Blanchard and Watson (1982), Meese (1986), Sargent and Wallace (1985), and Hamilton and Whiteman (1985).) Such models have been proposed as candidates for understanding the stock market, foreign exchange rates, and hyperinflations.
