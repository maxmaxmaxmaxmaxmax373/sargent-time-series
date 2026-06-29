# Sims's Application to Money and Income

Economists at the Federal Reserve Bank of St. Louis (Andersen and Jordan, 1968) have computed estimates of one-sided distributed lag regressions of (the log of) nominal income ($y_t$) against (the log of) money ($m_t$):

```{math}
:label: eq-133
y_t = \sum_{j=0}^\infty h_jm_{t-j} + \eta_t,
```

where $E \eta_t m_{t-j}=0$ for $j=0,1,2, \ldots$. Those economists recommend that the $h_j$ be taken seriously and be regarded as depicting the response of nominal income to exogenous impulses in the money supply. However, Keynesian economists have tended not to regard the $h_j$ as good estimates of the response pattern (or "dynamic multipliers") of nominal income to money. Their argument has two parts. First, in the kind of macroeconometric model the Keynesians have in mind, even were it true that money had been made to behave exogenously with respect to nominal income, the "final form" for money income has many additional right-hand-side variables not included in {eq}`eq-133`, e.g.,

```{math}
:label: eq-134
y_t = \sum_{j=0}^\infty v_jm_{t-j} + \sum_{j=0}^\infty w_jz_{t-j}+ \epsilon_t,
```

where $z_t$ is a vector of stochastic processes including government tax and expenditures parameters and $w_j$ is a vector conformable to $z_t$; the error term $\epsilon_t$ is a stationary stochastic process that obeys the orthogonality conditions $E \epsilon_t m_{t-j}=E \epsilon_t z_{t-j}=0$ for $j=0,\pm 1,\pm 2, \ldots$.

The strong condition that $\epsilon$ must be orthogonal to $m$ and $z$ at all leads and lags is the requirement that $m$ and $z$ be "strictly econometrically exogenous with respect to $y$" in relation {eq}`eq-134`. These orthogonality conditions characterize {eq}`eq-134` as a "final form" relationship. In {eq}`eq-134` the $v_j$ are the dynamic money multipliers that depict the average response of $y_t$ to a unit impulse in $m$, holding constant the $z$'s. Applying the law of iterated projections to {eq}`eq-134`, we obtain

$$
P[y_t|m_t,m_{t-1},\ldots] = \sum_{j=0}^\infty v_jm_{t-j} + \sum_{k=0}^\infty w_kP[z_{t-k}|m_t,m_{t-1},\ldots].
$$

Let $P[z_{t-k}|m_t,m_{t-1},\ldots] = \sum_{j=0}^\infty \alpha_{kj}m_{t-j}$. Then we have[^fn-sims-1]

$$
P[y_t|m_t,m_{t-1},\ldots] = \sum_{j=0}^\infty v_jm_{t-j} + \sum_{k=0}^\infty w_k\sum_{j=0}^\infty \alpha_{kj}m_{t-j}
$$

or

```{math}
:label: eq-135
y_t= \sum_{j=0}^\infty \left(v_j + \sum_{k=0}^\infty w_k\alpha_{kj}\right)m_{t-j} + \eta_t
```

where by the orthogonality principle we have $E \eta_t m_{t-j}=0$, $j=0,1,2,\ldots$. Now {eq}`eq-135` is identical with {eq}`eq-133`, so that the population $h_j$ of {eq}`eq-133` obey

$$
h_j = v_j + \sum_{k=0}^\infty w_k\alpha_{kj}.
$$

Therefore, the $h_j$ in general *do not* equal the money multipliers, the $v_j$. The $h_j$ are "mongrel" coefficients that do not indicate the typical average response of $y$ to exogenous impulses in $m$, everything else, namely the $z$'s, being held constant. For this reason, Keynesians would argue, estimating Equation {eq}`eq-133` is not a good way of estimating the dynamic multipliers, the $v_j$.

Now project both sides of {eq}`eq-134` against the entire sequence $\{m_{t-j}\}^\infty_{j=-\infty}$ to get

```{math}
:label: eq-136
y_t= \sum_{j=0}^\infty h_jm_{t-j} + \sum_{k=0}^\infty w_k \sum_{j=-\infty}^\infty \gamma_{kj}m_{t-j} + \xi_t
```

where $E\xi_tm_{t-j} = 0$ for all $j$ and

$$
P(z_{t-k}|\{m_{t-j}\}^\infty_{j=-\infty}) = \sum_{j=-\infty}^\infty \gamma_{kj}m_{t-j}
$$

where $\gamma_{kj}$ is a vector of coefficients. We can write {eq}`eq-136` as

$$
P(y_t|\{m_{t-j}\}^\infty_{j=-\infty}) = \sum_{j=-\infty}^\infty d_jm_{t-j}
$$

where

$$
d_j = \begin{cases}
h_j + \sum_{k=0}^\infty w_k\gamma_{kj}, & j\geq0 \\
\sum_{k=0}^\infty w_k\gamma_{kj}, & j<0
\end{cases}
$$

In general, so long as the processes $m_t$ and $z_t$ are correlated (as we had to assume to make the argument that the St. Louis $h_j$ are mongrel parameters), the $\gamma_{kj}$ and therefore the $d_j$ will not vanish for some $j<0$. That is because in general future $m$'s will help explain current and past $z_t$.[^fn-sims-2] Therefore, so long as the $w_k$ are not zero in the final form {eq}`eq-134`, i.e., so long as the $z$'s appear in the final form of $y_t$, the projection of $y_t$ on current, lagged, and future $m$'s is predicted to be two-sided. For this reason, a test of the null hypothesis that the projection of $y_t$ on the entire $\{m\}$ process is one-sided (i.e., it equals the projection of $y_t$ on current and past $m$'s alone) can be regarded as testing the null hypothesis that the $w_k$ in {eq}`eq-134` *are* zeros. But remember that the contention that the $w_k$ are not zero is what underlies the Keynesian objection against interpreting the St. Louis equation's $h_j$ as estimates of the dynamic money multipliers. So computing the two-sided projection

```{math}
:label: eq-137
y_t = \sum_{j=-\infty}^\infty \delta_jm_{t-j} + \hat{\eta}_t
```

where $E\hat{\eta}_tm_{t-j}=0$ for all $j$, and testing the null hypothesis that $\delta_j=0$ for all $j<0$ provides a means of testing the null hypothesis that the St. Louis equation is "properly specified" — i.e., that it is appropriate to set the $w_k$ equal to zero.

Using post-World War II U.S. data, Sims estimated {eq}`eq-137` and implemented the preceding test. He found that he could not reject with high confidence the hypothesis that future $m$'s bear zero coefficient in {eq}`eq-137`. In general, if the Keynesian objection to the St. Louis equation were correct, in large enough samples one would expect to reject the hypothesis tested by Sims. Sims's particular statistical results have provoked much controversy. Since his tests are subject to the usual kinds of type I and type II statistical errors, there is some room for disagreement about how far his reports go in confirming using the St. Louis equation to estimate money multipliers. Nevertheless, it should be recognized how much of a contribution Sims made in providing a formal statistical setting in which one could in principle subject to statistical testing the Keynesian claims made against the St. Louis approach. Before Sims's work, those claims were entirely *a priori* and, though they had been made repeatedly, had never been subjected to any empirical tests.

As it happens, the test implemented by Sims is also useful in discriminating against another hypothesis which has often been advanced to argue that the St. Louis equation {eq}`eq-133` is not a legitimate final form (i.e. does not have a disturbance that obeys the requirement that it be orthogonal to past, present, and future $m$'s). The argument is that the money supply fails to be exogenous in {eq}`eq-133` because the monetary authority has set $m$ via some sort of feedback rule on lagged $y$'s. For example, it is often asserted that the Federal Reserve "leans against the wind," increasing $m$ faster in recessions, and more slowly in a boom. If the Fed behaved this way, it could mean that the projection {eq}`eq-133` of $y$ on $m$ partly reflects this feedback from past $y$ to $m$ as well as the effect of $m$ on $y$. Furthermore, such behavior by the Fed would in general lead us to expect the projection of $y$ on the entire $m$ process to differ from the projection of $y$ on current and past $m$'s, so that the $\eta_t$ in {eq}`eq-133` would not obey the restrictions $E\eta_t m_{t-s} = 0$ for all $s$; i.e., {eq}`eq-133` would not be a final form.

Now Sims's theorems assure us that if the projection of $y_t$ on $\{ m_{t-k}\}^\infty_{-\infty}$ is one-sided on the present and past (as Sims was unable to reject), then there exists a representation (i.e., a model consistent with the data) of the form

$$
m_t=C^{11}(L)\epsilon_t, \qquad y_t = d(L)m_t + C^{22}(L)u_t
$$

where $Eu_t\epsilon_s = 0$ for all $t,s$, and $d(L),C^{11}(L),C^{22}(L)$ are one-sided on the present and past. This representation is one in which there is no feedback from $y$ to $m$. Thus, Sims's results are consistent with (but do not necessarily imply) the view that there was no systematic feedback from $y$ to $m$ in the sample period he studied.

Sims's work on money and income was important because it provided a valid framework for testing empirically some often-stated objections to interpreting St. Louis regressions as final form equations.

A cautionary counterpoint to this innovation-accounting program appears in
{doc}`A Difficulty in Interpreting Vector Autoregressions <36a_interpreting_vars>`: in a class
of rational expectations models the white noise a vector autoregression recovers is *not* the
white noise that is fundamental for the agents, so the impulse responses and variance
decompositions computed from the autoregression can misrepresent how the economy responds to
the surprises that actually move agents.

[^fn-sims-1]: This is a version of H. Theil's "omitted variable theorem." See Theil (1971, pp. 548–550).

[^fn-sims-2]: Unless $m_t$ is strictly exogenous with respect to the vector $z_t$ or, equivalently, the vector $z_t$ does not Granger cause $m_t$.
