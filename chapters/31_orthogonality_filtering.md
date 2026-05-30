# Preservation of Orthogonality Conditions under Filtering

Various rational expectations models often lead to restrictions that take the form of orthogonality conditions. Hayashi and Sims (1983) and Hansen and Sargent (1982a) have studied the effects that filtering an equation has on the validity of such orthogonality conditions. They were concerned with understanding the circumstances under which filtering to remove serial correlation from a set of disturbances could create a representation among filtered variables not satisfying the theoretical orthogonality conditions that characterized the original unfiltered representation.

Suppose that one is interested in the one-sided projection equation

```{math}
:label: eq-153
y_t=\sum_{j=0}^\infty h_j x_{t-j} + v_t
```

where $Ev_tx_{t-k} = 0$ for $j\geq0$. Suppose that $v_t$ is serially correlated, and has covariance generating function $g_v(z)$ which can be factored as

```{math}
:label: eq-154
g_v(z) = d(z)d(z^{-1})
```

It follows from {eq}`eq-154` that $v_t$ can be represented either as

```{math}
:label: eq-155
v_t = d(L)a_t, \qquad Ea_t^2=1
```

or

```{math}
:label: eq-156
v_t = d(L^{-1})\tilde{a}_t, \qquad E\tilde{a}_t^2=1.
```

where $d(L)=\sum_{j=0}^\infty d_jL^j$, where $a_t$ is a white noise satisfying $a_t=v_t-P[v_t|v_{t-1},v_{t-2},\ldots]$ and $\tilde{a}_t$ is a white noise satisfying $\tilde{a}_t=v_t-P[v_t|v_{t+1},v_{t+2},\ldots]$. Let $\gamma(L)=\sum_{j=0}^\infty \gamma_jL^j=d(L)^{-1}$.

Consider the effects on {eq}`eq-153` of "backward" filtering, that is, of operating with $\gamma(L)$ on both sides of equations {eq}`eq-153`. We obtain

```{math}
:label: eq-157
\gamma(L)y_t=\sum_{j=0}^\infty h_j\gamma(L)x_{t-j} + (\gamma_0v_t + \gamma_1v_{t-1} + \cdots).
```

Note that the residual in {eq}`eq-157` is not in general orthogonal to $x_{t-j}$ for $j\geq0$. For we have

```{math}
:label: eq-158
E(\gamma(L)v_t)x_{t-j} = E[\gamma_0v_t + \gamma_1v_{t-1} + \cdots]x_{t-j} = \sum_{k=0}^\infty\gamma_kEv_{t-k}x_{t-j}
```

Since {eq}`eq-153` imposes only that $Ex_{t-j}v_t=0$ for $j\geq0$, {eq}`eq-158` becomes

$$
E(\gamma(L)v_t)x_{t-j} = \sum_{k=j+1}^\infty\gamma_kEv_{t-k}x_{t-j}.
$$

The orthogonality conditions imposed in {eq}`eq-153` do not imply that $Ev_{t-k}x_{t-j}=0$ for $k>j$. Therefore, in general, the residual in {eq}`eq-157` is not orthogonal to the set $\{x_t,x_{t-1},\ldots\}$. For a variety of linear rational expectations models, this is important because it means that $\{x_{t-j}\}_{j\geq0}$ cannot be used via an instrumental variables technique to estimate $\{h_j\}_{j=0}^\infty$ in equation {eq}`eq-157`. More precisely, equation {eq}`eq-157` does *not* imply the restrictions

$$
E(\gamma(L)y_t)x_{t-k} = \sum_{j=0}^\infty h_jE(\gamma(L)x_{t-j})x_{t-k},\quad k\geq0,
$$

which would be imposed in such an instrumental variables estimator.

The idea of Hayashi and Sims was to "forward filter" {eq}`eq-153` to obtain

```{math}
:label: eq-159
\gamma(L^{-1})y_t=\sum_{j=0}^\infty h_j\gamma(L^{-1})x_{t-j} + (\gamma_0v_t + \gamma_1v_{t+1} + \cdots).
```

The orthogonality conditions in {eq}`eq-153` do imply that the residual in {eq}`eq-159` is orthogonal to $\{x_{t-j}\}_{j=0}^\infty$. For we have

$$
E(\gamma(L^{-1})v_t)x_{t-j}=\gamma_0Ev_tx_{t-j} + \gamma_1Ev_{t+1}x_{t-j} + \gamma_2Ev_{t+2}x_{t-j} + \cdots = 0,
$$

since $Ev_tx_{t-j}=0$ for $j\geq0$. Thus, the residual in {eq}`eq-159` is a white noise $(\gamma(L^{-1})v_t)=\tilde{a}_t$, and is orthogonal to $x_{t-j}$, $j\geq0$. Hayashi and Sims propose to estimate $\{h_j\}_{j=0}^\infty$ by imposing the orthogonality conditions

$$
E(\gamma(L^{-1})y_t)x_{t-k} = \sum_{j=0}^\infty h_jE(\gamma(L^{-1})x_{t-j})x_{t-k},\quad k\geq0.
$$

Estimators along these lines can be shown to be efficient as well as consistent. The key to retaining consistency is using a filter that preserves orthogonality of the residual to the potential instruments $\{x_{t-j}\}_{j=0}^\infty$. The key to obtaining efficiency is to find a filter that manages to make the residuals serially uncorrelated.
