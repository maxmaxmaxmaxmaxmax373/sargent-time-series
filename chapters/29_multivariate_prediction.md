# Multivariate Prediction Formulas

Continue to assume that $(x_t, y_t)$ is a jointly covariance stationary, strictly indeterministic process with a moving average representation

$$
\begin{bmatrix} x_t \\ y_t \end{bmatrix} =
\begin{bmatrix}
    C^{11}(L) & C^{12}(L) \\
    C^{21}(L) & C^{22}(L)
\end{bmatrix}
\begin{bmatrix} \epsilon_t \\ u_t \end{bmatrix}
= C(L)\begin{bmatrix} \epsilon_t \\ u_t \end{bmatrix}
$$

where $E\epsilon_t u_s=0$ for all $t,s$, $\{\epsilon_t,u_t\}$ are jointly fundamental for $(x_t,y_t)$, and where $C(L)^{-1}$ exists and is one-sided and convergent in nonnegative powers of $L$, so that $(x_t,y_t)$ has an autoregressive representation

$$
C(L)^{-1}\begin{bmatrix} x_t \\ y_t \end{bmatrix} =
\begin{bmatrix} \epsilon_t \\ u_t \end{bmatrix}
\quad\text{or}\quad A(L)\begin{bmatrix} x_t \\ y_t \end{bmatrix}=\begin{bmatrix} \epsilon_t \\ u_t \end{bmatrix}
$$

where $A(L)=C(L)^{-1}$. Paralleling our calculations in the univariate case it is easy to deduce that the projection of $(x_{t+1},y_{t+1})$ against $\{x_t,x_{t-1},\ldots,y_t,y_{t-1},\ldots\}$, call it $P_t\begin{bmatrix}x_{t+1}\\y_{t+1}\end{bmatrix}$, is

$$
P_t\begin{bmatrix}x_{t+1}\\y_{t+1}\end{bmatrix}=\left(\frac{C(L)}{L}\right)_+\begin{bmatrix}\epsilon_t\\u_t\end{bmatrix}=\left(\frac{C(L)}{L}\right)_+A(L)\begin{bmatrix}x_t\\y_t\end{bmatrix}
$$

More generally, we have

$$
P_t\begin{bmatrix}x_{t+j}\\y_{t+j}\end{bmatrix}=\left(\frac{C(L)}{L^j}\right)_+A(L)\begin{bmatrix}x_t\\y_t\end{bmatrix}
$$

These results extend in a natural way to $n$-dimensional stochastic processes. In particular, the $n$-variate version of Wold's theorem implies that if $\{y_t\}$ is an $n$-dimensional, jointly covariance stationary, strictly indeterministic stochastic process with mean zero, it has a moving average representation

```{math}
:label: eq-138
y_t=C(L)\epsilon_t
```

where $C(L)=C_0 + C_1L + \cdots$, $C_j$ being an $n \times n$ matrix and the $C_j$ being "square summable," where $\epsilon_t$ is an $(n \times 1)$-vector stochastic process, where the components $\epsilon_{it}$ are serially uncorrelated and mutually orthogonal (at all lags), $E\epsilon_{it}\epsilon'_{js}=0$ for all $t,s$ where $i \neq j$; and the $\epsilon_{it}$ are "jointly fundamental for $y_t$," i.e., for each $i$, $y_{it} - P(y_{it}|y_{t-1},y_{t-2},\ldots)$ is a linear combination of $\epsilon_{jt}$, $j=1, \ldots,n$. For the process {eq}`eq-138`, we have the prediction formula

$$
E_ty_{t+j} =\left[C(L)/L^j\right]_+\epsilon_t
$$

where $E_t(x) \equiv E[x|y_t,y_{t-1}, \ldots]$. Where $C(L)^{-1}$ exists, so that $y_t$ has a vector autoregressive representation, then we also have the formula

$$
E_ty_{t+j} =\left[C(L)/L^j\right]_+C(L)^{-1}y_t
$$

To take an example, let $R_{nt}$ be the rate of $n$-period bonds, and assume that $(R_{nt},R_{1t})$ has moving average representation

```{math}
:label: eq-139
R_{nt} = \alpha(L)\epsilon_t + \beta(L)u_t,\quad n>1,\qquad R_{1t} = \gamma(L)\epsilon_t + \delta(L)u_t
```

where all lag operators are one-sided on the present and past, and

$$
R_{nt} - P_{t-1}R_{nt}= \alpha_0\epsilon_t + \beta_0u_t,\qquad R_{1t} - P_{t-1}R_{1t} = \gamma_0\epsilon_t + \delta_0u_t
$$

The rational expectations theory of the term structure asserts that

$$
R_{nt} = \frac{1}{n}\left[R_{1t}+P_tR_{1,t+1}+\ldots+P_tR_{1,t+n-1}\right] = \frac{1}{n}\left[\gamma(L) + \frac{\gamma(L)}{L} + \ldots + \frac{\gamma(L)}{L^{n-1}}\right]_{+}\epsilon_t + \frac{1}{n}\left[\delta(L) + \frac{\delta(L)}{L} + \ldots + \frac{\delta(L)}{L^{n-1}}\right]_{+}u_t
$$

or

```{math}
:label: eq-140
R_{nt} = \frac{1}{n}\left[\left[\frac{1-L^{-n}}{1-L^{-1}}\right]\gamma(L)\right]_+\epsilon_t + \frac{1}{n}\left[\left[\frac{1-L^{-n}}{1-L^{-1}}\right]\delta(L)\right]_+u_t.
```

Thus, comparing {eq}`eq-139` with {eq}`eq-140`, it is seen that the rational expectations theory of the term structure imposes the following restrictions across the equations of the moving average representation of the $(R_{nt},R_{1t})$ process:

$$
\alpha(L) = \frac{1}{n}\left[\left[\frac{1-L^{-n}}{1-L^{-1}}\right]\gamma(L)\right]_+, \qquad \beta(L)= \frac{1}{n}\left[\left[\frac{1-L^{-n}}{1-L^{-1}}\right]\delta(L)\right]_+.
$$

These restrictions embody the content of the theory and are refutable. Hansen and Sargent (1981c) describe how to impose and test these restrictions econometrically.

Further applications of multivariate prediction formulas to rational expectations models occur in Whiteman (1983, 1984), Hansen and Sargent (1981c, 1981a), and Dagli and Taylor (1984).
