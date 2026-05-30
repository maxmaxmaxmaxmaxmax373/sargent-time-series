# Some Examples[^fn-ex-1]

#### First-Order-Markov:

Consider the first-order autoregressive process $(1-\lambda L)x_t = \epsilon_t$, $\epsilon_t$ white noise, $|\lambda| < 1$, $\epsilon_t = x_t - P[x_t|x_{t-1}, x_{t-2},\ldots]$; we can write $x_t = (1/(1-\lambda L))\epsilon_t$. We have

$$
\begin{aligned}
P_{t-1}x_t &= \left[L^{-1}(1 + \lambda L + \lambda^2 L^2 + \ldots)\right]_{+}(1-\lambda L)x_{t-1}\\
            &= ( \lambda + \lambda^2 L + \ldots)(1 - \lambda L)x_{t-1}\\
            &= \left(\frac{\lambda}{1 - \lambda L}\right)(1-\lambda L)x_{t-1} = \lambda x_{t-1}
\end{aligned}
$$

More generally,

$$
P_{t-k}x_t = \left[L^{-k}(1 + \lambda L + \ldots)\right]_{+}(1-\lambda L)x_{t-k} = \lambda^k x_{t-k}.
$$

Thus we have

$$
P^{t-k}x_{t+k} = \lambda^k x_t.
$$

#### First-Order Moving Average:

Suppose $x_t = (1 + \beta L)\epsilon_t$, $\epsilon_t$ white, $|\beta| < 1$, $\epsilon_t = x_t - P[x_t|x_{t-1}, x_{t-2},\ldots]$. Then we have

$$
P_{t-1}x_t = \left[L^{-1}(1 + \beta L)\right]_{+}\left(\frac{1}{1 + \beta L}\right)x_{t-1}, \qquad P_{t-1}x_t = \frac{\beta}{1 + \beta L}x_{t-1}.
$$

We also have that for $k \geq 2$,

$$
P_{t-k}x_t = \left[L^{-k}(1 + \beta L)\right]_{+}\left(\frac{1}{1 + \beta L}\right)x_{t-k} = 0,
$$

which can also be seen directly by projecting on $\{x_{t-k}, x_{t-k-1},\ldots\}$ both sides of $x_t = (1 + \beta L)\epsilon_t$.

#### First-Order Moving Average, Autoregressive:

Suppose we have

$$
x_t  = \left(\frac{1 + a L}{1 - \beta L}\right)\epsilon_t, \qquad \epsilon_t \text{ white},\quad |a|<1, |\beta| < 1,\\
\epsilon_t = x_t - P[x_t|x_{t-1}, x_{t-2},\ldots]
$$

We then have

$$
\begin{aligned}
P_{t-1}x_t &= \left(\frac{L^{-1}(1 + a L)}{(1 - \beta L)}\right)_{+}\left(\frac{1 - \beta L}{1 + a L}\right)x_{t-1}\\
            &= \left(\frac{L^{-1}}{1 - \beta L} + \frac{a}{1 - \beta L}\right)_{+}\left(\frac{1-\beta L}{1 + a L}\right)x_{t-1}\\
P_{t-1}x_t &= \left(\frac{\beta + a}{1 - \beta L}\right)\left(\frac{1 - \beta L}{1 + a L}\right)x_{t-1}, \qquad P_{t-1}x_t = \left\{\frac{a + \beta}{1 + a L}\right\}x_{t-1},
\end{aligned}
$$

which expresses the forecast of $x_t$ as a geometric distributed lag of past $x$'s. The first-order mixed moving average, autoregressive model for $x_t$ thus provides a rationalization for the familiar "adaptive expectations" model. As we let $\beta \to 1$ (from below, in order to assure that the roots condition $|\beta| < 1$ is met), $P_{t-1}x_t$ approaches

$$
P_{t-1}x_t = \{(1 + a)/(1 + a L)\}x_{t-1},
$$

which with $a < 0$ is equivalent with Cagan's (1956) adaptive expectations scheme

$$
P_{t-1}x_t = \{(1 - \lambda )/(1 + \lambda L)\}x_{t-1},
$$

with $a = - \lambda$. Notice that as $\beta \to 1$ (from below), we approach the situation in which $(1-L)x_t = (1 + a L)\epsilon_t$, so that the first difference of $x_t$ follows a first-order moving average. The parameter $a$ must be negative in order that $\lambda > 0$

For the general case in which $k \geq 1$, we have

$$
\begin{aligned}
P_{t-k}x_t &= \left(\frac{L^{-k}(1+aL)}{(1-\beta L)}\right)_{+}\left(\frac{1-\beta L}{1+aL}\right)x_{t-k}\\
            &=\left(\frac{\beta^k}{1-\beta L} + \frac{a \beta^{k-1}}{1-\beta L}\right)_{+}\left(\frac{1-\beta L}{1+aL}\right)x_{t-k}\\
            &=\left(\frac{\beta^k}{1-\beta L} + \frac{a \beta^{k-1}}{1-\beta L}\right)\left(\frac{1-\beta L}{1+aL}\right)x_{t-k} = \frac{\beta^{k-1}(\beta + a)}{(1+aL)} x_{t-k}.
\end{aligned}
$$

We can write this alternatively as

$$
P_t x_{t+k} = (\beta^{k-1}(\beta + a)/(1 + a L))x_t
$$

Notice that as $\beta \to 1$ (from below) we approach the situation in which

$$
P_t x_{t+k} = ((1 + a)/(1 + a L))x_t,
$$

so that the same forecast is made for all horizons $k \geq 1$. In this sense there is a well-defined concept of "permanent $x$." This was first pointed out in the economics literature by Muth (1960), who showed that the hypothesis of rational expectations in conjunction with the model for income $(1 - L)x_t = (1+aL)\epsilon_t$, provides a rationalization both for the concept of permanent income and the geometric distributed lag formula that Friedman had earlier used to estimate permanent income in empirical work.

[^fn-ex-1]: In these examples we continue to assume that $Ex_t = E\epsilon_t = 0$. Modifying the formulas to account for a nonzero mean of $x_t$ is trivial and involves adding constant terms to the formulas. Many more examples are worked out in Whittle (1983) and Nerlove (1967).
