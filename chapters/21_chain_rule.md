# The Chain Rule of Forecasting

The law of iterated projections implies a recursion relationship that is sometimes very useful in a forecasting context. The relationship is known as Wold's "chain rule of forecasting". It shows how projections $P_t x_{t+k}$ for all $k \geq 2$ can be calculated from knowledge of the form of $P_t x_{t+1}$ alone.[^fn-1]

Suppose that $\{x_t\}$ is a linearly indeterministic covariance stationary stochastic process for which

$$
P_t  x_{t+1} = \sum_{j=0}^\infty h_j x_{t-j}, \quad \sum_{j=0}^\infty h_j^2 < \infty
$$

It follows that

$$
P_{t+k} x_{t+k+1} = h_0 x_{t+k} + h_1 x_{t+k-1} + \ldots + h_k x_t + h_{k+1} x_{t-1} + \ldots.
$$

Projecting both sides of this equation on $(x_t, x_{t-1},\ldots)$ gives, via the law of iterated projections,

```{math}
:label: eq-92
P_{t+k} x_{t+k+1} = h_0 P_t x_{t+k} + h_1 P_t x_{t+k-1} + \ldots h_{k-1} P_t x_{t+1} + \sum_{i=0}^{\infty} h_{k+i} x_{t - i}.
```

This recursion relationship is the "chain rule of forecasting" which shows how to build up projections of $x_t$ arbitrarily far into the future from knowledge of the formula for the one-step-ahead projection alone.

To take an example, suppose that $\{x_t\}$ is a first-order Markov process so that

$$
P_t x_{t+1} = \lambda x_t, \quad |\lambda| < 1.
$$

From application of {eq}`eq-92` it follows that $P_t x_{t+j} = \lambda^j x_t$, $j \geq 1$.

[^fn-1]: Interesting applications at the chain rule of forecasting occur in Shiller (1972).
