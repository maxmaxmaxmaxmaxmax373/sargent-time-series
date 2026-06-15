"""
Figure 1 — Realizations of low-order stochastic difference equations.

Reproduces Figure 1 from Sargent, Macroeconomic Theory (1987), Chapter XI,
in a modern QuantEcon style:

  (a)  Y_t = 0.9 Y_{t-1} + eps_t              first-order stochastic (AR(1))
  (b)  Y_t = Y_{t-1} - 0.5 Y_{t-2} + eps_t    second-order stochastic (AR(2), complex roots)
  (c)  Y_t = Y_{t-1} - 0.5 Y_{t-2}            deterministic part, y_0 = y_1 = 1

eps_t ~ N(0, 1) i.i.d. white noise.

Output: ../figures/fig1_ar_realizations.png
"""

from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

# --- Output location -------------------------------------------------------
FIG_DIR = Path(__file__).resolve().parent.parent / "figures"
FIG_DIR.mkdir(exist_ok=True)

# --- Reproducibility -------------------------------------------------------
SEED = 42
rng = np.random.default_rng(SEED)


def simulate_ar(coeffs, n, shocks=None, init=None):
    """Recursively simulate y_t = sum_k coeffs[k] * y_{t-1-k} (+ shock_t).

    coeffs : list of AR coefficients [a_1, a_2, ...] for lags 1, 2, ...
    n      : number of periods to generate (t = 0, ..., n-1)
    shocks : length-n array of disturbances; if None, deterministic (zeros)
    init   : initial values for the first len(coeffs) periods; default zeros
    """
    p = len(coeffs)
    y = np.zeros(n)
    if shocks is None:
        shocks = np.zeros(n)
    if init is not None:
        y[:len(init)] = init
    start = 0 if init is None else len(init)
    for t in range(start, n):
        val = shocks[t]
        for k in range(p):
            if t - 1 - k >= 0:
                val += coeffs[k] * y[t - 1 - k]
        y[t] = val
    return y


# --- Panel (a): first-order stochastic, AR(1) ------------------------------
n_a = 151
eps_a = rng.standard_normal(n_a)
y_a = simulate_ar([0.9], n_a, shocks=eps_a)

# --- Panel (b): second-order stochastic, AR(2) with complex roots ----------
n_b = 151
eps_b = rng.standard_normal(n_b)
y_b = simulate_ar([1.0, -0.5], n_b, shocks=eps_b)

# --- Panel (c): deterministic part of the AR(2), y_0 = y_1 = 1 -------------
n_c = 26
y_c = simulate_ar([1.0, -0.5], n_c, shocks=None, init=[1.0, 1.0])

# --- Plot ------------------------------------------------------------------
plt.rcParams.update({
    "font.size": 11,
    "axes.spines.top": False,
    "axes.spines.right": False,
})

fig, axes = plt.subplots(3, 1, figsize=(9, 9))

# (a)
axes[0].plot(np.arange(n_a), y_a, color="C0", lw=1.0)
axes[0].axhline(0, color="black", lw=0.6, alpha=0.6)
axes[0].set_title(r"(a)  $Y_t = 0.9\,Y_{t-1} + \varepsilon_t$", loc="left")
axes[0].set_xlim(0, 150)

# (b)
axes[1].plot(np.arange(n_b), y_b, color="C1", lw=1.0)
axes[1].axhline(0, color="black", lw=0.6, alpha=0.6)
axes[1].set_title(r"(b)  $Y_t = Y_{t-1} - 0.5\,Y_{t-2} + \varepsilon_t$", loc="left")
axes[1].set_xlim(0, 150)

# (c)
axes[2].plot(np.arange(n_c), y_c, color="C2", lw=1.4, marker="o", ms=3)
axes[2].axhline(0, color="black", lw=0.6, alpha=0.6)
axes[2].set_title(r"(c)  $Y_t = Y_{t-1} - 0.5\,Y_{t-2}$ (deterministic, $y_0=y_1=1$)",
                  loc="left")
axes[2].set_xlim(0, 25)
axes[2].set_xlabel("$t$")

for ax in axes:
    ax.grid(True, alpha=0.3)

fig.tight_layout()

out_path = FIG_DIR / "fig1_ar_realizations.png"
fig.savefig(out_path, dpi=150, bbox_inches="tight")
print(f"saved {out_path}")
