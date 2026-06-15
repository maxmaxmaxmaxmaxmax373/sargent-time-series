"""
Figure 7 — Factoring the "Euler equation" of a first-order moving average.

Reproduces Figure 7 from Sargent, Macroeconomic Theory (1987), Chapter XI,
Section 16, in an enhanced QuantEcon style.

To find the Wold (fundamental) representation x_t = eps_t + d_1 eps_{t-1} of a
process whose covariogram vanishes for |tau| >= 2, one solves the "Euler
equation"

    lambda + 1/lambda = c(0)/c(1).

This figure plots the left-hand side f(lambda) = lambda + 1/lambda. Because
|lambda + 1/lambda| >= 2, a real solution requires |c(0)/c(1)| >= 2, i.e.
2|c(1)| <= c(0). Each admissible ratio is met by a reciprocal pair of roots,
lambda and 1/lambda: one inside the unit circle (the fundamental choice) and one
outside.

Output: ../figures/fig7_euler_factoring.png
"""

from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

FIG_DIR = Path(__file__).resolve().parent.parent / "figures"
FIG_DIR.mkdir(exist_ok=True)

LIM = 6.4

# --- The two branches of f(lambda) = lambda + 1/lambda ---------------------
pos = np.linspace(0.16, 6, 600)
neg = np.linspace(-6, -0.16, 600)

plt.rcParams.update({"font.size": 11})
fig, ax = plt.subplots(figsize=(8.5, 8.5))

# existence region |c(0)/c(1)| >= 2 (where a real d_1 exists)
ax.axhspan(2, LIM, color="C0", alpha=0.06)
ax.axhspan(-LIM, -2, color="C0", alpha=0.06)

# the curve
ax.plot(pos, pos + 1 / pos, color="C0", lw=2)
ax.plot(neg, neg + 1 / neg, color="C0", lw=2)

# thresholds |f| = 2
ax.axhline(2, color="0.5", ls="--", lw=1)
ax.axhline(-2, color="0.5", ls="--", lw=1)

# extrema: min 2 at lambda=1, max -2 at lambda=-1
ax.plot([1, -1], [2, -2], "o", color="C3", ms=6, zorder=5)
ax.annotate(r"$(1,\,2)$", xy=(1, 2), xytext=(1.25, 2.6), fontsize=10, color="C3")
ax.annotate(r"$(-1,\,-2)$", xy=(-1, -2), xytext=(-3.3, -2.7), fontsize=10, color="C3")

# worked example: c(0)/c(1) = 2.5  ->  lambda = 0.5 and lambda = 2 (reciprocals)
yex = 2.5
ax.plot([0, 6], [yex, yex], color="C2", lw=1, ls=":")
ax.plot([0.5, 2.0], [yex, yex], "o", color="C2", ms=6, zorder=5)
for lam in (0.5, 2.0):
    ax.plot([lam, lam], [0, yex], color="C2", lw=0.8, ls=":")
ax.annotate(r"$\lambda=\frac{1}{2}$" + "\n" + r"$(|\lambda|<1)$",
            xy=(0.5, yex), xytext=(0.62, 3.5), fontsize=9, color="C2",
            arrowprops=dict(arrowstyle="->", color="C2", lw=0.8))
ax.annotate(r"$\lambda=2$" + "\n" + r"$(|\lambda|>1)$",
            xy=(2.0, yex), xytext=(2.4, 3.5), fontsize=9, color="C2",
            arrowprops=dict(arrowstyle="->", color="C2", lw=0.8))
ax.text(6.05, yex, r"$\frac{c(0)}{c(1)}=2.5$", fontsize=9, color="C2",
        va="center")

# curve label
ax.text(5.3, 5.6, r"$\lambda + \dfrac{1}{\lambda}$", fontsize=14, color="C0")

# existence-region label
ax.text(-6.2, 4.6, r"$\left|\frac{c(0)}{c(1)}\right|\geq 2$:" + "\n"
        r"real $d_1$ exists", fontsize=9, color="C0", va="center")

# --- Axes through the origin ------------------------------------------------
ax.spines[["top", "right", "left", "bottom"]].set_visible(False)
ax.annotate("", xy=(LIM, 0), xytext=(-LIM, 0), arrowprops=dict(arrowstyle="->", lw=1))
ax.annotate("", xy=(0, LIM), xytext=(0, -LIM), arrowprops=dict(arrowstyle="->", lw=1))
ax.text(LIM, -0.5, r"$\lambda$", fontsize=14)
ax.text(0.2, LIM - 0.1, r"$\dfrac{c(0)}{c(1)}$", fontsize=14, va="top")

ticks = [-6, -4, -2, 2, 4, 6]
for x in ticks:
    ax.plot([x, x], [-0.12, 0.12], color="black", lw=0.8)
    ax.text(x, -0.5, f"${x}$", ha="center", fontsize=9)
for y in ticks:
    ax.plot([-0.12, 0.12], [y, y], color="black", lw=0.8)
    ax.text(-0.35, y, f"${y}$", va="center", ha="right", fontsize=9)

ax.set_xticks([])           # remove default frame ticks; keep the manual ones
ax.set_yticks([])
ax.set_xlim(-LIM - 0.3, LIM + 1.4)
ax.set_ylim(-LIM - 0.3, LIM + 0.3)
ax.set_aspect("equal")

fig.tight_layout()
out = FIG_DIR / "fig7_euler_factoring.png"
fig.savefig(out, dpi=150, bbox_inches="tight")
print(f"saved {out}")
