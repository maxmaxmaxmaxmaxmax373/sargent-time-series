"""
Chapter IX, Figure 4 — Determining the factorization roots of the firm's Euler
equation.

Factoring the characteristic polynomial of the employment Euler equation
requires solving

    -phi  =  b*lambda + 1/lambda ,      0 < b < 1,

for the two roots lambda_1 < lambda_2.  The right-hand side b*lambda + 1/lambda
is U-shaped, with a minimum value 2*sqrt(b) at lambda = 1/sqrt(b).  A horizontal
line at height -phi cuts it twice.

The figure shows that because -phi = f1/d + (1+b) > 1+b (whenever f1 > 0), and
because the line at height 1+b meets the curve exactly at lambda = 1 and
lambda = 1/b, the two roots straddle the unit circle and 1/b:

    lambda_1 < 1 < 1/b < lambda_2 .

The smaller root lambda_1 < 1 is the stable one used to build the firm's
optimal policy; the larger, lambda_2 > 1/b, is unstable and is handled by
solving forward.  Here b = 0.7 and -phi = 2.0.

Output: ../figures/fig-9-4.png
"""

from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

FIG_DIR = Path(__file__).resolve().parent.parent / "figures"
FIG_DIR.mkdir(exist_ok=True)

B = 0.7
NEG_PHI = 2.0


def main():
    lam = np.linspace(0.3, 3.0, 600)
    U = B * lam + 1.0 / lam

    # roots of b*lam + 1/lam = -phi  <=>  b lam^2 + phi lam + 1 = 0
    disc = NEG_PHI ** 2 - 4 * B
    l1 = (NEG_PHI - np.sqrt(disc)) / (2 * B)
    l2 = (NEG_PHI + np.sqrt(disc)) / (2 * B)
    lmin = 1 / np.sqrt(B)
    umin = 2 * np.sqrt(B)
    invb = 1 / B

    plt.rcParams.update({"font.size": 12})
    fig, ax = plt.subplots(figsize=(10, 7.5))

    ax.plot(lam, U, color="C0", lw=2.6, label=r"$b\lambda + 1/\lambda$")
    ax.plot(lam, B * lam, color="C2", lw=1.8, label=r"$b\lambda$")
    ax.plot(lam, 1.0 / lam, color="C1", lw=1.8, label=r"$1/\lambda$")

    # horizontal reference lines
    ax.axhline(NEG_PHI, color="C3", lw=1.6, ls="--")
    ax.axhline(1 + B, color="0.5", lw=1.2, ls=":")
    ax.text(0.33, NEG_PHI + 0.06, r"$-\phi$", color="C3", fontsize=13)
    ax.text(2.62, 1 + B - 0.16, r"$1+b$", color="0.4", fontsize=12)

    # intersection points lambda_1, lambda_2 (on the -phi line)
    for lr in (l1, l2):
        ax.plot(lr, NEG_PHI, "o", color="C3", ms=8, zorder=5)
        ax.vlines(lr, 0, NEG_PHI, color="C3", ls="--", lw=1.0)
    # minimum of the U-curve
    ax.plot(lmin, umin, "o", color="C0", ms=7, zorder=5)
    ax.vlines(lmin, 0, umin, color="C0", ls=":", lw=1.0)
    # where the 1+b line meets the curve: lambda = 1 and lambda = 1/b
    for lr in (1.0, invb):
        ax.plot(lr, 1 + B, "s", color="0.4", ms=6, zorder=5)
        ax.vlines(lr, 0, 1 + B, color="0.5", ls=":", lw=0.8)

    # x-axis tick labels at the special abscissae
    ax.set_xticks([l1, 1.0, lmin, invb, l2])
    ax.set_xticklabels(
        [fr"$\lambda_1$" "\n" fr"${l1:.2f}$", "$1$",
         r"$1/\sqrt{b}$", "$1/b$",
         fr"$\lambda_2$" "\n" fr"${l2:.2f}$"], fontsize=11)

    ax.annotate(r"$\lambda_1 < 1 < 1/b < \lambda_2$",
                xy=(1.55, 3.1), fontsize=14)
    ax.text(2.35, 0.55, f"$b = {B}$", fontsize=13,
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="0.7"))

    # curve labels near their right ends
    ax.text(2.5, B * 2.5 + 0.07, r"$b\lambda$", color="C2", fontsize=13)
    ax.text(2.62, 1 / 2.95 + 0.02, r"$1/\lambda$", color="C1", fontsize=13)
    ax.text(2.18, B * 2.18 + 1 / 2.18 + 0.08, r"$b\lambda + 1/\lambda$",
            color="C0", fontsize=13, rotation=24)

    ax.set_xlim(0.3, 3.0)
    ax.set_ylim(0, 3.5)
    ax.set_xlabel(r"$\lambda$", fontsize=14)
    ax.set_title(r"Solving $-\phi = b\lambda + 1/\lambda$ for the roots "
                 r"$\lambda_1,\lambda_2$", fontsize=13)
    ax.legend(loc="upper left", fontsize=11, framealpha=0.95)

    fig.tight_layout()
    out = FIG_DIR / "fig-9-4.png"
    fig.savefig(out, dpi=140)
    print("wrote", out)
    print(f"  lambda_1 = {l1:.4f}, lambda_2 = {l2:.4f}, "
          f"1/sqrt(b) = {lmin:.4f}, 1/b = {invb:.4f}, min = {umin:.4f}")


if __name__ == "__main__":
    main()
