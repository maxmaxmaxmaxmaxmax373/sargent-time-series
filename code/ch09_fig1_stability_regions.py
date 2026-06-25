"""
Chapter IX, Figure 1 — Stability regions of the second-order difference equation
in the (t1, t2) parameter plane.

For Y_t = t1 Y_{t-1} + t2 Y_{t-2}, the characteristic roots are

    lambda = ( t1 +/- sqrt(t1^2 + 4 t2) ) / 2 .

The (t1, t2) plane is carved into four behavior types by two curves:

  * the parabola  t1^2 + 4 t2 = 0  (i.e. t2 = -t1^2/4) separates real roots
    (above) from complex-conjugate roots (below);
  * the "stability triangle" with edges  t1 + t2 = 1  (right),
    t2 = 1 + t1  (left), and  t2 = -1  (bottom), vertices (0,1), (2,-1),
    (-2,-1), inside which both roots are < 1 in modulus.

This redraws Baumol's classic diagram (the source of the chapter's Figure 1)
by classifying a fine grid and shading the four regions, then overlaying the
boundary curves and the equations that define them.

Output: ../figures/fig-9-1.png
"""

from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm

FIG_DIR = Path(__file__).resolve().parent.parent / "figures"
FIG_DIR.mkdir(exist_ok=True)

# 0 = damped oscillation, 1 = damped monotone, 2 = explosive oscillation,
# 3 = explosive growth
COLORS = ["#bcd4e6", "#cde5cd", "#f6c9c0", "#f3e2a9"]


def classify(t1, t2):
    disc = t1 ** 2 + 4 * t2
    if disc < 0:                                   # complex pair
        r = np.sqrt(-t2)                           # modulus
        return 0 if r < 1 else 2
    sq = np.sqrt(disc)
    lam = np.array([(t1 + sq) / 2, (t1 - sq) / 2])
    dom = lam[np.argmax(np.abs(lam))]              # dominant root
    if np.abs(lam).max() < 1:
        return 1                                   # both inside: monotone-ish
    return 2 if dom < 0 else 3                      # explode: sign decides type


def main():
    t1 = np.linspace(-3, 3, 1200)
    t2 = np.linspace(-2.2, 1.7, 1000)
    T1, T2 = np.meshgrid(t1, t2)
    Z = np.vectorize(classify)(T1, T2)

    plt.rcParams.update({"font.size": 12})
    fig, ax = plt.subplots(figsize=(10, 8))

    cmap = ListedColormap(COLORS)
    norm = BoundaryNorm([-0.5, 0.5, 1.5, 2.5, 3.5], cmap.N)
    ax.pcolormesh(T1, T2, Z, cmap=cmap, norm=norm, shading="auto", alpha=0.85)

    # boundary curves -------------------------------------------------------
    x = np.linspace(-3, 3, 400)
    ax.plot(x, -x ** 2 / 4, "k-", lw=2)                     # parabola
    tri_x = np.linspace(-2, 2, 10)
    ax.plot([-2, 0], [-1, 1], "k-", lw=2)                   # t2 = 1 + t1
    ax.plot([0, 2], [1, -1], "k-", lw=2)                    # t1 + t2 = 1
    ax.plot([-2, 2], [-1, -1], "k-", lw=2)                  # t2 = -1

    # axes through the origin
    ax.axhline(0, color="0.4", lw=0.8)
    ax.axvline(0, color="0.4", lw=0.8)

    # vertices and apex
    for (px, py) in [(0, 1), (2, -1), (-2, -1)]:
        ax.plot(px, py, "ko", ms=4)
    ax.annotate("$(0,1)$", (0, 1), xytext=(8, 4),
                textcoords="offset points", fontsize=10)
    ax.annotate("$(2,-1)$", (2, -1), xytext=(8, 4),
                textcoords="offset points", fontsize=10)
    ax.annotate("$(-2,-1)$", (-2, -1), xytext=(-58, 4),
                textcoords="offset points", fontsize=10)

    # region labels (placed inside the correct region)
    ax.text(0.0, -0.62, "Damped\noscillations", ha="center", va="center",
            fontsize=12.5, weight="bold")
    ax.text(0.0, 0.46, "Damped,\nmonotone", ha="center", va="center",
            fontsize=10.5, weight="bold")
    ax.text(2.05, 1.15, "Explosive\ngrowth", ha="center", va="center",
            fontsize=12.5, weight="bold")
    ax.text(-2.05, 1.15, "Explosive\noscillations", ha="center", va="center",
            fontsize=12.5, weight="bold")
    ax.text(0.0, -1.75, "Explosive oscillations", ha="center", va="center",
            fontsize=12.5, weight="bold")

    # curve equations
    ax.text(1.62, -0.78, "$t_1+t_2=1$", rotation=-45, fontsize=11, color="0.15")
    ax.text(-2.35, -0.78, "$t_2=1+t_1$", rotation=45, fontsize=11, color="0.15")
    ax.text(0.92, -0.30, "$t_1^2+4t_2=0$", rotation=-24, fontsize=11,
            color="0.2")
    ax.text(-1.25, -0.92, "$t_2=-1$", fontsize=11, color="0.15")

    ax.set_xlim(-3, 3)
    ax.set_ylim(-2.2, 1.7)
    ax.set_xlabel("$t_1$  (sum of roots $\\lambda_1+\\lambda_2$)", fontsize=13)
    ax.set_ylabel("$t_2$  (minus the product $-\\lambda_1\\lambda_2$)", fontsize=13)
    ax.set_title("Behavior of $Y_t = t_1 Y_{t-1} + t_2 Y_{t-2}$ in the "
                 "$(t_1,t_2)$ plane", fontsize=13)

    fig.tight_layout()
    out = FIG_DIR / "fig-9-1.png"
    fig.savefig(out, dpi=140)
    print("wrote", out)


if __name__ == "__main__":
    main()
