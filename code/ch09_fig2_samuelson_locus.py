"""
Chapter IX, Figure 2 — Samuelson's multiplier-accelerator model as a path
through the (t1, t2) stability plane.

Samuelson's model reduces to  Y_t = (c+gamma) Y_{t-1} - gamma Y_{t-2} + alpha,
i.e.  t1 = c + gamma,  t2 = -gamma,  so that  t1 + t2 = c  for every gamma >= 0.
As the accelerator gamma rises from 0, the parameter point (t1, t2) slides DOWN
and to the RIGHT along the fixed line  t1 + t2 = c  (parallel to, and inside,
the stability boundary t1 + t2 = 1).

Reading where that line crosses the region boundaries of Figure 1 reads off,
at a glance, which values of gamma give monotone convergence, damped
oscillations, or explosive oscillations.  Here c = 0.8.

Output: ../figures/fig-9-2.png
"""

from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm

FIG_DIR = Path(__file__).resolve().parent.parent / "figures"
FIG_DIR.mkdir(exist_ok=True)

COLORS = ["#bcd4e6", "#cde5cd", "#f6c9c0", "#f3e2a9"]
C = 0.8                                            # marginal propensity to consume


def classify(t1, t2):
    disc = t1 ** 2 + 4 * t2
    if disc < 0:
        return 0 if np.sqrt(-t2) < 1 else 2
    sq = np.sqrt(disc)
    lam = np.array([(t1 + sq) / 2, (t1 - sq) / 2])
    if np.abs(lam).max() < 1:
        return 1
    return 2 if lam[np.argmax(np.abs(lam))] < 0 else 3


def main():
    t1 = np.linspace(-3, 3, 1000)
    t2 = np.linspace(-2.2, 1.7, 900)
    T1, T2 = np.meshgrid(t1, t2)
    Z = np.vectorize(classify)(T1, T2)

    plt.rcParams.update({"font.size": 12})
    fig, ax = plt.subplots(figsize=(10, 8))

    cmap = ListedColormap(COLORS)
    norm = BoundaryNorm([-0.5, 0.5, 1.5, 2.5, 3.5], cmap.N)
    ax.pcolormesh(T1, T2, Z, cmap=cmap, norm=norm, shading="auto", alpha=0.45)

    x = np.linspace(-3, 3, 400)
    ax.plot(x, -x ** 2 / 4, "k-", lw=1.5)
    ax.plot([-2, 0], [-1, 1], "k-", lw=1.5)
    ax.plot([0, 2], [1, -1], "k-", lw=1.5)
    ax.plot([-2, 2], [-1, -1], "k-", lw=1.5)
    ax.axhline(0, color="0.4", lw=0.8)
    ax.axvline(0, color="0.4", lw=0.8)

    # faint region labels
    ax.text(0.0, -0.6, "damped\noscillations", ha="center", va="center",
            fontsize=10, color="0.35")
    ax.text(0.0, 0.5, "damped,\nmonotone", ha="center", va="center",
            fontsize=9.5, color="0.35")
    ax.text(2.2, 1.2, "explosive\ngrowth", ha="center", va="center",
            fontsize=10, color="0.4")
    ax.text(-2.1, 1.2, "explosive\noscillations", ha="center", va="center",
            fontsize=10, color="0.4")

    # --- Samuelson locus  t1 + t2 = c  ------------------------------------
    gam = np.linspace(-0.9, 3.0, 200)
    lt1, lt2 = C + gam, -gam
    ax.plot(lt1, lt2, "--", color="C3", lw=2.4, zorder=5)

    # direction of increasing gamma (down-right)
    ax.annotate("", xy=(C + 2.4, -2.4), xytext=(C + 1.4, -1.4),
                arrowprops=dict(arrowstyle="-|>", color="C3", lw=2.4))
    ax.text(C + 2.05, -1.55, r"increasing $\gamma$", color="C3",
            rotation=-45, fontsize=11)

    # key points along the locus
    def mark(g, txt, dy=12, dx=8, ha="left"):
        p = (C + g, -g)
        ax.plot(*p, "o", color="k", ms=6, zorder=6)
        ax.annotate(txt, p, xytext=(dx, dy), textcoords="offset points",
                    fontsize=10, ha=ha,
                    bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="0.6",
                              alpha=0.9))
        return p

    mark(0.0, r"$\gamma=0$:  $(c,0)$, monotone")
    g_osc = 2 - 2 * np.sqrt(1 - C)                 # parabola crossing
    mark(g_osc - C, fr"$\gamma\approx{g_osc - C:.2f}$:  oscillations begin",
         dy=14, dx=10)
    mark(1.0, r"$\gamma=1$:  instability ($t_2=-1$)", dy=-22, dx=10)

    ax.text(-1.65, 0.95, r"$t_1+t_2=c$", color="C3", fontsize=13, rotation=-45)

    ax.set_xlim(-3, 3)
    ax.set_ylim(-2.2, 1.7)
    ax.set_xlabel("$t_1 = c+\\gamma$", fontsize=13)
    ax.set_ylabel("$t_2 = -\\gamma$", fontsize=13)
    ax.set_title("Samuelson's multiplier-accelerator: the locus "
                 "$t_1+t_2=c$ ($c=0.8$)", fontsize=13)

    fig.tight_layout()
    out = FIG_DIR / "fig-9-2.png"
    fig.savefig(out, dpi=140)
    print("wrote", out)
    print(f"oscillations begin at gamma = {g_osc - C:.3f}; "
          f"instability at gamma = 1")


if __name__ == "__main__":
    main()
