"""
Chapter IX, Figure 3 — Realizations of the second-order difference equation
Y_t = t1 Y_{t-1} + t2 Y_{t-2} for various (t1, t2).

Each panel iterates the homogeneous equation from the initial conditions
Y_0 = 0, Y_1 = 1 (so every path starts from the same unit "impulse") and is
annotated with its characteristic roots and the behavior type, letting the
reader check each against the region it occupies in Figure 1.  The eight cases
sweep the taxonomy: damped oscillations, damped monotone decay, an explosive
oscillation, and explosive growth.

Output: ../figures/fig-9-3.png
"""

from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

FIG_DIR = Path(__file__).resolve().parent.parent / "figures"
FIG_DIR.mkdir(exist_ok=True)

N = 26                                              # t = 0,1,...,25


def simulate(t1, t2, n=N, y0=0.0, y1=1.0):
    y = np.zeros(n)
    y[0], y[1] = y0, y1
    for t in range(2, n):
        y[t] = t1 * y[t - 1] + t2 * y[t - 2]
    return y


def describe(t1, t2):
    """Return (color, behavior label, root string)."""
    disc = t1 ** 2 + 4 * t2
    if disc < 0:
        r = np.sqrt(-t2)
        w = np.arccos(t1 / (2 * np.sqrt(-t2)))
        period = 2 * np.pi / w
        kind = "damped osc." if r < 1 else "explosive osc."
        col = "C0" if r < 1 else "C3"
        return col, kind, fr"$|\lambda|={r:.2f}$, period $\approx{period:.1f}$"
    sq = np.sqrt(disc)
    l1, l2 = (t1 + sq) / 2, (t1 - sq) / 2
    dom = l1 if abs(l1) >= abs(l2) else l2
    if abs(dom) < 1:
        kind, col = "damped monotone", "C2"
    else:
        kind = "explosive growth" if dom > 0 else "explosive osc."
        col = "C1" if dom > 0 else "C3"
    return col, kind, fr"$\lambda={l1:.2f},\,{l2:.2f}$"


# (t1, t2) cases, sweeping the taxonomy of Figure 1
CASES = [
    (1.0, -0.75),     # damped oscillation
    (0.0, -1.5),      # explosive oscillation
    (1.0, -0.5),      # damped oscillation
    (0.5, 0.3),       # damped, monotone (real roots)
    (-0.5, -0.3),     # damped oscillation (t1<0)
    (1.0, 0.5),       # explosive growth
    (1.2, -0.95),     # slowly damped oscillation (near unit circle)
    (1.5, -0.56),     # damped monotone, two positive real roots (0.8, 0.7)
]


def main():
    plt.rcParams.update({"font.size": 10})
    fig, axes = plt.subplots(4, 2, figsize=(12, 12))
    t = np.arange(N)

    for ax, (t1, t2) in zip(axes.ravel(), CASES):
        y = simulate(t1, t2)
        col, kind, roots = describe(t1, t2)
        ax.axhline(0, color="0.7", lw=0.7)
        ax.plot(t, y, "-", color=col, lw=1.6)
        ax.plot(t, y, ".", color=col, ms=4)
        s2 = f"+ {t2}" if t2 >= 0 else f"- {abs(t2)}"
        s1 = f"{t1}" if t1 != 1 else ""
        title = fr"$Y_t = {s1}Y_{{t-1}} {s2}Y_{{t-2}}$"
        if t1 == 0:
            title = fr"$Y_t = {s2}Y_{{t-2}}$".replace("+ ", "")
        ax.set_title(title, fontsize=12)
        ax.set_xlim(0, N - 1)
        ax.set_xlabel("$t$", fontsize=10)
        # annotation box: behavior + roots
        ax.text(0.97, 0.93, f"{kind}\n{roots}", transform=ax.transAxes,
                ha="right", va="top", fontsize=9,
                bbox=dict(boxstyle="round,pad=0.3", fc="white", ec=col,
                          alpha=0.9))

    fig.suptitle("Realizations of $Y_t = t_1 Y_{t-1} + t_2 Y_{t-2}$ "
                 "(initial conditions $Y_0=0,\\ Y_1=1$)", fontsize=14)
    fig.tight_layout(rect=(0, 0, 1, 0.975))
    out = FIG_DIR / "fig-9-3.png"
    fig.savefig(out, dpi=130)
    print("wrote", out)
    for (t1, t2) in CASES:
        print(f"  ({t1:>4},{t2:>5}): {describe(t1, t2)[1]}")


if __name__ == "__main__":
    main()
