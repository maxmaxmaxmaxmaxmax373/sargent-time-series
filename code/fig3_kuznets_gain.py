"""
Figure 3 — Gain function of Kuznets's two-step transformation.

Reproduces Figure 3 from Sargent, Macroeconomic Theory (1987), Chapter XI,
Section 9, in an enhanced QuantEcon style.

Kuznets first took a five-year moving average A(L), then a centered five-year
difference B(L). The squared gain applied to the spectrum of the input is

    G(omega) = 2 * (1/5)^2 * (1 - cos 5w)(1 - cos 10w) / (1 - cos w),
    omega in [0, pi].

Even if the input X_t is white noise (flat spectrum), G(omega) has a large
peak at a low frequency, so the transformed series y_t = A(L)B(L) X_t appears
to contain "long swings." These long swings are a statistical artifact of the
transformation, not a feature of the economic system. With annual data the
main peak corresponds to a cycle of about 20 years -- close to the Kuznets
cycle found in the data.

Output: ../figures/fig3_kuznets_gain.png
"""

from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

FIG_DIR = Path(__file__).resolve().parent.parent / "figures"
FIG_DIR.mkdir(exist_ok=True)


def gain(w):
    """Kuznets squared-gain G(omega); G(0) = 0 by the limit."""
    num = (1 - np.cos(5 * w)) * (1 - np.cos(10 * w))
    den = 1 - np.cos(w)
    with np.errstate(divide="ignore", invalid="ignore"):
        out = np.where(den == 0.0, 0.0, 2 * (1 / 5) ** 2 * num / den)
    return out


# --- Evaluate --------------------------------------------------------------
w = np.linspace(0, np.pi, 4000)
G = gain(w)

# Main peak (in the first lobe) and its implied cycle length for annual data
peak_idx = np.argmax(G)
w_peak = w[peak_idx]
G_peak = G[peak_idx]
period = 2 * np.pi / w_peak           # cycle length in periods (years, annual data)

# --- Plot ------------------------------------------------------------------
plt.rcParams.update({
    "font.size": 11,
    "axes.spines.top": False,
    "axes.spines.right": False,
})

fig, ax = plt.subplots(figsize=(9, 6))

# shade the main (spurious low-frequency) lobe up to the first zero at pi/5
lobe = w <= np.pi / 5
ax.fill_between(w[lobe], G[lobe], color="C0", alpha=0.15)

ax.plot(w, G, color="C0", lw=1.8)
ax.axhline(0, color="black", lw=0.6)

# mark the main peak and translate it to a cycle length
ax.plot([w_peak], [G_peak], "o", color="C3", ms=5)
ax.annotate(
    f"main peak at $\\omega \\approx {w_peak:.2f}$\n"
    f"$\\Rightarrow$ cycle $\\approx$ {period:.0f} years (annual data)\n"
    f"$\\approx$ the 20-year Kuznets cycle",
    xy=(w_peak, G_peak), xytext=(w_peak + 0.45, G_peak - 0.25),
    fontsize=10,
    arrowprops=dict(arrowstyle="->", color="C3", lw=1.2),
)

# note on the economic meaning
ax.text(np.pi / 5, G_peak * 0.62,
        'spurious "long swings":\nlow-frequency power created\nby the filter, not the data',
        fontsize=9, color="0.35", va="center")

# the transfer-function formula, as in the original figure
ax.text(0.97, 0.93,
        r"$G(\omega) = \left(\frac{1}{5}\right)^2\,"
        r"\frac{(1-\cos 5\omega)\cdot 2(1-\cos 10\omega)}{1-\cos\omega}$",
        transform=ax.transAxes, ha="right", va="top", fontsize=12,
        bbox=dict(facecolor="white", edgecolor="0.8", boxstyle="round,pad=0.4"))

# x-axis: zeros at multiples of pi/5
ticks = [np.pi * k / 5 for k in range(1, 6)]
labels = [r"$\frac{\pi}{5}$", r"$\frac{2\pi}{5}$", r"$\frac{3\pi}{5}$",
          r"$\frac{4\pi}{5}$", r"$\pi$"]
ax.set_xticks(ticks)
ax.set_xticklabels(labels, fontsize=12)

ax.set_xlim(0, np.pi)
ax.set_ylim(0, G_peak * 1.08)
ax.set_xlabel(r"frequency $\omega$")
ax.set_ylabel(r"gain $G(\omega)$")
ax.grid(True, alpha=0.3)

fig.tight_layout()
out_path = FIG_DIR / "fig3_kuznets_gain.png"
fig.savefig(out_path, dpi=150, bbox_inches="tight")
print(f"saved {out_path}")
print(f"main peak: omega={w_peak:.4f}, G={G_peak:.3f}, period={period:.2f}")
