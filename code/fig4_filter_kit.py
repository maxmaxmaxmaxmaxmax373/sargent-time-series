"""
Figure 4 — A small kit of frequency response functions h(e^{-i omega}).

Reproduces Figure 4 from Sargent, Macroeconomic Theory (1987), Chapter XI,
Section 10, regrouped by filter family (enhanced/teaching layout):

    fig4a_ma_filters.png        moving-average / differencing filters
    fig4b_ar_filters.png        autoregressive (inverse) filters
    fig4c_seasonal_filters.png  seasonal L^12 filters

For a filter h(L) = num(L) / den(L), the frequency response is obtained by
setting L = e^{-i omega}:

    h(e^{-i omega}) = sum_k num_k e^{-i omega k} / sum_k den_k e^{-i omega k}

We plot the amplitude |h(e^{-i omega})| and the phase arg h(e^{-i omega}) for
omega in [0, pi].
"""

from math import ceil
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

FIG_DIR = Path(__file__).resolve().parent.parent / "figures"
FIG_DIR.mkdir(exist_ok=True)

W = np.linspace(0, np.pi, 1000)


def evalpoly(coeffs, w):
    """Evaluate sum_k coeffs[k] e^{-i w k} (coeffs index = power of L)."""
    coeffs = np.asarray(coeffs, dtype=float)
    k = np.arange(coeffs.size)
    return (coeffs[:, None] * np.exp(-1j * np.outer(k, w))).sum(axis=0)


def response(num, den, w=W):
    """Frequency response h(e^{-i w}) = num(e^{-i w}) / den(e^{-i w})."""
    return evalpoly(num, w) / evalpoly(den, w)


def seasonal(coef0, coef12):
    """Coefficient array for  coef0 + coef12 * L^12."""
    c = np.zeros(13)
    c[0] = coef0
    c[12] = coef12
    return c


# --- Filter families (label is rendered as mathtext) -----------------------
ma_filters = [
    (r"$1 - L$",                 [1, -1],          [1]),
    (r"$1 + L$",                 [1, 1],           [1]),
    (r"$(1 - L)^2$",             [1, -2, 1],       [1]),
    (r"$1 + L + L^2$",           [1, 1, 1],        [1]),
    (r"$1 + L + L^2 + L^3$",     [1, 1, 1, 1],     [1]),
    (r"$1 - 0.5L$",              [1, -0.5],        [1]),
    (r"$1 - 0.9L$",              [1, -0.9],        [1]),
    (r"$1 + 0.5L$",              [1, 0.5],         [1]),
    (r"$1 - L + 0.8L^2$",        [1, -1, 0.8],     [1]),
    (r"$1 - 0.3L - 0.5L^2$",     [1, -0.3, -0.5],  [1]),
]

ar_filters = [
    (r"$(1 - 0.5L)^{-1}$",            [1], [1, -0.5]),
    (r"$(1 - 0.9L)^{-1}$",            [1], [1, -0.9]),
    (r"$(1 + 0.5L)^{-1}$",            [1], [1, 0.5]),
    (r"$(1 - L + 0.8L^2)^{-1}$",      [1], [1, -1, 0.8]),
    (r"$(1 - 0.3L - 0.8L^2)^{-1}$",   [1], [1, -0.3, -0.8]),
    (r"$(1 - 0.3L - 0.5L^2)^{-1}$",   [1], [1, -0.3, -0.5]),
]

seasonal_filters = [
    (r"$1 - 0.5L^{12}$",        seasonal(1, -0.5), [1]),
    (r"$1 - 0.9L^{12}$",        seasonal(1, -0.9), [1]),
    (r"$(1 - 0.5L^{12})^{-1}$", [1], seasonal(1, -0.5)),
    (r"$(1 - 0.9L^{12})^{-1}$", [1], seasonal(1, -0.9)),
]

XTICKS = [0, np.pi / 2, np.pi]
XLABELS = ["0", r"$\pi/2$", r"$\pi$"]
PHASE_TICKS = [-np.pi, 0, np.pi]
PHASE_LABELS = [r"$-\pi$", "0", r"$\pi$"]


def make_figure(filters, fname, suptitle, pairs_per_row=2):
    n = len(filters)
    nrows = ceil(n / pairs_per_row)
    ncols = 2 * pairs_per_row
    fig, axes = plt.subplots(nrows, ncols,
                             figsize=(3.0 * ncols, 2.3 * nrows),
                             squeeze=False)

    for i, (label, num, den) in enumerate(filters):
        r = i // pairs_per_row
        cp = i % pairs_per_row
        ax_a = axes[r][2 * cp]
        ax_p = axes[r][2 * cp + 1]

        h = response(num, den)
        ax_a.plot(W, np.abs(h), color="C0", lw=1.6)
        ax_p.plot(W, np.angle(h), color="C1", lw=1.6)

        ax_a.set_title(label, fontsize=11)
        ax_p.set_title("phase", fontsize=9, color="0.4")

        ax_a.set_ylim(bottom=0)
        ax_p.set_ylim(-np.pi - 0.2, np.pi + 0.2)
        ax_p.set_yticks(PHASE_TICKS)
        ax_p.set_yticklabels(PHASE_LABELS, fontsize=8)

        for ax in (ax_a, ax_p):
            ax.set_xlim(0, np.pi)
            ax.set_xticks(XTICKS)
            ax.set_xticklabels(XLABELS, fontsize=8)
            ax.axhline(0, color="black", lw=0.5, alpha=0.5)
            ax.grid(True, alpha=0.3)

        if cp == 0:
            ax_a.set_ylabel("amplitude", fontsize=9)

    # hide any unused panels
    for j in range(n, nrows * pairs_per_row):
        r = j // pairs_per_row
        cp = j % pairs_per_row
        axes[r][2 * cp].axis("off")
        axes[r][2 * cp + 1].axis("off")

    fig.suptitle(suptitle, fontsize=13, fontweight="bold")
    fig.supxlabel(r"frequency  $\omega$", fontsize=10)
    fig.tight_layout(rect=[0, 0.0, 1, 0.99])

    out = FIG_DIR / fname
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"saved {out}")


make_figure(ma_filters, "fig4a_ma_filters.png",
            "Figure 4a — Moving-average / differencing filters")
make_figure(ar_filters, "fig4b_ar_filters.png",
            "Figure 4b — Autoregressive (inverse) filters")
make_figure(seasonal_filters, "fig4c_seasonal_filters.png",
            "Figure 4c — Seasonal filters")
