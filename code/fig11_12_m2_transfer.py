"""
Figures 11 and 12 — Transfer function from M2 innovations to the geometric-lead
forecast y_t.

For each series, y_t = b(L) x_t and x_t = (1/a(L)) eps_t, so the transfer
function from the innovation eps_t to y_t is

    h(e^{-iw}) = B(e^{-iw}) / A(e^{-iw}),

with b(L) the geometric-lead filter (lambda = 0.9) and a(L) the AR(18)
polynomial. Figure 11 is the unadjusted case (b^{nsa}/a^{nsa}); Figure 12 is the
adjusted case (b^{sa}/a^{sa}). Each shows amplitude (log scale) and phase versus
frequency; the other case is overlaid faintly for comparison.

Outputs: ../figures/fig11_m2_transfer_nsa.png
         ../figures/fig12_m2_transfer_sa.png
"""

from pathlib import Path
import sys

import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
from m2_seasonal import (load_dlog_m2, fit_ar, lead_filter, transfer, LAMBDA,
                         SEASONAL_PERIODS)

FIG_DIR = Path(__file__).resolve().parent.parent / "figures"

sa, nsa = load_dlog_m2()
a_sa, _ = fit_ar(sa)
a_nsa, _ = fit_ar(nsa)
b_sa, b_nsa = lead_filter(a_sa), lead_filter(a_nsa)

omega = np.linspace(1e-3, np.pi, 2000)
h_nsa = transfer(a_nsa, b_nsa, omega)
h_sa = transfer(a_sa, b_sa, omega)

plt.rcParams.update({"font.size": 11, "axes.spines.top": False,
                     "axes.spines.right": False})


def make(focus, fname, fig_no):
    """focus = 'nsa' or 'sa'."""
    if focus == "nsa":
        hf, ho = h_nsa, h_sa
        flab, olab, title = "unadjusted (NSA)", "adjusted (SA)", "unadjusted M2 growth"
    else:
        hf, ho = h_sa, h_nsa
        flab, olab, title = "adjusted (SA)", "unadjusted (NSA)", "adjusted M2 growth"

    fig, (axA, axP) = plt.subplots(2, 1, figsize=(9, 7), sharex=True)

    for P in SEASONAL_PERIODS:
        for ax in (axA, axP):
            ax.axvline(2 * np.pi / P, color="0.8", ls=":", lw=0.8, zorder=0)

    # amplitude (log scale)
    axA.semilogy(omega, np.abs(ho), color="0.6", lw=1.0, ls="--", label=olab)
    axA.semilogy(omega, np.abs(hf), color="C3", lw=1.6, label=flab)
    axA.set_ylabel("amplitude of filter")
    axA.set_title(f"Figure {fig_no} — Transfer function from innovations to "
                  f"$y_t$, {title}\n"
                  rf"($h = b/a$, $\lambda={LAMBDA}$), 1959–present",
                  fontsize=11.5, fontweight="bold")
    axA.legend(frameon=False, fontsize=9, loc="upper right")
    axA.grid(True, which="both", alpha=0.2)

    # phase
    axP.plot(omega, np.angle(ho), color="0.6", lw=1.0, ls="--")
    axP.plot(omega, np.angle(hf), color="C0", lw=1.6)
    axP.set_ylabel("phase of filter")
    axP.set_xlabel(r"frequency  $\omega$   (period $=2\pi/\omega$ months)")
    axP.set_yticks([-np.pi, 0, np.pi])
    axP.set_yticklabels([r"$-\pi$", "0", r"$\pi$"])
    axP.grid(True, alpha=0.2)

    # period labels along the top axis
    ticks_P = [12, 6, 4, 3, 2.4, 2]
    axA.set_xticks([2 * np.pi / P for P in ticks_P])
    axA.set_xticklabels([f"{P:g}" for P in ticks_P])
    axP.set_xlim(0, np.pi)

    fig.tight_layout()
    out = FIG_DIR / fname
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"saved {out}")


make("nsa", "fig11_m2_transfer_nsa.png", 11)
make("sa", "fig12_m2_transfer_sa.png", 12)
