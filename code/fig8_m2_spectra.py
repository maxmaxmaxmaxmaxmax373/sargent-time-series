"""
Figure 8 — Spectral densities of seasonally adjusted and unadjusted M2.

Modernized reproduction of Figure 8 from Sargent, Macroeconomic Theory (1987),
Chapter XI, Section 33 (original: monthly M2, 1959:1-1986:2; here 1959-present).
An 18th-order autoregression is fit by OLS to log M2, seasonally adjusted
(M2SL) and not seasonally adjusted (M2NS); the figure plots the log of the
implied AR spectral density g(omega) = sigma^2 / |A(e^{-i omega})|^2.

The unadjusted series shows peaks at the seasonal frequencies (periods 12, 6,
4, 3, 2.4, 2 months); seasonal adjustment carves dips there.

(Sargent's figure caption reads "M1", but the text describes the construction
with M2; we follow the text and use M2.)

Output: ../figures/fig8_m2_spectra.png
"""

from pathlib import Path
import sys

import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
from m2_seasonal import (load_log_m2, fit_ar, ar_spectrum, P_AR,
                         SEASONAL_PERIODS)

FIG_DIR = Path(__file__).resolve().parent.parent / "figures"
FIG_DIR.mkdir(exist_ok=True)

sa, nsa = load_log_m2()
a_sa, s2_sa = fit_ar(sa)
a_nsa, s2_nsa = fit_ar(nsa)

omega = np.linspace(2 * np.pi / 20, np.pi, 1500)     # period 20 down to 2 months
log_g_sa = np.log(ar_spectrum(a_sa, s2_sa, omega))
log_g_nsa = np.log(ar_spectrum(a_nsa, s2_nsa, omega))

plt.rcParams.update({"font.size": 11, "axes.spines.top": False,
                     "axes.spines.right": False})
fig, ax = plt.subplots(figsize=(10, 6))

# seasonal frequencies (periods 12, 6, 4, 3, 2.4, 2 months)
for P in SEASONAL_PERIODS:
    ax.axvline(2 * np.pi / P, color="0.7", ls=":", lw=0.8, zorder=0)

ax.plot(omega, log_g_nsa, color="C1", lw=1.8, label="Unadjusted (NSA)")
ax.plot(omega, log_g_sa, color="C0", lw=1.8, label="Adjusted (SA)")
ax.fill_between(omega, log_g_sa, log_g_nsa, where=(log_g_nsa > log_g_sa),
                color="C1", alpha=0.10)

# x-axis labeled as period in months
ticks_P = [18, 12, 6, 4, 3, 2.4, 2]
ax.set_xticks([2 * np.pi / P for P in ticks_P])
ax.set_xticklabels([f"{P:g}" for P in ticks_P])
ax.set_xlim(2 * np.pi / 20, np.pi)
ax.set_xlabel("Period (months)")
ax.set_ylabel("Log spectral density")

# annotate the seasonal harmonics across the top
ymax = max(log_g_nsa.max(), log_g_sa.max())
ax.text(2 * np.pi / 12, ymax, "  seasonal harmonics: 12, 6, 4, 3, 2.4, 2 months",
        fontsize=8, color="0.4", va="bottom", ha="left")

ax.set_title(f"Figure 8 — AR({P_AR}) spectral densities of seasonally adjusted "
             "and unadjusted M2, 1959–present", fontsize=12, fontweight="bold")
ax.legend(loc="upper right", frameon=False)
ax.grid(True, alpha=0.25)

fig.tight_layout()
out = FIG_DIR / "fig8_m2_spectra.png"
fig.savefig(out, dpi=150, bbox_inches="tight")
print(f"saved {out}")
print(f"AR({P_AR}) sigma^2: SA={s2_sa:.3e}  NSA={s2_nsa:.3e}")
