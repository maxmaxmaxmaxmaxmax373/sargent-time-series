"""
Figure 9 — Autoregressive coefficients of seasonally adjusted and unadjusted M2.

Plots the OLS AR(18) coefficients a_1..a_18 fit to log M2, seasonally adjusted
(M2SL) and not seasonally adjusted (M2NS), 1959-present. Note the difference
near lag 12: the unadjusted series carries a seasonal autoregressive coefficient
that the adjusted series lacks.

Output: ../figures/fig9_m2_ar_coeffs.png
"""

from pathlib import Path
import sys

import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
from m2_seasonal import load_dlog_m2, fit_ar, P_AR

FIG_DIR = Path(__file__).resolve().parent.parent / "figures"

sa, nsa = load_dlog_m2()
a_sa, _ = fit_ar(sa)
a_nsa, _ = fit_ar(nsa)
lags = np.arange(1, P_AR + 1)

plt.rcParams.update({"font.size": 11, "axes.spines.top": False,
                     "axes.spines.right": False})
fig, ax = plt.subplots(figsize=(10, 5))

w = 0.4
ax.bar(lags - w / 2, a_nsa, width=w, color="C1", label="Unadjusted (NSA)")
ax.bar(lags + w / 2, a_sa, width=w, color="C0", label="Adjusted (SA)")
ax.axhline(0, color="black", lw=0.6)
ax.axvline(12, color="0.6", ls=":", lw=1)
ymax = ax.get_ylim()[1]
ax.annotate("seasonal lag 12", xy=(12, ymax), xytext=(12.4, 0.8 * ymax),
            fontsize=9, color="0.4")

ax.set_xticks(lags)
ax.set_xlabel("lag")
ax.set_ylabel("autoregressive coefficient $a_k$")
ax.set_title(f"Figure 9 — AR({P_AR}) coefficients of seasonally adjusted and "
             "unadjusted M2 growth, 1959–present", fontsize=12, fontweight="bold")
ax.legend(frameon=False)
ax.grid(True, axis="y", alpha=0.25)

fig.tight_layout()
out = FIG_DIR / "fig9_m2_ar_coeffs.png"
fig.savefig(out, dpi=150, bbox_inches="tight")
print(f"saved {out}")
print(f"a_12: NSA={a_nsa[11]:+.3f}  SA={a_sa[11]:+.3f}")
