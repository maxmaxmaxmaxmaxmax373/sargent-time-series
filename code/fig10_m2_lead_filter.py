"""
Figure 10 — Geometric-lead forecast-filter coefficients for M2.

The projection y_t = P_t sum_j lambda^j x_{t+j} (lambda = 0.9) equals b(L) x_t,
where b(L) is the Hansen-Sargent filter (formula (90)) implied by each series'
AR(18) representation. This figure plots b^{nsa}(L) and b^{sa}(L) versus lag,
using the AR(18) fits to unadjusted (M2NS) and adjusted (M2SL) log M2.

Note the difference near lag 12: the unadjusted forecast filter inherits the
seasonal autoregressive structure that the adjusted one lacks.

Output: ../figures/fig10_m2_lead_filter.png
"""

from pathlib import Path
import sys

import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
from m2_seasonal import load_dlog_m2, fit_ar, lead_filter, LAMBDA

FIG_DIR = Path(__file__).resolve().parent.parent / "figures"

sa, nsa = load_dlog_m2()
b_sa = lead_filter(fit_ar(sa)[0])
b_nsa = lead_filter(fit_ar(nsa)[0])
lags = np.arange(len(b_sa))                       # 0 .. 17

plt.rcParams.update({"font.size": 11, "axes.spines.top": False,
                     "axes.spines.right": False})
fig, ax = plt.subplots(figsize=(10, 5))

ax.plot(lags, b_nsa, "o-", color="C1", ms=4, lw=1.2, label="$b^{nsa}(L)$ (unadjusted)")
ax.plot(lags, b_sa, "s-", color="C0", ms=4, lw=1.2, label="$b^{sa}(L)$ (adjusted)")
ax.axhline(0, color="black", lw=0.6)
ax.axvline(12, color="0.6", ls=":", lw=1)
ax.annotate("seasonal lag 12", xy=(12, 0), xytext=(12.3, 0.6 * ax.get_ylim()[1]),
            fontsize=9, color="0.4")

ax.set_xticks(lags)
ax.set_xlabel("lag")
ax.set_ylabel("filter coefficient $b_j$")
ax.set_title(r"Figure 10 — Geometric-lead forecast filters $b(L)$ for M2 growth "
             rf"($\lambda={LAMBDA}$), 1959–present",
             fontsize=12, fontweight="bold")
ax.legend(frameon=False)
ax.grid(True, axis="y", alpha=0.25)

fig.tight_layout()
out = FIG_DIR / "fig10_m2_lead_filter.png"
fig.savefig(out, dpi=150, bbox_inches="tight")
print(f"saved {out}")
print(f"b_12: NSA={b_nsa[12]:+.4f}  SA={b_sa[12]:+.4f}")
