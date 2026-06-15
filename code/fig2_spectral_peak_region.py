"""
Figure 2 — Parameter regions for the second-order process
    Y_t = t1 * Y_{t-1} + t2 * Y_{t-2} + eps_t
in the (t1, t2) plane.

Enhanced (color-coded) version of the Jenkins & Watts (1969, p. 229) diagram
reproduced as Figure 2 in Sargent, Macroeconomic Theory (1987), Chapter XI.
This is an extension of the QuantEcon "Samuelson model" stability diagram
(https://python.quantecon.org/samuelson.html): the stationarity triangle and
the real/complex-root parabola, with the spectral peak/trough curves added.

Boundary curves
---------------
Stationarity triangle (vertices (-2,-1), (2,-1), (0,1)):
    t1 + t2 = 1            (right edge)
    t2 = 1 + t1           (left edge)
    t2 = -1               (bottom edge)
Real/complex root boundary (parabola):
    t1^2 + 4 t2 = 0        i.e.  t2 = -t1^2 / 4      (complex roots below)
Spectral peak/trough boundaries (Chapter XI, eqs 41-42):
    t1 (1 - t2) = -4 t2    i.e.  t2 = t1 / (t1 - 4)
    t1 (1 - t2) =  4 t2    i.e.  t2 = t1 / (t1 + 4)

Regions inside the triangle
---------------------------
  peak    : region 40 holds and t2 < 0     -> peak in spectrum
  trough  : region 40 holds and t2 > 0     -> trough in spectrum
  wing    : complex roots but not region 40 -> oscillatory covariogram,
                                               but NO spectral peak
  nopeak  : real roots, not region 40       -> no peak or trough

where region 40 is  |-t1 (1 - t2) / (4 t2)| < 1.

Output: ../figures/fig2_spectral_peak_region.png
"""

from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.patches import Patch

FIG_DIR = Path(__file__).resolve().parent.parent / "figures"
FIG_DIR.mkdir(exist_ok=True)

# --- Grid and region classification ----------------------------------------
t1 = np.linspace(-2.6, 2.6, 1300)
t2 = np.linspace(-1.4, 1.3, 1300)
T1, T2 = np.meshgrid(t1, t2)

in_tri = (T1 + T2 < 1) & (T2 - T1 < 1) & (T2 > -1)

with np.errstate(divide="ignore", invalid="ignore"):
    val = -T1 * (1 - T2) / (4 * T2)
in40 = np.abs(val) < 1                       # region (40): peak or trough
complex_roots = (T1**2 + 4 * T2 < 0)         # below parabola -> complex roots

# integer codes: 0 outside/explosive, 1 nopeak, 2 wing, 3 peak, 4 trough
code = np.zeros_like(T1, dtype=int)
code[in_tri & ~in40 & ~complex_roots] = 1    # real roots, no peak/trough
code[in_tri & ~in40 & complex_roots] = 2     # complex roots, no spectral peak
code[in_tri & in40 & (T2 < 0)] = 3           # peak
code[in_tri & in40 & (T2 > 0)] = 4           # trough

colors = ["white", "#ececec", "#a8ddb5", "#9ecae1", "#fdbe85"]
cmap = ListedColormap(colors)
norm = BoundaryNorm([-0.5, 0.5, 1.5, 2.5, 3.5, 4.5], cmap.N)

# --- Plot -------------------------------------------------------------------
plt.rcParams.update({"font.size": 11})
fig, ax = plt.subplots(figsize=(9.5, 8))

ax.pcolormesh(T1, T2, code, cmap=cmap, norm=norm, shading="auto")

# Boundary curves -----------------------------------------------------------
xt = np.linspace(-2, 2, 600)

# Stationarity triangle
ax.plot([0, 2], [1, -1], color="black", lw=1.6)      # right: t1 + t2 = 1
ax.plot([0, -2], [1, -1], color="black", lw=1.6)     # left:  t2 = 1 + t1
ax.plot([-2, 2], [-1, -1], color="black", lw=1.6)    # bottom: t2 = -1

# Parabola t1^2 + 4 t2 = 0
ax.plot(xt, -xt**2 / 4, color="#08519c", lw=1.8)

# Spectral peak/trough curves, clipped to the triangle
def clip_to_triangle(x, y):
    m = (x + y < 1) & (y - x < 1) & (y > -1)
    return x, np.where(m, y, np.nan)

c1 = xt / (xt - 4)        # t1(1 - t2) = -4 t2
c2 = xt / (xt + 4)        # t1(1 - t2) =  4 t2
ax.plot(*clip_to_triangle(xt, c1), color="#a63603", lw=1.6)
ax.plot(*clip_to_triangle(xt, c2), color="#a63603", lw=1.6)

# Axes through the origin ----------------------------------------------------
ax.axhline(0, color="black", lw=0.8)
ax.axvline(0, color="black", lw=0.8)
ax.annotate("", xy=(2.6, 0), xytext=(-2.6, 0),
            arrowprops=dict(arrowstyle="->", lw=0.8))
ax.annotate("", xy=(0, 1.3), xytext=(0, -1.4),
            arrowprops=dict(arrowstyle="->", lw=0.8))
ax.text(2.55, -0.12, r"$t_1$", fontsize=13)
ax.text(0.06, 1.24, r"$t_2$", fontsize=13)

tick_bbox = dict(facecolor="white", edgecolor="none", alpha=0.85, pad=0.5)
for x in (-2, -1, 1, 2):
    ax.plot([x, x], [-0.03, 0.03], color="black", lw=0.8)
    ax.text(x, -0.13, f"${x}$", ha="center", fontsize=9, bbox=tick_bbox)
for y in (-1, 1):
    ax.plot([-0.03, 0.03], [y, y], color="black", lw=0.8)
    ax.text(-0.17, y, f"${y}$", va="center", ha="right", fontsize=9, bbox=tick_bbox)

# Region labels --------------------------------------------------------------
ax.text(0, -0.55, "PEAK IN\nSPECTRUM", ha="center", va="center",
        fontsize=10, fontweight="bold")
ax.text(0, 0.42, "TROUGH IN\nSPECTRUM", ha="center", va="center", fontsize=9,
        fontweight="bold")

# Explosive (non-stationary) regions
ax.text(1.35, 0.7, "EXPLOSIVE\nGROWTH", ha="center", va="center",
        fontsize=9, color="0.35")
ax.text(-1.35, 0.7, "EXPLOSIVE\nOSCILLATIONS", ha="center", va="center",
        fontsize=9, color="0.35")
ax.text(0, -1.22, "EXPLOSIVE OSCILLATIONS", ha="center", va="center",
        fontsize=9, color="0.35")

# Boundary equation labels (white background so they read over lines/fills)
lbl_bbox = dict(facecolor="white", edgecolor="none", alpha=0.8, pad=0.6)
ax.text(0.62, 0.52, r"$t_1+t_2=1$", fontsize=9, rotation=-45, color="black",
        ha="center", va="center", bbox=lbl_bbox)
ax.text(-0.62, 0.52, r"$t_2=1+t_1$", fontsize=9, rotation=45, color="black",
        ha="center", va="center", bbox=lbl_bbox)
ax.text(-1.55, -0.93, r"$t_2=-1$", fontsize=9, bbox=lbl_bbox)
ax.text(1.30, -0.66, r"$t_1^2+4t_2=0$", fontsize=9, color="#08519c",
        rotation=-30, ha="center", va="center", bbox=lbl_bbox)
ax.text(0.86, 0.20, r"$t_1(1-t_2)=4t_2$", fontsize=8, color="#a63603",
        rotation=12, ha="center", va="center", bbox=lbl_bbox)
ax.text(-0.86, 0.20, r"$t_1(1-t_2)=-4t_2$", fontsize=8, color="#a63603",
        rotation=-12, ha="center", va="center", bbox=lbl_bbox)

# Legend ---------------------------------------------------------------------
legend_handles = [
    Patch(facecolor="#9ecae1", edgecolor="0.5", label="Peak in spectrum"),
    Patch(facecolor="#fdbe85", edgecolor="0.5", label="Trough in spectrum"),
    Patch(facecolor="#a8ddb5", edgecolor="0.5",
          label="Oscillatory covariogram, no spectral peak"),
    Patch(facecolor="#ececec", edgecolor="0.5", label="No peak or trough"),
]
ax.legend(handles=legend_handles, loc="lower center",
          bbox_to_anchor=(0.5, -0.16), ncol=2, frameon=False, fontsize=9)

ax.set_xlim(-2.6, 2.6)
ax.set_ylim(-1.4, 1.3)
ax.set_aspect("equal")
ax.axis("off")

fig.tight_layout()
out_path = FIG_DIR / "fig2_spectral_peak_region.png"
fig.savefig(out_path, dpi=150, bbox_inches="tight")
print(f"saved {out_path}")
