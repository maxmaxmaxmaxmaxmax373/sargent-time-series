"""
Figure 6 — Pairwise coherences among seven U.S. macro series, 1947-present.

A modernized reproduction of Figure 6 from Sargent, Macroeconomic Theory
(1987), Chapter XI, Section 11 (the original used 1948-1976). For each of the
21 unordered pairs of the seven variables, we plot the (magnitude-squared)
coherence

    coh(omega) = |g_xy(omega)|^2 / ( g_x(omega) g_y(omega) ),

estimated from smoothed auto- and cross-periodograms of the linearly detrended
series (smoothing is essential: the raw squared coherence of a single
realization is identically 1). The panels are arranged as the upper triangle of
the 7x7 coherence matrix.

x-axis: frequency omega in [0, pi], labeled as period 2*pi/omega = 8, 4, 2
quarters at omega = pi/4, pi/2, pi.
y-axis: each panel auto-scaled to its own [min, max], which are annotated.

Data: FRED (cached under ../data/fred_cache by fig5_business_cycle_series.py).
Output: ../figures/fig6_coherences.png
"""

import subprocess
import time
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

BASE = Path(__file__).resolve().parent.parent
FIG_DIR = BASE / "figures"
FIG_DIR.mkdir(exist_ok=True)
CACHE = BASE / "data" / "fred_cache"
CACHE.mkdir(parents=True, exist_ok=True)

START = pd.Timestamp("1947-01-01")


def fred(series_id):
    f = CACHE / f"{series_id}.csv"
    if not f.exists():
        url = f"https://fred.stlouisfed.org/graph/fredgraph.csv?id={series_id}"
        subprocess.run(["curl", "-s", "-m", "90", url, "-o", str(f)], check=True)
        time.sleep(1.0)
    df = pd.read_csv(f)
    df.columns = ["date", "value"]
    df["date"] = pd.to_datetime(df["date"])
    df["value"] = pd.to_numeric(df["value"], errors="coerce")
    return df.dropna().set_index("date")["value"]


def to_q(s):
    return s.resample("QS").mean()


def clip(s):
    return s[s.index >= START]


# --- Load data and build the seven transformed quarterly series ------------
gnp     = fred("GNPC96")
unrate  = fred("UNRATE")
baa     = fred("BAA")
m2      = fred("M2SL")
nber_m  = fred("M1444CUSM027SNBR")
cpi_sa  = fred("CPIAUCSL")
cpi_nsa = fred("CPIAUCNS")
ophnfb  = fred("OPHNFB")
aheman  = fred("AHEMAN")

anchor = pd.Timestamp("1959-01-01")
money = pd.concat([nber_m[nber_m.index < anchor] * (m2.loc[anchor] / nber_m.loc[anchor]),
                   m2]).sort_index()
money = money[~money.index.duplicated()]
real_money = (money / cpi_sa).dropna()
real_wage = (aheman / cpi_nsa).dropna()

series = {
    "Real GNP":                            np.log(clip(to_q(gnp))),
    "Unemployment Rate":                   clip(to_q(unrate)),
    "Interest Rate (Baa Bond Rate, NSA)":  clip(to_q(baa)),
    "Change in Real Money Supply":         clip(np.log(to_q(real_money)).diff() * 100),
    "Inflation Rate":                      clip(np.log(to_q(cpi_nsa)).diff() * 400),
    "Output per Worker Hour":              np.log(clip(to_q(ophnfb))),
    "Real Wage (NSA)":                     np.log(clip(to_q(real_wage))),
}

order = list(series.keys())

SHORT = {
    "Real GNP": "Real\nGNP",
    "Unemployment Rate": "Unemployment\nRate",
    "Interest Rate (Baa Bond Rate, NSA)": "Interest Rate\n(Baa, NSA)",
    "Change in Real Money Supply": "Change in\nReal Money",
    "Inflation Rate": "Inflation\nRate",
    "Output per Worker Hour": "Output per\nWorker Hour",
    "Real Wage (NSA)": "Real Wage\n(NSA)",
}


# --- Coherence estimation ---------------------------------------------------
def smooth(y, wlen):
    if wlen < 3:
        return y
    w = np.hanning(wlen)
    w /= w.sum()
    ypad = np.pad(y, wlen // 2, mode="reflect")
    out = np.convolve(ypad, w, mode="same")
    return out[wlen // 2: wlen // 2 + len(y)]


def detrend(x):
    t = np.arange(x.size)
    A = np.vstack([np.ones_like(t, dtype=float), t]).T
    coef, *_ = np.linalg.lstsq(A, x, rcond=None)
    return x - A @ coef


def coherence(sa, sb):
    df = pd.concat([sa, sb], axis=1, sort=True).dropna()
    x = detrend(df.iloc[:, 0].to_numpy(float))
    y = detrend(df.iloc[:, 1].to_numpy(float))
    n = x.size
    wlen = max(7, (n // 15) | 1)
    X, Y = np.fft.rfft(x), np.fft.rfft(y)
    Sxx = smooth((X * np.conj(X)).real, wlen)
    Syy = smooth((Y * np.conj(Y)).real, wlen)
    Sxy = smooth(X * np.conj(Y), wlen)            # complex cross-spectrum
    coh = np.abs(Sxy) ** 2 / (Sxx * Syy + 1e-18)
    freqs = np.linspace(0, np.pi, coh.size)
    return freqs, np.clip(coh, 0, 1)


# --- Plot: upper-triangular 6x6 grid ---------------------------------------
plt.rcParams.update({"font.size": 9})
nv = len(order)                       # 7
fig, axes = plt.subplots(nv - 1, nv - 1, figsize=(14, 12))
fig.subplots_adjust(left=0.13, right=0.99, top=0.88, bottom=0.07,
                    hspace=0.45, wspace=0.30)

XT = [np.pi / 4, np.pi / 2, np.pi]
XL = ["8", "4", "2"]

for i in range(nv - 1):               # row variable = order[i]
    for j in range(nv - 1):           # column variable = order[j + 1]
        ax = axes[i][j]
        if i <= j:
            a, b = order[i], order[j + 1]
            freqs, coh = coherence(series[a], series[b])
            ax.plot(freqs, coh, color="C0", lw=1.0)
            lo, hi = coh.min(), coh.max()
            pad = max(1e-3, 0.05 * (hi - lo))
            ax.set_ylim(lo - pad, hi + pad)
            ax.set_xlim(0, np.pi)
            ax.set_xticks(XT)
            ax.set_yticks([])
            ax.text(0.02, 0.97, f"{hi:.3f}", transform=ax.transAxes,
                    ha="left", va="top", fontsize=6.5, color="0.35")
            ax.text(0.02, 0.03, f"{lo:.3f}", transform=ax.transAxes,
                    ha="left", va="bottom", fontsize=6.5, color="0.35")
            if i == j:                 # bottom cell of this column -> x labels
                ax.set_xticklabels(XL, fontsize=7)
            else:
                ax.set_xticklabels([])
        else:
            ax.axis("off")

# Column headers (row 0 is fully populated)
for j in range(nv - 1):
    pos = axes[0][j].get_position()
    fig.text(pos.x0 + pos.width / 2, 0.90, SHORT[order[j + 1]],
             ha="center", va="bottom", fontsize=9, fontweight="bold")

# Row labels (left margin, aligned to each row)
for i in range(nv - 1):
    pos = axes[i][0].get_position()
    fig.text(0.015, pos.y0 + pos.height / 2, SHORT[order[i]],
             ha="left", va="center", fontsize=9, fontweight="bold")

fig.suptitle("Figure 6 — Pairwise coherences among seven U.S. quarterly series, "
             "1947–present", fontsize=13, fontweight="bold", y=0.96)
fig.text(0.5, 0.015,
         "Each panel: squared coherence vs frequency; x-axis is period in "
         "quarters (8, 4, 2 at $\\omega=\\pi/4,\\pi/2,\\pi$); "
         "panel min/max are annotated.",
         ha="center", fontsize=9, color="0.3")

out = FIG_DIR / "fig6_coherences.png"
fig.savefig(out, dpi=140, bbox_inches="tight")
print(f"saved {out}")

# quick numeric summary of low-frequency (business-cycle) coherence
print("\nLow-frequency (period > 8 quarters) mean squared coherence:")
for i in range(nv):
    for j in range(i + 1, nv):
        f, c = coherence(series[order[i]], series[order[j]])
        lowband = c[f < np.pi / 4]
        print(f"  {order[i][:18]:18s} x {order[j][:18]:18s}: {lowband.mean():.3f}")
