"""
Figure 5 — Seven U.S. macroeconomic time series and their log spectra.

A modernized reproduction of Figure 5 from Sargent, Macroeconomic Theory
(1987), Chapter XI, Section 11. The original covered 1947-1975; here we plot
the full span 1947-present, so the figure both replicates the original window
and extends it to the present.

For each variable: left panel = the realized series; right panel = the log of
an estimated spectral density (a smoothed periodogram of the linearly detrended
series, in the spirit of the QuantEcon "estspec" lecture).

Data: FRED (Federal Reserve Economic Data). The money supply before 1959 is
unavailable on FRED as a modern aggregate, so we splice the NBER macrohistory
series M1444CUSM027SNBR (Friedman-Schwartz / Rasche lineage, currency + demand
deposits + commercial-bank time deposits, ~M2) onto M2SL by growth-chaining at
the 1959-01 junction.

Output: ../figures/fig5_business_cycle_series.png
Cached input data: ../data/fred_cache/*.csv  (committed for reproducibility)
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


# --- FRED download (cached) -------------------------------------------------
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


START = pd.Timestamp("1947-01-01")


def to_quarterly(s):
    """Aggregate a monthly (or quarterly) series to quarter-start averages."""
    return s.resample("QS").mean()


def clip(s):
    """Restrict a series to the common sample 1947-present."""
    return s[s.index >= START]


# --- Download the raw series ------------------------------------------------
gnp      = fred("GNPC96")            # Real GNP, quarterly
unrate   = fred("UNRATE")            # Unemployment rate, monthly SA
baa      = fred("BAA")               # Moody's Baa corporate yield, monthly NSA
m2       = fred("M2SL")              # M2, monthly SA, 1959+
nber_m   = fred("M1444CUSM027SNBR")  # NBER ~M2, monthly, 1947-1969
cpi_sa   = fred("CPIAUCSL")          # CPI, monthly SA  (deflator for real money)
cpi_nsa  = fred("CPIAUCNS")          # CPI, monthly NSA (for inflation, real wage)
ophnfb   = fred("OPHNFB")            # Output per hour, nonfarm business, quarterly
aheman   = fred("AHEMAN")            # Avg hourly earnings, mfg, monthly NSA

# --- Splice a continuous nominal money series (monthly) --------------------
anchor = pd.Timestamp("1959-01-01")
scale = m2.loc[anchor] / nber_m.loc[anchor]      # growth-chain == ratio at anchor
money_pre = nber_m[nber_m.index < anchor] * scale
money = pd.concat([money_pre, m2]).sort_index()
money = money[~money.index.duplicated()]

# --- Build the seven variables ---------------------------------------------
real_money = (money / cpi_sa).dropna()
dlog_money = np.log(to_quarterly(real_money)).diff() * 100      # % per quarter

inflation = np.log(to_quarterly(cpi_nsa)).diff() * 400          # annualized %

real_wage = (aheman / cpi_nsa).dropna()                         # NSA real wage

variables = [
    # (name, display series, spectrum input, take log for spectrum?)
    ("Real GNP",                      clip(to_quarterly(gnp)),       clip(to_quarterly(gnp)),       True),
    ("Unemployment Rate",             clip(to_quarterly(unrate)),    clip(to_quarterly(unrate)),    False),
    ("Interest Rate (Baa Bond Rate, NSA)", clip(to_quarterly(baa)),  clip(to_quarterly(baa)),       False),
    ("Change in Real Money Supply",   clip(dlog_money),              clip(dlog_money),              False),
    ("Inflation Rate",                clip(inflation),               clip(inflation),               False),
    ("Output per Worker Hour",        clip(to_quarterly(ophnfb)),    clip(to_quarterly(ophnfb)),    True),
    ("Real Wage (NSA)",               clip(to_quarterly(real_wage)), clip(to_quarterly(real_wage)), True),
]


# --- Spectral estimation: smoothed periodogram of the detrended series ------
def smooth(y, wlen):
    if wlen < 3:
        return y
    w = np.hanning(wlen)
    w /= w.sum()
    ypad = np.pad(y, wlen // 2, mode="reflect")
    out = np.convolve(ypad, w, mode="same")
    return out[wlen // 2: wlen // 2 + len(y)]


def log_spectrum(series, use_log):
    x = series.dropna().to_numpy(dtype=float)
    if use_log:
        x = np.log(x)
    t = np.arange(x.size)
    A = np.vstack([np.ones_like(t, dtype=float), t]).T      # linear detrend
    coef, *_ = np.linalg.lstsq(A, x, rcond=None)
    x = x - A @ coef
    n = x.size
    I = np.abs(np.fft.rfft(x)) ** 2 / n                      # periodogram
    wlen = max(5, (n // 20) | 1)                             # odd smoothing window
    I = smooth(I, wlen)
    freqs = np.linspace(0, np.pi, I.size)
    return freqs, np.log(I + 1e-12)


# --- Plot -------------------------------------------------------------------
plt.rcParams.update({"font.size": 10, "axes.spines.top": False,
                     "axes.spines.right": False})

nvar = len(variables)
fig, axes = plt.subplots(nvar, 2, figsize=(11, 2.4 * nvar))

XT = [0, np.pi / 2, np.pi]
XL = ["0", r"$\pi/2$", r"$\pi$"]

for i, (name, disp, spec, use_log) in enumerate(variables):
    axL, axR = axes[i, 0], axes[i, 1]

    d = disp.dropna()
    axL.plot(d.index, d.to_numpy(), color="C0", lw=0.9)
    axL.set_title(name, fontsize=11, fontweight="bold")
    axL.grid(True, alpha=0.3)

    freqs, logI = log_spectrum(spec, use_log)
    axR.plot(freqs, logI, color="C3", lw=1.3)
    # seasonal frequencies for quarterly data: 4-quarter (pi/2) and 2-quarter (pi)
    for fs in (np.pi / 2, np.pi):
        axR.axvline(fs, color="0.6", ls=":", lw=0.8)
    axR.set_xticks(XT)
    axR.set_xticklabels(XL)
    axR.set_xlim(0, np.pi)
    axR.grid(True, alpha=0.3)

    if i == 0:
        axL.annotate("Variable", xy=(0.5, 1.28), xycoords="axes fraction",
                     ha="center", fontsize=12, fontweight="bold")
        axR.annotate("Log of Estimated Spectrum", xy=(0.5, 1.28),
                     xycoords="axes fraction", ha="center", fontsize=12,
                     fontweight="bold")
    if i == nvar - 1:
        axL.set_xlabel("year")
        axR.set_xlabel(r"frequency  $\omega$")

fig.suptitle("Figure 5 — Seven U.S. quarterly series, 1947–present, "
             "with log estimated spectra", fontsize=13, fontweight="bold",
             y=1.005)
fig.tight_layout(rect=[0, 0, 1, 0.99])

out = FIG_DIR / "fig5_business_cycle_series.png"
fig.savefig(out, dpi=140, bbox_inches="tight")
print(f"saved {out}")
