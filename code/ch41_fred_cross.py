"""
Chapter 41, Figure 5 — A moving cross-spectrum of the bill and long rates at the
seasonal band, on FRED data (Sargent 1968, Tables 1-2).

Sargent computed cross-spectra of the first-differenced bill rate against the
longer government rates and reported, at the twelve-month seasonal: coherences
that are "almost all very high", and phase shifts whose "predominance of negative
signs suggests that the longer rates generally lead the bill rate."  At period 12
his Table 1 gives a phase of -0.7978 rad for the 3-month-vs-10-year pair, i.e. a
lead of the bill over the ten-year of -0.7978 * 12/(2*pi) = -1.52 months (the
ten-year leads the bill by about a month and a half).

We reproduce this with the moving cross-spectrum of Chapter 41: demodulate the
first-differenced bill rate (TB3MS) and ten-year rate (GS10) at omega0 = 2*pi/12,
both NSA, smooth the local cross-products through time, and read off the moving
coherence and the moving phase lead (in months).

Output: ../figures/ch41_fred_cross.png
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

OMEGA0 = 2 * np.pi / 12.0
M = 18
ITERS = 2
SMOOTH = 18
SARGENT_LEAD = -0.7978 * 12 / (2 * np.pi)   # Table 1, 3mo-vs-10yr, period 12


def fred(series_id):
    f = CACHE / f"{series_id}.csv"
    if not f.exists():
        url = f"https://fred.stlouisfed.org/graph/fredgraph.csv?id={series_id}"
        subprocess.run(["curl", "-s", "-m", "90", url, "-o", str(f)], check=True)
        time.sleep(1.0)
    df = pd.read_csv(f)
    df.columns = ["date", "value"]
    df["date"] = pd.to_datetime(df["date"], errors="coerce")
    df["value"] = pd.to_numeric(df["value"], errors="coerce")
    return df.dropna().set_index("date")["value"]


def lowpass(z, m=M, iters=ITERS):
    kernel = np.ones(2 * m + 1) / (2 * m + 1)
    out = np.asarray(z, dtype=complex)
    for _ in range(iters):
        out = np.convolve(out, kernel, mode="same")
    return out


def tsmooth(z, m=SMOOTH):
    kernel = np.ones(2 * m + 1) / (2 * m + 1)
    return np.convolve(np.asarray(z, dtype=complex), kernel, mode="same")


def demod(x, omega0=OMEGA0, m=M):
    t = np.arange(len(x))
    return lowpass(x * np.exp(-1j * omega0 * t), m)


def main():
    bill = fred("TB3MS")
    g10 = fred("GS10")
    df = pd.concat({"bill": bill, "g10": g10}, axis=1).dropna().loc["1953":"1975"]
    dy = (df["bill"].diff().dropna() * 100.0)     # bill, basis points
    dx = (df["g10"].diff().dropna() * 100.0)      # ten-year, basis points
    idx = dy.index

    dyb, dxb = demod(dy.to_numpy()), demod(dx.to_numpy())
    Syx = tsmooth(dyb * np.conj(dxb))
    Syy = tsmooth(np.abs(dyb) ** 2).real
    Sxx = tsmooth(np.abs(dxb) ** 2).real
    coh = np.abs(Syx) ** 2 / (Syy * Sxx)
    lead = np.angle(Syx) / OMEGA0                 # lead of bill over 10yr, months

    edge = ITERS * M + SMOOTH
    sl = slice(edge, len(idx) - edge)

    plt.rcParams.update({"font.size": 10, "axes.spines.top": False,
                         "axes.spines.right": False})
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.4))

    ax1.plot(idx[sl], coh[sl], color="C0", lw=1.5)
    ax1.set_ylim(0, 1.05)
    ax1.set_title("(a) moving coherence, bill vs 10-year\nat the 12-month seasonal")
    ax1.set_xlabel("year")
    ax1.set_ylabel("coherence")
    ax1.axhline(np.median(coh[sl]), color="0.6", lw=0.8, ls="--")
    ax1.text(idx[sl][2], np.median(coh[sl]) - 0.08,
             f"median {np.median(coh[sl]):.2f}", fontsize=8, color="0.4")

    ax2.plot(idx[sl], lead[sl], color="C0", lw=1.5, label="demodulation estimate")
    ax2.axhline(0, color="0.7", lw=0.6)
    ax2.axhline(SARGENT_LEAD, color="C3", lw=1.0, ls="--",
                label=f"Sargent Table 1 ({SARGENT_LEAD:.2f} mo)")
    ax2.set_title("(b) moving phase lead of bill over 10-year\n(negative = 10-year leads)")
    ax2.set_xlabel("year")
    ax2.set_ylabel("lead (months)")
    ax2.legend(fontsize=8, loc="upper right")

    fig.suptitle(
        "Moving cross-spectrum at the seasonal band: the ten-year rate leads the "
        "bill rate (FRED, NSA, first differences)", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    out = FIG_DIR / "ch41_fred_cross.png"
    fig.savefig(out, dpi=130)
    print("wrote", out)
    print(f"coherence median {np.median(coh[sl]):.2f}; "
          f"lead median {np.median(lead[sl]):.2f} mo (Sargent {SARGENT_LEAD:.2f})")


if __name__ == "__main__":
    main()
