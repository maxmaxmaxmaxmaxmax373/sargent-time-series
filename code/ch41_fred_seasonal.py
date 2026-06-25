"""
Chapter 41, Figure 4 — Sargent (1968) replicated on FRED data:
the seasonal in 1950s interest rates, by complex demodulation.

Sargent studied monthly U.S. interest rates over 1951-1960, working with first
differences to remove the trend.  We reproduce his central seasonal findings with
the closest *not-seasonally-adjusted* monthly series on FRED:

    TB3MS  3-month Treasury bill, secondary market (from 1934), NSA  -> "bill"
    GS1    1-year Treasury constant maturity (from 1953-04),    NSA
    GS5    5-year Treasury constant maturity (from 1953-04),    NSA
    GS20   20-year Treasury constant maturity (from 1953-04),   NSA

(Interest-rate series on FRED are published NOT seasonally adjusted, which is
exactly what this exercise requires.  Sargent's two-year government and commercial
paper rates do not reach back to the 1950s on FRED, so we omit them.)

Three panels:
 (a) the smoothed-periodogram spectral density of the first-differenced bill rate
     over 1951-1972 shows sharp peaks at the seasonal periods 12, 6, 4, 3, 2.4, 2
     months -- Sargent's Figure 1a;
 (b) the moving seasonal component s_t = sum_k 2 Re[d_t(k*omega0) e^{i k omega0 t}]
     of the first-differenced bill rate (k = 1..4), basis points, over the whole
     sample -- the analog of Sargent's Figure 3a -- showing the bill-rate seasonal
     that was strong through the 1950s-1970s and faded after the mid-1980s;
 (c) the seasonal "share" (local variance of the seasonal band / local variance of
     the first difference) by maturity over 1957-1972, which declines with maturity
     -- Sargent's finding (b).

Output: ../figures/ch41_fred_seasonal.png
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
M = 18                      # Sargent's half-width
ITERS = 2                   # iterate twice -> Bartlett window
HARM = (1, 2, 3, 4)         # 12-month cycle and its first three harmonics


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


def demodulate(x, omega0, m=M):
    t = np.arange(len(x))
    return lowpass(x * np.exp(-1j * omega0 * t), m)


def moving_seasonal(x, m=M, harm=HARM):
    t = np.arange(len(x))
    s = np.zeros(len(x))
    for k in harm:
        d = demodulate(x, k * OMEGA0, m)
        s += 2.0 * np.real(d * np.exp(1j * k * OMEGA0 * t))
    return s


def msmooth(z, m):
    kernel = np.ones(2 * m + 1) / (2 * m + 1)
    return np.convolve(np.asarray(z, float), kernel, mode="same")


def smooth(y, wlen):
    if wlen < 3:
        return y
    w = np.hanning(wlen)
    w /= w.sum()
    ypad = np.pad(y, wlen // 2, mode="reflect")
    out = np.convolve(ypad, w, mode="same")
    return out[wlen // 2: wlen // 2 + len(y)]


def log_spectrum(x):
    x = np.asarray(x, float)
    x = x - x.mean()
    n = x.size
    I = np.abs(np.fft.rfft(x)) ** 2 / n
    I = smooth(I, max(5, (n // 25) | 1))
    freqs = np.linspace(0, np.pi, I.size)
    return freqs, np.log(I + 1e-12)


def diff_bp(series):
    """First difference in basis points."""
    return series.diff().dropna() * 100.0


def main():
    bill = fred("TB3MS")
    mats = {"bill": fred("TB3MS"), "1yr": fred("GS1"),
            "5yr": fred("GS5"), "20yr": fred("GS20")}

    plt.rcParams.update({"font.size": 10, "axes.spines.top": False,
                         "axes.spines.right": False})
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.4))

    # ---- (a) spectrum of the first-differenced bill rate, 1951-1972 --------
    dbill = diff_bp(bill).loc["1951":"1972"]
    freqs, logspec = log_spectrum(dbill.to_numpy())
    period = 2 * np.pi / np.clip(freqs, 1e-9, None)
    ax = axes[0]
    ax.plot(freqs, logspec, color="C0", lw=1.3)
    for P in [12, 6, 4, 3, 2.4, 2]:
        ax.axvline(2 * np.pi / P, color="C3", lw=0.7, ls=":")
    ax.set_xticks([2 * np.pi / P for P in [12, 6, 4, 3, 2]])
    ax.set_xticklabels(["12", "6", "4", "3", "2"])
    ax.set_xlabel("period (months)")
    ax.set_ylabel(r"$\log$ spectral density")
    ax.set_title("(a) bill rate $\\Delta$, 1951–1972:\nseasonal peaks (Sargent Fig. 1a)")

    # ---- (b) moving seasonal of the bill rate, Sargent era -----------------
    # The seasonal is plotted in basis points, so it must be read against the
    # low rate volatility of the 1950s-60s (just as Sargent's Fig. 3a was a
    # single short window).  Over the whole post-war sample the absolute seasonal
    # swings are dwarfed by the 1979-82 volatility, which is why we focus here on
    # the era Sargent studied; the *relative* strength of the seasonal across
    # eras is taken up in the text.
    db_full = diff_bp(bill)
    s = moving_seasonal(db_full.to_numpy())
    seas = pd.Series(s, index=db_full.index)
    R1 = pd.Series(2 * np.abs(demodulate(db_full.to_numpy(), OMEGA0)),
                   index=db_full.index)
    win_b = slice("1951", "1975")
    ax = axes[1]
    ax.plot(seas.loc[win_b].index, seas.loc[win_b].values, color="C0", lw=0.9,
            label="seasonal $s_t$")
    ax.plot(R1.loc[win_b].index, R1.loc[win_b].values, color="C3", lw=1.0,
            label=r"amplitude $\pm R_t$ (12-mo)")
    ax.plot(R1.loc[win_b].index, -R1.loc[win_b].values, color="C3", lw=1.0)
    ax.axhline(0, color="0.7", lw=0.6)
    ax.set_title("(b) moving seasonal of the bill rate, 1951–1975\n(basis points; Sargent Fig. 3a)")
    ax.set_xlabel("year")
    ax.set_ylabel("seasonal (bp)")
    ax.legend(fontsize=8, loc="upper left")

    # ---- (c) seasonal share by maturity, 1957-1967 -------------------------
    win = slice("1957", "1967")
    labels, shares = [], []
    for lab in ["bill", "1yr", "5yr", "20yr"]:
        dr = diff_bp(mats[lab])
        x = dr.to_numpy()
        seas = moving_seasonal(x, m=M)
        var_s = msmooth(seas ** 2, 30)
        var_t = msmooth((x - x.mean()) ** 2, 30)
        share = pd.Series(var_s / var_t, index=dr.index).loc[win]
        labels.append(lab)
        shares.append(share.median())
    ax = axes[2]
    ax.bar(labels, shares, color=["C0", "C1", "C2", "C3"])
    ax.set_ylim(0, max(shares) * 1.25)
    ax.set_title("(c) seasonal share by maturity, 1957–1967\n(declines with maturity; Sargent finding b)")
    ax.set_ylabel("seasonal var. / total var.")
    for i, v in enumerate(shares):
        ax.text(i, v + 0.01, f"{v:.2f}", ha="center", fontsize=9)

    fig.suptitle(
        "Sargent (1968) replicated on FRED: the seasonal in U.S. interest rates "
        "(not seasonally adjusted), by complex demodulation", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    out = FIG_DIR / "ch41_fred_seasonal.png"
    fig.savefig(out, dpi=130)
    print("wrote", out)
    print("seasonal shares:", dict(zip(labels, [round(v, 3) for v in shares])))


if __name__ == "__main__":
    main()
