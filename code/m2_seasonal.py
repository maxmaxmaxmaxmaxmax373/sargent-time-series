"""
Shared helpers for the seasonal-adjustment figures (8-12) of Chapter XI,
Section 33.

Fits an 18th-order autoregression (by OLS) to log M2, seasonally adjusted
(M2SL) and not seasonally adjusted (M2NS), monthly, 1959-present, and provides
the implied autoregressive spectrum. These are reused by

    fig8_m2_spectra.py        spectral densities (SA vs NSA)
    fig9_m2_ar_coeffs.py      autoregressive coefficients
    fig10_m2_lead_filter.py   geometric-lead forecast filter coefficients
    fig11_12_m2_transfer.py   transfer-function magnitude and phase

Data from FRED, cached under ../data/fred_cache/.
"""

import subprocess
import time
from pathlib import Path

import numpy as np
import pandas as pd

BASE = Path(__file__).resolve().parent.parent
CACHE = BASE / "data" / "fred_cache"
CACHE.mkdir(parents=True, exist_ok=True)

P_AR = 18                      # autoregression order (as in Sargent's example)
LAMBDA = 0.9                   # discount in the geometric-lead forecast


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


def load_log_m2():
    """Aligned log M2: (seasonally adjusted, not seasonally adjusted)."""
    sa = np.log(fred("M2SL")).rename("sa")
    nsa = np.log(fred("M2NS")).rename("nsa")
    df = pd.concat([sa, nsa], axis=1).dropna()
    return df["sa"].to_numpy(float), df["nsa"].to_numpy(float)


def load_dlog_m2():
    """Aligned monthly M2 growth, Delta log M2: (SA, NSA).

    Log M2 is very close to a random walk (the level AR(18) has sum of
    coefficients ~ 0.9998, i.e. a near-unit root). Fitting the AR to the
    *growth* rate removes that dominant low-frequency pole, so the seasonal
    structure -- the point of these figures -- is visible in the spectra and,
    especially, in the transfer functions (Figures 11-12).
    """
    sa, nsa = load_log_m2()
    return np.diff(sa) * 100.0, np.diff(nsa) * 100.0      # percent per month


def fit_ar(x, p=P_AR):
    """OLS AR(p) with intercept. Returns (a_1..a_p, sigma^2)."""
    x = np.asarray(x, float)
    n = len(x)
    Y = x[p:]
    cols = [np.ones(n - p)] + [x[p - k:n - k] for k in range(1, p + 1)]
    X = np.column_stack(cols)
    coef, *_ = np.linalg.lstsq(X, Y, rcond=None)
    resid = Y - X @ coef
    sigma2 = resid @ resid / (len(Y) - (p + 1))
    return coef[1:], sigma2


def ar_polynomial(a, omega):
    """A(e^{-i omega}) = 1 - sum_k a_k e^{-i omega k}."""
    k = np.arange(1, len(a) + 1)
    return 1 - (a[:, None] * np.exp(-1j * np.outer(k, omega))).sum(axis=0)


def ar_spectrum(a, sigma2, omega):
    """AR spectral density g(omega) = sigma^2 / |A(e^{-i omega})|^2."""
    return sigma2 / np.abs(ar_polynomial(a, omega)) ** 2


def lead_filter(a, lam=LAMBDA):
    """Geometric-lead forecast filter b(L) (Hansen-Sargent formula (90)).

    For the AR representation a(L) x_t = eps_t with
    a(L) = 1 - sum_{k=1}^r a_k L^k, the projection
    y_t = P_t sum_j lam^j x_{t+j} equals b(L) x_t, where
        b_0 = a(lam)^{-1},
        b_j = a(lam)^{-1} sum_{k=j+1}^r lam^{k-j} a_k,   j = 1..r-1.
    Returns b_0..b_{r-1}.
    """
    r = len(a)
    a_lam = 1 - sum(a[k - 1] * lam ** k for k in range(1, r + 1))
    b = np.zeros(r)
    b[0] = 1.0 / a_lam
    for j in range(1, r):
        b[j] = sum(lam ** (k - j) * a[k - 1] for k in range(j + 1, r + 1)) / a_lam
    return b


def transfer(a, b, omega):
    """Transfer function from innovations to y_t: h(e^{-iw}) = B(e^{-iw})/A(e^{-iw})."""
    j = np.arange(len(b))
    B = (b[:, None] * np.exp(-1j * np.outer(j, omega))).sum(axis=0)
    A = ar_polynomial(a, omega)
    return B / A


# Seasonal periods (months) for monthly data: harmonics of the annual cycle
SEASONAL_PERIODS = [12, 6, 4, 3, 2.4, 2]
