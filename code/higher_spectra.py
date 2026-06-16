"""
Self-contained higher-order spectral estimators for the nonlinear-representation
lab chapter (Chapter §40).

Replaces the external ``higher-spectrum`` package with a small, dependency-light
(NumPy-only) implementation of the *direct* bispectrum and bicoherence
estimators (Nikias and Petropulu 1993, "Higher-Order Spectra Analysis").

Conventions (chosen to match the rest of the book and the lab figures):
  * Frequencies are in cycles/sample; the returned ``waxis`` runs from -0.5 to
    just under +0.5 via ``fftshift(fftfreq(nfft))``.
  * ``bispectrumd`` returns the complex bispectrum on an ``nfft x nfft`` grid;
    ``bicoherence`` returns the real, [0,1]-valued, Kim-Powers normalized
    bicoherence on the same grid.
  * The signatures accept ``y`` shaped either (N,) or (N,1) so the lab code can
    call them exactly as it called the external package.

The estimator: remove the mean, split the series into overlapping records,
Hanning-window and DFT each record, then average the triple product
X(f1) X(f2) conj(X(f1+f2)) (indices taken mod nfft) over records. Segment
averaging supplies the variance reduction; no separate 2-D frequency smoothing
is applied.
"""

import numpy as np


def _records(y, nfft, nsamp, overlap):
    """Return an (n_records, nfft) array of windowed, mean-removed DFTs."""
    x = np.asarray(y, dtype=float).ravel()
    x = x - x.mean()
    step = max(1, int(round(nsamp * (1.0 - overlap / 100.0))))
    starts = range(0, len(x) - nsamp + 1, step)
    win = np.hanning(nsamp)
    win = win / np.sqrt(np.mean(win ** 2))      # unit-power window
    X = []
    for s in starts:
        seg = x[s:s + nsamp]
        X.append(np.fft.fft(seg * win, nfft) / nsamp)
    if not X:
        raise ValueError("series too short for the requested nsamp/overlap")
    return np.asarray(X)


def _index_grids(nfft):
    i = np.arange(nfft)
    I1, I2 = np.meshgrid(i, i, indexing="ij")
    return I1, I2, (I1 + I2) % nfft


def bispectrumd(y, nfft=256, nsamp=256, overlap=50, wind=None):
    """Direct-method bispectrum estimate.

    Parameters
    ----------
    y : array_like, shape (N,) or (N, 1)
    nfft, nsamp, overlap : int
        FFT length, segment length, percent overlap.
    wind : ignored
        Present only for call-signature compatibility; segment averaging
        already provides smoothing.

    Returns
    -------
    Bspec : (nfft, nfft) complex ndarray, fftshift-ed.
    waxis : (nfft,) ndarray of frequencies in cycles/sample.
    """
    X = _records(y, nfft, nsamp, overlap)
    I1, I2, I12 = _index_grids(nfft)
    B = np.zeros((nfft, nfft), dtype=complex)
    for Xr in X:
        B += Xr[I1] * Xr[I2] * np.conj(Xr[I12])
    B /= X.shape[0]
    waxis = np.fft.fftshift(np.fft.fftfreq(nfft))
    return np.fft.fftshift(B), waxis


def bicoherence(y, nfft=256, nsamp=256, overlap=50, wind=None):
    """Kim-Powers normalized bicoherence (real, in [0, 1]).

    b(f1,f2) = |E X(f1)X(f2)X*(f1+f2)|
               / sqrt( E|X(f1)X(f2)|^2 * E|X(f1+f2)|^2 ).

    Returns ``(bic, waxis)`` with the same axis convention as ``bispectrumd``.
    """
    X = _records(y, nfft, nsamp, overlap)
    I1, I2, I12 = _index_grids(nfft)
    num = np.zeros((nfft, nfft), dtype=complex)
    den1 = np.zeros((nfft, nfft))
    den2 = np.zeros((nfft, nfft))
    for Xr in X:
        tri = Xr[I1] * Xr[I2]
        num += tri * np.conj(Xr[I12])
        den1 += np.abs(tri) ** 2
        den2 += np.abs(Xr[I12]) ** 2
    denom = np.sqrt(den1 * den2)
    bic = np.abs(num) / np.where(denom > 0, denom, np.inf)
    waxis = np.fft.fftshift(np.fft.fftfreq(nfft))
    return np.fft.fftshift(bic), waxis


def acf_np(x, nlags):
    """Biased sample autocorrelation for lags 0..nlags (NumPy only)."""
    x = np.asarray(x, dtype=float)
    x = x - x.mean()
    n = len(x)
    var = np.dot(x, x) / n
    out = np.empty(nlags + 1)
    out[0] = 1.0
    for k in range(1, nlags + 1):
        out[k] = np.dot(x[:-k], x[k:]) / n / var
    return out
