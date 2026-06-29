"""
Chapter XI, Section 13 (Representation Theory) — Whittle's spectral
factorization by the cepstrum (Whittle, *Prediction and Regulation*, 2nd ed.,
U. of Minnesota Press, 1983, p. 26).

Given the spectral density g(e^{-i omega}) of an indeterministic covariance
stationary process, recover the *fundamental* Wold moving-average kernel
d(L) = 1 + d_1 L + d_2 L^2 + ... and the one-step-ahead prediction-error
variance sigma^2 = E eps_t^2, using nothing but log, FFT and inverse FFT.

This is a NumPy port of the MATLAB routine factor.m.  The factorization is

    g(e^{-i omega}) = sigma^2 |d(e^{-i omega})|^2 ,

and the construction reads it off the cepstrum c_k (the Fourier coefficients of
log g):  sigma^2 = exp(c_0) (the Kolmogorov formula), while the strictly
positive-power half of the cepstrum is log d(z), so d = exp of it.

Two examples:
  1. a NON-invertible MA(1), g = |1 + 2 e^{-i omega}|^2, whose fundamental
     factorization is (1 + 0.5 L) with sigma^2 = 4 -- the method picks the
     invertible factor automatically;
  2. an AR(2) with complex roots, whose Wold kernel is a damped sinusoid.

Output: ../figures/ch13_whittle_factor.png
"""

from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

FIG_DIR = Path(__file__).resolve().parent.parent / "figures"
FIG_DIR.mkdir(exist_ok=True)


def whittle_factor(s):
    """Whittle spectral factorization (Whittle 1983, p. 26).

    Parameters
    ----------
    s : array_like
        Spectral density g(e^{-i omega}) sampled on the equally spaced grid
        omega_j = 2*pi*j/N, j = 0, ..., N-1.  Must be real and positive.

    Returns
    -------
    d : ndarray
        Moving-average kernel d_0, d_1, d_2, ... of the fundamental Wold
        representation (d_0 == 1).
    sigma2 : float
        One-step-ahead prediction-error variance E eps_t^2 = exp(c_0).
    """
    s = np.asarray(s, dtype=float)
    N = s.size
    cep = np.fft.ifft(np.log(s))            # cepstrum c_k (Fourier coeffs of log g)
    sigma2 = float(np.exp(cep[0].real))     # Kolmogorov formula  sigma^2 = exp(c_0)
    half = N // 2
    p = np.zeros(N, dtype=complex)
    p[1:half + 1] = cep[1:half + 1]         # keep strictly positive powers => log d
    d = np.real(np.fft.ifft(np.exp(np.fft.fft(p))))
    return d, sigma2


def ar_impulse_response(a, n):
    """Impulse response d_j of x_t = sum_k a[k] x_{t-1-k} + eps_t."""
    d = np.zeros(n)
    d[0] = 1.0
    for j in range(1, n):
        d[j] = sum(a[k] * d[j - 1 - k] for k in range(len(a)) if j - 1 - k >= 0)
    return d


def main():
    N = 8192
    w = 2 * np.pi * np.arange(N) / N

    # --- Example 1: non-invertible MA(1)  g = |1 + 2 e^{-iw}|^2 ------------
    g1 = np.abs(1 + 2 * np.exp(-1j * w)) ** 2
    d1, s1 = whittle_factor(g1)

    # --- Example 2: AR(2) with complex roots ------------------------------
    a1, a2, sig = 0.6, -0.5, 2.5
    g2 = sig / np.abs(1 - a1 * np.exp(-1j * w) - a2 * np.exp(-2j * w)) ** 2
    d2, s2 = whittle_factor(g2)
    d2_true = ar_impulse_response([a1, a2], 16)

    plt.rcParams.update({"font.size": 11})
    fig, ax = plt.subplots(2, 2, figsize=(12, 8))
    ww = w[:N // 2 + 1]                                   # plot omega in [0, pi]

    # row 1, Example 1
    ax[0, 0].plot(ww, np.log(g1[:N // 2 + 1]), color="C0", lw=1.8)
    ax[0, 0].set_title(r"Example 1 spectrum: $g(e^{-i\omega})=|1+2e^{-i\omega}|^2$")
    ax[0, 0].set_xlabel(r"$\omega$")
    ax[0, 0].set_ylabel(r"$\log g$")
    ax[0, 0].set_xticks([0, np.pi / 2, np.pi])
    ax[0, 0].set_xticklabels(["0", r"$\pi/2$", r"$\pi$"])

    k = np.arange(6)
    ax[0, 1].stem(k, d1[:6], linefmt="C0-", markerfmt="C0o", basefmt="0.7")
    ax[0, 1].plot(k, [1, 0.5, 0, 0, 0, 0], "rx", ms=9, mew=2,
                  label="fundamental $(1+0.5L)$")
    ax[0, 1].set_title(r"recovered kernel $d_j$ — selects the invertible factor")
    ax[0, 1].set_xlabel("$j$")
    ax[0, 1].legend(fontsize=9)
    ax[0, 1].text(0.97, 0.80,
                  fr"$\sigma^2 = {s1:.3f}$" "\n(true $4$)",
                  transform=ax[0, 1].transAxes, ha="right", va="top",
                  bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="C0"))

    # row 2, Example 2
    ax[1, 0].plot(ww, np.log(g2[:N // 2 + 1]), color="C2", lw=1.8)
    ax[1, 0].set_title(r"Example 2 spectrum: AR(2), "
                       r"$x_t=0.6x_{t-1}-0.5x_{t-2}+\epsilon_t$")
    ax[1, 0].set_xlabel(r"$\omega$")
    ax[1, 0].set_ylabel(r"$\log g$")
    ax[1, 0].set_xticks([0, np.pi / 2, np.pi])
    ax[1, 0].set_xticklabels(["0", r"$\pi/2$", r"$\pi$"])

    k2 = np.arange(16)
    ax[1, 1].stem(k2, d2[:16], linefmt="C2-", markerfmt="C2o", basefmt="0.7")
    ax[1, 1].plot(k2, d2_true, "rx", ms=8, mew=1.8,
                  label="true impulse response")
    ax[1, 1].axhline(0, color="0.7", lw=0.7)
    ax[1, 1].set_title(r"recovered kernel $d_j$ — a damped sinusoid")
    ax[1, 1].set_xlabel("$j$")
    ax[1, 1].legend(fontsize=9)
    ax[1, 1].text(0.97, 0.80,
                  fr"$\sigma^2 = {s2:.3f}$" "\n(true $2.5$)",
                  transform=ax[1, 1].transAxes, ha="right", va="top",
                  bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="C2"))

    fig.suptitle("Whittle's spectral factorization: recovering the Wold kernel "
                 "$d(L)$ and $\\sigma^2$ from the spectral density", fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    out = FIG_DIR / "ch13_whittle_factor.png"
    fig.savefig(out, dpi=140)
    print("wrote", out)
    print(f"Ex1: d[:3] = {np.round(d1[:3], 4)}, sigma2 = {s1:.4f}")
    print(f"Ex2: d[:6] = {np.round(d2[:6], 4)}, sigma2 = {s2:.4f}")
    print(f"     true  = {np.round(d2_true[:6], 4)}")


if __name__ == "__main__":
    main()
