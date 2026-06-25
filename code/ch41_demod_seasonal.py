"""
Chapter 41, Figure 1 — Complex demodulation of a *moving* seasonal.

We build a synthetic monthly series whose 12-month seasonal component has an
amplitude and a phase that drift slowly through time,

    x_t = A_t cos(omega0 t + phi_t) + (low-frequency cycle) + (harmonic) + noise,

with omega0 = 2*pi/12.  Complex demodulation at omega0 recovers the *time-
varying* amplitude R_t = 2|d_t| and phase theta_t = arg(d_t) of the seasonal
band, where

    d_t = W(L)[ x_t e^{i omega0 t} ],   s'_t = W(L)[x_t cos omega0 t],
                                        s''_t = W(L)[x_t sin omega0 t],
    d_t = s'_t + i s''_t.

W(L) is the low-pass filter Sargent (1968) used: a centred (2m+1)-term moving
average iterated twice (a Bartlett/triangular window, Godfrey 1965 eq. 2.7).
Remodulation r_t = 2[s'_t cos omega0 t + s''_t sin omega0 t] reconstructs the
seasonal band component; x_t - r_t is the deseasonalized series.

Panels: (a) series with the true seasonal overlaid; (b) recovered amplitude
R_t vs the true A_t (the "moving spectrum" at the seasonal band); (c) recovered
phase vs the true -phi_t; (d) remodulated seasonal vs truth.

Output: ../figures/ch41_demod_seasonal.png
"""

from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

BASE = Path(__file__).resolve().parent.parent
FIG_DIR = BASE / "figures"
FIG_DIR.mkdir(exist_ok=True)

OMEGA0 = 2 * np.pi / 12.0          # 12-month seasonal frequency
M = 12                             # half-width of each moving average
ITERS = 2                          # iterate twice -> Bartlett window


def lowpass(z, m=M, iters=ITERS):
    """Centred (2m+1)-term moving average, iterated `iters` times."""
    kernel = np.ones(2 * m + 1) / (2 * m + 1)
    out = np.asarray(z, dtype=complex)
    for _ in range(iters):
        out = np.convolve(out, kernel, mode="same")
    return out


def demodulate(x, omega0=OMEGA0):
    """Complex demodulate of real series x at frequency omega0.

    We use the e^{-i omega0 t} sign convention (consistent with the spectrum
    g_x(e^{-i omega}) of earlier chapters); Sargent's footnote uses e^{+i omega0 t},
    which only conjugates d_t (it flips the sign of the phase, leaving the
    amplitude R_t and the remodulated band component unchanged).
    """
    t = np.arange(len(x))
    d = lowpass(x * np.exp(-1j * omega0 * t))
    R = 2.0 * np.abs(d)            # amplitude  (factor 2, Sargent fn. 9)
    theta = np.angle(d)            # phase
    r = 2.0 * np.real(d * np.exp(1j * omega0 * t))   # remodulated band
    return d, R, theta, r


def main():
    rng = np.random.default_rng(41)
    N = 240
    t = np.arange(N)

    # Slowly drifting amplitude (a smooth bump) and phase.
    A = 0.5 + 1.3 * np.exp(-((t - 120) / 55.0) ** 2)
    phi = 0.7 * np.sin(2 * np.pi * t / 240.0)
    seasonal = A * np.cos(OMEGA0 * t + phi)

    # Contaminants the narrow low-pass filter should reject:
    low_cycle = 0.6 * np.cos(2 * np.pi * t / 120.0)      # 10-year cycle
    harmonic = 0.4 * np.cos(2 * OMEGA0 * t + 0.3)        # 6-month harmonic
    noise = 0.35 * rng.standard_normal(N)
    x = seasonal + low_cycle + harmonic + noise

    d, R, theta, r = demodulate(x)

    edge = 2 * M  # samples near each end are contaminated by filter wrap
    valid = slice(edge, N - edge)

    fig, axes = plt.subplots(2, 2, figsize=(12, 7.5))

    ax = axes[0, 0]
    ax.plot(t, x, color="0.6", lw=0.8, label=r"$x_t$")
    ax.plot(t, seasonal, color="C3", lw=1.4, label="true seasonal")
    ax.plot(t, A, color="C0", lw=1.0, ls="--", label=r"$\pm A_t$")
    ax.plot(t, -A, color="C0", lw=1.0, ls="--")
    ax.set_title("(a) series and its moving seasonal")
    ax.set_xlabel("month $t$")
    ax.legend(fontsize=8, ncol=2, loc="upper right")

    ax = axes[0, 1]
    ax.plot(t, A, color="C3", lw=1.6, label=r"true $A_t$")
    ax.plot(t, R, color="C0", lw=1.4, label=r"recovered $R_t=2|d_t|$")
    ax.axvspan(0, edge, color="0.9")
    ax.axvspan(N - edge, N, color="0.9")
    ax.set_title("(b) moving amplitude (spectrum at the seasonal band)")
    ax.set_xlabel("month $t$")
    ax.legend(fontsize=9, loc="upper right")

    ax = axes[1, 0]
    ax.plot(t, phi, color="C3", lw=1.6, label=r"true $\phi_t$")
    ax.plot(t, theta, color="C0", lw=1.4, label=r"recovered $\theta_t$")
    ax.axvspan(0, edge, color="0.9")
    ax.axvspan(N - edge, N, color="0.9")
    ax.set_title("(c) moving phase")
    ax.set_xlabel("month $t$")
    ax.legend(fontsize=9, loc="upper right")

    ax = axes[1, 1]
    ax.plot(t[valid], seasonal[valid], color="C3", lw=1.6, label="true seasonal")
    ax.plot(t[valid], r[valid], color="C0", lw=1.1, ls="--",
            label=r"remodulated $r_t$")
    ax.plot(t[valid], (x - r)[valid], color="0.55", lw=0.8,
            label=r"deseasonalized $x_t-r_t$")
    ax.set_title("(d) remodulation recovers / removes the seasonal")
    ax.set_xlabel("month $t$")
    ax.legend(fontsize=8, loc="upper right")

    fig.suptitle(
        r"Complex demodulation at $\omega_0=2\pi/12$:"
        r" tracking a seasonal whose amplitude and phase drift in time",
        fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    out = FIG_DIR / "ch41_demod_seasonal.png"
    fig.savefig(out, dpi=130)
    print("wrote", out)

    # quick self-check on the interior
    err = np.max(np.abs((r - seasonal)[valid]))
    print(f"max remodulation error (interior): {err:.3f}")
    print(f"max amplitude error (interior): "
          f"{np.max(np.abs((R - A)[valid])):.3f}")


if __name__ == "__main__":
    main()
