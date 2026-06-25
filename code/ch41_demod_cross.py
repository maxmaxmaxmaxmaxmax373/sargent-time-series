"""
Chapter 41, Figure 3 — A *moving cross-spectrum* by complex demodulation.

Two monthly series share a 12-month seasonal but the phase lead of y over x and
the strength of their seasonal comovement both drift through time -- exactly the
situation Sargent (1968) studied for the interest-rate spread and the average
maturity of the public debt.  Demodulating both series at omega0 = 2*pi/12 gives
demodulates d^y_t, d^x_t.  Smoothing the local products through time,

    S_yx,t = <d^y_t conj(d^x_t)>,  S_yy,t = <|d^y_t|^2>,  S_xx,t = <|d^x_t|^2>,

yields a *time-varying* cross-spectrum at the seasonal band, hence

    coherence_t = |S_yx,t|^2 / (S_yy,t S_xx,t),
    phase_t     = arg(S_yx,t),     lead in months = phase_t / omega0,
    gain_t      = |S_yx,t| / S_xx,t.

(The extra time-smoothing <.> is essential: without it the squared coherence of
a single demodulate pair is identically 1, just as for the ordinary spectrum.)

Output: ../figures/ch41_demod_cross.png
"""

from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

BASE = Path(__file__).resolve().parent.parent
FIG_DIR = BASE / "figures"
FIG_DIR.mkdir(exist_ok=True)

OMEGA0 = 2 * np.pi / 12.0
M = 10
ITERS = 2
SMOOTH = 18          # extra time-smoothing half-width for the cross-products


def lowpass(z, m=M, iters=ITERS):
    kernel = np.ones(2 * m + 1) / (2 * m + 1)
    out = np.asarray(z, dtype=complex)
    for _ in range(iters):
        out = np.convolve(out, kernel, mode="same")
    return out


def tsmooth(z, m=SMOOTH):
    kernel = np.ones(2 * m + 1) / (2 * m + 1)
    return np.convolve(np.asarray(z, dtype=complex), kernel, mode="same")


def demod(x, omega0=OMEGA0):
    t = np.arange(len(x))
    return lowpass(x * np.exp(-1j * omega0 * t))


def main():
    rng = np.random.default_rng(123)
    N = 300
    t = np.arange(N)

    # True seasonal lead of y over x, drifting between roughly -1 and +2 months.
    lead = 0.5 + 1.5 * np.sin(2 * np.pi * t / 300.0)          # months
    common = np.cos(OMEGA0 * t)

    # Coherence is degraded by independent seasonal "noise" that is strongest in
    # the middle third of the sample.
    noise_amp = 0.2 + 2.2 * np.exp(-((t - 150) / 40.0) ** 2)
    ny = noise_amp * np.cos(OMEGA0 * t + 2 * np.pi * rng.random(N))
    nx = noise_amp * np.cos(OMEGA0 * t + 2 * np.pi * rng.random(N))

    y = np.cos(OMEGA0 * t) + ny + 0.25 * rng.standard_normal(N)
    x = np.cos(OMEGA0 * (t - lead)) + nx + 0.25 * rng.standard_normal(N)

    dy, dx = demod(y), demod(x)
    Syx = tsmooth(dy * np.conj(dx))
    Syy = tsmooth(np.abs(dy) ** 2).real
    Sxx = tsmooth(np.abs(dx) ** 2).real

    coh = np.abs(Syx) ** 2 / (Syy * Sxx)
    phase = np.angle(Syx)
    lead_est = phase / OMEGA0                                 # months

    edge = 2 * M + SMOOTH
    valid = slice(edge, N - edge)

    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.2))

    ax = axes[0]
    w = slice(60, 120)
    ax.plot(t[w], y[w], color="C0", lw=1.1, label=r"$y_t$")
    ax.plot(t[w], x[w], color="C3", lw=1.1, label=r"$x_t$")
    ax.set_title("(a) two seasonal series (excerpt)")
    ax.set_xlabel("month $t$")
    ax.legend(fontsize=9)

    ax = axes[1]
    ax.plot(t[valid], coh[valid], color="C0", lw=1.6)
    ax.axvspan(105, 195, color="0.92")
    ax.set_ylim(0, 1.05)
    ax.set_title("(b) moving coherence at the seasonal band")
    ax.set_xlabel("month $t$")
    ax.text(150, 0.06, "noisy\nmiddle", ha="center", fontsize=8, color="0.4")

    ax = axes[2]
    ax.plot(t[valid], lead[valid], color="C3", lw=1.6, label="true lead")
    ax.plot(t[valid], lead_est[valid], color="C0", lw=1.4, ls="--",
            label="estimated lead")
    ax.axhline(0, color="0.7", lw=0.8)
    ax.set_title("(c) moving phase lead of $y$ over $x$")
    ax.set_xlabel("month $t$")
    ax.set_ylabel("lead (months)")
    ax.legend(fontsize=9)

    fig.suptitle(
        "Complex demodulation gives a cross-spectrum that moves through time: "
        "coherence and lead at the 12-month band", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    out = FIG_DIR / "ch41_demod_cross.png"
    fig.savefig(out, dpi=130)
    print("wrote", out)


if __name__ == "__main__":
    main()
