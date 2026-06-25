"""
Chapter 41, Figure 2 — Locating a spectral peak with the demodulate phase.

Godfrey (1965, eqs. 2.9-2.12): a pure component x_t = a cos(omega1 t),
demodulated at a trial frequency omega0 and passed through an ideal low-pass
filter, yields a demodulate whose phase is *linear in time*,

    Phi_t = (omega0 - omega1) t,

with slope equal to the offset between the trial frequency and the true peak.
Demodulating just below, exactly on, and just above the peak therefore gives a
phase line that slopes down, is flat, and slopes up.  This turns peak-finding
into reading off the slope of a straight line — and shows how a demodulate that
is *not* centred on a peak still drifts in phase even for a perfectly
stationary sinusoid.

Output: ../figures/ch41_demod_peak.png
"""

from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

BASE = Path(__file__).resolve().parent.parent
FIG_DIR = BASE / "figures"
FIG_DIR.mkdir(exist_ok=True)

M = 12
ITERS = 2


def lowpass(z, m=M, iters=ITERS):
    kernel = np.ones(2 * m + 1) / (2 * m + 1)
    out = np.asarray(z, dtype=complex)
    for _ in range(iters):
        out = np.convolve(out, kernel, mode="same")
    return out


def demod_phase(x, omega0):
    t = np.arange(len(x))
    d = lowpass(x * np.exp(-1j * omega0 * t))
    return d


def main():
    rng = np.random.default_rng(7)
    N = 240
    t = np.arange(N)

    period1 = 12.0
    omega1 = 2 * np.pi / period1
    x = np.cos(omega1 * t) + 0.05 * rng.standard_normal(N)

    trials = [
        (2 * np.pi / 12.8, "below peak (period 12.8)", "C0"),
        (omega1, "on peak (period 12.0)", "C3"),
        (2 * np.pi / 11.2, "above peak (period 11.2)", "C2"),
    ]

    edge = 2 * M
    valid = slice(edge, N - edge)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.6))

    for omega0, label, c in trials:
        d = demod_phase(x, omega0)
        phase = np.unwrap(np.angle(d))
        ax1.plot(t[valid], phase[valid], color=c, lw=1.6, label=label)
        # theoretical slope (omega1 - omega0); align intercept to the data
        theo = (omega1 - omega0) * t
        theo = theo - theo[valid][0] + phase[valid][0]
        ax1.plot(t[valid], theo[valid], color=c, lw=0.9, ls=":")
        ax2.plot(t[valid], 2 * np.abs(d)[valid], color=c, lw=1.6, label=label)

    ax1.set_title("(a) demodulate phase is linear in time")
    ax1.set_xlabel("month $t$")
    ax1.set_ylabel(r"$\arg d_t$ (unwrapped, rad)")
    ax1.legend(fontsize=8, loc="upper right")
    ax1.text(0.02, 0.03,
             r"slope $=\omega_1-\omega_0$;  dotted = theory",
             transform=ax1.transAxes, fontsize=8)

    ax2.set_title("(b) demodulate amplitude $2|d_t|$")
    ax2.set_xlabel("month $t$")
    ax2.legend(fontsize=8, loc="upper right")

    fig.suptitle(
        "Finding a spectral peak: a stationary sinusoid demodulated off-peak "
        "drifts in phase at a constant rate", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    out = FIG_DIR / "ch41_demod_peak.png"
    fig.savefig(out, dpi=130)
    print("wrote", out)


if __name__ == "__main__":
    main()
