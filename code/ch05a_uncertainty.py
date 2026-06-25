"""
Chapter 05a, Figure — the time-frequency uncertainty principle for discrete-time
Fourier-transform pairs.

A square-summable sequence {c_k} and its Fourier transform f(omega) = sum_k c_k
e^{-i omega k} cannot both be sharply concentrated.  Treating |c_k|^2 and
|f(omega)|^2 as (energy-normalized) distributions over time index k and over
frequency omega in [-pi, pi], define the dispersions

    Dt^2 = sum_k (k - kbar)^2 |c_k|^2 / sum_k |c_k|^2,
    Dw^2 = (1/2pi) int (omega - wbar)^2 |f|^2 d omega / [(1/2pi) int |f|^2 d omega].

The discrete-time Heisenberg bound is  Dt * Dw >= (1/2)(1 - |f(pi)|^2 / ||f||^2),
which becomes the familiar  Dt * Dw >= 1/2  when the transform vanishes at the
Nyquist frequency.  Gaussian sequences nearly attain the bound; rectangular
windows sit far above it.

Panels:
 (a) |c_k| for three Gaussian sequences -- narrow, medium, wide in time;
 (b) the corresponding |f(omega)| -- wide, medium, narrow in frequency (the
     duality);
 (c) the uncertainty product Dt*Dw versus width, for the Gaussian family (which
     hugs the 1/2 floor) and the rectangular family (which does not).

Output: ../figures/ch05a_uncertainty.png
"""

from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

BASE = Path(__file__).resolve().parent.parent
FIG_DIR = BASE / "figures"
FIG_DIR.mkdir(exist_ok=True)

K = 400                                   # time half-range
k = np.arange(-K, K + 1)
NW = 8192
omega = np.linspace(-np.pi, np.pi, NW, endpoint=False)
DW = omega[1] - omega[0]
EXP = np.exp(-1j * np.outer(k, omega))    # (len k, len omega)


def ft(c):
    return c @ EXP


def dispersions(c):
    c = np.asarray(c, float)
    p = c ** 2
    p = p / p.sum()                       # |c_k|^2 normalized
    kbar = (k * p).sum()
    dt2 = ((k - kbar) ** 2 * p).sum()
    f = ft(c)
    pf = np.abs(f) ** 2
    norm = pf.sum() * DW / (2 * np.pi)    # == sum|c_k|^2 == 1 (Parseval)
    wbar = (omega * pf).sum() * DW / (2 * np.pi) / norm
    dw2 = ((omega - wbar) ** 2 * pf).sum() * DW / (2 * np.pi) / norm
    fpi2 = np.abs(f[0]) ** 2 / (pf.sum() * DW / (2 * np.pi))   # |f(-pi)|^2/||f||^2
    return np.sqrt(dt2), np.sqrt(dw2), f, fpi2


def gaussian(s):
    return np.exp(-(k ** 2) / (2.0 * s ** 2))


def rect(m):
    return (np.abs(k) <= m).astype(float)


def main():
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False,
                         "axes.spines.right": False})
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.4))

    # --- (a),(b) three Gaussian pairs --------------------------------------
    widths = [(2.0, "C0", "narrow in time"),
              (6.0, "C1", "medium"),
              (20.0, "C2", "wide in time")]
    for s, c, lab in widths:
        cg = gaussian(s)
        cg = cg / np.sqrt((cg ** 2).sum())
        f = ft(cg)
        fn = np.abs(f) / np.abs(f).max()
        axes[0].plot(k, cg / cg.max(), color=c, lw=1.4, label=f"s={s:g} ({lab})")
        axes[1].plot(omega, fn, color=c, lw=1.4, label=f"s={s:g}")

    axes[0].set_xlim(-70, 70)
    axes[0].set_title("(a) time domain $|c_k|$ (normalized)")
    axes[0].set_xlabel("time index $k$")
    axes[0].legend(fontsize=8)

    axes[1].set_xlim(-np.pi, np.pi)
    axes[1].set_xticks([-np.pi, 0, np.pi])
    axes[1].set_xticklabels([r"$-\pi$", "0", r"$\pi$"])
    axes[1].set_title(r"(b) frequency domain $|f(\omega)|$ (normalized)")
    axes[1].set_xlabel(r"frequency $\omega$")
    axes[1].legend(fontsize=8)

    # --- (c) uncertainty product versus width ------------------------------
    s_grid = np.linspace(1.0, 40.0, 40)
    g_prod = []
    for s in s_grid:
        dt, dw, _, _ = dispersions(gaussian(s))
        g_prod.append(dt * dw)
    m_grid = np.arange(1, 60)
    r_prod = []
    for m in m_grid:
        dt, dw, _, _ = dispersions(rect(m))
        r_prod.append(dt * dw)

    axes[2].plot(s_grid, g_prod, "o-", color="C0", ms=3, lw=1.0,
                 label="Gaussian sequences")
    axes[2].plot(m_grid / np.sqrt(3), r_prod, "s-", color="C3", ms=3, lw=1.0,
                 label="rectangular windows")
    axes[2].axhline(0.5, color="0.3", lw=1.2, ls="--",
                    label=r"floor $\Delta_t\Delta_\omega=\frac{1}{2}$")
    axes[2].set_ylim(0, 3.0)
    axes[2].set_title(r"(c) uncertainty product $\Delta_t\,\Delta_\omega$")
    axes[2].set_xlabel(r"time dispersion $\Delta_t$")
    axes[2].legend(fontsize=8, loc="upper left")

    fig.suptitle("The uncertainty principle: a sequence and its Fourier transform "
                 "cannot both be concentrated", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    out = FIG_DIR / "ch05a_uncertainty.png"
    fig.savefig(out, dpi=130)
    print("wrote", out)

    # report a few numbers used in the text
    for s in [2, 6, 20]:
        dt, dw, _, fpi2 = dispersions(gaussian(s))
        print(f"Gaussian s={s:2d}: Dt={dt:7.3f}  Dw={dw:6.3f}  "
              f"Dt*Dw={dt*dw:6.3f}  |f(pi)|^2/||f||^2={fpi2:.2e}")
    for m in [5, 20]:
        dt, dw, _, fpi2 = dispersions(rect(m))
        print(f"rect    m={m:2d}: Dt={dt:7.3f}  Dw={dw:6.3f}  Dt*Dw={dt*dw:6.3f}")


if __name__ == "__main__":
    main()
