"""
Chapter "Interpreting Vector Autoregressions" — the numerical example of
Hansen and Sargent (1991), "Two Difficulties in Interpreting Vector
Autoregressions," Section 1.

A single market for a good has price p_t and quantity q_t.  A representative
supplier and demander each solve a linear-quadratic dynamic problem; their Euler
equations, together with serially correlated supply and demand shocks
s_t = B_s(L) w_{st} and d_t = B_d(L) w_{dt}, pin down an equilibrium

    S(L) [q_t, p_t]' = R(L) [w_{dt}, w_{st}]'          (structural rep, (1.24))

driven by the two white noises (w_{dt}, w_{st}) that are *fundamental for agents*.
An econometrician who sees only (q_t, p_t) instead recovers the Wold / vector-
autoregression representation

    S(L) [q_t, p_t]' = R*(L) eps*_t                    (Wold rep, (1.25))

driven by innovations eps*_t that are fundamental for (q_t, p_t).  Because the
zeros of det R(z) lie inside the unit circle here, eps*_t differs from the agents'
shocks: the VAR innovations are NOT the structural shocks.

We solve the model with QuantEcon's `DLE` class, mapping it into the class of
linear-quadratic economies of Hansen and Sargent (1990) exactly as in the
chapter's appendix.  This is a NumPy/QuantEcon port of the authors' MATLAB
routine `twodiff1.m`.

The figure has three rows:
  Row 1 (structural shocks, rep (1.24)): response of (q, p) to a unit demand
        shock and to a unit supply shock.
  Row 2 (Wold innovations, rep (1.25)): response of (q, p) to the two fundamental
        innovations recovered by a VAR (Gram-Schmidt with q ordered first).
  Row 3 (innovations to innovations): response of the Wold innovations eps*_t to
        the structural shocks (w_{dt}, w_{st}); the demand shock is impounded
        almost contemporaneously, the supply shock as a distributed lag.

Output: ../figures/ch36a_two_difficulties.png
"""

from pathlib import Path

import numpy as np
from scipy.linalg import companion
from scipy.signal import convolve
import quantecon as qe
from quantecon import DLE
import matplotlib.pyplot as plt

FIG_DIR = Path(__file__).resolve().parent.parent / "figures"
FIG_DIR.mkdir(exist_ok=True)

T = 18
VAR_A = 0.5      # variance of the shock on channel 0 (the "price-moving" shock)
VAR_B = 4.0      # variance of the shock on channel 1


def build_economy():
    """Map the supply/demand model into a QuantEcon DLE (Hansen-Sargent 1990)."""
    beta = 1 / 1.05
    hd, hs, gs, gd = 1.0, 1.0, 10.0, 0.1
    g = np.array([1, 0.8, 0.6, 0.4, 0.2])              # a(L) coefficients
    # shock-process polynomials B_d(L), B_s(L)
    bd = convolve(convolve([1, 0.6], [1, 0.4]), [1, 0.2])
    bs = convolve(convolve([1, -0.8], [1, 0.4]), [1, 0.2])

    Deltah = companion(np.array([1, 0, 0, 0, 0]))
    Thetah = np.array([1, 0, 0, 0], dtype=float)[:, None]
    Gamma = np.array([1, 0, 0], dtype=float)[:, None]
    A22 = np.block([[companion(np.array([1, 0, 0, 0, 0])), np.zeros((4, 4))],
                    [np.zeros((4, 4)), companion(np.array([1, 0, 0, 0, 0]))]])
    C2 = np.array([[1, 0], [0, 0], [0, 0], [0, 0],
                   [0, 1], [0, 0], [0, 0], [0, 0]], dtype=float)

    f1 = np.sqrt(gd); f2 = np.sqrt(hd); f4 = -1 / np.sqrt(hd)
    f5 = np.sqrt(hs); f6 = 1 / np.sqrt(hs); f7 = np.sqrt(gs)
    Lambda = np.vstack([np.zeros(4), f1 * g[1:]])
    Phic = np.array([1, -f5, 0])[:, None]
    Phii = np.array([-1, 0, -f7])[:, None]
    Phig = np.array([[0, 0], [1, 0], [0, 1]], dtype=float)
    Pih = np.array([f2, f1 * g[0]])[:, None]
    Ud = np.vstack([np.zeros(8), np.append(np.zeros(4), bs / f4), np.zeros(8)])
    Ub = np.vstack([np.append(bd / f6, np.zeros(4)), np.zeros(8)])

    info = (A22, C2, Ub, Ud)
    tech = (Phic, Phig, Phii, Gamma, np.array([[1.0]]), np.array([[1.0]]))
    pref = (np.array([[beta]]), Lambda, Pih, Deltah, Thetah)
    return DLE(info, tech, pref)


def structural_irf(econ):
    """Response of (q, p) to the two structural shocks: rep (1.24)."""
    G = np.vstack([econ.Sc, econ.Mc])          # rows: quantity, price
    C, A0 = econ.C, econ.A0
    imp = [np.array([[VAR_A], [0]]), np.array([[0], [VAR_B]])]
    out = []
    for s in imp:
        out.append(np.array([(G @ np.linalg.matrix_power(A0, t) @ C @ s)[:, 0]
                             for t in range(T + 1)]))
    return out                                  # [to shock A, to shock B]


def wold_irf(econ, price_first=False):
    """Response of (q, p) to the fundamental (Wold) innovations: rep (1.25)."""
    G = np.vstack([econ.Mc, econ.Sc]) if price_first else np.vstack([econ.Sc, econ.Mc])
    lss = qe.LinearStateSpace(econ.A0, econ.C, G, 1e-10 * np.eye(2))
    kal = qe.Kalman(lss)
    ma = kal.stationary_coefficients(T + 1, "ma")
    cov = kal.stationary_innovation_covar()
    out = []
    for j, v in enumerate((VAR_A, VAR_B)):
        out.append(np.array([ma[t][:, j] * np.sqrt(cov[j, j] / v)
                             for t in range(T + 1)]))
    return out, cov


def innovation_irf(econ):
    """Response of the Wold innovations eps*_t to the structural shocks: rep (1.16)."""
    G = np.vstack([econ.Sc, econ.Mc])
    lss = qe.LinearStateSpace(econ.A0, econ.C, G, 1e-10 * np.eye(2))
    kal = qe.Kalman(lss)
    w_lss = kal.whitener_lss()
    wirf = w_lss.impulse_response(T)[1]
    out = []
    for j, v in enumerate((VAR_A, VAR_B)):
        out.append(np.array([wirf[t][:, j] * v for t in range(T + 1)]))
    return out


def main():
    econ = build_economy()
    s_irf = structural_irf(econ)
    w_irf, covstar = wold_irf(econ)
    i_irf = innovation_irf(econ)

    G = np.vstack([econ.Sc, econ.Mc])
    R0 = G @ econ.C @ econ.C.T @ G.T            # structural contemp covariance
    t = np.arange(T + 1)

    plt.rcParams.update({"font.size": 10})
    fig, ax = plt.subplots(3, 2, figsize=(12, 11))

    titles = [
        ("Structural shock (1.24): demand innovation $w_{dt}$",
         "Structural shock (1.24): supply innovation $w_{st}$"),
        ("Wold innovation (1.25): first fundamental shock $\\epsilon^*_{1t}$",
         "Wold innovation (1.25): second fundamental shock $\\epsilon^*_{2t}$"),
        ("Response of Wold innovations to the demand shock $w_{dt}$",
         "Response of Wold innovations to the supply shock $w_{st}$"),
    ]
    series = [(s_irf, ["quantity", "price"]),
              (w_irf, ["quantity", "price"]),
              (i_irf, ["innovation in $q$", "innovation in $p$"])]

    for r, (data, names) in enumerate(series):
        for c in range(2):
            a = ax[r, c]
            a.axhline(0, color="0.7", lw=0.7)
            a.plot(t, data[c][:, 0], color="C0", lw=1.8, label=names[0])
            a.plot(t, data[c][:, 1], color="C3", lw=1.8, ls="--", label=names[1])
            a.set_title(titles[r][c], fontsize=10)
            a.set_xlim(0, T)
            a.set_xlabel("lag")
            a.legend(fontsize=9, loc="best")

    fig.suptitle("Two difficulties in interpreting VARs: agents' shocks vs. the "
                 "innovations a VAR recovers", fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.975))
    out = FIG_DIR / "ch36a_two_difficulties.png"
    fig.savefig(out, dpi=130)
    print("wrote", out)
    print("structural contemp covariance R0 =\n", np.round(R0, 4))
    print("Wold innovation covariance R0* =\n", np.round(covstar, 4))
    print("eig(R0* - R0) =", np.round(np.linalg.eigvalsh(covstar - R0), 5),
          " (>= 0 => inequality (1.19))")


if __name__ == "__main__":
    main()
