"""
Chapter "A Difficulty in Interpreting Vector Autoregressions" — each market
participant as a price-taking discounted linear regulator (the "Big X, little x"
trick).

In the rational expectations equilibrium of the single-market example, the
representative supplier and the representative demander are both *price takers*:
each treats the equilibrium price process {p_t} as an exogenous stochastic
process it cannot influence, and chooses its own quantity to maximize a
discounted quadratic objective. This is exactly the logic of the QuantEcon
"dynamic Stackelberg" lecture, Section 50.7.1 ("Recursive formulation of a
follower's problem"): append the *exogenous* aggregate law of motion to the
agent's own state and solve an optimal linear regulator.

The equilibrium price process comes from the Hansen-Sargent `DLE` solution used
earlier in the chapter:

    X_{t+1} = A X_t + C w_{t+1},      p_t = G_p X_t ,

with A = econ.A0, C = econ.C, and G_p = econ.Mc (the DLE's price selector). Each
agent's composite state stacks this exogenous price-process state X_t ("Big X")
on top of its own lagged quantities ("little x"), and the optimal feedback rule
is the dynamic supply / demand curve of the previous subsection.

This script verifies that:
  * the supplier's regulator reproduces the dynamic-supply-curve own-state
    feedback  delta_s / beta  on q_{t-1};
  * the demander's regulator reproduces the dynamic-demand-curve own-state
    feedback coefficients gamma_{d,1..4} on q_{t-1}, ..., q_{t-4}.

It uses QuantEcon's `LQ` class. Run from the `code/` directory (it imports
`build_economy` from ch36a_two_difficulties.py).
"""

import sys
from pathlib import Path

import numpy as np
import quantecon as qe

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ch36a_two_difficulties import build_economy

BETA = 1 / 1.05
HS, GS = 1.0, 10.0
HD, GD = 1.0, 0.1
A_POLY = np.array([1, 0.8, 0.6, 0.4, 0.2])      # a(L) coefficients a_0,...,a_4


def supplier_regulator(A, C, Gp):
    """Supplier's discounted linear regulator facing exogenous price p = Gp X.

    State  hat_X = [X ; q_{t-1}],  control u = q_t.  Maximizes
        E sum beta^t { p_t q_t - (h_s/2) q_t^2 - (g_s/2)(q_t - q_{t-1})^2 } .
    Returns the optimal feedback rule q_t = -F hat_X (so coef on q_{t-1} = -F[-1]).
    """
    nx = A.shape[0]
    n = nx + 1
    Alq = np.zeros((n, n)); Alq[:nx, :nx] = A
    Blq = np.zeros((n, 1)); Blq[nx, 0] = 1.0
    Clq = np.zeros((n, 1)); Clq[:nx, 0] = np.asarray(C).ravel()[:nx]
    # cost = -[ p u - (h_s/2)u^2 - (g_s/2)(u - q_{t-1})^2 ]  (LQ minimizes)
    Q = np.array([[(HS + GS) / 2]])
    R = np.zeros((n, n)); R[nx, nx] = GS / 2
    N = np.zeros((1, n)); N[0, :nx] = -np.asarray(Gp).ravel() / 2; N[0, nx] = -GS / 2
    lq = qe.LQ(Q, R, Alq, Blq, C=Clq, N=N, beta=BETA)
    P, F, d = lq.stationary_values()
    return F.ravel()


def demander_regulator(A, C, Gp):
    """Demander's discounted linear regulator facing exogenous price p = Gp X.

    State  hat_X = [X ; q_{t-1}, q_{t-2}, q_{t-3}, q_{t-4}],  control u = q_t.
    Maximizes  E sum beta^t { -p_t q_t - (h_d/2)q_t^2 - (g_d/2)(a(L)q_t)^2 } .
    Returns q_t = -F hat_X (coefs on q_{t-1..4} = -F[nx:]).
    """
    nx = A.shape[0]
    n = nx + 4
    Alq = np.zeros((n, n)); Alq[:nx, :nx] = A
    # q-lag shift register: new [q_{t-1..4}] = [u, q_{t-1}, q_{t-2}, q_{t-3}]
    Alq[nx + 1, nx] = 1.0; Alq[nx + 2, nx + 1] = 1.0; Alq[nx + 3, nx + 2] = 1.0
    Blq = np.zeros((n, 1)); Blq[nx, 0] = 1.0
    Clq = np.zeros((n, 1)); Clq[:nx, 0] = np.asarray(C).ravel()[:nx]
    a0, v = A_POLY[0], A_POLY[1:5]
    Q = np.array([[(HD + GD * a0 ** 2) / 2]])
    R = np.zeros((n, n)); R[nx:nx + 4, nx:nx + 4] = (GD / 2) * np.outer(v, v)
    N = np.zeros((1, n)); N[0, :nx] = -np.asarray(Gp).ravel() / 2
    N[0, nx:nx + 4] = (GD * a0 / 2) * v
    lq = qe.LQ(Q, R, Alq, Blq, C=Clq, N=N, beta=BETA)
    P, F, d = lq.stationary_values()
    return F.ravel()


def dynamic_curve_coefficients():
    """The own-state feedback coefficients implied by the factored Euler equations."""
    delta_s = np.roots([1, -((1 + BETA) + HS / GS), BETA]).min()
    sup = delta_s / BETA
    # demander: stable roots delta_{d,i} and c_d(L)=prod(1 - delta/beta L)
    a = A_POLY
    p_abz = np.poly1d([a[k] * BETA ** k for k in range(5)])
    full = GD * (p_abz * np.poly1d(a[::-1])) + np.poly1d([HD, 0, 0, 0, 0])
    stab = sorted(np.roots(full.c), key=abs)[:4]
    cdL = np.array([1.0])
    for dd in stab:
        cdL = np.polynomial.polynomial.polymul(cdL, [1, -(dd / BETA)])
    gamma = -cdL[1:5].real
    return sup, gamma


def main():
    econ = build_economy()
    A, C, Gp = econ.A0, econ.C, econ.Mc      # equilibrium price process from the DLE
    print(f"DLE equilibrium state dimension: {A.shape[0]}")

    sup_target, gamma_target = dynamic_curve_coefficients()

    Fs = supplier_regulator(A, C, Gp)
    sup_lq = -Fs[-1]
    print("\nSUPPLIER (price-taking linear regulator):")
    print(f"  coef on q_(t-1) from LQ          = {sup_lq: .4f}")
    print(f"  delta_s/beta from supply curve   = {sup_target: .4f}  "
          f"-> match: {np.isclose(sup_lq, sup_target, atol=1e-3)}")

    Fd = demander_regulator(A, C, Gp)
    gamma_lq = -Fd[-4:]
    print("\nDEMANDER (price-taking linear regulator):")
    print(f"  coefs on q_(t-1..4) from LQ      = {np.round(gamma_lq, 4)}")
    print(f"  gamma_(d,1..4) from demand curve = {np.round(gamma_target, 4)}  "
          f"-> match: {np.allclose(gamma_lq, gamma_target, atol=1e-3)}")


if __name__ == "__main__":
    main()
