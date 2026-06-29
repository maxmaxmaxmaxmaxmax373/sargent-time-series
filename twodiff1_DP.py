# Python version of twodiff1.m

import numpy as np
from scipy.linalg import companion
from scipy.signal import convolve

import quantecon as qe
from quantecon import DLE

import matplotlib.pyplot as plt
# Base Parameters
beta = 1/1.05
hd = 1
hs = 1
gs = 10
gd = 0.1
a2 = np.sqrt(1/(2*hs))
g = np.array([1, 0.8, 0.6, 0.4, 0.2])
lam1 = 0.6
lam2 = 0.4
lam3 = 0.2
mu1 = -0.8
mu2 = 0.4
mu3 = 0.2
bd1 = np.array([1, lam1])
bd2 = np.array([1, lam2])
bd3 = np.array([1, lam3])
bs1 = np.array([1, mu1])
bs2 = np.array([1, mu2])
bs3 = np.array([1, mu3])

deltak = 1
thetak = 1

Deltah = companion(np.array([1, 0, 0, 0, 0]))
Thetah = np.array([1, 0, 0, 0], dtype=float)[:,np.newaxis]
Gamma = np.array([1, 0, 0], dtype=float)[:,np.newaxis]

A22 = np.block([[companion(np.array([1, 0, 0, 0, 0])), np.zeros((4, 4))],
                [np.zeros((4, 4)), companion(np.array([1, 0, 0, 0, 0]))]])
C2 = np.array([[1, 0], [0, 0], [0, 0], [0, 0], [0, 1], [0, 0], [0, 0], [0, 0]], dtype=float)

# Derived Parameters
f1 = np.sqrt(gd)
f2 = np.sqrt(hd)
f4 = -1/np.sqrt(hd)
f5 = np.sqrt(hs)
f6 = 1/np.sqrt(hs)
f7 = np.sqrt(gs)

Lambda = np.vstack([np.zeros(4), f1 * g[1:]])

Phic = np.array([1, -f5, 0])[:, np.newaxis]
Phii = np.array([-1, 0, -f7])[:, np.newaxis]
Phig = np.array([[0, 0],[1, 0],[0, 1]], dtype=float)

Pih = np.array([f2, f1*g[0]])[:,np.newaxis]

bb = convolve(bs1, bs2)
bs = convolve(bs3, bb)
bb = convolve(bd1, bd2)
bd = convolve(bb, bd3)

Ud = np.vstack([np.zeros(8), np.append(np.zeros(4), bs/f4), np.zeros(8)])
Ub = np.vstack([np.append(bd/f6, np.zeros(4)), np.zeros(8)])

beta = np.array([[beta]])
deltak = np.array([[deltak]])
thetak = np.array([[thetak]])

info = (A22, C2, Ub, Ud)
tech = (Phic, Phig, Phii, Gamma, deltak, thetak)
pref = (beta, Lambda, Pih, Deltah, Thetah)
econ = DLE(info, tech, pref)

# Figure 1
G = np.vstack([econ.Sc, econ.Mc])
C = econ.C
A0 = econ.A0

T = 18

var1 = 0.5
var2 = 4

shock1 = np.array([[var1], [0]])
shock2 = np.array([[0], [var2]])

shocks = [shock1, shock2]

irf_a = np.zeros((T+1, 2))
irf_b = np.zeros((T+1, 2))

for t in np.arange(T+1):
    irf_a[t, :] = (G@np.linalg.matrix_power(A0, t)@C@shocks[0])[:,0]
    irf_b[t, :] = (G@np.linalg.matrix_power(A0, t)@C@shocks[1])[:,0]

fig, ax = plt.subplots(1, 2)
handles = ax[0].plot(irf_a)
ax[0].set_title("1.A")
ax[0].set_xticks(np.arange(0, 20, 2))
ax[0].set_yticks(np.arange(-0.1, 0.6, 0.1))
ax[0].margins(x=0, y=0.05)
ax[0].legend(handles, ["Quantity", "Price"])

handles = ax[1].plot(irf_b)
ax[1].set_title("1.B")
ax[1].set_xticks(np.arange(0, 20, 2))
ax[1].set_yticks(np.arange(-0.25, 0.25, 0.05))
ax[1].margins(x=0, y=0.05)
ax[1].legend(handles, ["Quantity", "Price"])
fig.show()

# Figure 2
G = np.vstack([econ.Sc, econ.Mc])
H = 1e-8 * np.eye(2)
lss = qe.LinearStateSpace(A0, C, G, H)

hs_kal = qe.Kalman(lss)
w_lss = hs_kal.whitener_lss()
ma_coefs = hs_kal.stationary_coefficients(T+1, "ma")
cov = hs_kal.stationary_innovation_covar()

irf_a = np.zeros((T+1, 2))
irf_b = np.zeros((T+1, 2))

for t in range(T+1):
    irf_a[t, :] = ma_coefs[t][:, 0] * np.sqrt(cov[0, 0]/var1)
    irf_b[t, :] = ma_coefs[t][:, 1] * np.sqrt(cov[1, 1]/var2)

fig, ax = plt.subplots(1, 2)
handles = ax[0].plot(irf_a)
ax[0].set_title("2.A")
ax[0].set_xticks(np.arange(0, 20, 2))
ax[0].set_yticks(np.arange(-0.05, 0.3, 0.05))
ax[0].margins(x=0, y=0.05)
ax[0].legend(handles, ["Quantity", "Price"])

handles = ax[1].plot(irf_b)
ax[1].set_title("2.B")
ax[1].set_xticks(np.arange(0, 20, 2))
ax[1].set_yticks(np.arange(-0.1, 0.7, 0.1))
ax[1].margins(x=0, y=0.05)
ax[1].legend(handles, ["Quantity", "Price"])
fig.show()

# Figure 3
G = np.vstack([econ.Mc, econ.Sc])
H = 1e-8 * np.eye(2)
lss = qe.LinearStateSpace(A0, C, G, H)

hs_kal = qe.Kalman(lss)
w_lss = hs_kal.whitener_lss()
ma_coefs = hs_kal.stationary_coefficients(T+1, "ma")
cov = hs_kal.stationary_innovation_covar()

irf_a = np.zeros((T+1, 2))
irf_b = np.zeros((T+1, 2))

for t in range(T+1):
    irf_a[t, :] = ma_coefs[t][:, 0] * np.sqrt(cov[0, 0]/var2)
    irf_b[t, :] = ma_coefs[t][:, 1] * np.sqrt(cov[1, 1]/var1)

fig, ax = plt.subplots(1, 2)
handles = ax[0].plot(irf_a)
ax[0].set_title("3.A")
ax[0].set_xticks(np.arange(0, 20, 2))
ax[0].set_yticks(np.arange(-0.1, 0.7, 0.1))
ax[0].margins(x=0, y=0.05)
ax[0].legend(handles, ["Price", "Quantity"])

handles = ax[1].plot(irf_b)
ax[1].set_title("3.B")
ax[1].set_xticks(np.arange(0, 20, 2))
ax[1].set_yticks(np.arange(-0.05, 0.3, 0.05))
ax[1].margins(x=0, y=0.05)
ax[1].legend(handles, ["Price", "Quantity"])
fig.show()

# Figure 4
G = np.vstack([econ.Sc, econ.Mc])
H = 1e-8 * np.eye(2)
lss = qe.LinearStateSpace(A0, C, G, H)

hs_kal = qe.Kalman(lss)
w_lss = hs_kal.whitener_lss()
ma_coefs = hs_kal.stationary_coefficients(T+1, "ma")
cov = hs_kal.stationary_innovation_covar()
whitener_irf = w_lss.impulse_response(T)[1]

irf_a = np.zeros((T+1, 2))
irf_b = np.zeros((T+1, 2))

for t in range(T+1):
    irf_a[t, :] = whitener_irf[t][:, 0] * var1
    irf_b[t, :] = whitener_irf[t][:, 1] * var2

fig, ax = plt.subplots(1, 2)
handles = ax[0].plot(irf_a)
ax[0].set_title("4.A")
ax[0].set_xticks(np.arange(0, 20, 2))
ax[0].set_yticks(np.arange(-0.05, 0.5, 0.05))
ax[0].margins(x=0, y=0.05)
ax[0].legend(handles, ["Innov in q", "Innov in p"])

handles = ax[1].plot(irf_b)
ax[1].set_title("4.B")
ax[1].set_xticks(np.arange(0, 20, 2))
ax[1].set_yticks(np.arange(-0.25, 0.25, 0.05))
ax[1].margins(x=0, y=0.05)
ax[1].legend(handles, ["Innov in q", "Innov in p"])
fig.show()