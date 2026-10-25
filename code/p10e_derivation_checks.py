"""Independent checks of (i) the stationary Laplace-transform ODE and (ii) the reduction I_mm=m, I_ma=(1-mu m)/alpha, I_aa=(E lam-2mu+mu^2 m)/alpha^2,
by direct path simulation of the K=1 Hawkes process (no formula from the note is used in the simulation)."""
import numpy as np
from p10b_transient import grids
xq, wq = np.polynomial.legendre.leggauss(24); xq = 0.5*(xq+1); wq = 0.5*wq
mu, beta, eps = 1.0, 1.0, 0.2; alpha = (1-eps)*beta
rng = np.random.default_rng(99); S = 0.0; T = 0.0; burn = 5000; n_ev = 400000
us = np.array([0.05, 0.2, 0.5, 1.0, 2.0])
phi_acc = np.zeros(len(us)); A = np.zeros((2, 2)); lam_acc = 0.0
for i in range(n_ev):
    lam = mu + alpha*S; w = 0.0
    while True:
        w += rng.exponential(1/lam)
        if rng.random()*lam < mu + alpha*S*np.exp(-beta*w): break
    if i >= burn:
        s = w*xq; Sn = S*np.exp(-beta*s); l = mu + alpha*Sn
        phi_acc += w*(wq@np.exp(-np.outer(l, us)).T.T) if False else w*(np.exp(-np.outer(l, us))*wq[:, None]).sum(axis=0)
        a = np.stack([np.ones_like(Sn), Sn], 1)
        A += w*((wq/l)[:, None, None]*a[:, :, None]*a[:, None, :]).sum(axis=0)
        lam_acc += w*(wq@l); T += w
    S = S*np.exp(-beta*w) + 1.0
phi_mc = phi_acc/T; I_mc = A/T; Elam = lam_acc/T
g = grids(mu, beta, eps, N=60000, lo=1e-13)
phi_ode = np.interp(np.log(us), g['lx'], g['phi']); m = g['Cphi'][-1]
print(f"E lam: MC {Elam:.3f}  theory mu/eps = {mu/eps:.3f}")
print("phi(u)=E exp(-u lam):  u, MC, ODE")
for u, a_, b_ in zip(us, phi_mc, phi_ode): print(f"   {u:4.2f}  {a_:.4f}  {b_:.4f}")
Ith = np.array([[m, (1-mu*m)/alpha], [(1-mu*m)/alpha, (mu/eps-2*mu+mu**2*m)/alpha**2]])
print("I matrix MC:\n", I_mc, "\nI from reduction + exact m:\n", Ith)
