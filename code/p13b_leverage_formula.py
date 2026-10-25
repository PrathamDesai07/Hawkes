"""Derived leverage formulas (K=1) vs measured.  Theta = m E[lam]/... = mu m/eps;  r = mu/eps.
 point leverage:  l_pt <= (eta r /(c_F (Theta-1))) * (E lam/lam* - 2 + m lam*)   (from phi^T Iev^{-1} phi = r (E lam - 2 lam* + m lam*^2)/(m E lam - 1)),
 gap leverage:    xbar <= eta/((1-eta) c_F) * max_i tr(Iev^{-1} dF_i)  and  tr(Iev^{-1} dJ_i) <= r*E_i*(E lam/mu^2 + m)/(Theta-1) (E_i = Lambda*(tau_i) ~ Exp(1))."""
import numpy as np, math, sys
from p13b_engine import run
from p14_mech import Iev_matrix
from p10_fisher_k1_exact import exact_m
mu_s = 1.0
for eps, eta, n, R in ((0.4, 0.004, 15000, 4), (0.2, 0.001, 30000, 3), (0.2, 0.004, 30000, 3), (0.4, 0.016, 15000, 3)):
    Iev = Iev_matrix(mu_s, 1.0, eps); m = exact_m(mu_s, 1.0, eps); r = mu_s/eps; Elam = mu_s/eps; Th = m*Elam*mu_s/mu_s; Th = mu_s*m/eps
    out = [run([1.0], eps, lambda t: mu_s, mu_s, eta, n, 700+s, Iev, mu_s) for s in range(R)]
    cat = lambda k: np.concatenate([o[k] for o in out]); cF = cat('cF'); lev = cat('lev'); x = cat('x'); tr = cat('trIdF')
    cmed, c5 = np.median(cF), np.quantile(cF, 0.05)
    lp = lambda c: eta*r/(c*(Th-1))*(Elam/mu_s - 2 + m*mu_s)
    print(f"eps={eps} eta={eta} (eta/eps^2={eta/eps**2:.3f}, Theta={Th:.3f}): point leverage measured max {lev.max():.4f} (99.9% {np.quantile(lev,0.999):.4f}) | formula with c_F median {lp(cmed):.4f}, with c_F 5% {lp(c5):.4f}")
    print(f"      gap leverage xbar measured max {x.max():.3f} (99.9% {np.quantile(x,0.999):.3f}) | formula eta/((1-eta)c_F) max tr(I^-1 dF): median c_F {eta/((1-eta)*cmed)*tr.max():.3f}, 5% c_F {eta/((1-eta)*c5)*tr.max():.3f}  (max tr(I^-1 dF) = {tr.max():.2f})", flush=True)
