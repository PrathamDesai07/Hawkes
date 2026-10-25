"""Fully closed-form lower bound on lam_min(I) for the geometric grid (no simulation, no numerical eigenvalue of an unknown matrix):
   lam_min(I) >= (1/(2L)) [ LM - sqrt(E|a|^4) sqrt(E lam^p / L^p) ],     LM = (mu/(4 eps)) lam_G / (1 + 2 (mu/eps)^2 sum beta_k^-2),
   lam_G >= 1/ sum_ij 4 b_i b_j P_i P_j/(b_i+b_j)  (Cauchy inverse),  E Z^p <= C_p eps^-p (Lemma 1'),  E|a|^4 <= 2(1+(sum (beta_k/alpha_k)^2)^2 E Z^4),  E lam^p <= 2^(p-1)(mu^p + B^p E Z^p)."""
import numpy as np
from math import comb
def Cp_list(mu, b, B, P):
    C = [1.0]
    for p in range(1, P+1):
        a = {1: p*mu + comb(p, 2)*B}
        for i in range(2, p+1): a[i] = mu*comb(p, i) + B*(comb(p, i+1) if i+1 <= p else 0)
        C.append(sum(a[i]*C[p-i] for i in range(1, p+1))/(p*b))
    return C
def chain(K, r, eps, mu=1.0):
    beta = float(r)**np.arange(K); alpha = (1-eps)*beta/K; b, B = beta.min(), beta.max()
    P = np.ones(K)
    for i in range(K):
        for k in range(K):
            if k != i: P[i] *= (beta[i]+beta[k])/abs(beta[i]-beta[k])
    lG = 1.0/sum(4*beta[i]*beta[j]*P[i]*P[j]/(beta[i]+beta[j]) for i in range(K) for j in range(K))
    LM = (mu/(4*eps))*lG/(1+2*(mu/eps)**2*(beta**-2.0).sum())
    C = Cp_list(mu, b, B, 14); sk = ((beta/alpha)**2).sum()
    EA4 = 2*(1 + sk**2*C[4]/eps**4)
    best = 0.0; arg = None
    for p in range(2, 15):
        Elp = 2**(p-1)*(mu**p + B**p*C[p]/eps**p)
        for L in (mu/eps)*np.logspace(0, 40, 800):
            v = (LM - np.sqrt(EA4)*np.sqrt(Elp/L**p))/(2*L)
            if v > best: best, arg = v, (p, L*eps/mu)
    return LM, best, arg
mc = {(1, 2, .3): .1389, (1, 2, .1): .0654, (1, 2, .05): .0383, (2, 2, .3): .0337, (2, 2, .1): .0455, (2, 4, .1): .0677, (3, 2, .3): .0033, (3, 2, .1): .0045}
print(f"{'K':>2} {'r':>3} {'eps':>5} {'MC lmin(I)':>10} {'closed-form LM':>14} {'closed-form bound':>18} {'bound/MC':>10}  (p, L eps/mu)")
for (K, r, eps), m in mc.items():
    LM, bd, arg = chain(K, r, eps)
    print(f"{K:2d} {r:3d} {eps:5.2f} {m:10.4f} {LM:14.3e} {bd:18.3e} {bd/m:10.2e}  {arg}")
