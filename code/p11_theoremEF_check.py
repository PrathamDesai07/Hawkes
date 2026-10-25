"""Numerical consequences of Theorems E and F (any K).
E (k_K>1): E[1/W] <= gam/Delta, Delta = mu-c0-2kb*gam-(mu/2+C2)/z0 > 0;  E[1/lam] <= max(z0/mu, 1/(b-gam)) * gam/Delta  (optimise z0).
F (k_K<1): E[1/lam] >= c * eps^(k_K+eta): check the slope d ln m / d ln eps ~ k_K by Monte Carlo."""
import numpy as np
from p11_corrector_check import constants
from p2_identifiability import mc

def boundE(mu, beta, alpha):
    best = np.inf
    for z0 in np.logspace(0, 5, 400):
        c = constants(mu, beta, alpha, z0)
        D = mu - c['c0'] - 2*c['kb']*c['gam'] - (mu/2+c['C2'])/z0
        if D > 0:
            best = min(best, max(z0/mu, 1/(beta.min()-c['gam']))*c['gam']/D)
    return best

def kK(mu, beta, alpha):
    c = constants(mu, beta, alpha, 10.0); return 2*mu*c['nut']/c['n']

print("Theorem E: explicit bound vs Monte Carlo m (k_K>1)")
for beta, mu, eps in ((np.array([1., 2.]), 0.8, 0.03), (np.array([1., 2.]), 0.8, 0.01), (np.array([1., 2.]), 1.5, 0.03), (np.array([1., 2., 4.]), 2.5, 0.03)):
    K = len(beta); alpha = (1-eps)*beta/K
    I, _ = mc(mu, beta, eps, 250000, seed=2)
    print(f"  beta={beta.tolist()} mu={mu} eps={eps}: k_K={kK(mu,beta,alpha):.2f}  bound={boundE(mu,beta,alpha):.4f} (={boundE(mu,beta,alpha)/eps:.1f} eps)  MC m={I[0,0]:.4f} (={I[0,0]/eps:.2f} eps)", flush=True)

print("\nTheorem F: slope of ln m vs ln eps should approach k_K (<1)")
for beta, mu in ((np.array([1., 2.]), 0.4), (np.array([1., 3., 9.]), 1.0)):
    K = len(beta); ms = {}
    for eps in (0.1, 0.03, 0.01):
        alpha = (1-eps)*beta/K; I, _ = mc(mu, beta, eps, 200000 if eps > 0.01 else 500000, seed=4); ms[eps] = I[0, 0]
    es = sorted(ms, reverse=True)
    kk = [kK(mu, beta, (1-e)*beta/K) for e in es]
    sl = [np.log(ms[es[i]]/ms[es[i-1]])/np.log(es[i]/es[i-1]) for i in range(1, 3)]
    print(f"  beta={beta.tolist()} mu={mu}: k_K at eps={es}: {[round(x,3) for x in kk]};  m={[round(ms[e],4) for e in es]}; local slopes {[round(s,3) for s in sl]}", flush=True)
