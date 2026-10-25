"""Theorem C (general K): if 2mu > n^2 beta_max then E_pi[1/lam] <= 2 eps B /(b (2mu - n^2 B)), b=beta_min, B=beta_max.
Proof device: W = Z + mu/B, Z = sum (alpha_k/beta_k) S_k;  A log W = (mu - eps lam)/W + lam[log(1+n/W) - n/W]; log(1+x)>=x-x^2/2.
Checks: K=1 exact m (ODE) and K=2,3 stationary MC (fixed simulator) against the bound and the conjecture 2 eps/(2mu-1/nu1)."""
import numpy as np
from p10_fisher_k1_exact import exact_m
from p2_identifiability import mc

def bound(mu, beta, eps):
    n = 1-eps; B, b = beta.max(), beta.min(); den = 2*mu - n*n*B
    return np.nan if den <= 0 else 2*eps*B/(b*den)

print("K=1 (exact m): condition 2mu > n^2 beta  (k>1);  bound is tight in eps")
for mu in (1.0, 0.6):
    for eps in (0.3, 0.1, 0.01):
        beta = np.array([1.0]); m = exact_m(mu, 1.0, eps)
        print(f"  mu={mu} eps={eps}: m={m:.5f}  bound={bound(mu,beta,eps):.5f}  ratio m/bound={m/bound(mu,beta,eps):.3f}")
for beta, mu in ((np.array([1., 2.]), 1.5), (np.array([1., 2., 4.]), 2.5)):
    K = len(beta); print(f"\nK={K}, beta={beta.tolist()}, mu={mu}, alpha_k=(1-eps)beta_k/K, MC m=I_mumu")
    for eps in (0.1, 0.03):
        alpha = (1-eps)*beta/K; nu1 = (alpha/beta**2).sum(); conj = 2*eps/(2*mu-1/nu1)
        I, _ = mc(mu, beta, eps, 200000, seed=3)
        print(f"  eps={eps}: m(MC)={I[0,0]:.5f}  bound={bound(mu,beta,eps):.5f}  conj 2eps/(2mu-1/nu1)={conj:.5f}  m/eps={I[0,0]/eps:.3f}", flush=True)
