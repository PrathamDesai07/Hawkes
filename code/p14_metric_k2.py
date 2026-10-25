"""Problem 14: Fisher-metric equivalence constant kappa_I for K=2 (beta=(1,2), alpha_k=(1-eps)beta_k/2) over a drift box, I from stationary MC.
kappa_I = max over pairs of  lam_max/lam_min of  I(theta')^{-1/2} I(theta) I(theta')^{-1/2}.  MC noise ~5 % in I, so kappa_I is good to ~10-20 %."""
import numpy as np, sys, itertools
from scipy.linalg import eigh
from p2_identifiability import mc
eps0 = float(sys.argv[1]); n_ev = int(sys.argv[2]); beta = np.array([1., 2.]); mu0 = 1.5
pts = {(cm, ce): (mu0*cm, eps0*ce) for cm in (1.0, 1.5) for ce in (1.0, 1.5)}
Is = {}
for key, (mu, eps) in pts.items():
    Is[key], _ = mc(mu, beta, eps, n_ev, seed=11 + int(10*key[0]) + int(10*key[1])); print(f" computed I at mu={mu:.3f}, eps={eps:.3f}", flush=True)
def kap(keys):
    lo, hi = np.inf, 0.0
    for a, b in itertools.product(keys, keys):
        w = eigh(Is[a], Is[b], eigvals_only=True); lo = min(lo, w[0]); hi = max(hi, w[-1])
    return hi/lo, lo, hi
print(f"K=2, eps0={eps0}, mu0={mu0}: cond(I(theta_0)) = {np.linalg.cond(Is[(1.0,1.0)]):.1f}")
for name, keys in (("mu x1.5 (eps fixed)", [(1.0, 1.0), (1.5, 1.0)]), ("eps x1.5 (mu fixed)", [(1.0, 1.0), (1.0, 1.5)]), ("both x1.5", list(Is))):
    k, lo, hi = kap(keys); print(f"  {name}: kappa_I = {k:.2f}   generalized eigenvalue range {lo:.3f}-{hi:.3f}", flush=True)
