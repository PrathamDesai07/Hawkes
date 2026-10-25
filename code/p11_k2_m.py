"""K=2 stationary m=E[1/lam] (= I_mumu): conjecture m ~ 2 eps/(2 mu - 1/nu1), nu1 = sum alpha_k/beta_k^2 (n->1).
Also check the exact K-general identity  u^T I u = (1+eps)^2 m - (eps/mu)(1+2eps)  for u=(1,-eps*alpha/mu)."""
import numpy as np
from p2_identifiability import mc
mu = 1.0; beta = np.array([1., 2.])
print("K=2, beta=(1,2), alpha_k=(1-eps)beta_k/2, mu=1; nu1=0.75 -> predicted m/eps -> 2/(2-1/0.75)=3.0")
for eps in (0.1, 0.03):
    alpha = (1-eps)*beta/2
    I, _ = mc(mu, beta, eps, 250000, seed=5)
    m = I[0, 0]; u = np.concatenate([[1.0], -eps*alpha/mu]); qf = u@I@u
    print(f" eps={eps}: m={m:.5f}  m/eps={m/eps:.3f}   u^T I u (MC)={qf:.5f}   identity (1+e)^2 m-(e/mu)(1+2e)={(1+eps)**2*m-eps/mu*(1+2*eps):.5f}", flush=True)
