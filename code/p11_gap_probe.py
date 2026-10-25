"""K=2 gap-region probe for Problem 11: m/eps for beta=(1,2), mu=0.8 (1/nu1 < 2mu < beta_max); conjecture m/eps -> 2/(2mu-1/nu1)."""
import numpy as np
from p2_identifiability import mc
beta = np.array([1., 2.]); K = 2; mu = 0.8
for eps in (0.1, 0.03, 0.01):
    alpha = (1-eps)*beta/K; nu1 = (alpha/beta**2).sum()
    I, _ = mc(mu, beta, eps, 300000 if eps > 0.01 else 600000, seed=8)
    print(f"eps={eps}: m/eps={I[0,0]/eps:.2f}  conjectured limit 2/(2mu-1/nu1)={2/(2*mu-1/nu1):.2f}", flush=True)
