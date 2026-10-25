"""K>=2, k_K<1: efficient information for mu (alpha unknown) eff = 1/(I^-1)_mumu versus alpha-known m = I_mumu (stationary MC)."""
import numpy as np
from p2_identifiability import mc
for beta, mu in ((np.array([1., 2.]), 0.4), (np.array([1., 3., 9.]), 0.5)):
    K = len(beta)
    for eps in (0.1, 0.03, 0.01):
        alpha = (1-eps)*beta/K
        I, _ = mc(mu, beta, eps, 200000 if eps > 0.01 else 500000, seed=6)
        eff = 1/np.linalg.inv(I)[0, 0]; m = I[0, 0]
        print(f"beta={beta.tolist()} mu={mu} eps={eps}: m={m:.4f}  eff={eff:.4f}  eff/m={eff/m:.3f}   m/eps={m/eps:.1f} eff/eps={eff/eps:.1f}", flush=True)
