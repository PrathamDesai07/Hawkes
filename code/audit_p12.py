import numpy as np
from p2_identifiability import mc
for beta, mu, eps, n_ev in ((np.array([1., 2.]), 1.5, 0.05, 300000), (np.array([1., 2., 4.]), 1.5, 0.05, 300000)):
    K = len(beta); I, _ = mc(mu, beta, eps, n_ev, seed=21); ev, V = np.linalg.eigh(I)
    top = V[:, -1]; invb = np.concatenate([[0.0], 1/beta]); invb /= np.linalg.norm(invb)
    alp = (1-eps)*beta/K; ah = np.concatenate([[0.0], alp/np.linalg.norm(alp)])
    print(f"K={K} beta={beta.tolist()} eps={eps}: eigenvalues {np.round(ev,4).tolist()}; |cos(top, (0,1/beta))| = {abs(top@invb):.4f};  |cos(top,(0,alpha_hat))| = {abs(top@ah):.4f}; top component on mu = {abs(top[0]):.4f}", flush=True)
