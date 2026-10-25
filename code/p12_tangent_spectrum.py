"""Rigorous lower bounds for EVERY eigenvalue of I from the tangent minorant I >= N(lam0)  (Weyl monotonicity: lam_j(I) >= lam_j(N(lam0)) for each j and each lam0 > 0).
 Compared with the Monte-Carlo eigenvalues of Problem 12 (mu = 1.5, alpha_k = (1-eps) beta_k / K)."""
import numpy as np, mpmath as mp
from p2_tangent_bound import matrices
mp.mp.dps = 60

def eigs_N(M, M3, lam0):
    d = M.shape[0]; N = mp.matrix(d, d)
    for i in range(d):
        for j in range(d):
            N[i, j] = 2*mp.mpf(M[i, j])/lam0 - mp.mpf(M3[i, j])/lam0**2
    return sorted([float(x) for x in mp.eigsy(N, eigvals_only=True)])

cases = [("K=2 beta=(1,2) eps=0.1", [1., 2.], 0.1, [0.032, 0.065, 18.2]),
         ("K=2 beta=(1,2) eps=0.03", [1., 2.], 0.03, [0.015, 0.064, 57.2]),
         ("K=2 beta=(1,2) eps=0.01", [1., 2.], 0.01, [0.0053, 0.066, 171.0]),
         ("K=3 beta=(1,2,4) eps=0.1", [1., 2., 4.], 0.1, [0.0047, 0.0446, 0.1495, 19.1]),
         ("K=3 beta=(1,2,4) eps=0.03", [1., 2., 4.], 0.03, [0.0055, 0.0225, 0.1652, 60.9])]
mu = 1.5
for name, beta, eps, mc in cases:
    beta = np.array(beta); K = len(beta); alpha = (1-eps)*beta/K
    M, M3, El = matrices(mu, alpha, beta)
    best = np.full(K+1, -np.inf)
    for f in np.concatenate([np.linspace(0.3, 3.0, 28), np.linspace(3.5, 12.0, 18)]):
        ev = np.array(eigs_N(M, M3, mp.mpf(f*El)))
        best = np.maximum(best, ev)
    print(f"{name}\n   MC eigenvalues of I        : {mc}\n   rigorous lower bounds lam_j(N): {[float(f'{x:.4g}') for x in best]}\n   ratios bound/MC            : {[round(b/m, 2) for b, m in zip(best, mc)]}", flush=True)
