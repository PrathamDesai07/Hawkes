"""Problem 6 (D): the leverage-scale constant  C_hat = 1/(lam_bar * lam_min(S N S)),  S = diag(1, beta_k/lam_bar),  so that
   phi' N^-1 phi <= C_hat * lam_bar * (1 + sum_k (beta_k S_k / lam_bar)^2)   for every state S >= 0 (phi = (1,S)),
 with N = N(lam0) the rigorous tangent minorant of I (so I >= N, F_t >= t N/2 on the information event).  Also C_hat for I itself (Monte-Carlo I from p6_out3)."""
import sys, types, numpy as np
sys.modules.setdefault('mpmath', types.ModuleType('mpmath')); sys.modules['mpmath'].mp = types.SimpleNamespace(dps=60)
from p2_tangent_bound import matrices
def chat(mu, beta, eps, f):
    beta = np.asarray(beta, float); K = len(beta); alpha = (1-eps)*beta/K
    M, M3, El = matrices(mu, alpha, beta); l0 = f*El; N = 2*M/l0 - M3/l0**2; N = (N+N.T)/2
    lb = El; S = np.diag(np.concatenate([[1.0], beta/lb])); w = np.linalg.eigvalsh(S@N@S)
    return 1.0/(lb*w.min()) if w.min() > 0 else np.inf
print("C_hat(N(lam0 = f*lam_bar)) ;  columns f = 1.0 1.5 2.0 3.0 ; best (smallest) in last column")
for name, beta, mu in (("K=1 beta=1 mu=1", [1.0], 1.0), ("K=1 beta=1 mu=2", [1.0], 2.0), ("K=2 beta=1,2 mu=1.5", [1.0, 2.0], 1.5), ("K=3 beta=1,2,4 mu=2", [1.0, 2.0, 4.0], 2.0)):
    for eps in (0.3, 0.1, 0.03, 0.01):
        v = [chat(mu, beta, eps, f) for f in (1.0, 1.5, 2.0, 3.0)]
        print(f"{name:22s} eps={eps:<5} " + " ".join(f"{x:9.3f}" for x in v) + f"   best={min(v):.3f}   (k0={2*mu/beta[0]:.1f})")
