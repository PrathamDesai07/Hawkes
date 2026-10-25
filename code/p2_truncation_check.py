"""Rigorous: I = E[a a^T/lam] >= E[a a^T 1{lam<=L}]/L  (Loewner).  Empirical check of how tight
   B(L) = lam_min( E[a a^T 1{lam<=L}] ) / L   is, maximised over L = c*mu/eps, vs lam_min(I) (MC)."""
import numpy as np
from p2_window_bound import run
mu = 1.0
print(f"{'K':>2} {'r':>3} {'eps':>5} {'lmin(I)':>9} {'best B(L)':>10} {'B/lmin':>7} {'best c':>6}  lam_min(M) full 2nd-moment matrix")
for K, r in ((1, 2.0), (2, 2.0), (2, 4.0), (3, 2.0)):
    beta = r**np.arange(K)
    for eps in (0.3, 0.1, 0.05):
        A, _, Ss = run(mu, beta, eps, 80000, [0.1*eps], return_samples=True)
        alpha = (1-eps)*beta/K; a = np.hstack([np.ones((len(Ss), 1)), Ss]); lam = mu + Ss@alpha
        best = (0, 0)
        for c in (1, 1.5, 2, 3, 4, 6, 8, 12, 20):
            sel = lam <= c*mu/eps; M = (a[sel].T@a[sel])/len(Ss)
            b = np.linalg.eigvalsh(M)[0]/(c*mu/eps)
            if b > best[0]: best = (b, c)
        Mfull = a.T@a/len(Ss)
        lm = np.linalg.eigvalsh(A)[0]
        print(f"{K:2d} {r:3.0f} {eps:5.2f} {lm:9.4f} {best[0]:10.4f} {best[0]/lm:7.2f} {best[1]:6}  {np.linalg.eigvalsh(Mfull)[0]:.4f}", flush=True)
