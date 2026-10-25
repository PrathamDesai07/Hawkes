"""How loose is Lemma 3?  lam_min(I) >= (eps/mu) * inf_{|v|=1} (E|v^T a|)^2.  Compare with MC lam_min(I)."""
import numpy as np
from scipy.optimize import minimize
from p2_window_bound import run
mu = 1.0
print(f"{'K':>2} {'r':>3} {'eps':>5} {'lmin(I)':>9} {'inf (E|va|)^2':>14} {'Lemma3 bound':>13} {'bound/lmin':>10}")
for K, r in ((1, 2.0), (2, 2.0), (2, 4.0), (3, 2.0)):
    beta = r**np.arange(K)
    for eps in (0.3, 0.1, 0.05):
        A, _, Ss = run(mu, beta, eps, 60000, [0.1*eps], return_samples=True)
        a = np.hstack([np.ones((len(Ss), 1)), Ss])
        f = lambda v: np.mean(np.abs(a@(v/np.linalg.norm(v))))**2
        best = min(minimize(f, x0, method="Nelder-Mead", options={"xatol": 1e-9, "fatol": 1e-12, "maxiter": 4000}).fun
                   for x0 in np.random.default_rng(0).normal(size=(12, K+1)))
        lm = np.linalg.eigvalsh(A)[0]
        print(f"{K:2d} {r:3.0f} {eps:5.2f} {lm:9.4f} {best:14.5f} {eps/mu*best:13.5f} {eps/mu*best/lm:10.3f}", flush=True)
