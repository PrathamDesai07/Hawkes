"""Is the 1/eps burn-in cost an artefact of the prior weight rho0 or of the start point?  eps=0.02, stationary S."""
import math, sys
from p5_online_newton_k1 import run
mu_s, beta, eps, n = 1.0, 1.0, 0.02, 100000; reps = int(sys.argv[1]); warm = int(20/eps**2)
print(f"eps={eps}, n={n}, stationary start state; entries are E R_n - log n  (c(eps) ~ 0.5/eps = {0.5/eps:.0f} before)")
print(f"{'start (mu0, alpha0)':>24} | " + " ".join(f"rho0={r:>5}" for r in (0.2, 1, 5, 50)))
for th0 in ((2.0, 0.5), (1.0, 0.5), (2.0, 1-eps), (1.0, 0.9), (1.2, 1-eps)):
    row = []
    for rho0 in (0.2, 1, 5, 50):
        res = [run(mu_s, beta, eps, n, {n}, 7000+r, th0, rho0=rho0, warm=warm) for r in range(reps)]
        row.append(sum(x[n][0] for x in res)/reps - math.log(n))
    print(f"{str(th0):>24} | " + " ".join(f"{v:10.2f}" for v in row), flush=True)
