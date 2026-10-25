"""Box cost vs prior weight rho0 (K=1, eps=0.1, corner start): extra regret relative to the truth start."""
import math, sys
from p5_projected_exact import run
reps = int(sys.argv[1]); n = int(sys.argv[2]); mu_s, beta, eps = 1.0, 1.0, 0.1; al_s = (1-eps)*beta
for s in (5, 20, 100):
    lo = (mu_s/s, al_s/s); hi = (s*mu_s, beta*(1-eps/2))
    for rho0 in (1.0, 0.1, 0.01):
        res = {}
        for name, th0 in (("truth", (mu_s, al_s)), ("corner", (s*mu_s, al_s/s))):
            out = [run(1300+r, mu_s, beta, eps, n, th0, rho0, lo, hi, warm=300) for r in range(reps)]
            res[name] = sum(o[0] for o in out)/reps
        print(f"s={s:>3} rho0={rho0:5}: E R_n truth-start {res['truth']:7.2f}, corner-start {res['corner']:8.2f}, extra {res['corner']-res['truth']:8.2f}", flush=True)
