"""Empirical cost of the a-priori box (K=1, exact F-projection).  Box [mu*/s, s mu*] x [al*/s, beta(1-eps/2)], widened by s; start at the far corner (s mu*, al*/s) (worst case) and at the truth.
 Reports E R_n and E R_n - ln T_n at n events, and the last time the iterate is outside R(0.3)."""
import math, sys
from p5_projected_exact import run
reps = int(sys.argv[1]); n = int(sys.argv[2]); mu_s, beta = 1.0, 1.0
for eps in (0.3, 0.1):
    al_s = (1-eps)*beta
    for s in (2, 5, 20, 100):
        lo = (mu_s/s, al_s/s); hi = (s*mu_s, beta*(1-eps/2))
        for name, th0 in (("truth", (mu_s, al_s)), ("corner", (s*mu_s, al_s/s))):
            out = [run(900+r, mu_s, beta, eps, n, th0, 1.0, lo, hi, warm=int(3/eps**2)) for r in range(reps)]
            Rm = sum(o[0] for o in out)/reps; lT = sum(math.log(o[3]) for o in out)/reps
            print(f"eps={eps} s={s:>3} start={name:>6}: E R_n = {Rm:7.2f}, E R_n - ln T_n = {Rm-lT:7.2f}, c+ = {max(s, hi[1]/al_s):.2f}, c- = {1/s:.3f}", flush=True)
