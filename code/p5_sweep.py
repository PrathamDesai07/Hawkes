"""Plateau constant c(eps) = E R_n - log n at large n, for four start conditions."""
import math, sys
from p5_online_newton_k1 import run
mu_s, beta = 1.0, 1.0; reps = int(sys.argv[1]); 
print(f"{'eps':>6} {'n':>7} | {'far th, cold S':>15} {'far th, stat S':>15} {'true th, cold S':>16} {'true th, stat S':>16} | slope(1e4->n) far/stat")
for eps, n in ((0.3, 40000), (0.1, 40000), (0.05, 60000), (0.03, 100000), (0.02, 150000), (0.01, 400000)):
    cps = {n//4, n}
    cfg = [((2*mu_s, .5*beta), 0), ((2*mu_s, .5*beta), int(20/eps**2)), ((mu_s, (1-eps)*beta), 0), ((mu_s, (1-eps)*beta), int(20/eps**2))]
    vals = []; slopes = None
    for th0, warm in cfg:
        res = [run(mu_s, beta, eps, n, cps, 5000+r, th0, warm=warm) for r in range(reps)]
        Rn = sum(x[n][0] for x in res)/reps; vals.append(Rn - math.log(n))
        if th0[0] == 2*mu_s and warm: slopes = (Rn - sum(x[n//4][0] for x in res)/reps)/math.log(4)
    print(f"{eps:6.2f} {n:7d} | {vals[0]:15.2f} {vals[1]:15.2f} {vals[2]:16.2f} {vals[3]:16.2f} | {slopes:.2f}", flush=True)
