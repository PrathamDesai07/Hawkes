"""Robustness: interior starts and a wider box; exact F-norm projection vs diagonal clip."""
import math, sys
from p5_projected_exact import run
reps = int(sys.argv[1]); n = int(sys.argv[2]); mu_s, beta = 1.0, 1.0
print(f"{'eps':>5} {'box':>14} {'start':>14} {'mode':>6} | {'E R_n':>8} {'R-log T':>8} | {'active':>7} {'nonexp viol':>11}")
for eps in (0.1, 0.03):
    al_s = (1-eps)*beta
    for box_name, lo, hi in (("[.5,2]x[.5,hi]", (0.5, 0.5*al_s), (2.0, beta*(1-eps/2))), ("[.2,5]x[.1,hi]", (0.2, 0.1*al_s), (5.0, beta*(1-eps/2)))):
        for st_name, th0 in (("(1.4,0.8al*)", (1.4, 0.8*al_s)), ("(1.0,al*)=truth", (1.0, al_s)), ("(1.8,0.6al*)", (1.8, 0.6*al_s))):
            for mode in ("exact", "diag"):
                out = [run(900+r, mu_s, beta, eps, n, th0, 1.0, lo, hi, warm=int(3/eps**2), mode=mode) for r in range(reps)]
                m = lambda k: sum(o[k] for o in out)/reps
                print(f"{eps:5.2f} {box_name:>14} {st_name:>14} {mode:>6} | {m(0):8.2f} {m(0) - sum(math.log(o[3]) for o in out)/reps:8.2f} | {m(4):7.4f} {m(5):11.5f}", flush=True)
