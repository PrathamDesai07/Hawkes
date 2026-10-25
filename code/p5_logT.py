"""Test R_n ~ (d/2) log T_n + C with C eps-free.  Start at truth (theta0=theta*), weak prior rho0, moderately warm start."""
import math, sys
from p5_online_newton_k1 import run
mu_s, beta = 1.0, 1.0; reps = int(sys.argv[1])
print("E[R_n - log n]  and  E[R_n - log T_n]   (d/2 = 1);  rho0 = 1; n = checkpoints")
print(f"{'eps':>6} {'n':>7} {'R-log n':>9} {'R-log T':>9} {'se(R-log T)':>11} {'E e_n':>7}")
for eps, ns, warm in ((0.3, (1000, 10000, 40000), 200), (0.1, (1000, 10000, 40000), 2000), (0.03, (1000, 10000, 60000), 10000), (0.01, (10000, 60000), 30000)):
    n = ns[-1]; res = [run(mu_s, beta, eps, n, set(ns), 9000+r, (mu_s, (1-eps)*beta), rho0=1.0, warm=warm) for r in range(reps)]
    for c in ns:
        a = [x[c][0] - math.log(c) for x in res]; b = [x[c][0] - math.log(x[c][2]) for x in res]; e = [x[c][1] for x in res]
        mb = sum(b)/reps; se = (sum((v-mb)**2 for v in b)/(reps-1))**.5/reps**.5
        print(f"{eps:6.2f} {c:7d} {sum(a)/reps:9.2f} {mb:9.2f} {se:11.2f} {sum(e)/reps:7.2f}", flush=True)
