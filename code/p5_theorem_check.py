"""Evaluate the conditional bound  E R_n <= (P0 + A*D)/(1-B),  A=(1+eta) c+ r(x), B=(1+1/eta) ell/mu_min, r(x)=x/log(1+x), D=E log det(F_n/F_0),
with c+, x, ell taken as the worst values seen over all reps and all indices (empirical proxies for the a.s. assumptions).  Start at theta*, so P0=0."""
import math, sys
from p5_leverage import run
reps = int(sys.argv[1]); n = int(sys.argv[2]); mu_s, beta = 1.0, 1.0; mu_min = 0.5
print(f"mu_min={mu_min}, start at truth (P0=0), n={n}, reps={reps}")
print(f"{'eps':>5} {'rho0':>6} | {'ell_max':>8} {'x_max':>7} {'c+':>6} | {'B(eta*)':>8} {'eta*':>5} {'A':>6} | {'D':>6} {'bound':>8} | {'E R_n':>7} {'(d/2)logT':>9}")
for eps, rho0, warm in ((0.3, 5, 300), (0.1, 20, 3000), (0.1, 60, 3000), (0.03, 200, 20000), (0.03, 600, 20000)):
    out = [run(900+r, mu_s, beta, eps, n, (mu_s, (1-eps)*beta), float(rho0), warm, mu_min=mu_min) for r in range(reps)]
    ell = max(v[0] for o in out for v in o[0]); xm = max(v[1] for o in out for v in o[0]); cp = max(v[2] for o in out for v in o[0])
    D = sum(o[1] for o in out)/reps; R = sum(o[2] for o in out)/reps; best = None
    for eta in (0.25, 0.5, 0.75, 1, 1.5, 2, 3, 5, 10, 30):
        B = (1+1/eta)*ell/mu_min
        if B < 1:
            A = (1+eta)*cp*xm/math.log(1+xm); val = A*D/(1-B)
            if best is None or val < best[0]: best = (val, eta, B, A)
    T = n*eps/mu_s
    if best: print(f"{eps:5.2f} {rho0:6d} | {ell:8.3f} {xm:7.3f} {cp:6.2f} | {best[2]:8.2f} {best[1]:5} {best[3]:6.2f} | {D:6.2f} {best[0]:8.1f} | {R:7.2f} {math.log(T):9.2f}", flush=True)
    else: print(f"{eps:5.2f} {rho0:6d} | {ell:8.3f} {xm:7.3f} {cp:6.2f} |  no eta gives B<1 (B_min={(1+1/30)*ell/mu_min:.2f}) | D={D:.2f}  E R_n={R:.2f}", flush=True)
