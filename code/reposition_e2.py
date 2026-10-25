"""E2: K=1 online Fisher scoring (exact F-norm projection onto Theta_c = [.5mu*,2mu*]x[.5al*,beta(1-eps/2)], start at truth, rho0=1).
 Compare E[R_T - log T] at fixed T*gamma = 150 with  (1/2) log det I(theta*; eps)  (exact, stationary).
 If the log(1/eps) term of the regret constant is (1/2) log det I, then Q(eps) = E[R - log T] - (1/2) log det I is flat in eps."""
import math, sys
from p5_projected_exact import run
from p10_fisher_k1 import exact_m

def half_logdet(mu, eps, beta=1.0):
    al = (1 - eps)*beta; m = exact_m(mu, al, beta); A = m*mu/eps; return 0.5*math.log((A - 1)/al**2)

def block(mu, eps, reps, Tg=150.0):
    beta = 1.0; al_s = (1 - eps)*beta; lam_bar = mu/eps; n = int(lam_bar*Tg/(beta*eps))
    lo = (0.5*mu, 0.5*al_s); hi = (2*mu, beta*(1 - eps/2))
    out = [run(2000 + r, mu, beta, eps, n, (mu, al_s), 1.0, lo, hi, warm=int(3/eps**2)) for r in range(reps)]
    v = [o[0] - math.log(o[3]) for o in out]; mean = sum(v)/reps; se = (sum((x - mean)**2 for x in v)/(reps - 1))**.5/reps**.5
    return mean, se, n

if __name__ == "__main__":
    reps = int(sys.argv[1]); cases = eval(sys.argv[2])     # list of (mu, [eps...])
    for mu, epss in cases:
        print(f"mu={mu}  (k0 = 2 mu/beta = {2*mu:.1f}), T*gamma = 150, {reps} reps")
        print(f"   {'eps':>6} {'n events':>9} | {'E[R - log T]':>14} | {'1/2 log det I':>13} | {'Q = diff':>9}")
        for eps in epss:
            mean, se, n = block(mu, eps, reps); h = half_logdet(mu, eps)
            print(f"   {eps:6.3f} {n:9d} | {mean:7.2f} +- {se:4.2f} | {h:13.2f} | {mean - h:9.2f}", flush=True)
