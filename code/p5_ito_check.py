"""N4 redo: paired difference  D = sum_i [ s_i'F_{i-1}^-1 s_i - a_i'F_{i-1}^-1 a_i / lam*(tau_i)^2 ]  per path; mean and standard error.
 Also split off the first 50 events (heavy leverage) and use two prior strengths."""
import math, random, sys
from p5_online_newton_k1 import fisher_inc
from p5_exp_tails_check import gap
from p5_projected_exact import project_F

def path(seed, mu, beta, eps, n, rho0, lo, hi):
    al_s = (1-eps)*beta; lam_bar = mu/eps; rng = random.Random(seed); th = (1.4, 0.8*al_s)
    F = [[rho0/lam_bar, 0.0], [0.0, rho0*(lam_bar/beta)**2/lam_bar]]; s = 0.0
    for _ in range(int(3/eps**2)): s = s*math.exp(-beta*gap(rng, mu, al_s, beta, s)) + 1.0
    d_all = d_late = t1 = t2 = 0.0
    for i in range(n):
        w = gap(rng, mu, al_s, beta, s); e = math.exp(-beta*w); a1 = s*e; b1 = s*(1-e)/beta
        lam_st = mu + al_s*a1; s0_ = w - 1.0/lam_st; s1_ = b1 - a1/lam_st; det = F[0][0]*F[1][1] - F[0][1]**2
        q = lambda p, q_: (F[1][1]*p*p - 2*F[0][1]*p*q_ + F[0][0]*q_*q_)/det
        x = q(s0_, s1_) - q(1.0, a1)/lam_st**2; d_all += x; d_late += x if i >= 50 else 0.0; t1 += q(s0_, s1_); t2 += q(1.0, a1)/lam_st**2
        f0, f1, f2 = fisher_inc(th[0], th[1], beta, s, w); F[0][0] += f0; F[0][1] += f1; F[1][0] += f1; F[1][1] += f2
        lam_th = th[0] + th[1]*a1; g0 = w - 1.0/lam_th; g1 = b1 - a1/lam_th; det = F[0][0]*F[1][1] - F[0][1]**2
        d0 = (F[1][1]*g0 - F[0][1]*g1)/det; d1 = (-F[0][1]*g0 + F[0][0]*g1)/det
        th, _ = project_F(F, (th[0] - d0, th[1] - d1), lo, hi); s = a1 + 1.0
    return d_all, d_late, t1, t2

if __name__ == "__main__":
    reps = int(sys.argv[1]); mu, beta, eps = 1.0, 1.0, 0.1; al_s = (1-eps)*beta
    lo = (0.5*mu, 0.5*al_s); hi = (2*mu, beta*(1-eps/2))
    for rho0 in (1.0, 20.0):
        out = [path(3000+r, mu, beta, eps, 3000, rho0, lo, hi) for r in range(reps)]
        for name, k in (("all events", 0), ("events >= 50", 1)):
            v = [o[k] for o in out]; m = sum(v)/reps; se = (sum((x - m)**2 for x in v)/(reps-1))**.5/reps**.5
            print(f"rho0={rho0}: {name}: mean paired difference {m:+.3f} +- {se:.3f}   (E sum s'F^-1 s = {sum(o[2] for o in out)/reps:.2f}, E sum a'F^-1 a/lam*^2 = {sum(o[3] for o in out)/reps:.2f})", flush=True)
