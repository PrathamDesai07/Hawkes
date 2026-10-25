"""Check F_{i-1} - F_0 - J_T/c+ is PSD along the projected run (c+ = max(mu_hi/mu*, al_hi/al*)), exact F-norm projection, K=1."""
import math, random, sys
from p5_online_newton_k1 import fisher_inc
from p5_projected_exact import project_F

def one(seed, mu_s, beta, eps, n, th0, rho0, lo, hi, cplus, warm):
    rng = random.Random(seed); al_s = (1-eps)*beta; mu, al = th0; lam_bar = mu_s/eps
    F0 = [[rho0/lam_bar, 0.0], [0.0, rho0*(lam_bar/beta)**2/lam_bar]]; F = [r[:] for r in F0]; J = [0.0, 0.0, 0.0]; s = 0.0
    for _ in range(warm):
        lam0 = mu_s + al_s*s; w = 0.0
        while True:
            w += rng.expovariate(lam0)
            if rng.random()*lam0 <= mu_s + al_s*s*math.exp(-beta*w): break
        s = s*math.exp(-beta*w) + 1.0
    worst = float('inf')
    for i in range(n):
        lam0 = mu_s + al_s*s; w = 0.0
        while True:
            w += rng.expovariate(lam0)
            if rng.random()*lam0 <= mu_s + al_s*s*math.exp(-beta*w): break
        e = math.exp(-beta*w); a1 = s*e; b1 = s*(1-e)/beta
        f0, f1, f2 = fisher_inc(mu, al, beta, s, w); j0, j1, j2 = fisher_inc(mu_s, al_s, beta, s, w)
        F[0][0] += f0; F[0][1] += f1; F[1][0] += f1; F[1][1] += f2; J[0] += j0; J[1] += j1; J[2] += j2
        # PSD margin of D = F - F0 - J/c+ :  smallest eigenvalue, normalised by the largest eigenvalue of F
        a = F[0][0] - F0[0][0] - J[0]/cplus; b = F[0][1] - J[1]/cplus; c = F[1][1] - F0[1][1] - J[2]/cplus
        lm = 0.5*(a + c) - math.sqrt(0.25*(a - c)**2 + b*b); lmax = 0.5*(F[0][0] + F[1][1]) + math.sqrt(0.25*(F[0][0] - F[1][1])**2 + F[0][1]**2)
        worst = min(worst, lm/lmax)
        g0 = w - 1.0/(mu + al*a1); g1 = b1 - a1/(mu + al*a1)
        det = F[0][0]*F[1][1] - F[0][1]**2
        d0 = (F[1][1]*g0 - F[0][1]*g1)/det; d1 = (-F[0][1]*g0 + F[0][0]*g1)/det
        (mu, al), _ = project_F(F, (mu - d0, al - d1), lo, hi); s = a1 + 1.0
    return worst

if __name__ == "__main__":
    reps = int(sys.argv[1]); n = int(sys.argv[2]); mu_s, beta = 1.0, 1.0
    for eps in (0.3, 0.1, 0.03):
        al_s = (1-eps)*beta; lo = (0.5*mu_s, 0.5*al_s); hi = (2*mu_s, beta*(1-eps/2)); cp = max(hi[0]/mu_s, hi[1]/al_s)
        w = [one(60+r, mu_s, beta, eps, n, (1.4, 0.8*al_s), 1.0, lo, hi, cp, int(3/eps**2)) for r in range(reps)]
        print(f"eps={eps}: c+={cp:.3f}; min over runs and steps of lambda_min(F-F0-J/c+)/lambda_max(F) = {min(w):.3e}  (>= 0 means PSD holds)")
