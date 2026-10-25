"""Numerical check of the Problem-5 skeleton (K=1, no projection active).
 (A) expectation identity  E <g_i, th_i-th*> = E (th_i-th*)' dF_i (th_i-th*),  dF_i = Fisher increment at th_i over the gap
 (B) pathwise Newton identity  2*sum<g_i,v_i> = P0 + sum v_i'dF_i v_i + sum g_i'F_i^{-1}g_i - |v_{n+1}|^2_{F_n}   (v_i=th_i-th*; exact if no clipping)
 (C) bound   R_n^{Bregman} <= E sum<g,v> = E sum v'dF v <= P0 + E sum g'F^{-1}g."""
import math, random, sys
from p5_online_newton_k1 import fisher_inc, xq, wq

def one(seed, mu_s, beta, eps, n, th0, rho0, mu_min=1e-6, al_min=1e-6):
    rng = random.Random(seed); al_s = (1-eps)*beta; mu, al = th0; lam_bar = mu_s/eps
    F = [[rho0/lam_bar, 0.0], [0.0, rho0*(lam_bar/beta)**2/lam_bar]]
    P0 = (mu-mu_s)**2*F[0][0] + (al-al_s)**2*F[1][1]
    s = 0.0; R = Sg = Xp = G = 0.0; clipped = False
    for i in range(n):
        lam0 = mu_s + al_s*s; w = 0.0
        while True:
            w += rng.expovariate(lam0)
            if rng.random()*lam0 <= mu_s + al_s*s*math.exp(-beta*w): break
        e = math.exp(-beta*w); a1 = s*e; b1 = s*(1-e)/beta
        Lt = mu*w + al*b1; Ls = mu_s*w + al_s*b1; lg = 0.0
        for x, wt in zip(xq, wq):
            u = 0.5*w*(x+1); su = s*math.exp(-beta*u); ls = mu_s + al_s*su; lt = mu + al*su
            lg += 0.5*w*wt*ls*math.log(lt/ls)
        R += Lt - Ls - lg
        f0, f1, f2 = fisher_inc(mu, al, beta, s, w)
        v0, v1 = mu-mu_s, al-al_s
        g0 = w - 1.0/(mu + al*a1); g1 = b1 - a1/(mu + al*a1)
        Sg += g0*v0 + g1*v1; Xp += f0*v0*v0 + 2*f1*v0*v1 + f2*v1*v1
        F[0][0] += f0; F[0][1] += f1; F[1][0] += f1; F[1][1] += f2
        det = F[0][0]*F[1][1] - F[0][1]**2
        d0 = (F[1][1]*g0 - F[0][1]*g1)/det; d1 = (-F[0][1]*g0 + F[0][0]*g1)/det
        G += g0*d0 + g1*d1
        mu_n, al_n = mu - d0, al - d1
        if mu_n < mu_min or al_n < al_min or al_n > beta*(1-eps/2): clipped = True; break
        mu, al = mu_n, al_n; s = a1 + 1.0
    vT0, vT1 = mu-mu_s, al-al_s
    term = F[0][0]*vT0**2 + 2*F[0][1]*vT0*vT1 + F[1][1]*vT1**2
    resid = 2*Sg - (P0 + Xp + G - term)
    return clipped, R, Sg, Xp, G, P0, term, resid

if __name__ == "__main__":
    reps = int(sys.argv[1]); n = int(sys.argv[2]); mu_s, beta = 1.0, 1.0
    for eps, th0, rho0 in ((0.3, (2.0, 0.5), 1.0), (0.1, (2.0, 0.5), 1.0), (0.1, (1.0, 0.9), 1.0), (0.03, (1.0, 0.9), 1.0)):
        out = [one(100+r, mu_s, beta, eps, n, th0, rho0) for r in range(reps)]
        ok = [o for o in out if not o[0]]; m = len(ok)
        mean = lambda k: sum(o[k] for o in ok)/m
        se = lambda k: (sum((o[k]-mean(k))**2 for o in ok)/(m-1))**.5/m**.5
        print(f"eps={eps} th0={th0} n={n}: unclipped {m}/{reps}")
        print(f"   E R_breg {mean(1):8.2f}±{se(1):.2f} | E sum<g,v> {mean(2):8.2f}±{se(2):.2f} | E sum v'dF v {mean(3):8.2f}±{se(3):.2f} | E sum g'F^-1 g {mean(4):8.2f}±{se(4):.2f} | P0 {mean(5):.2f}")
        print(f"   bound P0+G = {mean(5)+mean(4):.2f}   max |pathwise identity residual| = {max(abs(o[7]) for o in ok):.2e}   frac runs with R <= P0+G: {sum(o[1] <= o[5]+o[4] for o in ok)/m:.2f}", flush=True)
