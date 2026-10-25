"""Exponential-tail ingredients (K=1).
 N1 exact identities:  (a) eps*Lam_T = mu*T + n(N_T - Lam_T) - n(S_T - S_0),  n=1-eps
                       (b) int_0^T (S_t - ES) dt = (S_0 - S_T)/gamma + (N_T - Lam_T)/gamma     (Poisson equation u(s)=(s-ES)/gamma)
 N2 lower tail of the integrated intensity: P(Lam_T <= (1-eta) mu T/eps) vs exp(-eta^2 mu eps T/(2 n^2))-type bound (after removing the boundary term)
 N3 empirical tail of an occupation functional P(int h <= T m/2): exponential in T*gamma?
 N4 Ito identity in expectation:  E sum s'F_{i-1}^{-1}s  =  E sum a'F_{i-1}^{-1}a/lam*(tau)^2"""
import math, random, sys
from p5_online_newton_k1 import fisher_inc
from p5_fisher_conc import stationary_law_moments
from p5_variance_check import occupation, trap
from p5_projected_exact import project_F

def gap(rng, mu, al, beta, s):
    lam0 = mu + al*s; w = 0.0
    while True:
        w += rng.expovariate(lam0)
        if rng.random()*lam0 <= mu + al*s*math.exp(-beta*w): return w

def n1_n2(reps, mu, beta, eps, T_target):
    al = (1-eps)*beta; n_ = 1 - eps; gam = beta*eps; ES = mu/gam; res = []; devs = []
    for r in range(reps):
        rng = random.Random(500+r); s = 0.0
        for _ in range(int(10/eps**2)): s = s*math.exp(-beta*gap(rng, mu, al, beta, s)) + 1.0   # stationary start
        s0 = s; t = 0.0; Lam = 0.0; N = 0; intS = 0.0
        while t < T_target:
            w = gap(rng, mu, al, beta, s); e = math.exp(-beta*w)
            Lam += mu*w + al*s*(1-e)/beta; intS += s*(1-e)/beta; t += w; N += 1; s = s*e + 1.0
        a = eps*Lam - (mu*t + n_*(N - Lam) - n_*(s - s0))
        b = (intS - ES*t) - ((s0 - s)/gam + (N - Lam)/gam)
        res.append((abs(a)/(mu*t), abs(b)/(ES*t), Lam/(mu*t/eps), t, s, s0))
    return res

def n4(reps, mu, beta, eps, n, lo, hi):
    al_s = (1-eps)*beta; tot1 = tot2 = 0.0
    lam_bar = mu/eps
    for r in range(reps):
        rng = random.Random(900+r); th = (1.4, 0.8*al_s); F = [[1.0/lam_bar, 0.0], [0.0, (lam_bar/beta)**2/lam_bar]]; s = 0.0
        for _ in range(int(3/eps**2)): s = s*math.exp(-beta*gap(rng, mu, al_s, beta, s)) + 1.0
        for i in range(n):
            w = gap(rng, mu, al_s, beta, s); e = math.exp(-beta*w); a1 = s*e; b1 = s*(1-e)/beta
            lam_st = mu + al_s*a1; s0_ = w - 1.0/lam_st; s1_ = b1 - a1/lam_st
            det = F[0][0]*F[1][1] - F[0][1]**2
            q = lambda p, q_: (F[1][1]*p*p - 2*F[0][1]*p*q_ + F[0][0]*q_*q_)/det
            tot1 += q(s0_, s1_); tot2 += q(1.0, a1)/lam_st**2
            f0, f1, f2 = fisher_inc(th[0], th[1], beta, s, w); F[0][0] += f0; F[0][1] += f1; F[1][0] += f1; F[1][1] += f2
            lam_th = th[0] + th[1]*a1; g0 = w - 1.0/lam_th; g1 = b1 - a1/lam_th; det = F[0][0]*F[1][1] - F[0][1]**2
            d0 = (F[1][1]*g0 - F[0][1]*g1)/det; d1 = (-F[0][1]*g0 + F[0][0]*g1)/det
            th, _ = project_F(F, (th[0] - d0, th[1] - d1), lo, hi); s = a1 + 1.0
    return tot1/reps, tot2/reps

if __name__ == "__main__":
    reps = int(sys.argv[1]); mu, beta = 1.0, 1.0
    print("N1: pathwise identities, max relative residual over runs")
    for eps in (0.3, 0.1):
        r = n1_n2(40, mu, beta, eps, 200/eps)
        print(f"   eps={eps}: (a) {max(x[0] for x in r):.2e}   (b) {max(x[1] for x in r):.2e}")
    print("N2: P(Lam_T <= (1-eta) mu T/eps)  vs  exp(-eta^2 mu eps T/(2 n^2)) (boundary term ignored) ; eta=0.5")
    for eps in (0.3, 0.1):
        for muepsT in (5, 10, 20):
            T = muepsT/(mu*eps); r = n1_n2(reps, mu, beta, eps, T); n_ = 1 - eps
            emp = sum(x[2] <= 0.5 for x in r)/reps; bd = math.exp(-0.25*muepsT/(2*n_**2))
            print(f"   eps={eps} mu*eps*T={muepsT:3d}: empirical {emp:.3f}   bound {bd:.3f}   (mean Lam_T/(mu T/eps) = {sum(x[2] for x in r)/reps:.3f})")
    print("N3: P(int h <= T m/2) for bump h on I2 = [c2-d,c2+d] (c2 = 1.8 ES, d = .15 ES, ramp d/2), eps=0.1  -- exponential in T*gamma?")
    eps = 0.1; ES, ES2 = stationary_law_moments(mu, beta, eps); gam = beta*eps; c = 1.8*ES; d = 0.15*ES
    h = lambda s: trap(s, c - d, c + d, d/2)
    vals = {}
    for Tg in (5, 10, 20, 40):
        rng = random.Random(77); v = [occupation(rng, mu, beta, eps, Tg/gam, int(20/eps**2), h) for _ in range(reps)]
        m = 0.047; vals[Tg] = sum(x <= Tg/gam*m/2 for x in v)/reps
        print(f"   T*gamma={Tg:3d}: P = {vals[Tg]:.3f}   -log P/(T gamma) = {(-math.log(vals[Tg])/Tg if vals[Tg] > 0 else float('nan')):.3f}", flush=True)
    print("N4: Ito identity in expectation (projected exact run, eps=0.1, n=5000)")
    eps = 0.1; al_s = (1-eps)*beta; a, b = n4(30, mu, beta, eps, 5000, (0.5*mu, 0.5*al_s), (2*mu, beta*(1-eps/2)))
    print(f"   E sum s'F^-1 s = {a:.3f}   vs   E sum a'F^-1 a/lam*^2 = {b:.3f}")
