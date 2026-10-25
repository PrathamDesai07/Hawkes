"""Check the three component steps of the conditional bound with iterates OFF the truth (v != 0).
 (1) information identity:   E sum s'F_{i-1}^-1 s   vs   E sum tr(F_{i-1}^-1 dJ_i)          (s = b - a/lam*, dJ_i = Fisher increment at theta*)
 (2) drift bound:            E sum h'F_{i-1}^-1 h   <=  (ell_max/mu_min) * E sum v'dF_i v    (h = g - s)
 (3) pathwise log-det:       sum tr(F_{i-1}^-1 dF_i)  <=  r(x_max) * log det(F_n/F_0),  r(x)=x/log(1+x)
 and the split  G' = E sum g'F_{i-1}^-1 g  vs  (1+eta) E sum s'..s + (1+1/eta) E sum h'..h."""
import math, random, sys
from p5_online_newton_k1 import fisher_inc, xq, wq

def one(seed, mu_s, beta, eps, n, th0, rho0, warm, mu_min):
    rng = random.Random(seed); al_s = (1-eps)*beta; mu, al = th0; lam_bar = mu_s/eps; al_max = beta*(1-eps/2)
    F = [[rho0/lam_bar, 0.0], [0.0, rho0*(lam_bar/beta)**2/lam_bar]]; s = 0.0
    for _ in range(warm):
        lam0 = mu_s + al_s*s; w = 0.0
        while True:
            w += rng.expovariate(lam0)
            if rng.random()*lam0 <= mu_s + al_s*s*math.exp(-beta*w): break
        s = s*math.exp(-beta*w) + 1.0
    ld0 = math.log(F[0][0]*F[1][1]); Gs = Gh = Gg = Xv = TJ = TF = 0.0; ellm = xm = 0.0; clip = 0
    for i in range(n):
        lam0 = mu_s + al_s*s; w = 0.0
        while True:
            w += rng.expovariate(lam0)
            if rng.random()*lam0 <= mu_s + al_s*s*math.exp(-beta*w): break
        e = math.exp(-beta*w); a1 = s*e; b1 = s*(1-e)/beta
        det = F[0][0]*F[1][1] - F[0][1]**2
        inv = lambda p, q: (F[1][1]*p*p - 2*F[0][1]*p*q + F[0][0]*q*q)/det
        for fr in (0, .25, .5, .75, 1.0, None):
            su = (s*math.exp(-beta*fr*w) if fr is not None else 0.0); ellm = max(ellm, inv(1.0, su)/(mu_s + al_s*su))
        f0, f1, f2 = fisher_inc(mu, al, beta, s, w); j0, j1, j2 = fisher_inc(mu_s, al_s, beta, s, w)
        m00 = (F[1][1]*f0 - F[0][1]*f1)/det; m01 = (F[1][1]*f1 - F[0][1]*f2)/det
        m10 = (-F[0][1]*f0 + F[0][0]*f1)/det; m11 = (-F[0][1]*f1 + F[0][0]*f2)/det
        tr = m00 + m11; dt = m00*m11 - m01*m10; xm = max(xm, 0.5*(tr + math.sqrt(max(tr*tr - 4*dt, 0.0)))); TF += tr
        lam_st = mu_s + al_s*a1; lam_th = mu + al*a1
        s0 = w - 1.0/lam_st; s1 = b1 - a1/lam_st                          # score at theta*
        g0 = w - 1.0/lam_th; g1 = b1 - a1/lam_th; h0, h1 = g0 - s0, g1 - s1
        quad = lambda p, q: inv(p, q) if False else (F[1][1]*p*p - 2*F[0][1]*p*q + F[0][0]*q*q)/det
        Gs += quad(s0, s1); Gh += quad(h0, h1); Gg += quad(g0, g1)
        TJ += (F[1][1]*j0 - 2*F[0][1]*j1 + F[0][0]*j2)/det
        v0, v1 = mu - mu_s, al - al_s; Xv += f0*v0*v0 + 2*f1*v0*v1 + f2*v1*v1
        F[0][0] += f0; F[0][1] += f1; F[1][0] += f1; F[1][1] += f2
        det2 = F[0][0]*F[1][1] - F[0][1]**2
        d0 = (F[1][1]*g0 - F[0][1]*g1)/det2; d1 = (-F[0][1]*g0 + F[0][0]*g1)/det2
        mu_n = max(mu - d0, mu_min); al_n = min(max(al - d1, 0.02), al_max); clip += (mu_n != mu - d0) or (al_n != al - d1)
        mu, al = mu_n, al_n; s = a1 + 1.0
    LD = math.log(F[0][0]*F[1][1] - F[0][1]**2) - ld0
    return Gs, TJ, Gh, Xv, Gg, TF, LD, ellm, xm, clip

if __name__ == "__main__":
    reps = int(sys.argv[1]); n = int(sys.argv[2]); mu_s, beta, mu_min = 1.0, 1.0, 0.5
    for eps, th0, rho0, warm in ((0.1, (1.2, 0.8), 60.0, 3000), (0.1, (1.2, 0.8), 20.0, 3000), (0.03, (1.2, 0.93), 300.0, 20000)):
        out = [one(400+r, mu_s, beta, eps, n, th0, rho0, warm, mu_min) for r in range(reps)]
        m = lambda k: sum(o[k] for o in out)/reps
        ell = max(o[7] for o in out); xm = max(o[8] for o in out); r_ = xm/math.log(1+xm)
        print(f"eps={eps} start={th0} rho0={rho0} n={n} reps={reps}  (steps that hit the projection: {sum(o[9] for o in out)/reps:.1f} per run)")
        print(f"  (1) E sum s'F^-1 s = {m(0):.3f}   vs   E sum tr(F^-1 dJ) = {m(1):.3f}")
        print(f"  (2) E sum h'F^-1 h = {m(2):.4f}   <=?  (ell_max/mu_min) E sum v'dF v = ({ell:.3f}/{mu_min})*{m(3):.3f} = {ell/mu_min*m(3):.4f}")
        print(f"  (3) max over runs of [sum tr(F^-1 dF) / logdet] = {max(o[5]/o[6] for o in out):.3f}   <=?  r(x_max={xm:.2f}) = {r_:.3f}")
        print(f"  split: E sum g'F^-1 g = {m(4):.3f};  (1+eta)E s-part + (1+1/eta)E h-part at eta=1: {2*m(0)+2*m(2):.3f}", flush=True)
