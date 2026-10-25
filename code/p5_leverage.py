"""Problem 7 / leverage quantities along the K=1 online Newton run (projected iteration).
 point leverage   ell_i = sup over gap of phi(u)' F_{i-1}^{-1} phi(u) / lam*(u)      (sampled at u = 0, tau/4, tau/2, 3tau/4, tau and u = infinity, i.e. phi=(1,0))
 gap leverage     x_i   = lambda_max( F_{i-1}^{-1} dF_i )
 ratio constant   c_i   = max(mu_i/mu*, al_i/al*)  (so lam_theta <= c_i lam*)
Reports per index window: mean and max over reps of ell, x, c, and the first index i0 after which B=(1+1/eta)*ell/mu_min < 1 for good (eta=1)."""
import math, random, sys
from p5_online_newton_k1 import fisher_inc, xq, wq

def run(seed, mu_s, beta, eps, n, th0, rho0, warm, mu_min=0.1, al_min=0.02):
    rng = random.Random(seed); al_s = (1-eps)*beta; eps_c = eps/2; mu, al = th0; lam_bar = mu_s/eps; al_max = beta*(1-eps_c)
    F = [[rho0/lam_bar, 0.0], [0.0, rho0*(lam_bar/beta)**2/lam_bar]]; s = 0.0
    for _ in range(warm):
        lam0 = mu_s + al_s*s; w = 0.0
        while True:
            w += rng.expovariate(lam0)
            if rng.random()*lam0 <= mu_s + al_s*s*math.exp(-beta*w): break
        s = s*math.exp(-beta*w) + 1.0
    ld0 = math.log(F[0][0]*F[1][1]); rec = []; R = 0.0                # (ell, x, c)
    for i in range(1, n+1):
        lam0 = mu_s + al_s*s; w = 0.0
        while True:
            w += rng.expovariate(lam0)
            if rng.random()*lam0 <= mu_s + al_s*s*math.exp(-beta*w): break
        e = math.exp(-beta*w); a1 = s*e; b1 = s*(1-e)/beta
        det = F[0][0]*F[1][1] - F[0][1]**2; ell = 0.0
        for fr in (0, .25, .5, .75, 1.0, None):
            su = (s*math.exp(-beta*fr*w) if fr is not None else 0.0); lam_s = mu_s + al_s*su
            q = (F[1][1] - 2*F[0][1]*su + F[0][0]*su*su)/det        # phi' F^-1 phi, phi=(1,su)
            ell = max(ell, q/lam_s)
        Lt = mu*w + al*b1; Ls = mu_s*w + al_s*b1; lg = 0.0
        for xx, wt in zip(xq, wq):
            u = 0.5*w*(xx+1); su = s*math.exp(-beta*u); ls_ = mu_s + al_s*su; lt = mu + al*su
            lg += 0.5*w*wt*ls_*math.log(lt/ls_)
        R += Lt - Ls - lg
        f0, f1, f2 = fisher_inc(mu, al, beta, s, w)
        # lambda_max(F^-1 dF) for 2x2
        m00 = (F[1][1]*f0 - F[0][1]*f1)/det; m01 = (F[1][1]*f1 - F[0][1]*f2)/det
        m10 = (-F[0][1]*f0 + F[0][0]*f1)/det; m11 = (-F[0][1]*f1 + F[0][0]*f2)/det
        tr = m00 + m11; dt = m00*m11 - m01*m10; x = 0.5*(tr + math.sqrt(max(tr*tr - 4*dt, 0.0)))
        rec.append((ell, x, max(mu/mu_s, al/al_s)))
        F[0][0] += f0; F[0][1] += f1; F[1][0] += f1; F[1][1] += f2
        g0 = w - 1.0/(mu + al*a1); g1 = b1 - a1/(mu + al*a1)
        det = F[0][0]*F[1][1] - F[0][1]**2
        d0 = (F[1][1]*g0 - F[0][1]*g1)/det; d1 = (-F[0][1]*g0 + F[0][0]*g1)/det
        mu = max(mu - d0, mu_min); al = min(max(al - d1, al_min), al_max); s = a1 + 1.0
    return rec, math.log(F[0][0]*F[1][1] - F[0][1]**2) - ld0, R

if __name__ == "__main__":
    reps = int(sys.argv[1]); n = int(sys.argv[2]); mu_s, beta = 1.0, 1.0
    wins = [(1, 10), (10, 100), (100, 1000), (1000, 10000), (10000, n)]
    for eps, th0, rho0, warm, label in ((0.1, (1.0, 0.9), 1.0, 3000, "start at truth, rho0=1"), (0.1, (1.0, 0.9), 20.0, 3000, "start at truth, rho0=20"),
                                         (0.03, (1.0, 0.97), 1.0, 20000, "start at truth, rho0=1"), (0.1, (2.0, 0.5), 1.0, 3000, "far start, rho0=1")):
        out = [run(300+r, mu_s, beta, eps, n, th0, rho0, warm) for r in range(reps)]
        print(f"eps={eps} {label}  (n={n}, reps={reps})   E[log det F_n/F_0] = {sum(o[1] for o in out)/reps:.2f}")
        print("   index window |   mean ell    max ell |   mean x     max x |  max c ")
        for lo, hi in wins:
            if lo >= n: continue
            hi = min(hi, n); sel = [[r[lo-1:hi] for r in [o[0]]][0] for o in out]
            ell = [v[0] for run_ in sel for v in run_]; xs = [v[1] for run_ in sel for v in run_]; cs = [v[2] for run_ in sel for v in run_]
            print(f"   [{lo:5d},{hi:5d}] | {sum(ell)/len(ell):9.4f} {max(ell):9.3f} | {sum(xs)/len(xs):8.4f} {max(xs):8.3f} | {max(cs):6.2f}")
        # burn-in for B<1 (eta=1 => B=2*ell/mu_min)
        i0s = []
        for o in out:
            last_bad = 0
            for i, r in enumerate(o[0], 1):
                if 2*r[0]/0.1 >= 1: last_bad = i
            i0s.append(last_bad + 1)
        i0s.sort(); print(f"   first index after which B=2*ell/mu_min<1 holds for the rest of the run: median {i0s[len(i0s)//2]}, max {i0s[-1]}", flush=True)
