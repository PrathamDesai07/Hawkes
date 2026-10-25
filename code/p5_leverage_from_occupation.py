"""Pathwise check (K=1): from J_T >= (O/L) Qm  to  point-leverage bound.
   ell_orc(T) = max_{s in [0,s_cap]} phi(s)' J_T^{-1} phi(s) / (mu+alpha s)           (oracle accumulated information, phi=(1,s))
   ell_bd(T)  = (L/O) * max_{s in [0,s_cap]} phi(s)' Qm^{-1} phi(s) / (mu+alpha s)      (consequence of the two-interval bound)
 Also the T0 required by the Chebyshev statement (q=0.1) next to the empirical time at which the 95% quantile of ell_orc drops below mu_min/4."""
import math, sys
from p5_fisher_conc import path, stationary_law_moments, lmin2

def phi_quad(M, s):          # (1,s) M^{-1} (1,s)', M = (a,b,c) symmetric 2x2
    a, b, c = M; det = a*c - b*b; return (c - 2*b*s + a*s*s)/det

if __name__ == "__main__":
    reps = int(sys.argv[1]); mu, beta = 1.0, 1.0
    for eps in (0.3, 0.1, 0.03):
        ES, ES2 = stationary_law_moments(mu, beta, eps); al = (1-eps)*beta; gam = beta*eps
        c1, c2, d = 0.5*ES, 1.8*ES, 0.15*ES; L = mu + al*(c2 + d); Qm = (1.0, 0.5*(c1+c2), 0.5*(c1*c1+c2*c2) - 2*d*d)
        s_cap = 4*ES; grid = [s_cap*k/200 for k in range(201)]
        Tgrid = [k/gam for k in (2, 5, 10, 20, 50)]
        res = [path(31+r, mu, beta, eps, Tgrid[-1], int(20/eps**2), Tgrid) for r in range(reps)]
        print(f"eps={eps}: s_cap=4ES={s_cap:.1f}, L={L:.1f}")
        print("   T*gamma | ell_orc: median  95%-q | ell_bd: median  95%-q | bd/orc (median) | any violation of ell_orc <= ell_bd?")
        for k, tg in enumerate(Tgrid):
            orc, bd, viol = [], [], 0
            for r in res:
                T, a, b, c, o1, o2 = r[k]; Jm = (a, b, c); O = min(o1, o2)
                if O <= 0: continue
                lo = max(phi_quad(Jm, s)/(mu + al*s) for s in grid); lb = (L/O)*max(phi_quad(Qm, s)/(mu + al*s) for s in grid)
                orc.append(lo); bd.append(lb); viol += lo > lb*(1+1e-9)
            orc.sort(); bd.sort(); n = len(orc)
            print(f"   {tg*gam:6.1f} | {orc[n//2]:9.4f} {orc[int(.95*n)]:9.4f} | {bd[n//2]:10.4f} {bd[int(.95*n)]:10.3f} | {bd[n//2]/orc[n//2]:6.1f} | {viol}/{n}", flush=True)
        # Chebyshev T0 (m from long-run plateau-with-ramps ~ pi-mass(I)/2 is a guess; use measured occupation fractions)
        m2 = 0.05; sig = math.sqrt(ES2 + 3*ES*ES); Lip = 1/(d/2); q = 0.1
        print(f"   Chebyshev T0*gamma (q=0.1, m~0.05, Lip=1/(d/2)): {16*Lip*sig/(gam*m2**1.5*q)*gam:.2e}")
