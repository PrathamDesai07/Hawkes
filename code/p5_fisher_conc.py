"""Problem 5/7 follow-up: concentration of the oracle accumulated Fisher information J_T (K=1), via occupation of two state intervals.
 Q-matrix (exact two-interval lower bound):  v'J_T v >= (O/L) * v' Qm v,  Qm = [[1,(c1+c2)/2],[(c1+c2)/2,(c1^2+c2^2)/2-2 d^2]],
 where O = min_j (time in I_j=[c_j-d,c_j+d]) and L = mu+alpha*(c_2+d)  (max intensity on the two intervals).
 Intervals: c1 = 0.5*ES, c2 = 1.8*ES, d = 0.15*ES, ES = mu/(beta*eps)."""
import math, random, sys
from p5_online_newton_k1 import fisher_inc

def occ(s, w, beta, lo, hi):
    """time within [0,w] during which s*exp(-beta*u) lies in [lo,hi]"""
    if s <= 0: return 0.0
    u_lo = 0.0 if s <= hi else math.log(s/hi)/beta                  # s e^{-bu} <= hi
    u_hi = w if s*math.exp(-beta*w) >= lo and False else (math.log(s/lo)/beta if s > lo else 0.0)   # s e^{-bu} >= lo
    a, b = max(u_lo, 0.0), min(u_hi, w)
    return max(b - a, 0.0)

def stationary_law_moments(mu, beta, eps):
    al = (1-eps)*beta; El = mu/eps; El2 = mu*(2*beta*mu + al*al)/(2*beta*eps*eps)
    ES = (El - mu)/al; ES2 = (El2 - 2*mu*El + mu*mu)/(al*al); return ES, ES2

def path(seed, mu, beta, eps, T_end, warm, Tgrid):
    rng = random.Random(seed); al = (1-eps)*beta; s = 0.0
    ES, ES2 = stationary_law_moments(mu, beta, eps)
    c1, c2, d = 0.5*ES, 1.8*ES, 0.15*ES
    for _ in range(warm):
        lam0 = mu + al*s; w = 0.0
        while True:
            w += rng.expovariate(lam0)
            if rng.random()*lam0 <= mu + al*s*math.exp(-beta*w): break
        s = s*math.exp(-beta*w) + 1.0
    t = 0.0; J = [0.0, 0.0, 0.0]; O = [0.0, 0.0]; out = []; gi = 0
    while gi < len(Tgrid):
        lam0 = mu + al*s; w = 0.0
        while True:
            w += rng.expovariate(lam0)
            if rng.random()*lam0 <= mu + al*s*math.exp(-beta*w): break
        i0, i1, i2 = fisher_inc(mu, al, beta, s, w); J[0] += i0; J[1] += i1; J[2] += i2
        O[0] += occ(s, w, beta, c1-d, c1+d); O[1] += occ(s, w, beta, c2-d, c2+d)
        t += w; s = s*math.exp(-beta*w) + 1.0
        while gi < len(Tgrid) and t >= Tgrid[gi]:
            out.append((t, J[0], J[1], J[2], O[0], O[1])); gi += 1
    return out

def lmin2(a, b, c):  # smallest eigenvalue of [[a,b],[b,c]]
    return 0.5*(a + c) - math.sqrt(0.25*(a - c)**2 + b*b)

if __name__ == "__main__":
    reps = int(sys.argv[1]); mu, beta = 1.0, 1.0
    for eps in (0.3, 0.1, 0.03):
        ES, ES2 = stationary_law_moments(mu, beta, eps); al = (1-eps)*beta
        c1, c2, d = 0.5*ES, 1.8*ES, 0.15*ES; L = mu + al*(c2 + d)
        Qm = (1.0, 0.5*(c1+c2), 0.5*(c1*c1 + c2*c2) - 2*d*d); lamQ = lmin2(Qm[0], Qm[1], Qm[2])
        gam = beta*eps; Tgrid = [k/gam for k in (0.5, 1, 2, 5, 10, 20, 50)]
        res = [path(11+r, mu, beta, eps, Tgrid[-1], int(20/eps**2), Tgrid) for r in range(reps)]
        # long-run occupation probabilities m_j (average of last grid point)
        m = [sum(r[-1][4+j]/r[-1][0] for r in res)/reps for j in (0, 1)]
        print(f"eps={eps}: ES={ES:.2f} intervals I1=[{c1-d:.2f},{c1+d:.2f}] I2=[{c2-d:.2f},{c2+d:.2f}]  pi-mass ~ ({m[0]:.3f},{m[1]:.3f})  L={L:.1f}  lam_min(Qm)={lamQ:.3f}")
        print("   T*gamma | T | lmin(J_T)/T: mean  5%-quantile  min | occupation min_j O_j/T: mean  5%-q | two-interval bound (O/L)lamQ /T: 5%-q | bound/actual(5%-q)")
        for k, tg in enumerate(Tgrid):
            T = res[0][k][0]; lm = sorted(lmin2(r[k][1], r[k][2], r[k][3])/r[k][0] for r in res)
            oc = sorted(min(r[k][4], r[k][5])/r[k][0] for r in res); q5 = lambda v: v[max(int(0.05*len(v)) - 0, 0)]
            bd = sorted((min(r[k][4], r[k][5])/L*lamQ)/r[k][0] for r in res)
            print(f"   {tg*gam:6.1f} | {T:8.1f} | {sum(lm)/len(lm):.5f} {q5(lm):.5f} {lm[0]:.5f} | {sum(oc)/len(oc):.4f} {q5(oc):.4f} | {q5(bd):.5f} | {q5(bd)/q5(lm):.2f}", flush=True)
