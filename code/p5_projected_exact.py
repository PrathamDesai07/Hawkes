"""K=1 online Fisher scoring with a compact box Theta_c and EXACT F-norm projection.
 Theta_c = [mu_lo, mu_hi] x [al_lo, al_hi]  (known a priori; contains theta*).  Ratio set: lam_theta/lam* in [c-, c+] holds by construction,
   c+ = max(mu_hi/mu*, al_hi/al*),  c- = min(mu_lo/mu*, al_lo/al*).
 Exact projection in the F-norm onto a box in R^2: the optimum of the convex QP over the box is among {y, edge minimisers, corners}.
 Checks: (i) non-expansiveness ||v_{i+1}||_{F_i} <= ||u_i||_{F_i} with u_i = theta_i - F_i^-1 g_i - theta*;
         (ii) E R_n  vs  P0 + E sum g'F_i^-1 g (unprojected-step quantity), (iii) fraction of steps where the projection is active, (iv) R_n - log T_n."""
import math, random, sys
from p5_online_newton_k1 import fisher_inc, xq, wq

def fnorm2(F, d0, d1): return F[0][0]*d0*d0 + 2*F[0][1]*d0*d1 + F[1][1]*d1*d1

def project_F(F, y, lo, hi):
    if lo[0] <= y[0] <= hi[0] and lo[1] <= y[1] <= hi[1]: return y, False
    cands = []
    for i in (0, 1):
        j = 1 - i
        for b in (lo[i], hi[i]):
            xj = y[j] - F[j][i]*(b - y[i])/F[j][j]; xj = min(max(xj, lo[j]), hi[j])
            x = [0, 0]; x[i] = b; x[j] = xj; cands.append(tuple(x))
    best = min(cands, key=lambda x: fnorm2(F, x[0]-y[0], x[1]-y[1]))
    return best, True

def run(seed, mu_s, beta, eps, n, th0, rho0, lo, hi, warm=0, mode="exact"):
    rng = random.Random(seed); al_s = (1-eps)*beta; mu, al = th0; lam_bar = mu_s/eps
    F = [[rho0/lam_bar, 0.0], [0.0, rho0*(lam_bar/beta)**2/lam_bar]]; s = 0.0
    for _ in range(warm):
        lam0 = mu_s + al_s*s; w = 0.0
        while True:
            w += rng.expovariate(lam0)
            if rng.random()*lam0 <= mu_s + al_s*s*math.exp(-beta*w): break
        s = s*math.exp(-beta*w) + 1.0
    P0 = fnorm2(F, mu-mu_s, al-al_s); R = G = 0.0; Tt = 0.0; active = 0; viol = 0; cmax = 0.0
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
        f0, f1, f2 = fisher_inc(mu, al, beta, s, w); F[0][0] += f0; F[0][1] += f1; F[1][0] += f1; F[1][1] += f2
        g0 = w - 1.0/(mu + al*a1); g1 = b1 - a1/(mu + al*a1)
        det = F[0][0]*F[1][1] - F[0][1]**2
        d0 = (F[1][1]*g0 - F[0][1]*g1)/det; d1 = (-F[0][1]*g0 + F[0][0]*g1)/det
        G += g0*d0 + g1*d1
        y = (mu - d0, al - d1)
        if mode == "exact": x, act = project_F(F, y, lo, hi)
        else: x = (min(max(y[0], lo[0]), hi[0]), min(max(y[1], lo[1]), hi[1])); act = x != y
        active += act
        lhs = fnorm2(F, x[0]-mu_s, x[1]-al_s); rhs = fnorm2(F, y[0]-mu_s, y[1]-al_s)
        viol += lhs > rhs*(1+1e-9)
        cmax = max(cmax, max(x[0]/mu_s, x[1]/al_s))
        mu, al = x; s = a1 + 1.0; Tt += w
    return R, G, P0, Tt, active/n, viol/n, cmax

if __name__ == "__main__":
    reps = int(sys.argv[1]); n = int(sys.argv[2]); mu_s, beta = 1.0, 1.0
    print("Theta_c = [0.5 mu*, 2 mu*] x [0.5 al*, beta(1-eps/2)];  start (2 mu*, 0.5 al*) = a corner of Theta_c, rho0 = 1")
    print(f"{'eps':>5} {'mode':>6} | {'E R_n':>7} {'R-log T':>8} | {'P0+E G':>8} {'bound ok':>8} | {'frac steps proj active':>22} | {'frac steps violating nonexp':>27}")
    for eps in (0.3, 0.1, 0.03):
        al_s = (1-eps)*beta; lo = (0.5*mu_s, 0.5*al_s); hi = (2*mu_s, beta*(1-eps/2))
        for mode in ("exact", "diag"):
            out = [run(500+r, mu_s, beta, eps, n, (2*mu_s, 0.5*al_s), 1.0, lo, hi, warm=int(3/eps**2), mode=mode) for r in range(reps)]
            m = lambda k: sum(o[k] for o in out)/reps
            Rm = m(0); bnd = m(2) + m(1)
            print(f"{eps:5.2f} {mode:>6} | {Rm:7.2f} {Rm - sum(math.log(o[3]) for o in out)/reps:8.2f} | {bnd:8.2f} {str(Rm <= bnd):>8} | {m(4):22.4f} | {m(5):27.5f}", flush=True)
