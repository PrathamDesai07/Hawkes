"""True-bootstrap ingredients (K=1, exact F-norm projection on Theta_c = [.5mu*,2mu*]x[.5al*,beta(1-eps/2)]):
 D_j = ||th_{j+1}-th*||^2_{F_j}.  Pathwise (projected) :  D_j <= D'_0 - X_j - M_j + G_j,   X=sum v'dF v,  G=sum g'F_i^-1 g,
 M=sum m_i,  m_i = -2 (a'v)/lam_th(tau) + 2 int_0^tau (phi'v) lam*/lam_th du   (zero-mean given the past),  <M> <= (4/c-) E X.
 Reports: (A) pathwise inequality margin, (B) E sum m^2 / E X, (C) quantiles of sup_j |M_j| vs Doob scale, (D) sup_j D_j vs log det,
          (E) last time the iterate is outside the ratio set R(delta), (F) slope of E R_n in log T."""
import math, random, sys
from p5_online_newton_k1 import fisher_inc, xq, wq
from p5_projected_exact import project_F, fnorm2

def run(seed, mu_s, beta, eps, n, th0, rho0, lo, hi, warm, deltas, cps):
    rng = random.Random(seed); al_s = (1-eps)*beta; mu, al = th0; lam_bar = mu_s/eps
    F = [[rho0/lam_bar, 0.0], [0.0, rho0*(lam_bar/beta)**2/lam_bar]]; F0det = F[0][0]*F[1][1]; s = 0.0
    for _ in range(warm):
        lam0 = mu_s + al_s*s; w = 0.0
        while True:
            w += rng.expovariate(lam0)
            if rng.random()*lam0 <= mu_s + al_s*s*math.exp(-beta*w): break
        s = s*math.exp(-beta*w) + 1.0
    D0 = fnorm2(F, mu-mu_s, al-al_s); X = M = G = R = Tt = 0.0; supD = 0.0; supM = 0.0; m2 = 0.0; marg = -1e300
    last_out = {d: 0.0 for d in deltas}; snap = {}
    for i in range(1, n+1):
        lam0 = mu_s + al_s*s; w = 0.0
        while True:
            w += rng.expovariate(lam0)
            if rng.random()*lam0 <= mu_s + al_s*s*math.exp(-beta*w): break
        e = math.exp(-beta*w); a1 = s*e; b1 = s*(1-e)/beta; v0, v1 = mu - mu_s, al - al_s
        Lt = mu*w + al*b1; Ls = mu_s*w + al_s*b1; lg = 0.0; mi_int = 0.0
        for x, wt in zip(xq, wq):
            u = 0.5*w*(x+1); su = s*math.exp(-beta*u); ls = mu_s + al_s*su; lt = mu + al*su; pv = v0 + v1*su
            lg += 0.5*w*wt*ls*math.log(lt/ls); mi_int += 0.5*w*wt*pv*ls/lt
        R += Lt - Ls - lg
        lam_th = mu + al*a1; m_i = -2*(v0 + v1*a1)/lam_th + 2*mi_int; M += m_i; m2 += m_i*m_i
        f0, f1, f2 = fisher_inc(mu, al, beta, s, w); X += f0*v0*v0 + 2*f1*v0*v1 + f2*v1*v1
        F[0][0] += f0; F[0][1] += f1; F[1][0] += f1; F[1][1] += f2
        g0 = w - 1.0/lam_th; g1 = b1 - a1/lam_th; det = F[0][0]*F[1][1] - F[0][1]**2
        d0 = (F[1][1]*g0 - F[0][1]*g1)/det; d1 = (-F[0][1]*g0 + F[0][0]*g1)/det; G += g0*d0 + g1*d1
        (mu, al), _ = project_F(F, (mu - d0, al - d1), lo, hi); s = a1 + 1.0; Tt += w
        Dj = fnorm2(F, mu-mu_s, al-al_s); supD = max(supD, Dj); supM = max(supM, abs(M)); marg = max(marg, Dj - (D0 - X - M + G))
        dev = max(abs(mu/mu_s - 1), abs(al/al_s - 1))
        for d in deltas:
            if dev > d: last_out[d] = Tt
        if i in cps: snap[i] = (Tt, R, supD, supM, math.log(F[0][0]*F[1][1]/F0det))
    return dict(D0=D0, X=X, G=G, M2=m2, supD=supD, supM=supM, marg=marg, last=last_out, snap=snap, T=Tt, R=R)

if __name__ == "__main__":
    reps = int(sys.argv[1]); mu_s, beta = 1.0, 1.0
    for eps, n in ((0.3, 40000), (0.1, 80000), (0.03, 200000)):
        al_s = (1-eps)*beta; lo = (0.5*mu_s, 0.5*al_s); hi = (2*mu_s, beta*(1-eps/2)); deltas = (0.5, 0.3, 0.2, 0.1)
        cps = {n//16, n//8, n//4, n//2, n}
        out = [run(1100+r, mu_s, beta, eps, n, (1.4, 0.8*al_s), 1.0, lo, hi, int(3/eps**2), deltas, cps) for r in range(reps)]
        mean = lambda f: sum(f(o) for o in out)/reps; q = lambda v, p: sorted(v)[min(int(p*len(v)), len(v)-1)]
        print(f"eps={eps} n={n} reps={reps} (d=2, T_end ~ {mean(lambda o: o['T']):.0f}, T_end*gamma ~ {mean(lambda o: o['T'])*beta*eps:.0f})")
        print(f"  (A) max over runs/steps of [D_j - (D0 - X - M + G)] = {max(o['marg'] for o in out):.2e}   (<= 0 means pathwise inequality holds; {'ok' if max(o['marg'] for o in out) <= 1e-8 else 'VIOLATED'})")
        print(f"  (B) E sum m_i^2 / E X = {mean(lambda o: o['M2'])/mean(lambda o: o['X']):.3f}   (theory <= 4/c- = {4/0.5:.1f}, c-=.5 box value)")
        print(f"  (C) sup_j |M_j|: median {q([o['supM'] for o in out], .5):.2f}, 95% {q([o['supM'] for o in out], .95):.2f};  Doob scale sqrt(16 E X/(c- q)) at q=.05: {math.sqrt(16*mean(lambda o: o['X'])/(0.5*0.05)):.1f}")
        print(f"  (D) sup_j D_j: median {q([o['supD'] for o in out], .5):.2f}, 95% {q([o['supD'] for o in out], .95):.2f}, max {max(o['supD'] for o in out):.2f};  E X = {mean(lambda o: o['X']):.2f}, E G = {mean(lambda o: o['G']):.2f}, E logdet(F_n/F_0) = {mean(lambda o: o['snap'][n][4]):.2f}")
        for d in deltas:
            ls = [o['last'][d] for o in out]
            print(f"  (E) last time outside R(delta={d}): median T*gamma = {q(ls, .5)*beta*eps:6.1f}, 90% = {q(ls, .9)*beta*eps:6.1f}, max = {max(ls)*beta*eps:6.1f}   (fraction of runs never outside: {sum(x == 0 for x in ls)/reps:.2f})")
        ks = sorted(cps); 
        Rm = [mean(lambda o, k=k: o['snap'][k][1]) for k in ks]; Tm = [mean(lambda o, k=k: o['snap'][k][0]) for k in ks]
        sl = [(Rm[i+1]-Rm[i])/math.log(Tm[i+1]/Tm[i]) for i in range(len(ks)-1)]
        print(f"  (F) E R_n = {[round(r,2) for r in Rm]};  slope dE R/dlog T between checkpoints = {[round(x,2) for x in sl]}  (target d/2 = 1)", flush=True)
