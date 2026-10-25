"""E5 ingredients (K=1, projected online Newton as in p5_bootstrap_check).  Per gap i (F = F_{i-1}, v = th_i - th*):
 Q_i = s'F^-1 s (s = b - a/lam*),  q_i = tr(F^-1 DJ_i) = E[Q_i|past],  lbar_i = sup_{y in [0,s_i]} phi'F^-1phi/lam*,  E_i = Lambda*(tau_i) ~ Exp(1),
 pathwise claim (P1): Q_i <= (2 lbar_i/mu*)(1+E_i^2);   (P2): h'F^-1h <= (lbar_i/mu_lo)(a'v)^2/(lam* lam_th);
 identity (P3): (a'v)^2/(lam* lam_th)(tau) = v'dF v + nu_i,  nu_i = increment of the martingale with integrand w_A=(phi'v)^2/(lam_th lam*), |w_A|<=kap1^2 c+ ."""
import math, random, sys
from p5_online_newton_k1 import fisher_inc, xq, wq
from p5_projected_exact import project_F

def run(seed, mu_s, beta, eps, n, th0, rho0, lo, hi, warm):
    rng = random.Random(seed); al_s = (1-eps)*beta; mu, al = th0; lam_bar = mu_s/eps
    F = [[rho0/lam_bar, 0.0], [0.0, rho0*(lam_bar/beta)**2/lam_bar]]; s = 0.0
    for _ in range(warm):
        lam0 = mu_s + al_s*s; w = 0.0
        while True:
            w += rng.expovariate(lam0)
            if rng.random()*lam0 <= mu_s + al_s*s*math.exp(-beta*w): break
        s = s*math.exp(-beta*w) + 1.0
    viol1 = 0; viol2 = 0; SQ = Sq = SH = SX = SNA = SG = 0.0; bAmax = 0.0; lb_hist = []; worstP1 = 0.0
    Qc = []; qc = []; Hc = []; NAc = []; Xc = []
    for i in range(1, n+1):
        lam0 = mu_s + al_s*s; w = 0.0
        while True:
            w += rng.expovariate(lam0)
            if rng.random()*lam0 <= mu_s + al_s*s*math.exp(-beta*w): break
        e = math.exp(-beta*w); a1 = s*e; b1 = s*(1-e)/beta; v0, v1 = mu - mu_s, al - al_s
        Ls = mu_s*w + al_s*b1                                   # Lambda*(tau) = E_i
        det = F[0][0]*F[1][1] - F[0][1]**2
        Fi = ((F[1][1]/det, -F[0][1]/det), (-F[0][1]/det, F[0][0]/det))
        quad = lambda x0, x1: Fi[0][0]*x0*x0 + 2*Fi[0][1]*x0*x1 + Fi[1][1]*x1*x1
        lbar = max(quad(1.0, y)/(mu_s + al_s*y) for y in [s*k/60 for k in range(61)])
        s0 = w - 1.0/(mu_s + al_s*a1); s1 = b1 - a1/(mu_s + al_s*a1)         # s = b - a/lam*  (b=(w,b1), a=(1,a1))
        Q = quad(s0, s1)
        j0, j1, j2 = fisher_inc(mu_s, al_s, beta, s, w); q = Fi[0][0]*j0 + 2*Fi[0][1]*j1 + Fi[1][1]*j2
        lam_th = mu + al*a1; lam_s = mu_s + al_s*a1
        c = (2*lbar/mu_s)*(1 + Ls*Ls)
        if Q > c*(1+1e-9): viol1 += 1
        worstP1 = max(worstP1, Q/c)
        av = v0 + v1*a1; g0 = w - 1.0/lam_th; g1 = b1 - a1/lam_th
        h0 = (1.0/lam_s - 1.0/lam_th)*0 ; # h = g - s : (a/lam* - a/lam_th) with sign: g=b-a/lam_th, s=b-a/lam*  => h = a/lam* - a/lam_th
        h0 = 1.0/lam_s - 1.0/lam_th; h1 = a1/lam_s - a1/lam_th
        H = quad(h0, h1); Hb = (lbar/lo[0])*av*av/(lam_s*lam_th)
        if H > Hb*(1+1e-9): viol2 += 1
        f0, f1, f2 = fisher_inc(mu, al, beta, s, w); Xi = f0*v0*v0 + 2*f1*v0*v1 + f2*v1*v1
        nu = av*av/(lam_s*lam_th) - Xi   # increment of N^A (E[nu|past]=0 under theta*) ; note int (phi'v)^2/lam_th du = Xi
        SQ += Q; Sq += q; SH += H; SX += Xi; SNA += nu
        F[0][0] += f0; F[0][1] += f1; F[1][0] += f1; F[1][1] += f2
        det = F[0][0]*F[1][1] - F[0][1]**2
        d0 = (F[1][1]*g0 - F[0][1]*g1)/det; d1 = (-F[0][1]*g0 + F[0][0]*g1)/det; SG += g0*d0 + g1*d1
        (mu, al), _ = project_F(F, (mu - d0, al - d1), lo, hi); s = a1 + 1.0
    return dict(viol1=viol1, viol2=viol2, SQ=SQ, Sq=Sq, SH=SH, SX=SX, SNA=SNA, SG=SG, worstP1=worstP1, ld=math.log(det/(rho0*rho0*(lam_bar/beta)**2/lam_bar**2)))

if __name__ == "__main__":
    reps = int(sys.argv[1]); mu_s, beta = 1.0, 1.0
    for eps, n in ((0.3, 40000), (0.1, 80000), (0.03, 120000)):
        al_s = (1-eps)*beta; lo = (0.5*mu_s, 0.5*al_s); hi = (2*mu_s, beta*(1-eps/2))
        out = [run(2100+r, mu_s, beta, eps, n, (1.4, 0.8*al_s), 1.0, lo, hi, int(3/eps**2)) for r in range(reps)]
        mean = lambda f: sum(f(o) for o in out)/reps; q = lambda v, p: sorted(v)[min(int(p*len(v)), len(v)-1)]
        print(f"eps={eps} n={n} reps={reps}")
        print(f"  (P1) violations of Q<=(2 lbar/mu)(1+E^2): {sum(o['viol1'] for o in out)} of {reps*n};  max Q/bound = {max(o['worstP1'] for o in out):.3f}")
        print(f"  (P2) violations of h'F^-1h <= (lbar/mu_lo)(a'v)^2/(lam* lam_th): {sum(o['viol2'] for o in out)}")
        r = [o['SQ']/o['Sq'] for o in out]
        print(f"  sum Q / sum q: median {q(r,.5):.3f}, 5% {q(r,.05):.3f}, 95% {q(r,.95):.3f}, max {max(r):.3f}   (E5: <= 1+eps')")
        print(f"  E sum q = {mean(lambda o:o['Sq']):.2f}, E logdet(F_n/F_0) = {mean(lambda o:o['ld']):.2f}, E G = {mean(lambda o:o['SG']):.2f}, E sum Q = {mean(lambda o:o['SQ']):.2f}, E sum H = {mean(lambda o:o['SH']):.4f}")
        z = [(o['SNA'], o['SX']) for o in out]
        print(f"  N^A_n / X_n: median {q([a/b for a,b in z],.5):.3f}, 95% {q([a/b for a,b in z],.95):.3f}, max {max(a/b for a,b in z):.3f}  (claim: N^A <= 0.72 X + b_A x)", flush=True)
