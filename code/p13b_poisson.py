"""Poisson-equation concentration (K=1): for f with |f(S)-f(S')| <= L|S-S'|, the fluctuation of int_0^T f(S_s) ds around T*pi(f) is a martingale with
 <M>_T <= (L/gamma)^2 int lam ds, plus boundary terms |u(S_0)|+|u(S_T)| <= (L/gamma)(S_0+S_T+E_pi S);  gamma = eps*beta.
Check: sd over windows of (int_window f - T pi(f)) vs  (L/gamma) sqrt(E Lambda_T)  for f = lam (L=alpha) and f = g = 1/max(lam, z) (L_g = alpha/z^2), z = c1/eps."""
import numpy as np, math, random, sys
eps = float(sys.argv[1]); c1 = float(sys.argv[2]); nev = int(sys.argv[3]); mu = 1.0; beta = 1.0; al = (1-eps)*beta; r = mu/eps; gam = beta*eps; z = c1/eps
xg, wg = np.polynomial.legendre.leggauss(12)
rng = random.Random(11); S = 0.0; Tt = 0.0; ts = []; Ilam = []; Ig = []
for i in range(nev):
    lam0 = mu + al*S; w = 0.0
    while True:
        w += rng.expovariate(lam0)
        if rng.random()*lam0 <= mu + al*S*math.exp(-beta*w): break
    u = w*(xg+1)/2; wt = w*wg/2; lam = mu + al*S*np.exp(-beta*u)
    Ilam.append(np.sum(wt*lam)); Ig.append(np.sum(wt/np.maximum(lam, z))); ts.append(Tt); S = S*math.exp(-beta*w) + 1.0; Tt += w
ts = np.array(ts); cl = np.concatenate([[0], np.cumsum(Ilam)]); cg = np.concatenate([[0], np.cumsum(Ig)]); Ttot = Tt
mean_lam = cl[-1]/Ttot; mean_g = cg[-1]/Ttot
print(f"eps={eps}, z=c1/eps={z:.2f}: pi(lam)~{mean_lam:.3f} (theory mu/eps={mu/eps:.3f}), pi(g)~{mean_g:.5f}")
print(f"{'T0*eps':>7} {'sd(int lam)':>12} {'bound (L/g)sqrt(ELam)+bdry':>28} {'sd(int g)':>11} {'bound':>11}")
sig_S = math.sqrt(max(np.var(np.array([0.0])), 0.0))
for C in (4, 16, 64):
    T0 = C/eps; vl = []; vg = []
    for e_ in range(int(0.05*nev), nev, max(1, nev//3000)):
        s_ = np.searchsorted(ts, ts[e_] - T0)
        if s_ >= e_: continue
        vl.append((cl[e_]-cl[s_]) - T0*mean_lam); vg.append((cg[e_]-cg[s_]) - T0*mean_g)
    L1 = al; Lg = al/z**2; ES = (mu/eps - mu)/al; ELam = mu/eps*T0
    bd = lambda L: (L/gam)*(math.sqrt(ELam) + 2*(ES*1.5))        # martingale sd + generous boundary allowance (S_0+S_T+ES ~ 3 ES)
    print(f"{C:7d} {np.std(vl):12.3f} {bd(L1):28.3f} {np.std(vg):11.5f} {bd(Lg):11.5f}")
