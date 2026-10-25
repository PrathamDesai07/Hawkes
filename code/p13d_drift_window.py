"""Drift in mu only (alpha fixed).
 (1) coupling contraction under time-varying mu_t: ordered pair (S^max >= S^min) driven by the same Poisson measure, E d_l(S^max_t,S^min_t) should equal l'e^{Mt}D_0 <= e^{-gamma t} d_0 for ANY mu_t.
 (2) window information under drift, K=1: c_win(T0=C/eps) relative to the FIXED reference metric Iev(mu_ref) and relative to the LOCAL stationary information Iev(mu_t)."""
import numpy as np, math, random, sys
from p13b_engine import quad_info
from p14_mech import Iev_matrix
from scipy.linalg import expm
# ---- (1)
beta = np.array([1.0, 2.0]); eps = 0.4; al = (1-eps)*beta/2; mu_s = 1.5; A = 0.4; P = 60.0     # fast drift to stress the claim
fn = lambda t: mu_s*(1 + A*math.sin(2*math.pi*t/P)); mu_max = mu_s*(1+A)
M = -np.diag(beta) + np.outer(np.ones(2), al); g = -max(np.linalg.eigvals(M).real); om = al/(beta-g)
D0 = np.array([3.0, 2.0]); acc = {t: 0.0 for t in (2, 5, 10, 20)}; R = 1500
for r_ in range(R):
    rng = random.Random(r_); Smin = np.array([1.0, 1.0]); Smax = Smin + D0; t = 0.0; targets = sorted(acc)
    # thinning for the upper process, mark u in [0, lam_bound]; each process fires iff u <= its intensity
    ti = 0
    while ti < len(targets):
        lam_b = mu_max + al@Smax; w = rng.expovariate(lam_b); tn = t + w
        if tn > targets[ti]:
            d_ = targets[ti]-t; Smax = Smax*np.exp(-beta*d_); Smin = Smin*np.exp(-beta*d_); t = targets[ti]; acc[targets[ti]] += om@(Smax-Smin); ti += 1; continue
        Smax = Smax*np.exp(-beta*w); Smin = Smin*np.exp(-beta*w); t = tn; u = rng.random()*lam_b
        if u <= fn(t) + al@Smax: Smax = Smax + 1
        if u <= fn(t) + al@Smin: Smin = Smin + 1
d0 = om@D0
print(f"(1) K=2 beta=(1,2), eps={eps}, mu_t = {mu_s}(1+{A} sin(2 pi t/{P})) (fast drift), gamma = {g:.4f}")
for t in sorted(acc): print(f"   t={t:3d}: E d_l = {acc[t]/R:.4f}; exact l'e^(Mt)D0 = {om@(expm(M*t)@D0):.4f}; bound d0 e^(-gamma t) = {d0*math.exp(-g*t):.4f}")
# ---- (2) K=1 drift windows
mu = 1.0; beta1 = np.array([1.0]); eps = 0.4; al1 = (1-eps)*np.array([1.0]); r = 1/eps; A = 0.3; P = 1200.0
fn1 = lambda t: mu*(1 + A*np.sin(2*np.pi*np.asarray(t)/P)); Iref = Iev_matrix(mu, 1.0, eps); Ihr = np.linalg.cholesky(np.linalg.inv(Iref))
rng = random.Random(41); S = np.zeros(1); Tt = 0.0; rows = []
for i in range(220000):
    lam0 = mu*(1+A) + al1@S; w = 0.0
    while True:
        w += rng.expovariate(lam0)
        if rng.random()*lam0 <= float(fn1(Tt+w)) + al1@(S*np.exp(-beta1*w)): break
    rows.append((Tt, w, quad_info(fn1, al1, beta1, S, w, Tt), float(fn1(Tt)))); S = S*np.exp(-beta1*w) + 1.0; Tt += w
ts = np.array([x[0] for x in rows]); Js = np.array([x[2] for x in rows]); cs = np.concatenate([np.zeros((1, 2, 2)), np.cumsum(Js, axis=0)], axis=0); mus = np.array([x[3] for x in rows])
print("(2) K=1 under sinusoidal mu drift (A=0.3, P=1200): c_win(T0=C/eps) relative to the FIXED reference metric (mu_ref=1) | relative to the LOCAL stationary information")
for C in (16, 64):
    T0 = C/eps; vf = []; vl = []
    for e_ in range(int(0.1*len(rows)), len(rows), max(1, len(rows)//2000)):
        s_ = np.searchsorted(ts, ts[e_]-T0)
        if s_ >= e_: continue
        Jw = (cs[e_]-cs[s_])/T0; mloc = mus[e_]; rl = mloc/eps; Il = Iev_matrix(mloc, 1.0, eps); Ihl = np.linalg.cholesky(np.linalg.inv(Il))
        vf.append(np.min(np.linalg.eigvalsh(Ihr.T@(Jw/r)@Ihr))); vl.append(np.min(np.linalg.eigvalsh(Ihl.T@(Jw/rl)@Ihl)))
    vf, vl = np.array(vf), np.array(vl)
    print(f"   C={C}: fixed-metric  median {np.median(vf):.3f}, 1% {np.quantile(vf,0.01):.3f}, min {vf.min():.3f} | local-metric median {np.median(vl):.3f}, 1% {np.quantile(vl,0.01):.3f}, min {vl.min():.3f}   (stationary K=1 reference: 0.66/0.26 at C=16, 0.84/0.55 at C=64)")
