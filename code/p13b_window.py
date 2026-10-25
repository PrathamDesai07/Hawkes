"""A4 lemma ingredients, K=1, constant parameters.
 (i)  quadratic-form reduction: for u=(u_mu,u_al), p = u_mu - mu u_al/alpha, q = u_al/alpha:  u^T J_w u = p^2 X + 2 p q T + q^2 L  (X=int 1/lam, L=int lam, T=window length) -> checked against quadrature
 (ii) window concentration: c_win(T0) = lambda_min( Iev^{-1/2} [J_w/(T0 r)] Iev^{-1/2} ) over windows of length T0 = C/eps (J_w = sum of oracle gap informations dJ over whole gaps inside the window)
 (iii) discounted-sum inequality:  eta F_{i-1} >= (e^{-eta m}/c+) eta * J(last m events)  (checked with the realised F and c+)"""
import numpy as np, math, random, sys
from p13b_engine import quad_info
from p14_mech import Iev_matrix
eps = float(sys.argv[1]); nev = int(sys.argv[2]); mu = 1.0; beta = 1.0; al = (1-eps)*beta; r = mu/eps
rng = random.Random(5); S = 0.0; Tt = 0.0; rows = []
fn = lambda t: mu
for i in range(nev):
    lam0 = mu + al*S; w = 0.0
    while True:
        w += rng.expovariate(lam0)
        if rng.random()*lam0 <= mu + al*S*math.exp(-beta*w): break
    dJ = quad_info(fn, np.array([al]), np.array([beta]), np.array([S]), w, Tt)
    # X = int 1/lam, L = int lam over the gap
    u = w*np.polynomial.legendre.leggauss(8)[0]/2 + w/2; wt = w*np.polynomial.legendre.leggauss(8)[1]/2; lam = mu + al*S*np.exp(-beta*u)
    rows.append((Tt, w, S, dJ, np.sum(wt/lam), np.sum(wt*lam)))
    S = S*math.exp(-beta*w) + 1.0; Tt += w
Iev = Iev_matrix(mu, beta, eps); Ih = np.linalg.cholesky(np.linalg.inv(Iev))
# (i) identity on the whole record
Jt = sum(x[3] for x in rows); X = sum(x[4] for x in rows); L = sum(x[5] for x in rows); T = Tt
u_ = np.array([0.7, -0.3]); p = u_[0] - mu*u_[1]/al; q = u_[1]/al
print(f"eps={eps}: (i) u'J u (quadrature) = {u_@Jt@u_:.6f}  vs  p^2 X + 2pq T + q^2 L = {p*p*X + 2*p*q*T + q*q*L:.6f}")
ts = np.array([x[0] for x in rows]); cumJ = np.cumsum([x[3] for x in rows], axis=0) if False else None
Js = np.array([x[3] for x in rows]); cs = np.concatenate([np.zeros((1, 2, 2)), np.cumsum(Js, axis=0)], axis=0)
print("(ii) c_win(T0): quantiles of lambda_min of whitened window information per unit time / its stationary value, windows of length T0 = C/eps")
for C in (1, 2, 4, 8, 16, 32, 64):
    T0 = C/eps; vals = []
    for e_ in range(int(0.1*nev), nev, max(1, nev//2500)):
        te = ts[e_]; s_ = np.searchsorted(ts, te - T0)
        if s_ >= e_: continue
        Jw = cs[e_] - cs[s_]; M = Ih.T@(Jw/(T0*r))@Ih; vals.append(np.min(np.linalg.eigvalsh(M)))
    vals = np.array(vals); print(f"   C={C:3d} (T0={T0:7.1f}, m_events~{T0*r:7.0f}): c_win median {np.median(vals):.3f}, 10% {np.quantile(vals,0.1):.3f}, 1% {np.quantile(vals,0.01):.3f}, min {vals.min():.3f}   (windows {len(vals)})")
