"""K=2 (beta=(1,2), mu=1.5) stationary: (1) inverse moments E[lam^-p]/eps^p, p = 1, 1.2 (needs p < k_K); (2) window-information concentration c_win(T0 = C/eps) (whitened by an independent MC estimate of I);
(3) conditioning of the correlation matrix of the information inner product on the basis (1, lam, f) (f = fast direction, sum w_k/beta_k = 0), vs the K=1 value 1 - Theta^{-1/2}."""
import numpy as np, math, random, sys
from p13b_engine import quad_info
from p2_identifiability import mc
eps = float(sys.argv[1]); nev = int(sys.argv[2]); beta = np.array([1.0, 2.0]); K = 2; mu = 1.5; al = (1-eps)*beta/K; r = mu/eps
I_ind, _ = mc(mu, beta, eps, 300000, seed=17); Iev = I_ind/r; Ih = np.linalg.cholesky(np.linalg.inv(Iev))
nu1 = (al/beta**2).sum(); kK = 2*mu*nu1/((al/beta).sum())
rng = random.Random(23); S = np.zeros(K); Tt = 0.0; rows = []; xg, wg = np.polynomial.legendre.leggauss(8); fn = lambda t: mu
inv1 = inv12 = 0.0
for i in range(nev):
    lam0 = mu + al@S; w = 0.0
    while True:
        w += rng.expovariate(lam0)
        if rng.random()*lam0 <= mu + al@(S*np.exp(-beta*w)): break
    u = w*(xg+1)/2; wt = w*wg/2; Su = S[None, :]*np.exp(-np.outer(u, beta)); lam = mu + Su@al
    inv1 += np.sum(wt/lam); inv12 += np.sum(wt/lam**1.2)
    rows.append((Tt, quad_info(fn, al, beta, S, w, Tt))); S = S*np.exp(-beta*w) + 1.0; Tt += w
print(f"K=2 eps={eps}, mu={mu}, k_K={kK:.2f}: stationary-run time {Tt:.0f}, event rate {nev/Tt:.2f} (theory {r:.2f})")
print(f" (1) inverse moments: E[1/lam] = {inv1/Tt:.5f} = {inv1/Tt/eps:.3f} eps ;  E[lam^-1.2] = {inv12/Tt:.5f} = {inv12/Tt/eps**1.2:.3f} eps^1.2")
ts = np.array([x[0] for x in rows]); Js = np.array([x[1] for x in rows]); cs = np.concatenate([np.zeros((1, 3, 3)), np.cumsum(Js, axis=0)], axis=0)
print(" (2) c_win(T0=C/eps): median / 10% / 1% / min of lambda_min of the whitened window information per unit time")
for C in (4, 16, 64):
    T0 = C/eps; vals = []
    for e_ in range(int(0.1*nev), nev, max(1, nev//2000)):
        s_ = np.searchsorted(ts, ts[e_] - T0)
        if s_ >= e_: continue
        M = Ih.T@((cs[e_]-cs[s_])/(T0*r))@Ih; vals.append(np.min(np.linalg.eigvalsh(M)))
    vals = np.array(vals); print(f"   C={C:3d}: {np.median(vals):.3f} / {np.quantile(vals,0.1):.3f} / {np.quantile(vals,0.01):.3f} / {vals.min():.3f}  ({len(vals)} windows)")
# (3) correlation matrix of the information inner product on the basis b0=(1,0,0) [constant], b1=(mu,alpha) [lambda], b2=(0,w) [fast]
w = np.array([1/beta[1], -1/beta[0]]); w = w/np.linalg.norm(w)       # sum w_k/beta_k = 0
B = np.array([[1.0, 0, 0], [mu, al[0], al[1]], [0, w[0], w[1]]]); Gm = B@I_ind@B.T; Cr = Gm/np.sqrt(np.outer(np.diag(Gm), np.diag(Gm)))
Th = mu*(I_ind[0, 0])/eps
print(f" (3) correlation matrix of (1, lam, f):\n{np.round(Cr,3)}\n     lambda_min = {np.linalg.eigvalsh(Cr)[0]:.3f};  K=1 analogue 1-Theta^(-1/2) with Theta = mu m/eps = {Th:.3f}: {1-Th**-0.5:.3f}")
