"""General-K tracking engine with CAUSAL sampling (mu_hat at the last event before each grid time).
 truth: beta (K,), alpha_k = (1-eps) beta_k w_k (w=equal n-shares by default), mu_t = mu_fn(t) (piecewise-constant on dt=0.1 for rw).
 estimators:  'adaptive' : theta=(mu, alpha_1..alpha_K) unknown, F <- (1-eta)F + dF_i (quadrature), theta -= F^{-1} g_i
              'known'    : only mu estimated (alpha true), F <- (1-eta)F + int du/lam_hat, mu -= g/F"""
import numpy as np, math, random
xg, wg = np.polynomial.legendre.leggauss(8); xg = 0.5*(xg+1); wg = 0.5*wg

def info_inc(th, beta, S, w, nsub=3):
    K = len(beta); F = np.zeros((K+1, K+1)); h = w/nsub
    for j in range(nsub):
        u = h*(j + xg); Su = S[None, :]*np.exp(-np.outer(u, beta)); phi = np.hstack([np.ones((len(u), 1)), Su]); lam = phi@th
        F += (h*wg/lam*phi.T)@phi
    return F

def run(est, beta, eps, mu_fn, mu_max, eta, T, seed, mu_ref, w_share=None, grid=5.0, alpha_true=None):
    rng = random.Random(seed); beta = np.asarray(beta, float); K = len(beta)
    ws = np.ones(K)/K if w_share is None else np.asarray(w_share)
    al_s = (1-eps)*ws*beta if alpha_true is None else np.asarray(alpha_true)
    th = np.concatenate([[mu_ref], al_s]); r = mu_ref/eps; S = np.zeros(K); Tt = 0.0
    nest = K+1 if est == 'adaptive' else 1
    F = np.eye(nest)*(1.0/eta)*0.05; ev = []; tmu = []
    n_cap = 1 - eps/2
    while Tt < T:
        lam0 = mu_max + al_s@S; w = 0.0
        while True:
            w += rng.expovariate(lam0)
            if rng.random()*lam0 <= mu_fn(Tt+w) + al_s@(S*np.exp(-beta*w)): break
        e = np.exp(-beta*w); a = np.concatenate([[1.0], S*e]); b = np.concatenate([[w], S*(1-e)/beta])
        g = b - a/(a@th)
        dF = info_inc(th, beta, S, w)
        if est == 'adaptive':
            F = (1-eta)*F + dF; th = th - np.linalg.solve(F, g)
            th[0] = max(th[0], 0.1); th[1:] = np.maximum(th[1:], 0.01); nh = (th[1:]/beta).sum()
            if nh > n_cap: th[1:] *= n_cap/nh
        else:
            F = (1-eta)*F + dF[0, 0]; th[0] = max(th[0] - g[0]/F[0, 0], 0.1)
        S = S*e + 1.0; Tt += w; ev.append((Tt, th[0]))
    ev = np.array(ev); tg = np.arange(1, int(T/grid)+1)*grid
    idx = np.clip(np.searchsorted(ev[:, 0], tg, side='right') - 1, 0, len(ev)-1)
    return tg, ev[idx, 1]
