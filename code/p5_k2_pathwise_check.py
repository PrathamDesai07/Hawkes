"""Pathwise checks for K >= 2 with a polytope F-projection (audit):
   (A) Lemma C1 per step:  D_i <= D_{i-1} - X_i - m_i + G_i,  D_i = ||th_{i+1}-th*||^2_{F_i}, X_i = v'dF v, m_i = -2 a'v/lam_th(tau) + 2 int (phi'v) lam*/lam_th du,
       G_i = g'F_i^{-1} g;  (B) R_i <= X_i pathwise;  (C) Lemma E5a: s'F_{i-1}^{-1}s <= (2 lbar/mu*)(1+E^2) with lbar = sup_u phi'F^-1 phi/lam*;
   (D) lemma E5d identity: (a'v)^2/(lam* lam_th)(tau) = X_i + nu_i with E[nu_i|past] = 0 (checked through the mean of nu_i);
   (E) trace bound used in the (corrected) Lemma E5c:  tr(F^-1 dF) <= lbar E/(c^- mu*)  (the version with c^+ in the first draft was wrong)."""
import numpy as np, math, random, sys
from scipy.optimize import minimize
xg, wg = np.polynomial.legendre.leggauss(12); xg = 0.5*(xg+1); wg = 0.5*wg

def project(F, y, lo, hi, beta, ncap):
    if np.all(y >= lo) and np.all(y <= hi) and (y[1:]/beta).sum() <= ncap: return y
    cons = [{'type': 'ineq', 'fun': lambda x: ncap - (x[1:]/beta).sum(), 'jac': lambda x: np.concatenate([[0.0], -1.0/beta])}]
    res = minimize(lambda x: (x-y)@F@(x-y), np.clip(y, lo, hi), jac=lambda x: 2*F@(x-y), bounds=list(zip(lo, hi)), constraints=cons, method='SLSQP', options={'ftol': 1e-14, 'maxiter': 300})
    return res.x

def run(seed, beta, eps, mu_s, n, th0, rho0, lo, hi):
    rng = random.Random(seed); beta = np.asarray(beta, float); K = len(beta)
    al_s = (1-eps)*beta/K; ths = np.concatenate([[mu_s], al_s]); ncap = 1 - eps/2; lam_bar = mu_s/eps
    F = np.diag(np.concatenate([[rho0/lam_bar], rho0*(lam_bar/beta)**2/lam_bar])); th = np.array(th0, float); S = np.zeros(K)
    cmin = min(lo[0]/mu_s, (lo[1:]/al_s).min()); cmax = max(hi[0]/mu_s, (hi[1:]/al_s).max())
    def gap(S):
        lam0 = ths[0] + al_s@S; w = 0.0
        while True:
            w += rng.expovariate(lam0)
            if rng.random()*lam0 <= ths[0] + al_s@(S*np.exp(-beta*w)): return w
    for _ in range(int(3/eps**2)):
        w = gap(S); S = S*np.exp(-beta*w) + 1.0
    D_prev = (th-ths)@F@(th-ths); worst = {'C1': -1e9, 'R': -1e9, 'E5a': -1e9, 'E5c': -1e9}; nus = []; nproj = 0
    for i in range(n):
        w = gap(S); e = np.exp(-beta*w); a = np.concatenate([[1.0], S*e]); b = np.concatenate([[w], S*(1-e)/beta]); v = th - ths
        u = w*xg; phi = np.hstack([np.ones((len(u), 1)), S[None, :]*np.exp(-np.outer(u, beta))]); lt = phi@th; ls = phi@ths; pv = phi@v
        Xi = (w*wg*pv**2/lt).sum(); dF = (w*wg/lt*phi.T)@phi
        mi = -2*(a@v)/(a@th) + 2*(w*wg*pv*ls/lt).sum()
        Ri = (w*wg*ls*(lt/ls - 1 - np.log(lt/ls))).sum()
        # E5a / E5c quantities (with F = F_{i-1})
        Fi_prev = np.linalg.inv(F); su = np.linspace(0, S.max()*1.0, 1)  # placeholder
        uu = np.linspace(0, w, 400); phg = np.hstack([np.ones((400, 1)), S[None, :]*np.exp(-np.outer(uu, beta))]); lam_g = phg@ths
        lbar = (np.einsum('ij,jk,ik->i', phg, Fi_prev, phg)/lam_g).max()
        # sup over the whole possible gap (u in [0,infty)) is attained at u=0 or u->inf; add the quiet limit phi=(1,0)
        lbar = max(lbar, (Fi_prev[0, 0])/ths[0])
        Ei = (w*wg*ls).sum(); g = b - a/(a@th); s = b - a/(a@ths)
        Q = s@Fi_prev@s; E5a = Q - (2*lbar/mu_s)*(1 + Ei**2)
        trdF = np.trace(Fi_prev@dF); E5c = trdF - lbar*Ei/(cmin*mu_s)
        nu = (a@v)**2/((a@ths)*(a@th)) - Xi
        F = F + dF; y = th - np.linalg.solve(F, g); th_new = project(F, y, lo, hi, beta, ncap); nproj += (th_new != y).any()
        G = g@np.linalg.solve(F, g); D = (th_new-ths)@F@(th_new-ths)
        worst['C1'] = max(worst['C1'], D - (D_prev - Xi - mi + G)); worst['R'] = max(worst['R'], Ri - Xi)
        worst['E5a'] = max(worst['E5a'], E5a); worst['E5c'] = max(worst['E5c'], E5c); nus.append(nu)
        th = th_new; D_prev = D; S = S*e + 1.0
    return worst, np.mean(nus), np.std(nus)/math.sqrt(len(nus)), nproj/n

if __name__ == "__main__":
    reps = int(sys.argv[1]); n = int(sys.argv[2]); beta = [1.0, 2.0]; mu_s = 1.5
    for eps in (0.4, 0.2):
        al_s = (1-eps)*np.array(beta)/2
        lo = np.concatenate([[0.5*mu_s], 0.4*al_s]); hi = np.concatenate([[2*mu_s], 2.0*al_s]); th0 = np.concatenate([[1.4*mu_s], 0.7*al_s])
        out = [run(11+r, beta, eps, mu_s, n, th0, 1.0, lo, hi) for r in range(reps)]
        w = {k: max(o[0][k] for o in out) for k in out[0][0]}
        print(f"K=2 eps={eps}: max over steps/reps of  (C1) D_i-(D_prev-X-m+G) = {w['C1']:.2e};  (B) R_i-X_i = {w['R']:.2e};  (E5a) Q-bound = {w['E5a']:.2e};  (E5c, c^- version) tr(F^-1 dF)-bound = {w['E5c']:.2e};  mean nu = {np.mean([o[1] for o in out]):.2e} +- {np.mean([o[2] for o in out]):.1e};  projection active {np.mean([o[3] for o in out]):.3f}", flush=True)
