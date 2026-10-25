"""General-K projected online Newton (eta = 0, i.e. Problem 5's recursion) with an F-norm projection onto the polytope
   Theta_c = {mu in [mu_lo,mu_hi], alpha_k in [al_lo_k, al_hi_k], sum alpha_k/beta_k <= 1 - eps/2}.
 Reports E R_n at checkpoints and the slope dE R_n / d ln T_n (target (K+1)/2).  R_n = sum_i int_gap lam* psi(lam_theta/lam*) du.
 The F-projection is a convex QP; solved by SLSQP only when the Newton point is infeasible."""
import numpy as np, math, random, sys
from scipy.optimize import minimize
xg, wg = np.polynomial.legendre.leggauss(10); xg = 0.5*(xg+1); wg = 0.5*wg

def project(F, y, lo, hi, beta, ncap):
    ok = np.all(y >= lo) and np.all(y <= hi) and (y[1:]/beta).sum() <= ncap
    if ok: return y, False
    cons = [{'type': 'ineq', 'fun': lambda x: ncap - (x[1:]/beta).sum(), 'jac': lambda x: np.concatenate([[0.0], -1.0/beta])}]
    res = minimize(lambda x: (x-y)@F@(x-y), np.clip(y, lo, hi), jac=lambda x: 2*F@(x-y), bounds=list(zip(lo, hi)), constraints=cons, method='SLSQP', options={'ftol': 1e-12, 'maxiter': 200})
    return res.x, True

def run(seed, beta, eps, mu_s, n, th0, rho0, lo, hi, warm, cps, nsub=2):
    rng = random.Random(seed); beta = np.asarray(beta, float); K = len(beta)
    al_s = (1-eps)*beta/K; ths = np.concatenate([[mu_s], al_s]); ncap = 1 - eps/2
    lam_bar = mu_s/eps
    F = np.diag(np.concatenate([[rho0/lam_bar], rho0*lam_bar/beta**2/ 1.0 * (1.0/ lam_bar) * lam_bar])) if False else np.diag(np.concatenate([[rho0/lam_bar], rho0*(lam_bar/beta)**2/lam_bar]))
    th = np.array(th0, float); S = np.zeros(K); R = 0.0; Tt = 0.0; snap = {}
    def gap(S):
        lam0 = ths[0] + al_s@S; w = 0.0
        while True:
            w += rng.expovariate(lam0)
            if rng.random()*lam0 <= ths[0] + al_s@(S*np.exp(-beta*w)): return w
    for _ in range(warm):
        w = gap(S); S = S*np.exp(-beta*w) + 1.0
    nproj = 0
    for i in range(1, n+1):
        w = gap(S); e = np.exp(-beta*w)
        a = np.concatenate([[1.0], S*e]); b = np.concatenate([[w], S*(1-e)/beta])
        dF = np.zeros((K+1, K+1)); h = w/nsub; lg = 0.0
        for j in range(nsub):
            u = h*(j + xg); phi = np.hstack([np.ones((len(u), 1)), S[None, :]*np.exp(-np.outer(u, beta))])
            lt = phi@th; ls = phi@ths
            dF += (h*wg/lt*phi.T)@phi; lg += (h*wg*ls*np.log(lt/ls)).sum()
        R += (b@th - b@ths) - lg
        g = b - a/(a@th); F = F + dF
        y = th - np.linalg.solve(F, g)
        th, act = project(F, y, lo, hi, beta, ncap); nproj += act
        S = S*e + 1.0; Tt += w
        if i in cps: snap[i] = (Tt, R, np.linalg.slogdet(F)[1])
    return snap, nproj/n

if __name__ == "__main__":
    reps = int(sys.argv[1]); n = int(sys.argv[2]); beta = [1.0, 2.0]; mu_s = 1.5
    for eps in (0.4, 0.3):
        K = len(beta); al_s = (1-eps)*np.array(beta)/K
        lo = np.concatenate([[0.5*mu_s], 0.4*al_s]); hi = np.concatenate([[2*mu_s], 2.0*al_s])
        cps = {n//16, n//8, n//4, n//2, n}
        th0 = np.concatenate([[1.3*mu_s], 0.8*al_s])
        outs = [run(3000+r, beta, eps, mu_s, n, th0, 1.0, lo, hi, int(3/eps**2), cps) for r in range(reps)]
        ks = sorted(cps); Rm = [np.mean([o[0][k][1] for o in outs]) for k in ks]; Tm = [np.mean([o[0][k][0] for o in outs]) for k in ks]
        ld = [np.mean([o[0][k][2] for o in outs]) for k in ks]
        sl = [(Rm[i+1]-Rm[i])/math.log(Tm[i+1]/Tm[i]) for i in range(len(ks)-1)]
        sl2 = [(Rm[i+1]-Rm[i])/(ld[i+1]-ld[i]) for i in range(len(ks)-1)]
        print(f"K=2 beta={beta} mu*={mu_s} eps={eps}: E R_n = {[round(x,2) for x in Rm]}  T = {[round(x) for x in Tm]}")
        print(f"   slope dER/dlnT = {[round(x,2) for x in sl]} (target (K+1)/2 = {(K+1)/2});  dER/dlogdet F = {[round(x,2) for x in sl2]} (target 0.5);  projection active in {np.mean([o[1] for o in outs]):.3f} of steps", flush=True)
