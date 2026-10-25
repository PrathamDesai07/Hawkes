"""Problem 3, general K. A=-diag(beta)+1 alpha^T, -gamma = slowest eigenvalue, l = left eigenvector (>0), c=l^T 1 (normalised to 1).
(a) Z=l^T S obeys L Z = -gamma Z + c*mu exactly; moment recursion (upper bound)
    p*gamma*M_p <= p*c*mu*M_{p-1} + sum_{j>=2} C(p,j) c^j (mu*M_{p-j} + kappa*M_{p-j+1}),  kappa=max_k alpha_k/l_k  (lam <= mu + kappa Z).
(b) Monotone coupling (same Poisson random measure, thinning): E[l^T|S1_t - S2_t|] <= exp(-gamma t) l^T|s1-s2|, equality if ordered."""
import numpy as np
from math import comb

def setup(alpha, beta):
    K = len(beta); A = -np.diag(beta) + np.outer(np.ones(K), alpha)
    w, V = np.linalg.eig(A.T); i = np.argmax(w.real); gamma = -w[i].real
    l = np.abs(V[:, i].real); l = l/l.sum()
    kappa = np.max(alpha/l)
    return gamma, l, kappa

def zmoments_bound(mu, gamma, kappa, P):
    M = [1.0]
    for p in range(1, P+1):
        r = p*mu*M[p-1] + sum(comb(p, j)*(mu*M[p-j] + kappa*M[p-j+1]) for j in range(2, p+1))
        M.append(r/(p*gamma))
    return M

def sim_stationary_Z(mu, alpha, beta, l, T, P, seed=0):
    xq, wq = np.polynomial.legendre.leggauss(16); xq = .5*(xq+1); wq = .5*wq
    rng = np.random.default_rng(seed); S = np.zeros(len(beta)); t = 0.0; acc = np.zeros(P+1)
    for _ in range(int(T)):
        lam = mu + alpha@S
        w = 0.0
        while True:
            w += rng.exponential(1/lam)
            if rng.random()*lam < mu + alpha@(S*np.exp(-beta*w)): break
        s = w*xq; Zs = (S[None, :]*np.exp(-np.outer(s, beta)))@l
        for p in range(1, P+1): acc[p] += w*np.sum(wq*Zs**p)
        S = S*np.exp(-beta*w) + 1.0; t += w
    return acc/t

def coupled(mu, alpha, beta, s1, s2, t_grid, n_rep, seed=0):
    rng = np.random.default_rng(seed); tmax = t_grid[-1]; out = np.zeros(len(t_grid))
    for _ in range(n_rep):
        S1, S2 = s1.copy(), s2.copy(); t = 0.0; gi = 0; rec = np.zeros(len(t_grid))
        while gi < len(t_grid):
            L1, L2 = mu + alpha@S1, mu + alpha@S2; Lm = max(L1, L2)
            w = rng.exponential(1/Lm)
            while gi < len(t_grid) and t + w >= t_grid[gi]:   # record decayed state at grid time
                d = t_grid[gi]-t; rec[gi] = (l*np.abs(S1*np.exp(-beta*d)-S2*np.exp(-beta*d))).sum(); gi += 1
            if gi >= len(t_grid): break
            e = np.exp(-beta*w); S1 = S1*e; S2 = S2*e; t += w
            u = rng.random()*Lm
            if u <= mu + alpha@S1: S1 = S1 + 1.0
            if u <= mu + alpha@S2: S2 = S2 + 1.0
        out += rec
    return out/n_rep

if __name__ == "__main__":
    mu = 1.0
    for beta in (np.array([1.0]), np.array([1.0, 4.0])):
        K = len(beta)
        for eps in (0.2, 0.05):
            alpha = (1-eps)*beta/K
            g, l, kap = setup(alpha, beta)
            M = zmoments_bound(mu, g, kap, 4)
            sm = sim_stationary_Z(mu, alpha, beta, l, 3e4/eps**0.5, 4)
            print(f"beta={beta} eps={eps}: gamma={g:.4f}  E Z^p bound {[f'{x:.3g}' for x in M[1:]]}  MC {[f'{x:.3g}' for x in sm[1:]]}  (E Z = mu/gamma = {mu/g:.3g})")
    print("\nCoupling contraction, beta=[1,4], eps=0.2")
    beta = np.array([1.0, 4.0]); eps = .2; alpha = (1-eps)*beta/2; g, l, kap = setup(alpha, beta)
    tg = np.array([0.5, 1, 2, 4, 6]); print("gamma", g)
    for name, s1, s2 in (("ordered", np.array([3., 2.]), np.array([0., 0.])), ("unordered", np.array([4., 0.]), np.array([0., 3.]))):
        d0 = (l*np.abs(s1-s2)).sum(); r = coupled(mu, alpha, beta, s1, s2, tg, 4000)
        print(name, "E d_l(t):", np.round(r, 4), " bound d0*exp(-g t):", np.round(d0*np.exp(-g*tg), 4))
