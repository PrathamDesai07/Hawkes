"""Problem 2: identifiability / conditioning of the Fisher information I(theta*)=E[a a^T/lam], a=(1,S_1..S_K).
Monte Carlo (exact thinning, Gauss-Legendre quadrature of a a^T/lam within each gap) + explicit Cauchy-matrix lower bound.
Bound: lam_min(I) >= rho*p*min(1,s^2)/(L*tr C(rho+L)^{-1}), L=2mu/eps, p=P(lam<=L, S_k>=s for all k),
C(r)_ij = 1/(g_i+g_j), g_0=r/2, g_k=beta_k+r/2, tr C^{-1} = sum_i 2 g_i prod_{j!=i}((g_i+g_j)/(g_i-g_j))^2."""
import sys, numpy as np
xq, wq = np.polynomial.legendre.leggauss(24); xq = 0.5*(xq+1); wq = 0.5*wq

def mc(mu, beta, eps, n_ev, seed=1, s_lo=0.5):
    K = len(beta); alpha = (1-eps)*beta/K          # n = sum alpha/beta = 1-eps
    rng = np.random.default_rng(seed)
    S = np.zeros(K); lam = mu; T = 0.0
    A = np.zeros((K+1, K+1)); pw = 0.0; Lcut = 2*mu/eps
    for _ in range(n_ev):
        lam = mu + alpha@S
        w = 0.0
        while True:                                  # thinning; lam decays between events
            w += rng.exponential(1/lam)
            Sw = S*np.exp(-beta*w)
            if rng.random()*lam < mu + alpha@Sw: break
        s = w*xq                                     # quadrature nodes in the gap
        Sn = S[None, :]*np.exp(-np.outer(s, beta))   # (Q,K)
        a = np.hstack([np.ones((len(s), 1)), Sn]); l = mu + Sn@alpha
        A += (w*wq/l*a.T)@a
        if lam <= Lcut and S.min() >= s_lo: pw += w
        T += w; S = Sw + 1.0
    return A/T, pw/T

def cauchy_trinv(beta, r):
    g = np.concatenate([[r/2], beta + r/2]); t = 0.0
    for i in range(len(g)):
        o = np.delete(g, i); t += 2*g[i]*np.prod(((g[i]+o)/(g[i]-o))**2)
    return t

def bound(mu, beta, eps, p, s_lo=0.5):
    L = 2*mu/eps; best = 0
    for rho in np.logspace(-3, 2, 60):
        best = max(best, rho*p*min(1, s_lo**2)/(L*cauchy_trinv(beta, rho+L)))
    return best

if __name__ == "__main__":
    n_ev = int(sys.argv[1]) if len(sys.argv) > 1 else 150000
    mu = 1.0
    print(f"{'K':>2} {'r':>3} {'eps':>5} {'lam_min':>9} {'lmin/eps':>9} {'lam_max':>9} {'cond':>9} {'cond*eps^2':>10} {'p':>5} {'bound':>9}")
    for K in (2, 3):
        for r in (2.0, 4.0):
            beta = r**np.arange(K)
            for eps in (0.3, 0.1, 0.05):
                I, p = mc(mu, beta, eps, n_ev)
                e = np.linalg.eigvalsh(I)
                print(f"{K:2d} {r:3.0f} {eps:5.2f} {e[0]:9.4f} {e[0]/eps:9.3f} {e[-1]:9.3f} {e[-1]/e[0]:9.1f} {e[-1]/e[0]*eps**2:10.2f} {p:5.2f} {bound(mu,beta,eps,p):9.2e}", flush=True)
