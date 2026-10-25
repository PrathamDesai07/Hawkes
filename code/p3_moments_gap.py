"""Problem 3: uniform moments and spectral gap as n -> 1.
K=1 stationary moments from  p*beta*eps*M_p = p*beta*mu*M_{p-1} + sum_{j>=2} C(p,j) alpha^j M_{p-j+1}
(derived from the generator L f = -beta(lam-mu) f' + lam (f(lam+alpha)-f(lam)) applied to lam^p).
Gap: L lam^j = -j*beta*eps*lam^j + lower order, so eigenvalues are -j*beta*eps (K=1).
General K: mean dynamics dE[S]=A E[S], A=-diag(beta)+1 alpha^T; slowest rate gamma = -max Re eig(A) ~ eps/sum(alpha_k/beta_k^2)."""
import numpy as np
from math import comb

def moments(mu, a, b, P):
    eps = 1 - a/b; M = [1.0]
    for p in range(1, P+1):
        rhs = p*b*mu*M[p-1] + sum(comb(p, j)*a**j*M[p-j+1] for j in range(2, p+1)) if p >= 2 else b*mu
        if p >= 2:   # j=p+... term M_{p-j+1} for j=p is M_1: all indices <= p-1, fine
            pass
        M.append(rhs/(p*b*eps))
    return M

def sim_moments(mu, a, b, T, P, seed=0):
    rng = np.random.default_rng(seed); lam = mu; t = 0.0; acc = np.zeros(P+1)
    while t < T:
        c = lam - mu; w = 0.0
        while True:
            w += rng.exponential(1/lam)
            if rng.random()*lam < mu + c*np.exp(-b*w): break
        # int_0^w (mu + c e^{-b s})^p ds, binomial expansion
        for p in range(1, P+1):
            acc[p] += sum(comb(p, j)*mu**(p-j)*c**j*(1-np.exp(-j*b*w))/(j*b) if j else mu**p*w for j in range(p+1))
        lam = mu + c*np.exp(-b*w) + a; t += w
    return acc/t

def gap_general(mu_unused, alpha, beta):
    K = len(beta); A = -np.diag(beta) + np.outer(np.ones(K), alpha)
    ev = np.linalg.eigvals(A); return -ev.real.max()

if __name__ == "__main__":
    mu, b = 1.0, 1.0
    print("K=1 moments: recursion vs simulation, and M_p*(eps/mu)^p (should stay bounded)")
    for eps in (0.4, 0.2, 0.1):
        a = (1-eps)*b; M = moments(mu, a, b, 4); s = sim_moments(mu, a, b, 4e4*eps*5, 4)
        print(f"eps={eps}: rec", [f"{x:.4g}" for x in M[1:]], "sim", [f"{x:.4g}" for x in s[1:]],
              "norm", [f"{M[p]*(eps/mu)**p:.3f}" for p in range(1, 5)])
    print("\nnormalised moments M_p (eps/mu)^p as eps->0 (K=1):")
    for eps in (0.1, 0.01, 0.001, 1e-4):
        a = (1-eps)*b; M = moments(mu, a, b, 6)
        print(f"eps={eps:g}", [f"{M[p]*(eps/mu)**p:.3f}" for p in range(1, 7)])
    print("\nspectral gap (slowest mean-relaxation rate) vs eps/sum(alpha/beta^2):")
    for beta in (np.array([1.0]), np.array([1.0, 4.0]), np.array([0.5, 2.0, 8.0])):
        K = len(beta)
        for eps in (0.2, 0.05, 0.01):
            alpha = (1-eps)*beta/K
            g = gap_general(None, alpha, beta); pred = eps/np.sum(alpha/beta**2)
            print(f"beta={beta} eps={eps}: gap={g:.5f} pred={pred:.5f} gap/eps={g/eps:.3f}")
