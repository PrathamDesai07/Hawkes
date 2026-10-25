"""High-precision (mpmath, 60 digits) version of the tangent bound: the moment linear systems are solved in multiprecision rational-like arithmetic, so that the only
 floating-point content is the 60-digit evaluation (no double-precision linear algebra).  Also the information level for the K=2 example used in K2_regret_chain.md."""
import numpy as np, itertools, mpmath as mp
mp.mp.dps = 60

def moments_mp(mu, alpha, beta):
    K = len(beta); mu = mp.mpf(mu); al = [mp.mpf(a) for a in alpha]; be = [mp.mpf(b) for b in beta]
    A = mp.matrix(K, K)
    for i in range(K):
        for j in range(K): A[i, j] = -(be[i] if i == j else 0) + al[j]
    m = mp.lu_solve(-A, mp.matrix([mu]*K)); El = mu + sum(al[k]*m[k] for k in range(K))
    n = K*K; Mx = mp.zeros(n, n); rhs = mp.zeros(n, 1)
    for i in range(K):
        for j in range(K):
            r = i*K + j; Mx[r, r] += be[i] + be[j]; rhs[r] = mu*(m[i] + m[j]) + El
            for k in range(K): Mx[r, k*K + i] -= al[k]; Mx[r, k*K + j] -= al[k]
    Qv = mp.lu_solve(Mx, rhs); Q = [[Qv[i*K + j] for j in range(K)] for i in range(K)]
    N3 = K**3; M3x = mp.zeros(N3, N3); r3 = mp.zeros(N3, 1); idx = lambda i, j, l: (i*K + j)*K + l
    for i, j, l in itertools.product(range(K), repeat=3):
        r = idx(i, j, l); M3x[r, r] += be[i] + be[j] + be[l]
        for (p, q) in ((i, j), (i, l), (j, l)):
            r3[r] += mu*Q[p][q]
            for k in range(K): M3x[r, idx(k, p, q)] -= al[k]
        for p in (i, j, l): r3[r] += mu*m[p] + sum(al[k]*Q[k][p] for k in range(K))
        r3[r] += El
    Tv = mp.lu_solve(M3x, r3); T = [[[Tv[idx(i, j, l)] for l in range(K)] for j in range(K)] for i in range(K)]
    return mu, al, m, Q, T, El

def tangent(mu, alpha, beta, lam0_factors):
    mu, al, m, Q, T, El = moments_mp(mu, alpha, beta); K = len(beta); d = K + 1
    M = mp.zeros(d, d); M[0, 0] = 1
    for i in range(K):
        M[0, i+1] = M[i+1, 0] = m[i]
        for j in range(K): M[i+1, j+1] = Q[i][j]
    M3 = mu*M
    for k in range(K):
        Ek = mp.zeros(d, d); Ek[0, 0] = m[k]
        for i in range(K): Ek[0, i+1] = Ek[i+1, 0] = Q[k][i]
        for i in range(K):
            for j in range(K): Ek[i+1, j+1] = T[k][i][j]
        M3 += al[k]*Ek
    best = None
    for f in lam0_factors:
        lam0 = f*El; N = 2*M/lam0 - M3/lam0**2; ev = mp.eigsy(N, eigvals_only=True); v = min(ev)
        if best is None or v > best[0]: best = (v, f)
    return best, El

if __name__ == "__main__":
    fac = [mp.mpf(x) for x in np.concatenate([np.linspace(0.5, 3.0, 26), np.linspace(3.5, 8.0, 10)])]
    for name, K, beta, mu, eps in (("K=1 r=2", 1, [1.0], 1.0, 0.3), ("K=2 r=2", 2, [1.0, 2.0], 1.0, 0.1), ("K=3 r=2", 3, [1.0, 2.0, 4.0], 1.0, 0.1),
                                   ("K=2 beta=(1,2) mu=1.5 eps=0.4 (K2 note example)", 2, [1.0, 2.0], 1.5, 0.4),
                                   ("K=2 beta=(1,2) mu=1.5 eps=0.3", 2, [1.0, 2.0], 1.5, 0.3)):
        beta = np.array(beta); alpha = (1-eps)*beta/K
        (v, f), El = tangent(mu, alpha, beta, fac)
        print(f"{name}: tangent lower bound on lmin(I) = {float(v):.6g}  (lam0/E lam = {float(f):.2f}, 60-digit arithmetic)", flush=True)
