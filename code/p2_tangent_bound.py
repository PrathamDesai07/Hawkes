"""Problem 2, a much sharper rigorous lower bound on I = E[a a'/lam] from a TANGENT-LINE minorant of the convex function 1/lam:
      1/lam >= 2/lam0 - lam/lam0^2   for every lam, lam0 > 0     (convexity),
  hence   I  >=  N(lam0) := (2/lam0) M - (1/lam0^2) M3,     M = E[a a'],  M3 = E[lam a a'],
  both EXACT (linear systems for the moments of S of order <= 3, from E[A f] = 0), so   lam_min(I) >= sup_{lam0} lam_min(N(lam0)).
  Moments: A(S_i S_j S_l) = -(b_i+b_j+b_l) S_i S_j S_l + lam [ S_iS_j + S_iS_l + S_jS_l + S_i + S_j + S_l + 1 ] with lam = mu + sum_k alpha_k S_k.
  Compared with: MC lam_min(I) (from the notes) and the previous closed-form chain.  mpmath is used for the eigenvalues (cancellation at small eps)."""
import numpy as np, itertools, mpmath as mp
mp.mp.dps = 60

def moments(mu, alpha, beta):
    K = len(beta)
    A = -np.diag(beta) + np.outer(np.ones(K), alpha)
    m = np.linalg.solve(-A, mu*np.ones(K)); El = mu + alpha@m
    # second moments Q_ij : (b_i+b_j) Q_ij = E[lam S_i] + E[lam S_j] + E lam,  E[lam S_i] = mu m_i + sum_k alpha_k Q_ki
    n = K*K; Mx = np.zeros((n, n)); rhs = np.zeros(n)
    for i in range(K):
        for j in range(K):
            r = i*K + j; Mx[r, r] += beta[i] + beta[j]; rhs[r] = mu*(m[i] + m[j]) + El
            for k in range(K): Mx[r, k*K + i] -= alpha[k]; Mx[r, k*K + j] -= alpha[k]
    Q = np.linalg.solve(Mx, rhs).reshape(K, K)
    # third moments T_ijl
    N3 = K**3; M3x = np.zeros((N3, N3)); r3 = np.zeros(N3); idx = lambda i, j, l: (i*K + j)*K + l
    for i, j, l in itertools.product(range(K), repeat=3):
        r = idx(i, j, l); M3x[r, r] += beta[i] + beta[j] + beta[l]
        # E[lam S_iS_j] = mu Q_ij + sum_k alpha_k T_kij  (and the two other pairs), E[lam S_i] = mu m_i + sum_k alpha_k Q_ki
        for (p, q) in ((i, j), (i, l), (j, l)):
            r3[r] += mu*Q[p, q]
            for k in range(K): M3x[r, idx(k, p, q)] -= alpha[k]
        for p in (i, j, l):
            r3[r] += mu*m[p] + alpha@Q[:, p]
        r3[r] += El
    T = np.linalg.solve(M3x, r3).reshape(K, K, K)
    return m, Q, T, El

def matrices(mu, alpha, beta):
    K = len(beta); m, Q, T, El = moments(mu, alpha, beta); d = K + 1
    M = np.zeros((d, d)); M[0, 0] = 1; M[0, 1:] = M[1:, 0] = m; M[1:, 1:] = Q
    # M3 = E[lam a a'] = mu M + sum_k alpha_k E[S_k a a']
    M3 = mu*M.copy()
    for k in range(K):
        Ek = np.zeros((d, d)); Ek[0, 0] = m[k]; Ek[0, 1:] = Ek[1:, 0] = Q[k, :]; Ek[1:, 1:] = T[k]
        M3 += alpha[k]*Ek
    return M, M3, El

def lam_min_N(M, M3, lam0):
    d = M.shape[0]
    Nmat = mp.matrix(d, d)
    for i in range(d):
        for j in range(d):
            Nmat[i, j] = 2*mp.mpf(M[i, j])/lam0 - mp.mpf(M3[i, j])/lam0**2
    ev = mp.eigsy(Nmat, eigvals_only=True)
    return min(ev)

def best(mu, alpha, beta, eps):
    M, M3, El = matrices(mu, alpha, beta)
    bestv, arg = -1, None
    for f in np.concatenate([np.linspace(0.2, 3.0, 57), np.linspace(3.2, 10.0, 18)]):
        lam0 = mp.mpf(f*El); v = lam_min_N(M, M3, lam0)
        if v > bestv: bestv, arg = v, f
    return float(bestv), arg

if __name__ == "__main__":
    mc = {(1, 2, .3): .1389, (1, 2, .1): .0654, (1, 2, .05): .0383, (2, 2, .3): .0337, (2, 2, .1): .0455, (2, 2, .05): .0435,
          (2, 4, .1): .0677, (3, 2, .3): .0033, (3, 2, .1): .0045, (3, 2, .05): .0051}
    print(f"{'K':>2} {'r':>3} {'eps':>5} {'MC lmin(I)':>10} {'tangent bound':>14} {'bound/MC':>9}  best lam0/E[lam]")
    for (K, r, eps), mval in mc.items():
        beta = float(r)**np.arange(K); alpha = (1-eps)*beta/K
        b, f = best(1.0, alpha, beta, eps)
        print(f"{K:2d} {r:3d} {eps:5.2f} {mval:10.4f} {b:14.5f} {b/mval:9.3f}  {f:.2f}", flush=True)
