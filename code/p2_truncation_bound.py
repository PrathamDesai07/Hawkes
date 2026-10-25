"""Problem 2: closed-form-ish bound  lam_min(I) >= (1/L) [ lam_min(M) - sqrt(E|a|^4) * sqrt(E lam^p / L^p) ],  M = E[a a^T].
Rigorous pieces: (i) I >= E[a a^T 1{lam<=L}]/L; (ii) E[(v'a)^2 1{lam>L}] <= sqrt(E (v'a)^4) sqrt(P(lam>L)) (Cauchy-Schwarz),
P(lam>L) <= E lam^p / L^p; (iii) M from the generator (linear system below); (iv) moments of Z=l'S from Lemma B (p3_generalK).
M: E S = m solves A m = -mu*1.  E[S_i S_j] = Q_ij solves (beta_i+beta_j) Q_ij = E[lam S_i] + E[lam S_j] + E[lam],
E[lam S_i] = mu m_i + sum_k alpha_k Q_ki,  E lam = mu + alpha'm."""
import numpy as np
from p3_generalK import setup, zmoments_bound

def second_moment(mu, alpha, beta):
    K = len(beta); A = -np.diag(beta) + np.outer(np.ones(K), alpha)
    m = np.linalg.solve(-A, mu*np.ones(K)); El = mu + alpha@m
    # unknown Q (K x K, symmetric). Equation for each (i,j): (b_i+b_j) Q_ij - sum_k alpha_k (Q_ki + Q_kj) = mu (m_i + m_j) + El
    n = K*K; Mx = np.zeros((n, n)); rhs = np.zeros(n)
    for i in range(K):
        for j in range(K):
            r = i*K + j; Mx[r, r] += beta[i] + beta[j]; rhs[r] = mu*(m[i] + m[j]) + El
            for k in range(K): Mx[r, k*K + i] -= alpha[k]; Mx[r, k*K + j] -= alpha[k]
    Q = np.linalg.solve(Mx, rhs).reshape(K, K)
    M = np.zeros((K+1, K+1)); M[0, 0] = 1; M[0, 1:] = M[1:, 0] = m; M[1:, 1:] = Q
    return M

def bound(mu, alpha, beta, eps, P=(2, 4, 6, 8, 10, 12)):
    K = len(beta); g, l, kap = setup(alpha, beta); lmin = l.min()
    Mz = zmoments_bound(mu, g, kap, max(P))
    lamM = np.linalg.eigvalsh(second_moment(mu, alpha, beta))[0]
    EA4 = 2*(1 + K**2*Mz[4]/lmin**4)                       # |a|^4 <= 2(1 + (sum S_k^2)^2), S_k <= Z/l_min
    best = (0, None)
    for p in P:
        Elp = 2**(p-1)*(mu**p + kap**p*Mz[p])              # lam <= mu + kap Z
        for L in (mu/eps)*np.logspace(0, 12, 400):
            v = (lamM - np.sqrt(EA4)*np.sqrt(Elp/L**p))/L
            if v > best[0]: best = (v, (p, L*eps/mu))
    return best, lamM

if __name__ == "__main__":
    mc = {(1, 2, .3): .1389, (1, 2, .1): .0654, (1, 2, .05): .0383, (2, 2, .3): .0337, (2, 2, .1): .0455, (2, 2, .05): .0435,
          (2, 4, .3): .0522, (2, 4, .1): .0677, (2, 4, .05): .0635, (3, 2, .3): .0033, (3, 2, .1): .0045, (3, 2, .05): .0051}
    print(f"{'K':>2} {'r':>3} {'eps':>5} {'MC lmin(I)':>10} {'lmin(M)':>8} {'bound':>10} {'bound/MC':>9}  (p, L*eps/mu)")
    for (K, r, eps), m in mc.items():
        beta = float(r)**np.arange(K); alpha = (1-eps)*beta/K
        (b, par), lm = bound(1.0, alpha, beta, eps)
        print(f"{K:2d} {r:3d} {eps:5.2f} {m:10.4f} {lm:8.4f} {b:10.3e} {b/m:9.2e}  {par}")
