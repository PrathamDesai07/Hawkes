"""Problem 3: moments and mixing of the Markov state S=(S_1..S_K) as n->1.
Z = sum_k (alpha_k/beta_k) S_k satisfies  A Z = mu - eps*lam  (generator), lam = mu + alpha.S.
(1) exact stationary 2nd moments (closed linear system) vs the proved bound on E Z^2;
(2) exact mean relaxation rate gamma (root of sum alpha_k/(beta_k-gamma)=1) vs eps*beta_min and eps/sum(alpha/beta^2);
(3) ordered-coupling check: E[Delta Z_t] follows exp(M t), M=-diag(beta)+1 alpha^T."""
import numpy as np

def setup(mu, beta, eps, w=None):
    K = len(beta); w = np.ones(K)/K if w is None else w
    return (1-eps)*w*beta, K            # alpha_k = (1-eps) w_k beta_k, n = 1-eps

def moments(mu, beta, alpha):
    K = len(beta); n = (alpha/beta).sum()
    M = -np.diag(beta) + np.outer(np.ones(K), alpha)
    m = np.linalg.solve(-M, mu*np.ones(K))                 # E S
    # second moments: -(b_j+b_k) M2_jk + sum_m a_m (M2_mj + M2_mk) = -[mu(m_j+m_k+1) + sum a_m m_m]
    I = np.eye(K); L = -(np.kron(np.diag(beta), I) + np.kron(I, np.diag(beta)))
    L += np.kron(np.outer(np.ones(K), alpha), I) * 0 
    # build explicitly (index (j,k) -> j*K+k)
    A = np.zeros((K*K, K*K)); rhs = np.zeros(K*K)
    for j in range(K):
        for k in range(K):
            r = j*K+k; A[r, r] -= beta[j]+beta[k]
            for mm in range(K):
                A[r, mm*K+k] += alpha[mm]; A[r, mm*K+j] += alpha[mm]
            rhs[r] = -(mu*(m[j]+m[k]+1) + alpha@m)
    M2 = np.linalg.solve(A, rhs).reshape(K, K)
    return m, M2

def gamma_root(alpha, beta):
    # smallest decay rate: largest real eigenvalue of M is -gamma
    K = len(beta); M = -np.diag(beta) + np.outer(np.ones(K), alpha)
    return -max(np.linalg.eigvals(M).real)

if __name__ == "__main__":
    mu = 1.0
    print("(1) E Z^2 vs bound  (mu/(2 b_min eps^2)) [ (2mu + n^2 b_max)/b_min + n^2 eps ]")
    print(f"{'beta':>12} {'eps':>6} {'E lam':>8} {'E lam^2':>10} {'eps^2 E lam^2':>13} {'E Z^2':>10} {'bound':>10} {'ratio':>6}")
    for beta in (np.array([1.]), np.array([1., 2.]), np.array([1., 4.]), np.array([1., 2., 4.]), np.array([1., 10.])):
        for eps in (0.3, 0.1, 0.03, 0.01, 0.003):
            alpha, K = setup(mu, beta, eps); n = 1-eps
            m, M2 = moments(mu, beta, alpha)
            El = mu + alpha@m; El2 = mu**2 + 2*mu*alpha@m + alpha@M2@alpha
            c = alpha/beta; EZ2 = c@M2@c
            bd = mu/(2*beta.min()*eps**2)*((2*mu+n**2*beta.max())/beta.min() + n**2*eps)
            print(f"{str(beta.tolist()):>12} {eps:6.3f} {El:8.2f} {El2:10.2f} {eps**2*El2:13.3f} {EZ2:10.2f} {bd:10.2f} {EZ2/bd:6.2f}")
    print("\n(2) relaxation rate gamma of E S (exact), vs eps*b_min and eps/sum(alpha/beta^2)")
    print(f"{'beta':>12} {'eps':>6} {'gamma':>9} {'gamma/eps':>10} {'eps*b_min':>10} {'eps/sum(a/b^2)':>15}")
    for beta in (np.array([1.]), np.array([1., 2.]), np.array([1., 4.]), np.array([1., 10.]), np.array([0.3, 1., 3.])):
        for eps in (0.3, 0.1, 0.01):
            alpha, K = setup(mu, beta, eps); g = gamma_root(alpha, beta)
            print(f"{str(beta.tolist()):>12} {eps:6.3f} {g:9.5f} {g/eps:10.4f} {eps*beta.min():10.5f} {eps/(alpha/beta**2).sum():15.5f}")
