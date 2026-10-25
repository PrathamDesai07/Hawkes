"""Problem 12: (a) grid dependence of the fast-block eigenvalues (n=1, alpha_k/beta_k=1/K);
(b) K=3 eigenvector overlaps: which eigenvector is the ridge, which is the shape mode, at eps=0.1 and 0.03.
Ridge eigenvalue prediction  lambda_ridge ~ eps/(mu (2 mu nu1 - 1)),  nu1=sum alpha_k/beta_k^2."""
import numpy as np
from p12_fastblock import fast_block
from p2_identifiability import mc

print("(a) smallest / largest fast eigenvalue, geometric grid beta_k=r^(k-1):")
for r in (1.5, 2.0, 4.0):
    row = []
    for K in (2, 3, 4, 5):
        beta = r**np.arange(K); ev = fast_block(beta, beta/K); row.append(f"K={K}: {ev[0]:.2e}/{ev[-1]:.3f}")
    print(f"  r={r}:  " + "   ".join(row), flush=True)

beta = np.array([1., 2., 4.]); K = 3; mu = 1.5
print("\n(b) K=3, beta=(1,2,4), mu=1.5")
for eps, n_ev in ((0.1, 250000), (0.03, 350000)):
    alpha = (1-eps)*beta/K; nu1 = (alpha/beta**2).sum()
    I, _ = mc(mu, beta, eps, n_ev, seed=3); ev, V = np.linalg.eigh(I)
    ridge = np.concatenate([[1.0], -eps*alpha/mu]); ridge /= np.linalg.norm(ridge)
    # predicted fast eigenvectors (in alpha-space, mu-component 0), from the corrected block at this alpha
    u = 1/beta; u /= np.linalg.norm(u); Q, _ = np.linalg.qr(np.column_stack([u, np.eye(K)[:, :K-1]])); B = Q[:, 1:]
    a_vec = lambda om: 1/(beta*(beta+1j*om)); D = lambda om: (alpha/(beta*(beta+1j*om))).sum()
    from scipy.integrate import quad
    M = np.zeros((K-1, K-1))
    for i in range(K-1):
        for j in range(K-1):
            M[i, j] = quad(lambda om: ((np.conj(a_vec(om)@B[:, i])*(a_vec(om)@B[:, j]))/abs(D(om))**2).real, 0, np.inf, limit=400)[0]/np.pi
    fe, fv = np.linalg.eigh((M+M.T)/2); fvecs = [np.concatenate([[0.0], B@fv]) for fv in fv.T] if False else [np.concatenate([[0.0], B@fv]) for fv in fv_ ] if False else None
    fvecs = [np.concatenate([[0.0], B@fvv]) for fvv in fv.T]
    print(f" eps={eps}: eigenvalues {np.round(ev,4).tolist()}; predicted ridge {eps/(mu*(2*mu*nu1-1)):.4f}, predicted fast {np.round(fe,4).tolist()}")
    for i in range(K+1):
        ov = [abs(V[:, i]@ridge)] + [abs(V[:, i]@f) for f in fvecs]
        print(f"    eigvec {i} (lambda={ev[i]:.4f}): overlap with ridge {ov[0]:.2f}, with predicted fast modes {np.round(ov[1:],2).tolist()}", flush=True)
