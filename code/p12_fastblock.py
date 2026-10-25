"""Corrected fast-block prediction for the middle eigenvalues of I (n=1):
 for w with sum w_k/beta_k=0:  w^T R w ~ (1/2pi) int |a_w(om)|^2/|D(om)|^2 d om,
 a_w=sum_k w_k/(beta_k (beta_k+i om)),  D(om)=sum_k alpha_k/(beta_k (beta_k+i om))   [1/(1-phi_hat)=1/(i om D)].
Eigenvalues of the quadratic form restricted to the orthogonal complement of u=(1/beta_k) (|w|=1)."""
import numpy as np
from scipy.integrate import quad

def fast_block(beta, alpha):
    K = len(beta); u = 1/beta; u = u/np.linalg.norm(u)
    Q, _ = np.linalg.qr(np.column_stack([u, np.eye(K)[:, :K-1]])); B = Q[:, 1:]    # (K, K-1)
    def entry(i, j):
        f = lambda om: ((np.conj(a_vec(om)@B[:, i]) * (a_vec(om)@B[:, j]))/abs(D(om))**2).real
        v = quad(f, 0, np.inf, limit=400)[0]*2/(2*np.pi)               # even in omega
        return v
    a_vec = lambda om: 1/(beta*(beta+1j*om)); D = lambda om: (alpha/(beta*(beta+1j*om))).sum()
    M = np.array([[entry(i, j) for j in range(K-1)] for i in range(K-1)])
    return np.linalg.eigvalsh((M+M.T)/2)

if __name__ == "__main__":
    for beta in (np.array([1., 2.]), np.array([1., 2., 4.]), np.array([1., 4.]), np.array([1., 3., 9.])):
        K = len(beta); alpha = beta/K     # n = 1
        print(f"beta={beta.tolist()}: corrected fast eigenvalues {np.round(fast_block(beta, alpha),4).tolist()}")
