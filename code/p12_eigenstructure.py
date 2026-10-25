"""Problem 12: eigenstructure of the Fisher information I(theta*) in (mu, alpha_1..alpha_K) coordinates, K=2,3.
Predictions (stationary, n=1-eps):
  ridge (smallest) direction  v_r ~ (1, -eps*alpha/mu): mu up, all alpha_k scaled by 1-eps*delta/mu  (keeps E lam = mu/eps fixed, kernel shape fixed);
     Rayleigh upper bound  v^T I v/|v|^2 = ((1+eps)^2 m - (eps/mu)(1+2eps))/(1+eps^2|alpha|^2/mu^2), m=I_mumu
  top direction ~ (0, alpha_hat) with eigenvalue >= (mu/eps-2mu)/|alpha|^2;
  K-1 middle directions: fast modes, w with sum_k w_k/beta_k = 0, eigenvalues ~ those of P G P, G_jk=1/(beta_j+beta_k)."""
import numpy as np, sys
from p2_identifiability import mc

def fast_gram_eigs(beta):
    K = len(beta); G = 1/(beta[:, None]+beta[None, :])
    u = 1/beta; u = u/np.linalg.norm(u)           # slow direction in w-space: sum w_k/beta_k
    # basis of the orthogonal complement of u
    Q, _ = np.linalg.qr(np.column_stack([u, np.eye(K)[:, :K-1]])); B = Q[:, 1:]
    return np.linalg.eigvalsh(B.T@G@B)

def analyse(beta, mu, eps, n_ev, seed):
    K = len(beta); alpha = (1-eps)*beta/K
    I, _ = mc(mu, beta, eps, n_ev, seed=seed); m = I[0, 0]
    ev, V = np.linalg.eigh(I)
    vr = np.concatenate([[1.0], -eps*alpha/mu]); vr /= np.linalg.norm(vr)
    cos_r = abs(V[:, 0]@vr)
    va = np.concatenate([[0.0], alpha/np.linalg.norm(alpha)]); cos_t = abs(V[:, -1]@va)
    ray = ((1+eps)**2*m - eps/mu*(1+2*eps))/(1+eps**2*(alpha@alpha)/mu**2)
    lam_top_lb = (mu/eps-2*mu)/(alpha@alpha)
    print(f" beta={beta.tolist()} mu={mu} eps={eps}: eigenvalues {np.round(ev,4).tolist()}")
    print(f"    smallest: |cos(eigvec, (1,-eps alpha/mu))| = {cos_r:.4f};  Rayleigh bound {ray:.4f} vs lambda_min {ev[0]:.4f};  lambda_min/eps = {ev[0]/eps:.2f}")
    print(f"    largest : |cos(eigvec, (0,alpha_hat))| = {cos_t:.4f};  lambda_max {ev[-1]:.2f} >= {lam_top_lb:.2f} (rigorous);  lambda_max*eps = {ev[-1]*eps:.3f}")
    print(f"    middle  : {np.round(ev[1:-1],4).tolist()}   predicted fast-Gram eigenvalues {np.round(fast_gram_eigs(beta),4).tolist()}", flush=True)

if __name__ == "__main__":
    which = sys.argv[1]
    if which == "K2":
        for eps, n in ((0.1, 250000), (0.03, 350000), (0.01, 700000)): analyse(np.array([1., 2.]), 1.5, eps, n, 3)
    else:
        for eps, n in ((0.1, 250000), (0.03, 350000), (0.01, 700000)): analyse(np.array([1., 2., 4.]), 1.5, eps, n, 3)
