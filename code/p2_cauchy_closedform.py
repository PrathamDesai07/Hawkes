"""Closed-form pieces for Problem 2.  (a) rigorous: Cov_pi(S) >= (E lam / 4) G,  G_ij = 1/(beta_i+beta_j)  [spectral argument: |1-phi_hat|<=1+n<=2];
 (b) lam_min(M) >= lam_min(C)/(1+2|m|^2), m_k = mu/(eps beta_k) exactly;  (c) lam_min(G) >= 1/sum_ij |(G^-1)_ij|, (G^-1)_ij explicit (Cauchy inverse),
 and the geometric-grid closed form  |(G^-1)_ij| <= 4 beta_i beta_j P_i P_j/(beta_i+beta_j),  P_i = prod_{k!=i} (beta_i+beta_k)/|beta_i-beta_k| <= exp(4 r/(r-1)^2).
 Checks each against the exact M (linear system of p2_truncation_bound.second_moment)."""
import numpy as np
from p2_truncation_bound import second_moment

def cauchy_inv_bound(beta):
    K = len(beta); P = np.ones(K)
    for i in range(K):
        for k in range(K):
            if k != i: P[i] *= (beta[i]+beta[k])/abs(beta[i]-beta[k])
    Ginv = np.zeros((K, K))
    for i in range(K):
        for j in range(K):
            # exact Cauchy inverse entry magnitude: 4 b_i b_j P_i P_j/(b_i+b_j)
            Ginv[i, j] = 4*beta[i]*beta[j]*P[i]*P[j]/(beta[i]+beta[j])
    return Ginv, P

print(f"{'K':>2} {'r':>4} {'eps':>5} | {'lmin(G)':>9} {'1/sum|Ginv|':>11} | {'lmin(C) exact':>13} {'(mu/4eps)lmin(G)':>17} {'ok':>3} | {'lmin(M)':>8} {'lmin(C)/(1+2|m|^2)':>19} {'ok':>3}")
mu = 1.0
for K, r in ((1, 2), (2, 2), (2, 4), (3, 2), (3, 4), (4, 4)):
    beta = float(r)**np.arange(K)
    for eps in (0.3, 0.1, 0.03):
        alpha = (1-eps)*beta/K; M = second_moment(mu, alpha, beta)
        m = M[0, 1:]; C = M[1:, 1:] - np.outer(m, m)
        assert np.allclose(m, mu/(eps*beta), rtol=1e-8), (m, mu/(eps*beta))
        G = 1.0/(beta[:, None] + beta[None, :]); lG = np.linalg.eigvalsh(G)[0]
        Ginv, P = cauchy_inv_bound(beta); lbG = 1.0/Ginv.sum()
        # verify |Ginv exact| entries
        Ge = np.abs(np.linalg.inv(G)); assert np.allclose(Ge, Ginv, rtol=1e-6, atol=1e-9*Ginv.max()), "Cauchy inverse formula"
        lC = np.linalg.eigvalsh(C)[0]; lbC = (mu/(4*eps))*lG
        lM = np.linalg.eigvalsh(M)[0]; lbM = lC/(1+2*(m@m)); lbM2 = lbC/(1+2*(m@m))
        okC = lC >= lbC*(1-1e-9); okM = lM >= lbM*(1-1e-9) and lM >= lbM2*(1-1e-9)
        print(f"{K:2d} {r:4d} {eps:5.2f} | {lG:9.2e} {lbG:11.2e} | {lC:13.3e} {lbC:17.3e} {str(okC):>3} | {lM:8.4f} {lbM2:19.3e} {str(okM):>3}")
