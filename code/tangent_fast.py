"""Fast tangent bound: O(K^4) moments (Kronecker-sum diagonalisation) + ONE-DIMENSIONAL concave maximisation over t = 1/lam0.
 N(t) = 2 t M - t^2 M3 is operator-concave in t (M3 = E[lam a a'] >= 0), and lam_min is concave and monotone, so t -> lam_min(N(t)) is concave: Brent search instead of a grid."""
import numpy as np, time
from scipy.optimize import minimize_scalar
from complexity_bench import moments_fast

def matrices_fast(mu, alpha, beta):
    K = len(beta); m, Q, T, El = moments_fast(mu, alpha, beta); d = K + 1
    M = np.zeros((d, d)); M[0, 0] = 1; M[0, 1:] = M[1:, 0] = m; M[1:, 1:] = Q
    M3 = mu*M.copy()
    for k in range(K):
        Ek = np.zeros((d, d)); Ek[0, 0] = m[k]; Ek[0, 1:] = Ek[1:, 0] = Q[k, :]; Ek[1:, 1:] = T[k]
        M3 += alpha[k]*Ek
    return M, M3, El

def tangent_fast(mu, alpha, beta):
    M, M3, El = matrices_fast(mu, alpha, beta)
    f = lambda t: -np.linalg.eigvalsh(2*t*M - t*t*M3)[0]
    res = minimize_scalar(f, bounds=(0.05/El, 5.0/El), method='bounded', options={'xatol': 1e-10/El})
    return -res.fun, 1.0/(res.x*El), res.nfev

if __name__ == "__main__":
    for K, r, eps, ref in ((1, 2, .3, 0.0791027), (2, 2, .1, 0.0239641), (3, 2, .1, 0.00285116), (2, 2, .01, 0.00260), (3, 2, .01, 0.00268)):
        beta = float(r)**np.arange(K); alpha = (1-eps)*beta/K
        t0 = time.perf_counter(); v, f, nf = tangent_fast(1.0, alpha, beta); dt = time.perf_counter() - t0
        print(f"K={K} r={r} eps={eps}: fast tangent bound {v:.6g} (multiprecision-grid reference {ref}), lam0/E lam = {f:.2f}, {nf} eigen-decompositions, {1000*dt:.1f} ms")
    for K in (10, 20, 40):
        beta = 1.2**np.arange(K); alpha = 0.9*beta/K; t0 = time.perf_counter(); v, f, nf = tangent_fast(1.0, alpha, beta)
        print(f"K={K}: bound {v:.4g}, {nf} evaluations, {1000*(time.perf_counter()-t0):.0f} ms total")
