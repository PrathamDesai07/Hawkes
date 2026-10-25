"""Problem 6, K>=2: rigorous lower bounds on every eigenvalue of I(theta*) from the tangent minorant (double precision, numpy),
and hence on det I and on (1/2) ln det I, as functions of eps.  Reuses code/p2_tangent_bound.py (moments); mpmath is not available in this Python,
so eigenvalues are in double precision (cancellation limits the reliable range to eps >= ~0.003)."""
import sys, numpy as np
sys.modules.setdefault('mpmath', type(sys)('mpmath'))   # p2_tangent_bound imports mpmath at module level; only moments/matrices are used here
import types
mp = sys.modules['mpmath']; mp.mp = types.SimpleNamespace(dps=60)
from p2_tangent_bound import matrices

def lower_bounds(mu, beta, eps):
    beta = np.asarray(beta, float); K = len(beta); alpha = (1-eps)*beta/K
    M, M3, El = matrices(mu, alpha, beta); best = np.full(K+1, -np.inf)
    for f in np.concatenate([np.linspace(0.3, 3.0, 55), np.linspace(3.2, 12.0, 23)]):
        l0 = f*El; N = 2*M/l0 - M3/l0**2; ev = np.sort(np.linalg.eigvalsh((N+N.T)/2)); best = np.maximum(best, ev)
    return best, El

if __name__ == "__main__":
    cases = [("K=1 beta=1 mu=1", [1.0], 1.0), ("K=2 beta=1,2 mu=1.5", [1.0, 2.0], 1.5), ("K=3 beta=1,2,4 mu=2", [1.0, 2.0, 4.0], 2.0)]
    for name, beta, mu in cases:
        print(name)
        print(f"  {'eps':>6s} {'lower bounds on lam_j(I)':>46s} {'ln det lower':>13s} {'eps*top':>8s} {'bottom/eps':>10s}")
        for eps in (0.3, 0.1, 0.03, 0.01):
            lb, El = lower_bounds(mu, beta, eps)
            print(f"  {eps:6.3f} {str([float(f'{x:.4g}') for x in lb]):>46s} {np.log(np.prod(np.maximum(lb,1e-300))):13.3f} {lb[-1]*eps:8.3f} {lb[0]/eps:10.4f}")
