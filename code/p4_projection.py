"""Problem 4: projection onto Theta_eps = {mu >= mu_min, alpha >= 0, sum_k alpha_k/beta_k <= 1-eps}
in the diagonal metric D = diag(d_mu, d_1..d_K).
mu is separable: mu = max(y_mu, mu_min).  For alpha: min sum d_k (a_k-y_k)^2 s.t. a>=0, w.a <= r (w_k=1/beta_k).
KKT: a_k = max(0, y_k - tau*w_k/d_k), tau >= 0 chosen so that w.a = r when the clipped y violates the budget.
g(tau) = sum_k w_k max(0, y_k - tau*w_k/d_k) is piecewise linear, decreasing; breakpoints t_k = y_k d_k / w_k.
Sort breakpoints (O(K log K)), or use a median-of-breakpoints search (O(K) expected)."""
import numpy as np, time
from scipy.optimize import minimize

def project_diag(y_mu, y_a, d_mu, d_a, beta, mu_min, eps):
    w = 1.0/beta; r = 1.0 - eps
    mu = max(y_mu, mu_min)
    a = np.maximum(y_a, 0.0)
    if w @ a <= r: return mu, a
    pos = y_a > 0
    ys, ws, ds = y_a[pos], w[pos], d_a[pos]
    t = ys*ds/ws                              # tau at which coordinate k reaches zero
    o = np.argsort(t); t, ys, ws, ds = t[o], ys[o], ws[o], ds[o]
    # active set = coordinates with t_k > tau, i.e. a suffix j.. of the ascending order
    Sy = np.cumsum((ws*ys)[::-1])[::-1]       # sum_{k>=j} w_k y_k
    Sw = np.cumsum((ws*ws/ds)[::-1])[::-1]    # sum_{k>=j} w_k^2/d_k
    tau_j = (Sy - r)/Sw
    left = np.concatenate([[0.0], t[:-1]])    # tau must satisfy t_{j-1} <= tau <= t_j
    ok = (tau_j >= left - 1e-12) & (tau_j <= t + 1e-12)
    tau = tau_j[np.argmax(ok)]
    return mu, np.maximum(y_a - tau*w/d_a, 0.0)

def project_generic(y, P, beta, mu_min, eps):
    """reference: min (th-y)^T P (th-y) over Theta_eps, SLSQP (any PSD metric P)."""
    K = len(beta); w = 1.0/beta
    f = lambda th: (th-y)@P@(th-y); g = lambda th: 2*P@(th-y)
    cons = [{"type": "ineq", "fun": lambda th: (1-eps) - w@th[1:], "jac": lambda th: np.concatenate([[0.], -w])}]
    bnds = [(mu_min, None)] + [(0, None)]*K
    x0 = np.concatenate([[max(y[0], mu_min)], np.maximum(y[1:], 0)*0.0 + 1e-3])
    return minimize(f, x0, jac=g, bounds=bnds, constraints=cons, method="SLSQP", options={"ftol": 1e-15, "maxiter": 500}).x

if __name__ == "__main__":
    rng = np.random.default_rng(0)
    print("correctness vs SLSQP (diagonal metric); max |diff| over random instances")
    worst = 0
    for trial in range(300):
        K = int(rng.integers(1, 8)); beta = np.exp(rng.normal(0, 1, K)); eps = rng.uniform(0.01, 0.5); mu_min = rng.uniform(0, .5)
        y = np.concatenate([[rng.normal(0.5, 1)], rng.normal(0.4, 1, K)*beta])
        d = np.exp(rng.normal(0, 1.5, K+1))
        mu, a = project_diag(y[0], y[1:], d[0], d[1:], beta, mu_min, eps)
        ref = project_generic(y, np.diag(d), beta, mu_min, eps)
        worst = max(worst, np.abs(np.concatenate([[mu], a]) - ref).max())
        assert (1.0/beta)@a <= 1-eps+1e-10 and a.min() >= 0 and mu >= mu_min-1e-12
    print("  worst abs diff:", worst)
    print("timing (K, ms per projection, violating input)")
    for K in (10, 100, 1000, 10000, 100000):
        beta = np.exp(rng.normal(0, 1, K)); y = np.abs(rng.normal(1, 1, K))*beta; d = np.ones(K)
        t0 = time.perf_counter()
        for _ in range(5): project_diag(1.0, y, 1.0, d, beta, .1, .05)
        print(f"  K={K:6d}: {(time.perf_counter()-t0)/5*1e3:.3f} ms")
