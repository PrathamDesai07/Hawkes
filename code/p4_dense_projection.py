"""Exact projection in a DENSE metric F onto Theta = {lo <= x <= hi, c'x <= r}   (x = (mu, alpha_1..alpha_K), c = (0, 1/beta), r = 1 - eps).
 min 1/2 (x-y)'F(x-y)  s.t. box and one budget constraint.   Dual in tau >= 0 for the budget: x(tau) = argmin_box 1/2 (x-y)'F(x-y) + tau c'x,
 g(tau) = c'x(tau) is continuous and non-increasing (convex QP), so tau is found by bracketing + bisection; each inner problem is a box QP solved by a
 finite active-set (Bertsekas projected-Newton / free-set) method.  Cost O(K^3) per inner solve (dense F), about 40-60 inner solves; typically
 the budget is inactive and one box-QP solve (or none, when y is feasible) suffices.
 Test: (i) vs brute-force enumeration of active sets (d <= 5) and vs SLSQP; (ii) Pythagorean/non-expansive property ||Pi(y)-x||_F <= ||y-x||_F for all feasible x."""
import numpy as np, itertools, time
from scipy.optimize import minimize

def boxqp(F, q, lo, hi, maxit=None):
    """minimise 1/2 x'Fx - q'x on lo <= x <= hi (F SPD) = bounded least squares ||L'(x - F^-1 q)||^2, F = L L'.
    Solved by BVLS (Stark-Parker), a finite exact active-set method."""
    from scipy.optimize import lsq_linear
    L = np.linalg.cholesky(F); z = np.linalg.solve(F, q)
    sol = lsq_linear(L.T, L.T@z, bounds=(lo, hi), method='bvls', tol=1e-15, lsq_solver='exact', max_iter=1000)
    return sol.x

def project_dense(F, y, lo, hi, c, r):
    y = np.asarray(y, float)
    if np.all(y >= lo) and np.all(y <= hi) and c@y <= r: return y
    Fy = F@y
    x0 = boxqp(F, Fy, lo, hi)
    if c@x0 <= r + 1e-14: return x0
    tlo, thi = 0.0, 1.0
    while c@boxqp(F, Fy - thi*c, lo, hi) > r: thi *= 2.0
    for _ in range(80):
        tm = 0.5*(tlo + thi); xm = boxqp(F, Fy - tm*c, lo, hi)
        if c@xm > r: tlo = tm
        else: thi = tm
    return boxqp(F, Fy - thi*c, lo, hi)

def brute(F, y, lo, hi, c, r):
    """reference: enumerate for each coordinate {free, at lo, at hi} and budget {inactive, active}; keep the best feasible KKT point."""
    d = len(y); best = None
    for pat in itertools.product((0, 1, 2), repeat=d):
        for act in (0, 1):
            x = np.zeros(d); fixed = np.array([p != 0 for p in pat]); free = ~fixed
            x[fixed] = np.where(np.array(pat)[fixed] == 1, lo[fixed], hi[fixed])
            if act == 0:
                if free.any():
                    x[free] = np.linalg.solve(F[np.ix_(free, free)], (F@y)[free] - F[np.ix_(free, fixed)]@x[fixed])
            else:
                # free vars with budget multiplier: solve [F_ff c_f; c_f' 0][x_f; tau] = [(Fy)_f - F_fx x_x ; r - c_x x_x]
                nf = free.sum()
                if nf == 0:
                    if abs(c@x - r) > 1e-9: continue
                else:
                    Mx = np.zeros((nf+1, nf+1)); Mx[:nf, :nf] = F[np.ix_(free, free)]; Mx[:nf, nf] = c[free]; Mx[nf, :nf] = c[free]
                    rhs = np.concatenate([(F@y)[free] - F[np.ix_(free, fixed)]@x[fixed], [r - c[fixed]@x[fixed]]])
                    try: sol = np.linalg.solve(Mx, rhs)
                    except np.linalg.LinAlgError: continue
                    x[free] = sol[:nf]
            if np.all(x >= lo - 1e-12) and np.all(x <= hi + 1e-12) and c@x <= r + 1e-12:
                val = 0.5*(x-y)@F@(x-y)
                if best is None or val < best[0]: best = (val, x.copy())
    return best[1]

if __name__ == "__main__":
    rng = np.random.default_rng(7); worst_b = worst_s = 0.0; worst_ne = 0.0; cnt_active = 0; cnt = 0
    for trial in range(400):
        K = int(rng.integers(1, 5)); d = K + 1
        beta = np.sort(rng.uniform(0.5, 6.0, K)); eps = rng.uniform(0.02, 0.5)
        A = rng.normal(size=(d, d)); F = A@A.T + np.diag(rng.uniform(0.01, 2.0, d))*rng.choice([1e-2, 1, 1e2])
        lo = np.concatenate([[rng.uniform(0.05, 0.5)], np.zeros(K) + rng.uniform(0, 0.05, K)]); hi = np.concatenate([[rng.uniform(1.5, 5)], rng.uniform(0.3, 3, K)])
        c = np.concatenate([[0.0], 1.0/beta]); r = 1 - eps
        if lo[1:]@c[1:] > r: continue
        y = np.concatenate([[rng.uniform(-1, 6)], rng.uniform(-1, 4, K)])*rng.choice([1, 3])
        x = project_dense(F, y, lo, hi, c, r); xb = brute(F, y, lo, hi, c, r)
        res = minimize(lambda z: (z-y)@F@(z-y), np.clip(y, lo, hi), jac=lambda z: 2*F@(z-y), bounds=list(zip(lo, hi)),
                       constraints=[{'type': 'ineq', 'fun': lambda z: r - c@z, 'jac': lambda z: -c}], method='SLSQP', options={'ftol': 1e-14, 'maxiter': 500})
        obj = lambda z: 0.5*(z-y)@F@(z-y)
        worst_b = max(worst_b, (obj(x) - obj(xb))/max(1e-12, abs(obj(xb))+1e-12)); worst_s = max(worst_s, (obj(x) - obj(res.x))/max(1e-12, abs(obj(res.x))+1e-12))
        assert np.all(x >= lo - 1e-9) and np.all(x <= hi + 1e-9) and c@x <= r + 1e-9, "infeasible"
        cnt += 1; cnt_active += not (np.all(y >= lo) and np.all(y <= hi) and c@y <= r)
        # non-expansiveness / Pythagorean: for random feasible z, ||x - z||_F <= ||y - z||_F
        for _ in range(5):
            z = lo + rng.uniform(0, 1, d)*(hi - lo)
            if c@z > r: z = lo.copy()
            if np.all(z >= lo) and c@z <= r:
                worst_ne = max(worst_ne, np.sqrt((x-z)@F@(x-z)) - np.sqrt((y-z)@F@(y-z)))
    print(f"{cnt} instances ({cnt_active} with an active projection): worst relative objective excess vs brute-force enumeration {worst_b:.2e}, vs SLSQP {worst_s:.2e}")
    print(f"max over (instance, feasible z) of ||Pi(y)-z||_F - ||y-z||_F = {worst_ne:.2e}  (<= 0 up to rounding means non-expansive)")
    # timing, dense F
    for K in (2, 5, 10, 20, 50):
        d = K + 1; beta = np.geomspace(1, 8, K); A = rng.normal(size=(d, d)); F = A@A.T + np.eye(d)
        lo = np.concatenate([[0.1], np.zeros(K)]); hi = np.concatenate([[5.0], np.full(K, 3.0)]); c = np.concatenate([[0.0], 1.0/beta]); r = 0.9
        y = np.concatenate([[3.0], np.full(K, 2.0)]); t0 = time.time(); reps = 20
        for _ in range(reps): project_dense(F, y, lo, hi, c, r)
        print(f"K={K:3d}: {1000*(time.time()-t0)/reps:8.2f} ms per projection (budget active, dense F)")
