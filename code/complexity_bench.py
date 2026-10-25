"""Complexity improvements, implemented and measured:
 (1) third/second moments of S by Kronecker-sum diagonalisation, O(K^4) time and O(K^3) space, instead of the dense K^3 x K^3 solve (O(K^9) time, O(K^6) space);
 (2) projection onto {box, budget}: safeguarded Newton on the budget multiplier (breakpoint-aware), a handful of inner box-QP solves instead of ~55 bisection steps;
 (3) online Newton step: maintain F^{-1} by q rank-one Sherman-Morrison updates, O(q d^2) per event, instead of a dense solve, O(d^3), with an accuracy check."""
import numpy as np, itertools, time, math
from scipy.optimize import lsq_linear

# ---------------- (1) moments ----------------
def moments_dense(mu, alpha, beta):
    from p2_tangent_bound import moments
    return moments(mu, alpha, beta)

def moments_fast(mu, alpha, beta):
    K = len(beta); A = np.diag(beta) - np.outer(np.ones(K), alpha)          # = -(drift matrix)
    lam_, V = np.linalg.eig(A); lam_ = lam_.real; V = V.real; Vi = np.linalg.inv(V)
    m = mu/(beta*(1 - (alpha/beta).sum()))                                    # E S_k = mu/(eps beta_k)
    El = mu + alpha@m
    # Q: A Q + Q A^T = R,  R_ij = mu(m_i+m_j) + El
    R2 = mu*(m[:, None] + m[None, :]) + El
    Qt = (Vi@R2@Vi.T)/(lam_[:, None] + lam_[None, :]); Q = V@Qt@V.T
    # T: (A+A+A) T = R3,  R3_ijl = sum_pairs mu Q_pq + sum_p (mu m_p + alpha.Q[:,p]) + El
    P = mu*m + alpha@Q                                                        # vector, P_p = mu m_p + sum_k alpha_k Q_kp
    R3 = (mu*Q[:, :, None] + mu*Q[:, None, :] + mu*Q[None, :, :]) * 0         # placeholder shape
    R3 = mu*(Q[:, :, None] + Q[:, None, :] + Q[None, :, :]) + (P[:, None, None] + P[None, :, None] + P[None, None, :]) + El
    # transform each mode by Vi
    X = np.einsum('ai,ijl->ajl', Vi, R3); X = np.einsum('bj,ajl->abl', Vi, X); X = np.einsum('cl,abl->abc', Vi, X)
    X = X/(lam_[:, None, None] + lam_[None, :, None] + lam_[None, None, :])
    T = np.einsum('ia,abc->ibc', V, X); T = np.einsum('jb,ibc->ijc', V, T); T = np.einsum('lc,ijc->ijl', V, T)
    return m, Q, T, El

print("(1) second/third moments of S:  dense K^3 x K^3 solve vs Kronecker-sum diagonalisation")
for K in (2, 3, 4, 6):
    beta = 1.5**np.arange(K); eps = 0.1; alpha = (1 - eps)*beta/K; mu = 1.0
    t0 = time.perf_counter(); m1, Q1, T1, E1 = moments_dense(mu, alpha, beta); td = time.perf_counter() - t0
    t0 = time.perf_counter(); m2, Q2, T2, E2 = moments_fast(mu, alpha, beta); tf = time.perf_counter() - t0
    err = max(np.abs(T1/T2 - 1).max(), np.abs(Q1/Q2 - 1).max())
    print(f"   K={K}: dense {1000*td:9.2f} ms, fast {1000*tf:7.2f} ms (x{td/tf:6.1f}); max rel diff {err:.1e}")
for K in (10, 20, 40):
    beta = 1.2**np.arange(K); alpha = 0.9*beta/K; t0 = time.perf_counter(); moments_fast(1.0, alpha, beta); print(f"   K={K}: fast only {1000*(time.perf_counter()-t0):9.1f} ms (dense would need a {K**3}x{K**3} matrix = {K**6*8/1e9:.2f} GB)")

# ---------------- (2) projection ----------------
def boxqp(F, q, lo, hi):
    L = np.linalg.cholesky(F); z = np.linalg.solve(F, q)
    sol = lsq_linear(L.T, L.T@z, bounds=(lo, hi), method='bvls', tol=1e-15, lsq_solver='exact', max_iter=1000)
    return sol.x

def project_bisect(F, y, lo, hi, c, r):
    calls = [0]
    def bq(q): calls[0] += 1; return boxqp(F, q, lo, hi)
    if np.all(y >= lo) and np.all(y <= hi) and c@y <= r: return y, 0
    Fy = F@y; x0 = bq(Fy)
    if c@x0 <= r + 1e-14: return x0, calls[0]
    tlo, thi = 0.0, 1.0
    while c@bq(Fy - thi*c) > r: thi *= 2.0
    for _ in range(80):
        tm = 0.5*(tlo + thi); xm = bq(Fy - tm*c)
        if c@xm > r: tlo = tm
        else: thi = tm
    return bq(Fy - thi*c), calls[0]

def project_newton(F, y, lo, hi, c, r, tol=1e-13):
    calls = [0]
    def bq(q): calls[0] += 1; return boxqp(F, q, lo, hi)
    if np.all(y >= lo) and np.all(y <= hi) and c@y <= r: return y, 0
    Fy = F@y; x = bq(Fy)
    if c@x <= r + 1e-14: return x, calls[0]
    tlo, thi, tau = 0.0, None, 0.0
    for _ in range(60):
        g = c@x - r                                                       # > 0 at tau (left of the root)
        free = (x > lo + 1e-12) & (x < hi - 1e-12)
        if free.any():
            slope = c[free]@np.linalg.solve(F[np.ix_(free, free)], c[free])  # -dg/dtau on the current piece
        else: slope = 0.0
        tlo = tau if g > 0 else tlo
        if g > 0: tlo = tau
        else: thi = tau if thi is None else min(thi, tau)
        step = tau + g/slope if slope > 1e-300 else (tau*2 + 1.0)
        if thi is not None and not (tlo < step < thi): step = 0.5*(tlo + thi)
        tau = step; x = bq(Fy - tau*c)
        if abs(c@x - r) < tol: break
    return x, calls[0]

print("\n(2) dense-metric projection: BVLS calls per projection (budget active)")
rng = np.random.default_rng(3); ratios = []; errs = []
for K in (2, 5, 10, 20):
    d = K + 1; cb = []; cn = []; tb = tn = 0.0
    for _ in range(30):
        beta = np.geomspace(1, 8, K); Am = rng.normal(size=(d, d)); F = Am@Am.T + np.eye(d)
        lo = np.concatenate([[0.1], np.zeros(K)]); hi = np.concatenate([[5.0], np.full(K, 3.0)]); c = np.concatenate([[0.0], 1.0/beta]); r = 0.9
        y = np.concatenate([[3.0], rng.uniform(0.5, 3.0, K)])
        t0 = time.perf_counter(); xb, nb = project_bisect(F, y, lo, hi, c, r); tb += time.perf_counter() - t0
        t0 = time.perf_counter(); xn, nn = project_newton(F, y, lo, hi, c, r); tn += time.perf_counter() - t0
        cb.append(nb); cn.append(nn); errs.append(np.sqrt((xb - xn)@F@(xb - xn)))
    print(f"   K={K:2d}: bisection {np.mean(cb):5.1f} inner solves ({1000*tb/30:7.1f} ms), safeguarded Newton {np.mean(cn):4.1f} ({1000*tn/30:6.1f} ms): speed-up x{tb/tn:5.1f}")
print(f"   max F-norm difference between the two projections: {max(errs):.2e}")

# ---------------- (3) inverse maintenance ----------------
print("\n(3) per-event Newton direction: dense solve vs Sherman-Morrison maintenance of F^{-1} (q = 12 rank-one terms per event)")
q = 12
for d in (3, 11, 51, 101):
    rngb = np.random.default_rng(d); F = np.eye(d)*5.0; Finv = np.linalg.inv(F); nev = 2000
    phis = [np.abs(rngb.normal(size=(q, d))) for _ in range(50)]; ws = rngb.uniform(0.01, 0.1, q)
    t0 = time.perf_counter(); Fd = F.copy()
    for e in range(nev):
        P = phis[e % 50]; Fd = Fd + (P.T*ws)@P; g = rngb.normal(size=d); dirn = np.linalg.solve(Fd, g)
    t_dense = time.perf_counter() - t0
    t0 = time.perf_counter(); Fi = Finv.copy(); Fm = F.copy()
    for e in range(nev):
        P = phis[e % 50]
        for j in range(q):
            u = Fi@P[j]; Fi -= (ws[j]/(1 + ws[j]*(P[j]@u)))*np.outer(u, u)
        Fm = Fm + (P.T*ws)@P; g = rngb.normal(size=d); dirn2 = Fi@g
    t_sm = time.perf_counter() - t0
    relerr = np.linalg.norm(Fi@Fm - np.eye(d))/math.sqrt(d)
    print(f"   d={d:3d}: dense {1e6*t_dense/nev:8.1f} us/event, Sherman-Morrison {1e6*t_sm/nev:8.1f} us/event (x{t_dense/t_sm:4.2f}); ||F^-1 F - I|| after {nev} events = {relerr:.1e}")
