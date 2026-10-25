"""Independent check of the exact moments used by the tangent bound (second and third moments of S) and of M3 = E[lam a a'] against a long simulation,
   plus the tangent bound for smaller eps (vs the previous closed-form chain)."""
import numpy as np, random, math, sys
from p2_tangent_bound import moments, matrices, best

def simulate_moments(mu, alpha, beta, nev, seed):
    rng = random.Random(seed); K = len(beta); S = np.zeros(K); t = 0.0
    acc3 = np.zeros((K, K, K)); acc2 = np.zeros((K, K)); acc1 = np.zeros(K); T = 0.0
    for _ in range(int(3000)):  # burn-in
        lam0 = mu + alpha@S; w = 0.0
        while True:
            w += rng.expovariate(lam0)
            if rng.random()*lam0 <= mu + alpha@(S*np.exp(-beta*w)): break
        S = S*np.exp(-beta*w) + 1.0
    xs = np.linspace(0, 1, 9)
    for _ in range(nev):
        lam0 = mu + alpha@S; w = 0.0
        while True:
            w += rng.expovariate(lam0)
            if rng.random()*lam0 <= mu + alpha@(S*np.exp(-beta*w)): break
        # time average of S_i S_j S_l over the gap (Simpson on exact exponentials)
        n = 33; u = np.linspace(0, w, n); wts = np.ones(n); wts[1:-1:2] = 4; wts[2:-1:2] = 2; wts *= (w/(n-1))/3
        Su = S[None, :]*np.exp(-np.outer(u, beta))
        acc1 += wts@Su; acc2 += np.einsum('t,ti,tj->ij', wts, Su, Su); acc3 += np.einsum('t,ti,tj,tl->ijl', wts, Su, Su, Su); T += w
        S = S*np.exp(-beta*w) + 1.0
    return acc1/T, acc2/T, acc3/T

if __name__ == "__main__":
    mu = 1.0; beta = np.array([1.0, 2.0]); eps = 0.3; alpha = (1-eps)*beta/2
    m, Q, T3, El = moments(mu, alpha, beta)
    m_s, Q_s, T_s = simulate_moments(mu, alpha, beta, 400000, 5)
    print("max rel diff  m:", np.max(np.abs(m_s/m-1)), " Q:", np.max(np.abs(Q_s/Q-1)), " T3:", np.max(np.abs(T_s/T3-1)), "(simulation, 4e5 events; expect ~1e-2)")
    print("\nTangent bound for smaller eps (mu=1, K, r):")
    for K, r in ((1, 2), (2, 2), (3, 2)):
        beta = float(r)**np.arange(K)
        for eps in (0.3, 0.1, 0.03, 0.01):
            alpha = (1-eps)*beta/K; b, f = best(1.0, alpha, beta, eps)
            print(f"K={K} r={r} eps={eps}: tangent lower bound on lmin(I) = {b:.5f}   (bound/eps = {b/eps:.3f}),  best lam0/E lam = {f:.2f}", flush=True)
