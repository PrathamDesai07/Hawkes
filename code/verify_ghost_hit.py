"""End-to-end check of the coupling step (any K).  A process started at time -L from state s equals  E (empty-start immigrants, rate mu)  plus 'ghost' immigrants with intensity
   g_s(tau) = sum_k alpha_k s_k exp(-beta_k tau), tau = t + L,  each with an independent cluster.  The records on [0,T] of two such processes (s ~ pi_theta, s' ~ pi_theta') coincide whenever
   neither ghost family has an event in [L, L+T] (in tau-coordinates).  We estimate  P(ghost family hits [0,T])  by direct simulation for s ~ pi_theta (stationary states from a long run)
   and compare with the proof's bound  C' mu e^{-gamma L}  (survival bound q(s)<=C' gamma e^{-gamma s}, E Z = mu/gamma)  and with the first-moment (expected hits) value.
   Also the exact identity  zeta' = -gamma zeta - r(|Y|)  with the analytic derivative."""
import numpy as np, math, random
from scipy.optimize import brentq
from scipy.integrate import solve_ivp
from verify_survival import setup, flow

def stationary_states(beta, alpha, mu, n, seed=0, spacing=None):
    rng = np.random.default_rng(seed); K = len(beta); S = np.zeros(K); t = 0.0; out = []; nxt = None
    g = brentq(lambda x: np.sum(alpha/(beta - x)) - 1.0, 0.0, beta.min()*(1 - 1e-12)); sp = spacing or 2/g; burn = 30/g; nxt = burn
    while len(out) < n:
        lam = mu + alpha@S; w = 0.0
        while True:
            w += rng.exponential(1/lam)
            if rng.random()*lam < mu + alpha@(S*np.exp(-beta*w)): break
        while nxt < t + w and len(out) < n: out.append(S*np.exp(-beta*(nxt - t))); nxt += sp
        S = S*np.exp(-beta*w) + 1.0; t += w
    return np.array(out)

def cluster_has_event_in(rng, beta, alpha, t0, lo, hi):
    """BFS cluster rooted at time t0 (root event included); True if any event time in [lo, hi]"""
    nn = float(np.sum(alpha/beta)); pk = np.cumsum((alpha/beta)/nn); stack = [t0]
    while stack:
        t = stack.pop()
        if lo <= t <= hi: return True
        if t > hi: continue          # descendants are later than t > hi: cannot hit
        k = rng.poisson(nn)
        for _ in range(k):
            j = min(int(np.searchsorted(pk, rng.random())), len(beta) - 1); stack.append(t + rng.exponential(1/beta[j]))
    return False

def ghost_hit(rng, beta, alpha, s, L, T):
    a = alpha*s; tot = L + T; w = a*(1 - np.exp(-beta*tot))/beta; N = rng.poisson(w.sum())
    for _ in range(N):
        k = rng.choice(len(beta), p=w/w.sum()); u = rng.random()      # truncated exponential on [0, tot] with rate beta_k
        tau = -math.log(1 - u*(1 - math.exp(-beta[k]*tot)))/beta[k]
        if cluster_has_event_in(rng, beta, alpha, tau, L, L + T): return True
    return False

if __name__ == "__main__":
    # exact derivative check of the zeta identity
    be, al, g = setup([1.0, 2.0], [0.5, 0.5], 0.05); S = np.linspace(0, 40, 400); Y = flow(be, al, S)
    Yp = -be*Y + al[None, :]*(1 - np.exp(-Y.sum(1)))[:, None]; zeta = (Y/(be - g)).sum(1); zp = (Yp/(be - g)).sum(1); absY = Y.sum(1)
    print("zeta' = -gamma zeta - r(|Y|), exact derivative: max abs diff =", np.max(np.abs(zp - (-g*zeta - (absY - 1 + np.exp(-absY))))))
    beta = np.array([1.0, 2.0]); eps = 0.05; mu = 1.5; w = np.array([.5, .5]); alpha = (1 - eps)*beta*w
    be, al, g = setup(beta, w, eps); T = 1.0/eps; states = stationary_states(beta, alpha, mu, 4000, seed=3)
    Zs = states@(alpha/(beta - g)); print(f"K=2 beta=(1,2), eps={eps}, mu={mu}, gamma={g:.4f}, T={T}: E Z = {Zs.mean():.2f} (theory mu/gamma = {mu/g:.2f}); {len(states)} stationary states")
    S_grid = np.linspace(0, 12/g, 4000); Y = flow(be, al, S_grid); qf = 1 - np.exp(-Y.sum(1)); rng = np.random.default_rng(5)
    print(f"{'K0':>3} {'L':>7} | {'P(ghost hits) MC':>16} | {'sum_k ghost-mass x q (first-moment-in-clusters)':>47} | proof-type bound C mu e^-gamma L (C=2.28 observed const)")
    for K0 in (1, 2, 3, 4):
        L = K0/g; hits = sum(ghost_hit(rng, be, al, s, L, T) for s in states[:3000]) / 3000
        # exact expected number of touching ghost immigrants:  int g_s(tau) q(L - tau) dtau  (+ window part), averaged over states
        taus = np.linspace(0, L + T, 3000); dt = taus[1] - taus[0]; exp_touch = 0.0
        for s in states[:600]:
            gs = (al*s)[None, :]*np.exp(-np.outer(taus, be)); gs = gs.sum(1)
            qq = np.where(taus < L, np.interp(L - taus, S_grid, qf), 1.0); exp_touch += np.sum(gs*qq)*dt
        exp_touch /= 600
        print(f"{K0:3d} {L:7.1f} | {hits:16.4f} | {exp_touch:47.4f} | {2.28*mu*math.exp(-K0):.4f}", flush=True)
