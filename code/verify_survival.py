"""Survival probability q(s) = P(a cluster rooted at 0 has an event at time > s) for the exponential-dictionary Hawkes cluster.
 Flow:  Y_k' = -beta_k Y_k + alpha_k (1 - exp(-|Y|)),  Y_k(0) = alpha_k/beta_k,   q(s) = 1 - exp(-|Y(s)|)       [derived from the offspring Poisson structure]
 Slow coordinate  zeta = sum_k Y_k/(beta_k - gamma):  zeta' = -gamma zeta - r(|Y|),  r(x) = x - 1 + e^-x,  gamma: sum alpha_k/(beta_k-gamma) = 1.
 Checks: (1) flow vs Monte Carlo cluster simulation, (2) the identity for zeta', (3) Bernoulli bound zeta <= zeta0 e^{-gamma s}/(1+(c zeta0/gamma)(1-e^{-gamma s})),
         (4) q(s) <= C' gamma e^{-gamma s} for s >= 1/gamma with the explicit C', versus the first-moment bound (Sum alpha_k/gamma) e^{-gamma s}."""
import numpy as np, math, random
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

def setup(beta, w, eps):
    beta = np.array(beta, float); w = np.array(w, float); alpha = (1 - eps)*beta*w/w.sum()
    g = brentq(lambda x: np.sum(alpha/(beta - x)) - 1.0, 0.0, beta.min()*(1 - 1e-12)); return beta, alpha, g

def flow(beta, alpha, S):
    f = lambda s, Y: -beta*Y + alpha*(1 - math.exp(-Y.sum()))
    sol = solve_ivp(f, [0, S[-1]], alpha/beta, t_eval=S, rtol=1e-11, atol=1e-14, method="LSODA"); return sol.y.T

def mc_survival(beta, alpha, s_grid, n, seed=0):
    rng = random.Random(seed); nn = float(np.sum(alpha/beta)); pk = (alpha/beta)/nn; cnt = np.zeros(len(s_grid))
    for _ in range(n):
        stack = [0.0]; tmax = 0.0
        while stack:
            t0 = stack.pop()
            # number of offspring ~ Poisson(nn)
            L = nn; k = 0; p = math.exp(-L); u = rng.random(); cdf = p
            while u > cdf: k += 1; p *= L/k; cdf += p
            for _ in range(k):
                j = np.searchsorted(np.cumsum(pk), rng.random()); t1 = t0 + rng.expovariate(beta[min(j, len(beta)-1)]); stack.append(t1); tmax = max(tmax, t1)
        cnt += (tmax > np.array(s_grid))
    return cnt/n

if __name__ == "__main__":
    for name, beta, w in (("K=1 beta=1", [1.0], [1.0]), ("K=2 beta=(1,2)", [1.0, 2.0], [0.5, 0.5]), ("K=3 beta=(1,2,4)", [1.0, 2.0, 4.0], [0.4, 0.3, 0.3])):
        for eps in (0.2, 0.05):
            be, al, g = setup(beta, w, eps); b, B = be.min(), be.max(); om_t = 1/(be - g)
            S = np.linspace(0, 6/g, 3001); Y = flow(be, al, S); absY = Y.sum(1); q = 1 - np.exp(-absY); zeta = (Y*om_t).sum(1)
            # (2) zeta' identity by finite differences
            dz = np.gradient(zeta, S); rhs = -g*zeta - (absY - 1 + np.exp(-absY)); id_err = np.max(np.abs(dz - rhs))/np.max(np.abs(dz))
            # (3) Bernoulli bound; R0, c, C'
            R0 = (B - g)/(b - g); c_r = (R0 - 1 + math.exp(-R0))/R0**2; c = c_r*(b - g)**2; z0 = zeta[0]
            bern = z0*np.exp(-g*S)/(1 + (c*z0/g)*(1 - np.exp(-g*S))); ok3 = np.all(zeta <= bern*(1 + 1e-6))
            Cp = (B - g)/(c*(1 - math.exp(-1))); m = S >= 1/g
            ratio = q[m]/(g*np.exp(-g*S[m])); fm = (al.sum()/g)*np.exp(-g*S[m])
            print(f"{name}, eps={eps}: gamma={g:.4f}, q(0)={q[0]:.3f}; zeta' identity max rel err {id_err:.1e}; Bernoulli bound holds: {ok3}; "
                  f"q/(gamma e^-gamma s) for s>=1/gamma: max {ratio.max():.2f} (proof constant C' = {Cp:.1f}); first-moment bound / q at s=2/gamma: {(al.sum()/g*np.exp(-2))/q[np.searchsorted(S, 2/g)]:.0f}x")
    print("\nMonte Carlo check of the flow (cluster simulation, 60000 clusters):")
    for name, beta, w, eps in (("K=1 beta=1", [1.0], [1.0], 0.1), ("K=2 beta=(1,2)", [1.0, 2.0], [0.5, 0.5], 0.1)):
        be, al, g = setup(beta, w, eps); sg = [0.5, 2, 5, 10, 20, 40]
        qm = mc_survival(be, al, sg, 60000); Y = flow(be, al, np.array(sg, float)); qf = 1 - np.exp(-Y.sum(1))
        print(f"  {name}, eps={eps}: s={sg}\n      MC   {np.round(qm, 4).tolist()}\n      flow {np.round(qf, 4).tolist()}", flush=True)
