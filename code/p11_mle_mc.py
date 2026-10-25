"""Problem 11: does the MLE attain 1/I_T? K=1, cold start, T=c/eps. Pure-python thinning (cumulative candidate times).
alpha-known: root of sum 1/(mu+alpha S_i) = T.  Joint (mu,alpha): L-BFGS-B on exact log-likelihood.
Compare Var(mu_hat) with 1/M_T (alpha known) and [I_T^{-1}]_{mu mu} (joint), I_T from the exact transient formula."""
import numpy as np, math, random, sys
from scipy.optimize import brentq, minimize
from p10b_transient import grids, info_matrix

def simulate(mu, beta, alpha, T, rng):
    S = 0.0; t = 0.0; Ss = []; ts = []
    while True:
        lam = mu + alpha*S; w = 0.0
        while True:
            w += rng.expovariate(lam)
            if t + w > T: return np.array(Ss), np.array(ts)
            if rng.random()*lam < mu + alpha*S*math.exp(-beta*w): break
        t += w; S = S*math.exp(-beta*w) + 1.0; Ss.append(S - 1.0); ts.append(t)   # S just before the event

def fit(Ss, ts, T, beta, mu0, a0, alpha_known):
    J = ((1-np.exp(-beta*(T-ts)))/beta).sum()                   # int_0^T S dt
    if alpha_known:
        f = lambda m: (1/(m + a0*Ss)).sum() - T
        if f(1e-6) < 0: return 1e-6, a0
        return brentq(f, 1e-6, 100.0), a0
    def nll(p):
        m, a = p; l = m + a*Ss
        return -(np.log(l).sum() - m*T - a*J), np.array([-(1/l).sum() + T, -(Ss/l).sum() + J])
    r = minimize(nll, [mu0, a0], jac=True, method="L-BFGS-B", bounds=[(1e-4, 30), (1e-4, 1.2*beta)])
    return r.x[0], r.x[1]

def run(mu, beta, T, c, R, seed):
    eps = c/T; alpha = (1-eps)*beta; rng = random.Random(seed)
    g = grids(mu, beta, eps, N=80000, lo=1e-14); I, M, L = info_matrix(mu, beta, eps, T, g)
    th_known, th_joint = [], []
    for _ in range(R):
        Ss, ts = simulate(mu, beta, alpha, T, rng)
        th_known.append(fit(Ss, ts, T, beta, mu, alpha, True)[0]); th_joint.append(fit(Ss, ts, T, beta, mu, alpha, False)[0])
    th_known = np.array(th_known); th_joint = np.array(th_joint)
    rob = lambda x: (np.percentile(x, 75)-np.percentile(x, 25))/1.349
    print(f" mu={mu} T={T} eps={eps:.4f} R={R}")
    print(f"   alpha known: mean mu_hat {th_known.mean():.3f}  sd {th_known.std():.3f} (robust {rob(th_known):.3f})   theory 1/sqrt(M_T) = {1/math.sqrt(M):.3f}")
    print(f"   joint      : mean mu_hat {th_joint.mean():.3f}  sd {th_joint.std():.3f} (robust {rob(th_joint):.3f})   theory sqrt(I^-1_mm) = {math.sqrt(np.linalg.inv(I)[0,0]):.3f}", flush=True)

if __name__ == "__main__":
    R = int(sys.argv[1]) if len(sys.argv) > 1 else 150
    for mu in (1.0, 0.3):
        for T in (100, 400): run(mu, 1.0, T, 2.0, R if T == 100 else R//2, seed=int(10*mu)+T)
