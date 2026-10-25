"""Is the MLE of mu (alpha known, K=1, cold start, eps=c/T) mixed normal? Observed information J=int_0^T dt/lam_t is random.
Check: CV of J; z_obs=(mu_hat-mu)*sqrt(J) vs z_fix=(mu_hat-mu)*sqrt(E J)."""
import numpy as np, math, random
from scipy.optimize import brentq
from scipy.stats import kurtosis
from p10b_transient import grids, info_matrix

def path(mu, beta, alpha, T, rng):
    S = 0.0; t = 0.0; Ss = []; J = 0.0
    while True:
        lam = mu + alpha*S; w = 0.0
        while True:
            w += rng.expovariate(lam)
            if t + w > T: w2 = T - t; break
            if rng.random()*lam < mu + alpha*S*math.exp(-beta*w): w2 = w; break
        Se = S*math.exp(-beta*w2)
        J += (w2 + math.log((mu+alpha*Se)/(mu+alpha*S))/beta)/mu           # int dt/lam over the piece
        if t + w > T: return np.array(Ss), J
        t += w; S = Se + 1.0; Ss.append(S - 1.0)

def run(mu, T, c, R, seed):
    beta = 1.0; eps = c/T; alpha = (1-eps)*beta; rng = random.Random(seed)
    g = grids(mu, beta, eps, N=80000, lo=1e-14); M = info_matrix(mu, beta, eps, T, g)[1]
    mh, Js = [], []
    for _ in range(R):
        Ss, J = path(mu, beta, alpha, T, rng)
        f = lambda m: (1/(m + alpha*Ss)).sum() - T
        mh.append(brentq(f, 1e-6, 100.0) if f(1e-6) > 0 else 1e-6); Js.append(J)
    mh = np.array(mh); Js = np.array(Js); e = mh - mu
    zo = e*np.sqrt(Js); zf = e*math.sqrt(M)
    print(f" mu={mu} T={T} eps={eps:.3f} R={R}: E J (MC)={Js.mean():.2f}  exact M_T={M:.2f}  CV(J)={Js.std()/Js.mean():.2f}")
    print(f"    z_obs=(mu_hat-mu)sqrt(J): mean {zo.mean():+.2f} sd {zo.std():.2f} excess kurtosis {kurtosis(zo):+.2f}   |  z_fix=(mu_hat-mu)sqrt(M_T): mean {zf.mean():+.2f} sd {zf.std():.2f} excess kurtosis {kurtosis(zf):+.2f}", flush=True)

for mu in (1.0, 0.3): run(mu, 100, 2.0, 400, seed=int(10*mu)+1)
