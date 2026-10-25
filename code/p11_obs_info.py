"""Problem 11 (K = 1, alpha known, eps = c/T, beta = 1): law of the OBSERVED information I_T = int dt/lam_t and of the MLE error.
 Start: lam_0 ~ Gamma(k, scale a/eps) (the Gamma approximation of the stationary law), S_0 = (lam_0 - mu)/alpha clipped at 0.
 Reports, for each T: quantiles of I_T, P(I_T <= x) for fixed x, quantiles of |mu_hat - mu|, and the studentised z = (mu_hat - mu) sqrt(I_T) (LAMN check)."""
import math, random, sys, time
import numpy as np

def one_path(mu, eps, T, seed, beta=1.0, c_target=None):
    rng = random.Random(seed); al = (1 - eps)*beta; a = al*al/(2*beta); k = 2*beta*mu/(al*al)
    lam0 = rng.gammavariate(k, a/eps); s = max(lam0 - mu, 0.0)/al
    t = 0.0; info = 0.0; S_ev = []
    exp = math.exp; log = math.log
    while True:
        lam_s = mu + al*s; w = 0.0; ok = False
        while True:
            w += rng.expovariate(lam_s)
            if t + w >= T: break
            if rng.random()*lam_s <= mu + al*s*exp(-beta*w): ok = True; break
        if not ok:   # no event before T: integrate to T
            ww = T - t; info += (ww + log((mu + al*s*exp(-beta*ww))/(mu + al*s))/beta)/mu; break
        info += (w + log((mu + al*s*exp(-beta*w))/(mu + al*s))/beta)/mu
        s_minus = s*exp(-beta*w); S_ev.append(s_minus); s = s_minus + 1.0; t += w
    S_ev = np.array(S_ev)
    # MLE: sum 1/(m + al S_i) = T   (score decreasing in m); bisection
    lo, hi = 1e-9, 50.0*mu
    for _ in range(80):
        mid = 0.5*(lo + hi)
        if np.sum(1.0/(mid + al*S_ev)) > T: lo = mid
        else: hi = mid
    return info, 0.5*(lo + hi), len(S_ev)

if __name__ == "__main__":
    mu = float(sys.argv[1]); c = float(sys.argv[2]); Ts = [float(x) for x in sys.argv[3].split(',')]; R = [int(x) for x in sys.argv[4].split(',')]
    for T, Rr in zip(Ts, R):
        eps = c/T; t0 = time.time(); res = [one_path(mu, eps, T, 1000 + i) for i in range(Rr)]
        I = np.array([r[0] for r in res]); mh = np.array([r[1] for r in res]); ne = np.mean([r[2] for r in res])
        err = np.abs(mh - mu); z = (mh - mu)*np.sqrt(I)
        q = lambda v, p: np.quantile(v, p)
        print(f"mu={mu} (k={2*mu:.1f}) c={c} T={T:g} eps={eps:.4f} ({Rr} paths, {ne:.0f} events/path, {time.time()-t0:.0f}s)", flush=True)
        print(f"   I_T quantiles 5/25/50/75/95%: {q(I,.05):.2f} {q(I,.25):.2f} {q(I,.5):.2f} {q(I,.75):.2f} {q(I,.95):.2f};  mean {I.mean():.2f};  P(I_T <= 5) = {np.mean(I <= 5):.2f}, P(I_T <= 10) = {np.mean(I <= 10):.2f}", flush=True)
        print(f"   |mu_hat - mu| quantiles 25/50/75/90%: {q(err,.25):.3f} {q(err,.5):.3f} {q(err,.75):.3f} {q(err,.9):.3f};  P(|err| > 0.5 mu) = {np.mean(err > 0.5*mu):.2f};  z = (mu_hat-mu) sqrt(I_T): mean {z.mean():+.2f}, sd {z.std():.2f}, P(|z|>2) = {np.mean(np.abs(z) > 2):.2f}", flush=True)
