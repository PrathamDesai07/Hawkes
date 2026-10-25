"""Problem 19: validity of the Ogata-thinning sampler used in p19_sim.cs (same algorithm, pure Python, seeded).
Time-rescaling theorem: with the exact compensator Lambda(t_{i-1}, t_i) = int lam du, the increments E_i must be i.i.d. Exp(1).
Tests: Kolmogorov-Smirnov distance to Exp(1), lag-1 autocorrelation, mean and variance; mean intensity vs mu/eps; E[1/lam] vs the Theorem A value."""
import math, random, numpy as np
from scipy import stats
def run(mu, alpha, beta, nev, seed):
    rng = random.Random(seed); s = 0.0; E = []; T = 0.0; invlam_int = 0.0
    for _ in range(nev):
        lam0 = mu + alpha*s; w = 0.0
        while True:                                               # Ogata thinning: lam(u) is non-increasing between events, lam0 dominates
            w += rng.expovariate(lam0)
            if rng.random()*lam0 <= mu + alpha*s*math.exp(-beta*w): break
        e = math.exp(-beta*w); E.append(mu*w + alpha*s*(1-e)/beta)       # exact compensator over the gap
        invlam_int += (w + math.log((mu + alpha*s*e)/(mu + alpha*s))/beta)/mu    # int du/lam over the gap
        T += w; s = s*e + 1.0
    return np.array(E), T, invlam_int
print(f"{'case':26s} {'KS D':>7s} {'p':>6s} {'mean E':>7s} {'var E':>7s} {'rho_1':>7s} {'lam_bar est':>11s} {'mu/eps':>7s} {'E[1/lam] est':>12s}")
for (mu, beta, eps) in ((1.0, 1.0, 0.3), (1.0, 1.0, 0.1), (1.0, 1.0, 0.03), (2.0, 1.0, 0.1), (0.5, 1.0, 0.1)):
    alpha = (1-eps)*beta; nev = 300000
    E, T, inv = run(mu, alpha, beta, nev, 7)
    ks = stats.kstest(E[1000:], 'expon'); rho = np.corrcoef(E[1001:], E[1000:-1])[0, 1]
    print(f"mu={mu} eps={eps:<5}          {ks.statistic:7.4f} {ks.pvalue:6.3f} {E.mean():7.3f} {E.var():7.3f} {rho:+7.3f} {nev/T:11.3f} {mu/eps:7.3f} {inv/T:12.4f}")
