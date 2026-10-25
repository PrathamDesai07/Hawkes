"""Average realised gain matrix G_i = F_i^{-1} dF_i of the forgetting recursion (stationary truth, start at truth), alpha unknown.
If G were deterministic and equal to eta*I the error recursion would be a pure eta-filter in every direction.
Reports eigenvalues of mean(G)/eta and the average of the *product* structure E[G e e^T G-type] is not needed here."""
import math, random, numpy as np, sys
from p5_online_newton_k1 import fisher_inc
def run(mu_s, beta, eps, eta, n, seed, rho0=5.0):
    rng = random.Random(seed); al_s = (1-eps)*beta; mu, al = mu_s, al_s; lam_bar = mu_s/eps
    F = np.array([[rho0/lam_bar, 0.0], [0.0, rho0*(lam_bar/beta)**2/lam_bar]]); s = 0.0
    Gs = np.zeros((2, 2)); cnt = 0; trF = []
    for i in range(n):
        lam0 = mu_s + al_s*s; w = 0.0
        while True:
            w += rng.expovariate(lam0)
            if rng.random()*lam0 <= mu_s + al_s*s*math.exp(-beta*w): break
        e = math.exp(-beta*w); a1 = s*e; b1 = s*(1-e)/beta
        f0, f1, f2 = fisher_inc(mu, al, beta, s, w); dF = np.array([[f0, f1], [f1, f2]])
        F = (1-eta)*F + dF
        g = np.array([w - 1.0/(mu + al*a1), b1 - a1/(mu + al*a1)])
        step = np.linalg.solve(F, g); mu = max(mu - step[0], 0.1); al = min(max(al - step[1], 0.02), beta*(1-eps/2))
        if i > 0.2*n: Gs += np.linalg.solve(F, dF); cnt += 1
        s = a1 + 1.0
    return Gs/cnt
if __name__ == "__main__":
    for eps, eta in ((0.4, 0.004), (0.2, 0.001)):
        n = int(40000*(0.4/eps)); Gm = np.mean([run(1.0, 1.0, eps, eta, n, 50+k)/eta for k in range(6)], axis=0)
        ev = np.linalg.eigvals(Gm)
        print(f"eps={eps}, eta={eta}: mean(G)/eta =\n{np.round(Gm,3)}\n   eigenvalues {np.round(np.sort(ev.real),3)}", flush=True)
