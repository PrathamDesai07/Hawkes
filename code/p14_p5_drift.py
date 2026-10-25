"""Problem 14: Problem 5's recursive-Newton estimator (copied from p5_online_newton_k1.run: same update, same clipping, same prior) on DRIFTING data,
with an optional forgetting factor eta:  F <- (1-eta) F + dF.  eta = 0 is Problem 5's algorithm verbatim.
Truth: mu_t = mu_s (1 + A sin(2 pi t / P)); alpha, beta fixed; both mu and alpha are estimated (alpha UNKNOWN).
Gap sampling: thinning with bound lam0 = mu_max + alpha_s*s; when A = 0, mu_max = mu_s and the RNG call sequence equals Problem 5's."""
import math, random, numpy as np
from p5_online_newton_k1 import fisher_inc, run as p5_run

def run_drift(mu_s, beta, eps, seed, th0, T=None, n_max=None, rho0=5.0, eta=0.0, A=0.0, P=1200.0,
              mu_min=0.1, al_min=0.02, eps_c=None):
    rng = random.Random(seed); al_s = (1-eps)*beta; eps_c = eps/2 if eps_c is None else eps_c
    mu, al = th0; lam_bar = mu_s/eps; mu_max = mu_s*(1+abs(A))
    F = [[rho0/lam_bar, 0.0], [0.0, rho0*(lam_bar/beta)**2/lam_bar]]
    J = [0.0, 0.0, 0.0]; s = 0.0; Tt = 0.0; al_max = beta*(1-eps_c); ev = []; i = 0; e_n = None
    while True:
        if (n_max is not None and i >= n_max) or (T is not None and Tt >= T): break
        i += 1
        lam0 = mu_max + al_s*s; w = 0.0
        while True:
            w += rng.expovariate(lam0)
            mt = mu_s*(1 + A*math.sin(2*math.pi*(Tt+w)/P)) if A != 0 else mu_s
            if rng.random()*lam0 <= mt + al_s*s*math.exp(-beta*w): break
        e = math.exp(-beta*w); a1 = s*e; b1 = s*(1-e)/beta
        i0, i1, i2 = fisher_inc(mu_s, al_s, beta, s, w); J[0] += i0; J[1] += i1; J[2] += i2
        f0, f1, f2 = fisher_inc(mu, al, beta, s, w)
        if eta:
            F[0][0] *= (1-eta); F[0][1] *= (1-eta); F[1][0] *= (1-eta); F[1][1] *= (1-eta)
        F[0][0] += f0; F[0][1] += f1; F[1][0] += f1; F[1][1] += f2
        g0 = w - 1.0/(mu + al*a1); g1 = b1 - a1/(mu + al*a1)
        det = F[0][0]*F[1][1] - F[0][1]**2
        d0 = (F[1][1]*g0 - F[0][1]*g1)/det; d1 = (-F[0][1]*g0 + F[0][0]*g1)/det
        mu = max(mu - d0, mu_min); al = min(max(al - d1, al_min), al_max)
        s = a1 + 1.0; Tt += w
        mt_end = mu_s*(1 + A*math.sin(2*math.pi*Tt/P)) if A != 0 else mu_s
        ev.append((Tt, mu, al, mt_end))
        x0, x1 = mu - mu_s, al - al_s
        e_n = J[0]*x0*x0 + 2*J[1]*x0*x1 + J[2]*x1*x1
    return np.array(ev), e_n

if __name__ == "__main__":
    print("Fidelity: normalised error e_n of the copy (no drift, eta=0) vs Problem 5's run() at the same seed")
    for eps, seed in ((0.3, 1000), (0.1, 1001), (0.03, 1002)):
        ref = p5_run(1.0, 1.0, eps, 2000, {2000}, seed, (2.0, 0.5))[2000][1]
        ev, e_n = run_drift(1.0, 1.0, eps, seed, (2.0, 0.5), n_max=2000)
        print(f"  eps={eps}: Problem 5 e_n = {ref:.10f}   copy e_n = {e_n:.10f}   difference {abs(ref-e_n):.2e}")
