"""Problem 8, T1 and T2 at the level of a single gap (K=1): for a fixed start state s and fixed theta, sample the gap exactly by inversion of the compensator
   Lambda*(u) = mu* u + alpha* s (1-e^{-beta u})/beta, and compare the means of  r = l(theta) - l(theta*)  (realised log-loss regret) and  R = int_0^tau lam* psi(lam_theta/lam*) du
   (pseudo-regret); also test  Lemma 8.1 on a grid  and  Lemma 8.2  gap by gap (violations)."""
import numpy as np
rng = np.random.default_rng(1)
def gl(n=40): x, w = np.polynomial.legendre.leggauss(n); return 0.5*(x+1), 0.5*w
xg, wg = gl()
def sample(mu_s, al_s, beta, s, N):
    E = rng.exponential(size=N); tau = np.empty(N)
    # invert Lambda*(tau) = E by bisection (Lambda* is increasing)
    lo = np.zeros(N); hi = E/mu_s + 1.0
    for _ in range(80):
        mid = 0.5*(lo+hi); L = mu_s*mid + al_s*s*(1-np.exp(-beta*mid))/beta
        big = L > E; hi = np.where(big, mid, hi); lo = np.where(big, lo, mid)
    return 0.5*(lo+hi)
def onestep(mu_s, al_s, beta, s, th, N=400000):
    tau = sample(mu_s, al_s, beta, s, N); mu, al = th
    lam_s = lambda u: mu_s + al_s*s*np.exp(-beta*u); lam_t = lambda u: mu + al*s*np.exp(-beta*u)
    Ls = mu_s*tau + al_s*s*(1-np.exp(-beta*tau))/beta; Lt = mu*tau + al*s*(1-np.exp(-beta*tau))/beta
    r = (Lt - Ls) - np.log(lam_t(tau)/lam_s(tau))
    u = tau[:, None]*xg[None, :]; w = tau[:, None]*wg[None, :]
    ls = lam_s(u); lt = lam_t(u)
    R = (w*(lt - ls - ls*np.log(lt/ls))).sum(1); Q = (w*(lt-ls)**2/ls).sum(1)
    rat = lt/ls; cp = np.maximum(1, rat.max(1)); cm = np.minimum(1, rat.min(1))
    lo_ok = R >= Q/(2*cp**2)*(1-1e-9); hi_ok = R <= Q/(2*cm**2)*(1+1e-9); lo2 = R >= Q/(2*cp)*(1-1e-9)
    return r.mean(), r.std()/np.sqrt(N), R.mean(), R.std()/np.sqrt(N), (r-R).mean(), (r-R).std()/np.sqrt(N), (~lo_ok).sum()+(~hi_ok).sum(), (~lo2).sum()
print(f"{'mu*,al*,beta':>16s} {'s':>5s} {'theta':>12s} {'E r':>8s} {'E R':>8s} {'E(r-R)':>10s} {'se':>7s} {'viol (8.2)':>10s} {'viol R>=Q/2c+':>13s}")
for (mu_s, al_s, beta) in ((1.0, 0.9, 1.0), (1.0, 0.7, 1.0), (1.5, 1.2, 1.0)):
    for s in (0.0, 1.0, 8.0):
        for th in ((mu_s, al_s), (1.5*mu_s, 0.6*al_s), (0.6*mu_s, 1.05*al_s), (2.5*mu_s, 0.3*al_s)):
            if th[1] >= beta: continue
            a, b, c, d, e, f, v, v2 = onestep(mu_s, al_s, beta, s, th)
            print(f"{str((mu_s,al_s,beta)):>16s} {s:5.1f} {str((round(th[0],2),round(th[1],2))):>12s} {a:8.4f} {c:8.4f} {e:10.5f} {f:7.5f} {v:10d} {v2:13d}")
# Lemma 8.1 on a grid
r = np.exp(np.linspace(-8, 8, 400001)); psi = r-1-np.log(r)
rp = np.maximum(1, r); rm = np.minimum(1, r)
ok = [(psi >= (r-1)**2/(2*rp**2)*(1-1e-12)).all(), (psi <= (r-1)**2/(2*rm**2)*(1+1e-12)).all(), (psi >= (r-1)**2/(2*rp)*(1-1e-12)).all(), (psi <= (r-1)**2/r*(1+1e-12)).all()]
print("Lemma 8.1 on r in [e^-8, e^8], 4e5 points: lower(sq), upper(sq), lower(linear), upper (r-1)^2/r all hold:", ok)
