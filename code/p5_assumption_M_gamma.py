"""Assumption M via the Gamma comparator (Theorem T, Step 7): pi(I) >= Gamma(I) - C eps^eta.  Trapezoid h: support [c-d,c+d], ramps width d/2, plateau [c-d/2,c+d/2].
S = (lam-mu)/alpha;  lam ~ Gamma(k, scale a/eps) (comparator), a = alpha^2/(2 beta), k = 2 beta mu/alpha^2.  Report E_Gamma h and the plateau mass for each eps, vs the simulated m (0.15, 0.05)."""
import numpy as np
from scipy import stats, integrate
def trap(s, c, d):
    x = np.abs(s-c); return np.clip((d - x)/(d/2), 0, 1)
for mu, beta in ((1.0, 1.0), (0.6, 1.0)):
    for eps in (0.3, 0.1, 0.03, 0.01, 0.0):
        al = (1-eps)*beta; a = al*al/(2*beta); k = 2*beta*mu/al**2
        scale = a/max(eps, 1e-12) if eps > 0 else None
        ES = mu/(beta*eps) if eps > 0 else None
        row = []
        for c in (0.5, 1.8):
            if eps > 0:
                g = stats.gamma(k, scale=scale)
                f = lambda s: trap(s, c*ES, 0.15*ES)*g.pdf(mu + al*s)*al
                Eh = integrate.quad(f, (c-0.15)*ES, (c+0.15)*ES, points=[(c-.075)*ES, (c+.075)*ES], limit=200)[0]
                pl = g.cdf(mu + al*(c+.075)*ES) - g.cdf(mu + al*(c-.075)*ES)
            else:  # limit: Gamma(k0, mean 1) in units of the mean
                k0 = 2*mu/beta; g = stats.gamma(k0, scale=1/k0)
                Eh = integrate.quad(lambda x: trap(x, c, .15)*g.pdf(x), c-.15, c+.15, points=[c-.075, c+.075])[0]
                pl = g.cdf(c+.075) - g.cdf(c-.075)
            row.append((round(Eh, 4), round(pl, 4)))
        print(f"mu={mu} k0={2*mu/beta:.1f} eps={eps}:  I1 (E h, plateau mass) = {row[0]},  I2 = {row[1]}")
