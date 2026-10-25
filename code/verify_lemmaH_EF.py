"""Independent check of Lemma H, Theorem E and Theorem F (problem11_dichotomy.md, Sec 4b). Own implementation of the generator.
 A f(S) = sum_k (-beta_k S_k) df/dS_k + lam (f(S+1) - f(S)),  f = log W + psi/W^2,  W = Z + z,  Z = omega.S,  psi = kappa.S."""
import numpy as np, math, sys
from scipy.optimize import brentq

def setup(beta, alpha, mu):
    beta = np.asarray(beta, float); alpha = np.asarray(alpha, float); n = float(np.sum(alpha/beta)); b, B = beta.min(), beta.max()
    g = brentq(lambda x: np.sum(alpha/(beta - x)) - 1.0, 0.0, b*(1 - 1e-12) if False else b - 1e-12*b) if n < 1 else 0.0
    om = alpha/(beta - g); nu = float(np.sum(om/beta)); c0 = n/(2*nu); kap = (c0*om - alpha/2)/beta; kbar = float(np.max(np.abs(kap)/om))
    return dict(beta=beta, alpha=alpha, mu=mu, n=n, b=b, B=B, gam=g, om=om, nu=nu, c0=c0, kap=kap, kbar=kbar, kK=mu/c0, eps=1 - n)

def consts(P, z):
    Bz = max(P['mu']/z, P['B'] - P['gam']); C2 = 2*P['kbar']*P['mu'] + 3*Bz*P['kbar']; return Bz, C2, Bz/3 + C2

def Af(P, S, z):
    """exact generator of f_z at states S (m x K)"""
    be, al, om, ka, mu = P['beta'], P['alpha'], P['om'], P['kap'], P['mu']
    Z = S@om; psi = S@ka; W = Z + z; lam = mu + S@al
    f = lambda Zv, psiv: np.log(Zv + z) + psiv/(Zv + z)**2
    grad_coef = om[None, :]/W[:, None] + ka[None, :]/W[:, None]**2 - 2*psi[:, None]*om[None, :]/W[:, None]**3     # df/dS_k
    drift = np.sum(-be[None, :]*S*grad_coef, axis=1)
    jump = lam*(f(Z + 1.0, psi + ka.sum()) - f(Z, psi))             # sum(kappa)=0 so psi does not jump
    return drift + jump, W

def lemmaH(P, z, m=200000, rng=np.random.default_rng(1)):
    K = len(P['beta']); scale = P['mu']/max(P['gam'], 1e-9)
    S = np.exp(rng.uniform(np.log(1e-3), np.log(5*scale), size=(m, K)))*rng.uniform(0, 1, size=(m, K))
    A, W = Af(P, S, z); g = P['gam']; mu, c0, kb = P['mu'], P['c0'], P['kbar']; Bz, C2, C4 = consts(P, z)
    lo = (mu - c0 - 2*kb*g)/W - g - (mu/2 + C2)/W**2
    hi = (mu - c0 + g*z + 2*kb*g)/W - g + (c0*z + C4)/W**2
    return float(np.min(A - lo)), float(np.min(hi - A)), float(np.min(((A - lo)/(np.abs(lo) + 1e-12)))), float(np.min((hi - A)/(np.abs(hi) + 1e-12)))

if __name__ == "__main__":
    print("Lemma H: minimum slack (A f - lower) and (upper - A f) over 2e5 random states per (parameter set, z).  Negative = violated.")
    sets = [("K=1 beta=1, mu=1.5 (kK>1)", [1.0], 0.99, 1.5, 1), ("K=1 beta=1, mu=0.3 (kK<1)", [1.0], 0.99, 0.3, 1),
            ("K=2 beta=(1,2) mu=1.5", [1.0, 2.0], 0.95, 1.5, 2), ("K=2 beta=(1,2) mu=0.4", [1.0, 2.0], 0.95, 0.4, 2),
            ("K=3 beta=(1,3,9) mu=1", [1.0, 3.0, 9.0], 0.98, 1.0, 3), ("K=3 beta=(1,2,4) mu=2.5", [1.0, 2.0, 4.0], 0.9, 2.5, 3)]
    worst = 1e9
    for name, beta, n, mu, K in sets:
        beta = np.array(beta); alpha = n*beta/K            # sum alpha/beta = n
        P = setup(beta, alpha, mu); row = []
        for z in (0.1, 1.0, 10.0, 100.0, 1000.0):
            lo, hi, rl, rh = lemmaH(P, z); row.append((z, lo, hi)); worst = min(worst, lo, hi)
        print(f"  {name:28s} eps={P['eps']:.2f} gamma={P['gam']:.4f} k_K={P['kK']:.2f} kbar={P['kbar']:.2f}: " + "; ".join(f"z={z:g}: {lo:+.1e}/{hi:+.1e}" for z, lo, hi in row))
    print(f"  overall minimum slack = {worst:.3e}")
