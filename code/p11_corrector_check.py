"""Problem 11: pathwise check of the corrected Lyapunov function f = log W + psi/W^2 (any K).
 gamma: root of sum alpha_k/(beta_k-gamma)=1;  omega_k=alpha_k/(beta_k-gamma) (sum=1);  Z=omega.S, W=Z+z0;
 nu~=sum omega_k/beta_k;  c0=n/(2 nu~);  kappa_k=(-alpha_k/2+c0 omega_k)/beta_k (sum=0);  psi=kappa.S.
 Claimed:  LB <= A f <= UB  with
  LB=(mu-c0-2kb*gam)/W-gam-(mu/2+C2)/W^2,    UB=(mu-c0+gam*z0+2kb*gam)/W-gam+(c0*z0+C4)/W^2,
  kb=max|kappa_k|/omega_k, B=max(mu/z0, beta_max-gam), C2=2kb*mu+3B*kb, C4=B/3+2kb*mu+3B*kb."""
import numpy as np
from scipy.optimize import brentq

def constants(mu, beta, alpha, z0):
    n = (alpha/beta).sum(); eps = 1-n
    g = brentq(lambda x: (alpha/(beta-x)).sum()-1, 0, beta.min()*(1-1e-12))
    om = alpha/(beta-g); nut = (om/beta).sum(); c0 = n/(2*nut)
    kap = (-alpha/2 + c0*om)/beta; kb = np.max(np.abs(kap)/om)
    B = max(mu/z0, beta.max()-g)
    return dict(n=n, eps=eps, gam=g, om=om, nut=nut, c0=c0, kap=kap, kb=kb, B=B,
                C2=2*kb*mu+3*B*kb, C4=B/3+2*kb*mu+3*B*kb, z0=z0, mu=mu, beta=beta, alpha=alpha)

def Af(S, c):
    mu, beta, alpha, om, kap, z0 = c['mu'], c['beta'], c['alpha'], c['om'], c['kap'], c['z0']
    lam = mu + alpha@S; Z = om@S; W = Z+z0; psi = kap@S
    f = lambda Z_, psi_: np.log(Z_+z0) + psi_/(Z_+z0)**2
    drift_Z = -(om*beta)@S; drift_psi = -(kap*beta)@S
    # d/dt of f along the flow:
    df = drift_Z/W + (drift_psi/W**2 - 2*psi*drift_Z/W**3)
    jump = lam*(f(Z+om.sum(), psi+kap.sum()) - f(Z, psi))
    return df + jump, W

if __name__ == "__main__":
    rng = np.random.default_rng(0); worst = [np.inf, -np.inf]
    for (beta, mu, eps, z0) in ((np.array([1.]), 0.3, 0.01, 20.), (np.array([1., 2.]), 0.8, 0.01, 30.),
                                (np.array([1., 2.]), 0.4, 0.05, 50.), (np.array([1., 3., 9.]), 1.0, 0.02, 40.),
                                (np.array([0.5, 1., 4.]), 0.7, 0.1, 25.)):
        K = len(beta); w = rng.dirichlet(np.ones(K)); alpha = (1-eps)*w*beta
        c = constants(mu, beta, alpha, z0); minslack = np.inf; maxslack = np.inf
        for _ in range(200000):
            S = rng.exponential(1, K)*10**rng.uniform(-3, 4, K)*rng.random(K)
            a, W = Af(S, c)
            LB = (mu-c['c0']-2*c['kb']*c['gam'])/W - c['gam'] - (mu/2+c['C2'])/W**2
            UB = (mu-c['c0']+c['gam']*z0+2*c['kb']*c['gam'])/W - c['gam'] + (c['c0']*z0+c['C4'])/W**2
            minslack = min(minslack, a-LB); maxslack = min(maxslack, UB-a)
        print(f"beta={beta.tolist()} mu={mu} eps={eps}: gamma={c['gam']:.4f} c0={c['c0']:.4f} 2mu*nu~/n={2*mu*c['nut']/c['n']:.3f} "
              f"| min(Af-LB)={minslack:.2e}  min(UB-Af)={maxslack:.2e}  (both must be >= 0)")
