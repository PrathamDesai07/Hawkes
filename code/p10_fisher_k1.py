"""Problem 10: Fisher information of K=1 Hawkes, lambda = mu + alpha*S, dS=-beta S dt + dN.
(Restored 25 Oct 2026: the original file had disappeared from this folder while other scripts, e.g. p1_check.py, import it.
 Same interface as before: exact_m(mu, alpha, beta), simulate(mu, alpha, beta, T, seed).)
Stationary Laplace transform phi(u)=E exp(-u lam) solves
  phi'(u) = -u*beta*mu*phi / (u*beta + exp(-u*alpha) - 1).
Then m=E[1/lam]=int_0^inf phi, E lam = mu/(1-n), and
  I_mm=m, I_ma=(1-mu*m)/alpha, I_aa=(E lam-2mu+mu^2 m)/alpha^2, det=(m*E lam-1)/alpha^2.
Gamma (CIR) approx: shape k=2*beta*mu/alpha^2, scale alpha^2/(2(beta-alpha)).
Cross-checked 25 Oct against Herrera-Marin arXiv 2609.38268 Thm 4.1 (same entries after the reparametrisation (log mu, log(1-rho)))."""
import numpy as np
from scipy.integrate import solve_ivp


def exact_m(mu, a, b):
    """m = E[1/lam] at stationarity, for alpha=a, beta=b."""
    U = 60.0/mu
    s = solve_ivp(lambda u, y: [-u*b*mu/(u*b + np.expm1(-u*a))*1.0, np.exp(y[0])],
                  [1e-12, U], [0.0, 0.0], rtol=1e-10, atol=1e-13, method="LSODA")
    return s.y[1, -1]


def simulate(mu, a, b, T, seed=0):
    """Time average of 1/lam over [0,T] by exact thinning (waiting time accumulated over rejections)."""
    rng = np.random.default_rng(seed)
    lam = mu; t = 0.0; I = 0.0
    while t < T:
        c = lam - mu
        w = 0.0
        while True:
            w += rng.exponential(1/lam)
            if rng.random() < (mu + c*np.exp(-b*w))/lam: break
        I += (w + np.log((mu + c*np.exp(-b*w))/(mu + c))/b)/mu
        lam = mu + c*np.exp(-b*w) + a
        t += w
    return I/t
