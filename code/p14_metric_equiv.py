"""Problem 14: uniform equivalence of the Fisher metrics along a drift box (K=1, exact I from m=E[1/lam]).
kappa_I(box) = max over theta, theta' in box of  lam_max / lam_min of  I(theta')^{-1/2} I(theta) I(theta')^{-1/2}.
Drift box: mu in [mu0, c_mu*mu0], eps in [eps0, c_eps*eps0], beta=1 fixed.  Compare with the condition number of a single I (~eps^-2)."""
import numpy as np
from scipy.linalg import eigh
from p10_fisher_k1_exact import exact_m

def Imat(mu, eps, beta=1.0):
    a = (1-eps)*beta; m = exact_m(mu, beta, eps); El = mu/eps
    return np.array([[m, (1-mu*m)/a], [(1-mu*m)/a, (El-2*mu+mu**2*m)/a**2]])

def kappa(mu0, c_mu, eps0, c_eps, n=5):
    pts = [(mu0*c_mu**(i/(n-1)), eps0*c_eps**(j/(n-1))) for i in range(n) for j in range(n)]
    Is = [Imat(m_, e_) for m_, e_ in pts]; best = (1.0, 1.0)
    lo, hi = np.inf, 0.0
    for A in Is:
        for B in Is:
            w = eigh(A, B, eigvals_only=True); lo = min(lo, w[0]); hi = max(hi, w[-1])
    return hi/lo, lo, hi

if __name__ == "__main__":
    print("kappa_I over drift boxes (mu0=1, beta=1):")
    print(f"{'eps0':>6} {'mu range':>9} {'eps range':>10} {'kappa_I':>9} {'gen.eig range':>16} {'cond(I(theta0))':>16}")
    for eps0 in (0.2, 0.1, 0.05):
        I0 = Imat(1.0, eps0); c0 = np.linalg.cond(I0)
        for c_mu, c_eps in ((1.5, 1.0), (1.0, 1.5), (1.5, 1.5), (1.0, 2.0), (2.0, 2.0)):
            k, lo, hi = kappa(1.0, c_mu, eps0, c_eps)
            print(f"{eps0:6.2f} {'x'+str(c_mu):>9} {'x'+str(c_eps):>10} {k:9.2f} {lo:7.3f}-{hi:7.3f} {c0:16.1f}")
