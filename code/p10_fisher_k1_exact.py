"""Problem 10: K=1 Fisher information, exact m=E[1/lam] and near-critical asymptotics.
phi(u)=E e^{-u lam}=exp(-beta*mu*int_0^u G),  G(s)=s/(eps*beta*s + h(alpha*s)),  h(x)=expm1(-x)+x,
m=int_0^inf phi(u)du.  I_mm=m, I_ma=(1-mu m)/alpha, I_aa=(mu/eps-2mu+mu^2 m)/alpha^2, det=(m*mu/eps-1)/alpha^2.
Regimes (k=2*mu/beta at criticality):  k>1: m~2 eps/(2mu-beta);  k=1: m~(2/beta)eps*ln(1/eps);
k<1: m~C eps^k,  C=(2/beta)^k int_0^inf u^{-k} F(u)du, F(u)=exp(-beta*mu int_0^u [g0(s)-2/(beta^2 s)]ds), g0(s)=s/(beta s+e^{-beta s}-1)."""
import numpy as np
from scipy.integrate import solve_ivp, quad

def h(x):
    x = np.asarray(x, float)
    return np.where(x < 0.05, x*x*(0.5 - x/6 + x*x/24 - x**3/120 + x**4/720), np.expm1(-x)+x)

def exact_m(mu, beta, eps):
    alpha = (1-eps)*beta
    G = lambda s: s/(eps*beta*s + float(h(alpha*s)))
    def rhs(t, y):
        u = np.exp(t); return [u*G(u), np.exp(-beta*mu*y[0])*u]
    t0, t1 = np.log(1e-9), np.log(80/mu)
    y0 = [1e-9*G(1e-9), 1e-9]
    s = solve_ivp(rhs, [t0, t1], y0, method="LSODA", rtol=1e-11, atol=1e-16)
    return s.y[1, -1]

def C_below(mu, beta):
    k = 2*mu/beta
    g0 = lambda s: s/(beta*s + float(h(beta*s))) if False else s/(beta*s - beta*s + float(h(beta*s)))   # eps=0: denominator h(beta s)
    def corr(u):   # int_0^u [g0 - 2/(beta^2 s)] ds
        f = lambda s: (g0(s) - 2/(beta**2*s)) if s > 1e-8 else 2/(3*beta)
        return quad(f, 0, u, limit=200)[0]
    integrand = lambda u: u**(-k)*np.exp(-beta*mu*corr(u))
    I = quad(integrand, 0, 1, limit=200)[0] + quad(integrand, 1, 400, limit=400)[0]
    return (2/beta)**k*I

def fisher(mu, beta, eps, m):
    a = (1-eps)*beta; El = mu/eps
    I = np.array([[m, (1-mu*m)/a], [(1-mu*m)/a, (El-2*mu+mu**2*m)/a**2]])
    return I

if __name__ == "__main__":
    beta = 1.0
    for mu, name in ((1.0, "k=2>1"), (0.5, "k=1"), (0.3, "k=0.6<1")):
        k = 2*mu/beta
        print(f"\n== mu={mu}, beta={beta}  ({name}) ==")
        pred = {True: None}
        if k < 1: C = C_below(mu, beta); print(f"predicted C = {C:.5f}")
        print(f"{'eps':>8} {'m':>11} {'m/eps':>9} {'m/eps^k':>9} {'dlog m/dlog eps':>16} {'lam_min(I)':>11} {'det*beta^2':>10}")
        prev = None
        for eps in (1e-1, 3e-2, 1e-2, 3e-3, 1e-3, 3e-4, 1e-4):
            m = exact_m(mu, beta, eps); I = fisher(mu, beta, eps, m)
            sl = "" if prev is None else f"{np.log(m/prev[1])/np.log(eps/prev[0]):16.3f}"
            print(f"{eps:8.0e} {m:11.5e} {m/eps:9.4f} {m/eps**k:9.4f} {sl:>16} {np.linalg.eigvalsh(I)[0]:11.4e} {np.linalg.det(I)*beta**2*(1-eps)**2:10.4f}")
            prev = (eps, m)
        if k > 1: print(f"predicted m/eps -> 2/(2mu-beta) = {2/(2*mu-beta):.4f};  det*beta^2 -> beta/(2mu-beta)... = {beta/(2*mu-beta):.4f}")
