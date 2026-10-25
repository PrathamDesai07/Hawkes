"""Numerical checks of the lemmas used in the rigorous K=1 stationary asymptotics (problem10 note, Theorem A)."""
import numpy as np
from p10_fisher_k1_exact import h
from p10b_transient import grids

# (L1) h(x) between x^2/2 (1-x/3) and x^2/2 ; N(x)=x^2/2-h(x) >= 0.132 x^2 for x>=1
x = np.concatenate([np.linspace(1e-6, 3, 20000), np.linspace(3, 200, 20000)])
hx = h(x)
print("L1 h<=x^2/2:", bool((hx <= x*x/2+1e-15).all()), "  h>=x^2/2(1-x/3):", bool((hx >= x*x/2*(1-x/3)-1e-15).all()))
xm = x[x >= 1]; print("   min N(x)/x^2 on x>=1 =", float(((xm*xm/2-h(xm))/xm**2).min()), "(claim >= 0.132)")

# (L2) phi <= Gamma LT P, D_eps >= 0, F decreasing, D_eps >= 0.13/beta for s>=1/alpha
for mu, beta, eps in ((1.0, 1.0, 0.1), (0.5, 1.0, 0.01), (0.3, 1.0, 1e-3), (2.0, 3.0, 0.2)):
    g = grids(mu, beta, eps, N=40000, lo=1e-13, hi=60/mu)
    X, phi, al = g['x'], g['phi'], g['alpha']
    a_ = al**2/(2*beta); k_e = 2*beta*mu/al**2
    P = (1+a_*X/eps)**(-k_e); F = phi/P
    s = X; D = (al**2*s**2/2 - h(al*s))/((eps*beta*s+h(al*s))*(eps*beta+al**2*s/2))
    big = s >= 1/al
    print(f"L2 mu={mu} beta={beta} eps={eps}: phi<=P {bool((phi<=P*(1+1e-6)).all())}, F nonincreasing {bool((np.diff(F)<=1e-9).all())}, "
          f"D>=0 {bool((D>=-1e-12).all())}, min D (s>=1/alpha) = {D[big].min():.4f} vs 0.132/alpha={0.132/al:.4f}")

# (L3) k=1 constant: slope of m/eps against ln(1/eps) should approach 2/beta
beta, mu = 1.0, 0.5
ms = {}
for e in (1e-4, 1e-5, 1e-6, 1e-7, 1e-8):
    g = grids(mu, beta, e, N=120000, lo=1e-18, hi=80); ms[e] = g['Cphi'][-1]
es = sorted(ms, reverse=True)
print("L3 k=1: eps, m/eps, local slope d(m/eps)/d ln(1/eps) (predicted -> 2/beta=2):")
for i, e in enumerate(es):
    sl = "" if i == 0 else f"{(ms[e]/e-ms[es[i-1]]/es[i-1])/np.log(es[i-1]/e):.4f}"
    print(f"   {e:8.0e} {ms[e]/e:10.4f} {sl}")
