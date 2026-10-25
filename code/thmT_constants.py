"""Theorem T with explicit constants (K=1, k0=2mu/beta>1).  Dimensionless variable x = alpha*eps*t.  r = eps/(1-eps), k = k0/(1-eps)^2, a t = ((1-eps)/(2 eps)) x.
 d(y) = alpha*D_eps(-i y/alpha) = N/(E1 (r - i y/2)),  N=(1-cos y-y^2/2)+i(y-sin y),  E1 = c(-i y/alpha) = -(1-cos y)-i((1+r)y-sin y).
 z(x) = i (k/2) int_0^x d,  F~ = exp(z),  phi = G F~,  G = (1 - i b x)^-k with b=(1-eps)/(2 eps).  phi'/phi = L' with L'(x) = (k/2) x / E1(x)  (per dx).
Checks (A) |d|<=19, (B) |L'| majorants, (C) |phi| majorants, then the certified H^1 bound Phi(eps) vs the true ||phi-G||_H1, and TV bound Phi/(2 sqrt 2)."""
import numpy as np, math, sys
from scipy import special, integrate
def d_fun(y, r):
    y = np.asarray(y, float); out = np.empty(y.shape, complex)
    sm = y < 0.2; ys = y[sm]
    re = -ys**4/24 + ys**6/720 - ys**8/40320 + ys**10/3628800
    im = ys**3/6 - ys**5/120 + ys**7/5040 - ys**9/362880
    c1 = ys**2/2 - ys**4/24 + ys**6/720 - ys**8/40320 + ys**10/3628800     # 1-cos
    im1 = r*ys + im                                                         # (1+r)y - sin y
    out[sm] = (re + 1j*im)/((-c1 - 1j*im1)*(r - 0.5j*ys))
    yl = y[~sm]; N = (1-np.cos(yl)-yl**2/2) + 1j*(yl-np.sin(yl)); E1 = -(1-np.cos(yl)) - 1j*((1+r)*yl-np.sin(yl))
    out[~sm] = N/(E1*(r-0.5j*yl)); return out
def E1_fun(y, r):
    return -(1-np.cos(y)) - 1j*((1+r)*y-np.sin(y))
def grid(Xmax):
    return np.concatenate([np.geomspace(1e-13, 1.0, 8000)[:-1], np.arange(1.0, 2000.0, 0.04), np.arange(2000.0, Xmax, 0.2)])
def true_h1(k0, eps, Xmax=1.0e5):
    r = eps/(1-eps); k = k0/(1-eps)**2; b = (1-eps)/(2*eps); alpha = 1-eps      # beta = 1
    x = grid(Xmax); d = d_fun(x, r)
    cum = np.concatenate([[0], np.cumsum(0.5*(d[1:]+d[:-1])*np.diff(x))])
    z = 1j*(k/2)*cum; G = (1-1j*b*x)**(-k); Gp = 1j*k*b*(1-1j*b*x)**(-k-1)    # dG/dx
    Dl = G*np.expm1(z); Dlp = Gp*np.expm1(z) + G*(1j*(k/2)*d)*np.exp(z)
    sc = alpha*eps                                                             # t = x/sc
    n2 = 2*np.trapezoid(np.abs(Dl)**2, x)/sc; d2 = 2*np.trapezoid(np.abs(Dlp)**2, x)*sc  # both signs of t
    return math.sqrt(n2 + d2), x, d, z, G
def certified(k0, eps):
    r = eps/(1-eps); k = k0/(1-eps)**2; alpha = 1-eps; a = alpha/2; M1 = 19*k*alpha/2          # beta = 1: a = alpha^2/2 /beta... a=alpha^2/(2 beta)
    a = alpha**2/2; M1 = 19*k*alpha/2
    t1 = 1/(M1*eps); X0 = max(5.0, 6.0/(k0-1)); ts = X0/(alpha*eps); q = k0/(2*(1+3/X0)); assert q > 0.5
    at1 = a*t1; amp = math.e*at1**(-k)                                         # |phi| <= amp on [t1,t*]
    # head
    J = (a**-3*special.beta(1.5, k-1.5)) if k > 1.5 else 2*integrate.quad(lambda t: t*t*(1+a*a*t*t)**(-k), 0, t1, limit=400)[0]
    I0 = (1/a)*special.beta(0.5, k-0.5)
    head1 = (math.e-1)**2*M1**2*eps**2*J
    head2 = M1**2*eps**2*(k*(math.e-1)+math.e)**2*I0
    # tail pieces of phi
    Iphi = 2*(amp**2*(ts-t1) + amp**2*math.exp(2*k0/X0)*ts/(2*q-1))
    t1b = 1/(alpha*eps)
    Iphip = 2*((1.1*k)**2*amp**2*(1/t1 - 1/t1b) + (3.2*k*alpha*eps)**2*(amp**2*(ts-t1b) + amp**2*math.exp(2*k0/X0)*ts/(2*q-1)))
    IG = 2*a**(-2*k)*t1**(1-2*k)/(2*k-1); IGp = 2*k*k*a*a*a**(-2*k-2)*t1**(-2*k-1)/(2*k+1)
    N2 = head1 + 2*Iphi + 2*IG; D2 = head2 + 2*Iphip + 2*IGp
    return math.sqrt(N2 + D2), dict(head=math.sqrt(head1+head2), phi_tail=math.sqrt(2*Iphi+2*Iphip), G_tail=math.sqrt(2*IG+2*IGp), t1=t1, X0=X0, q=q)
if __name__ == "__main__":
    k0 = float(sys.argv[1]) if len(sys.argv) > 1 else 2.0
    # (A) sup |d| over y, r
    ys = np.concatenate([np.geomspace(1e-8, 1, 400), np.linspace(1, 400, 400000)])
    print("(A) sup_y |d(y)| over eps in {.9,.5,.1,1e-2,1e-4,1e-6}:", [round(float(np.max(np.abs(d_fun(ys, e/(1-e))))), 4) for e in (.9, .5, .1, 1e-2, 1e-4, 1e-6)], " (claim <= 19)")
    for eps in (0.1, 0.03, 0.01, 1e-3, 1e-4, 1e-5, 1e-6, 1e-8):
        th, x, d, z, G = true_h1(k0, eps); ph = G*np.exp(z)
        cb, parts = certified(k0, eps)
        r = eps/(1-eps); k = k0/(1-eps)**2; alpha = 1-eps; a = alpha**2/2; ts = x/(alpha*eps)
        Lp = (k/2)*x/E1_fun(x, r)                                              # d log phi / dx
        # (C): |phi| majorant on t>=t1
        t1 = parts['t1']; X0 = parts['X0']; q = parts['q']; M = 19*k*alpha/2; amp = math.e*(a*t1)**(-k); tstar = X0/(alpha*eps)
        maj = np.where(ts <= tstar, amp, amp*math.exp(k0/X0)*(tstar/np.maximum(ts, 1e-300))**q); mask = ts >= t1
        okC = bool(np.all(np.abs(ph[mask]) <= maj[mask]*(1+1e-9)))
        print(f"eps={eps:g}: true H1 = {th:.3e} (TV bound {th/2.828:.3e}); certified Phi = {cb:.3e} (TV bound {cb/2.828:.3e}); parts head={parts['head']:.2e} phi_tail={parts['phi_tail']:.2e} G_tail={parts['G_tail']:.2e}; (C) |phi|<=majorant: {okC}", flush=True)
