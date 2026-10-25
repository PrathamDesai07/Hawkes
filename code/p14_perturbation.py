"""Problem 14, amplitude effect: SECOND-ORDER perturbation theory (sympy, harmonic algebra) for the mean error  e' = -Gamma(t) e - A w cos(w t),  time t,
   h := eta*r0 (gain rate), w := 2 pi/P.
   fixed    : Gamma/h = Phi(t) = I(mu_t)/I_ref = 1 + a1 A s + a2 A^2 s^2 + ...,        s = sin(w t),  a1 = kappa_I = dln I/dln mu,  a2 = I''(1)/(2 I(1))
   adaptive : Gamma/h = Phi/psi,   psi' = h(-(1 + A s) psi + Phi)   (psi = lagged information in units of I_ref; the window rate eta r_t is modulated by the event rate)
 Result: closed forms  R = |E_1/e_1|^2 = 1 + C2 A^2 + O(A^4)  for the fundamental error^2 relative to the constant-gain filter, then numerical comparison with the exact ODE."""
import sympy as sp, numpy as np, math, sys
h, w, a1, a2, A = sp.symbols('h w a1 a2 A', positive=True)
I = sp.I

class Sig:
    """signal sum_n c_n e^{i n w t}, coefficients are sympy expressions (real signal: c_{-n} = conj(c_n))"""
    def __init__(self, d=None): self.d = {k: sp.simplify(v) for k, v in (d or {}).items() if v != 0}
    def __add__(self, o):
        r = dict(self.d)
        for k, v in o.d.items(): r[k] = r.get(k, 0) + v
        return Sig(r)
    def __sub__(self, o): return self + o.scale(-1)
    def scale(self, c): return Sig({k: c*v for k, v in self.d.items()})
    def __mul__(self, o):
        r = {}
        for k1, v1 in self.d.items():
            for k2, v2 in o.d.items(): r[k1 + k2] = r.get(k1 + k2, 0) + v1*v2
        return Sig(r)
    def solve_lin(self, rate):   # solve x' + rate*x = self  (rate scalar expression) harmonic by harmonic
        return Sig({k: v/(rate + I*k*w) for k, v in self.d.items()})
    def get(self, n): return self.d.get(n, 0)

s = Sig({1: 1/(2*I), -1: -1/(2*I)})            # sin(wt)
cosw = Sig({1: sp.Rational(1, 2), -1: sp.Rational(1, 2)})
one = Sig({0: 1})

def fundamental_ratio(Gamma_over_h):
    """Gamma_over_h: list [G0, G1, G2] of Sig (orders A^0, A^1, A^2) of Gamma/h.  Solve e = A e1 + A^2 e2 + A^3 e3."""
    G1, G2 = Gamma_over_h[1], Gamma_over_h[2]
    e1 = cosw.scale(-w).solve_lin(h)                                   # e1' + h e1 = -w cos
    e2 = (G1*e1).scale(-h).solve_lin(h)                                # e2' + h e2 = -h G1 e1
    e3 = ((G1*e2) + (G2*e1)).scale(-h).solve_lin(h)                    # e3' + h e3 = -h (G1 e2 + G2 e1)
    return sp.simplify(e3.get(1)/e1.get(1))                            # relative correction of the fundamental amplitude at O(A^2)

# fixed preconditioner
G1f = s.scale(a1); G2f = (s*s).scale(a2)
rel_fixed = fundamental_ratio([one, G1f, G2f])
C2_fixed = sp.simplify(2*sp.re(sp.simplify(rel_fixed)))
# adaptive: psi = 1 + A psi1 + A^2 psi2
psi1 = (s.scale(a1 - 1)).scale(h).solve_lin(h)                          # psi1' + h psi1 = h (a1 - 1) s
psi2 = ((s*psi1).scale(-1) + (s*s).scale(a2)).scale(h).solve_lin(h)     # psi2' + h psi2 = h(-s psi1 + a2 s^2)
# Gamma/h = Phi/psi = 1 + A(a1 s - psi1) + A^2 (a2 s^2 - a1 s psi1 + psi1^2 - psi2)
G1a = s.scale(a1) - psi1
G2a = (s*s).scale(a2) - (s*psi1).scale(a1) + (psi1*psi1) - psi2
rel_adapt = fundamental_ratio([one, G1a, G2a])
C2_adapt = sp.simplify(2*sp.re(sp.simplify(rel_adapt)))
print("C2_fixed   =", sp.simplify(C2_fixed))
print("C2_adaptive=", sp.simplify(C2_adapt))

# numerical evaluation with the exact information function
sys.path.insert(0, '.')
from p10_fisher_k1 import exact_m
def info_functions(eps, kind):
    alpha = 1 - eps; Ef = lambda mu: mu/eps
    def Iscalar(mu):                                               # alpha known: information about mu per unit time
        return exact_m(mu, alpha, 1.0)
    def Ieff(mu):                                                  # alpha unknown: efficient information
        m = exact_m(mu, alpha, 1.0); return (m - eps/mu)/(1 - 2*eps + eps*mu*m)
    return Iscalar if kind == 'known' else Ieff
def derivs(Ifun, hh=2e-3):
    I0 = Ifun(1.0); Ip = (Ifun(1 + hh) - Ifun(1 - hh))/(2*hh); Ipp = (Ifun(1 + hh) - 2*I0 + Ifun(1 - hh))/hh**2
    return Ip/I0, Ipp/(2*I0)

def exact_ratio(a1v, a2v, Av, hv, wv, kind):
    """exact periodic solution via long ODE integration of the cubic-truncated... (use the actual Phi(x) = 1 + a1 x + a2 x^2 + higher via Taylor to 2nd order)"""
    from scipy.integrate import solve_ivp
    Phi = lambda t: 1 + a1v*Av*math.sin(wv*t) + a2v*(Av*math.sin(wv*t))**2
    if kind == 'fixed':
        rhs = lambda t, y: [-hv*Phi(t)*y[0] - Av*wv*math.cos(wv*t)]; y0 = [0.0]
    else:
        rhs = lambda t, y: [-hv*Phi(t)/y[1]*y[0] - Av*wv*math.cos(wv*t), hv*(-(1 + Av*math.sin(wv*t))*y[1] + Phi(t))]; y0 = [0.0, 1.0]
    Pp = 2*math.pi/wv; ncyc = 40
    sol = solve_ivp(rhs, [0, ncyc*Pp], y0, t_eval=np.linspace(0, ncyc*Pp, ncyc*400 + 1), rtol=1e-10, atol=1e-13, method='LSODA')
    t = sol.t; e = sol.y[0]; sel = t > 0.5*ncyc*Pp; tt = t[sel]; yy = e[sel]
    X = np.column_stack([np.ones_like(tt)] + [f(k*wv*tt) for k in (1, 2, 3) for f in (np.sin, np.cos)])
    c, *_ = np.linalg.lstsq(X, yy, rcond=None); fund = (c[1]**2 + c[2]**2)/2
    pred = (Av*Av/2)*wv*wv/(hv*hv + wv*wv)
    return fund/pred

if __name__ == "__main__":
    eps = float(sys.argv[1]) if len(sys.argv) > 1 else 0.4; eta = float(sys.argv[2]) if len(sys.argv) > 2 else 0.004; P = 1200.0
    hv = eta/eps; wv = 2*math.pi/P
    for kind, label in (('known', 'alpha known (scalar)'), ('eff', 'alpha unknown (efficient information, scalar projection)')):
        a1v, a2v = derivs(info_functions(eps, kind))
        subs = {h: hv, w: wv, a1: a1v, a2: a2v}
        cf = float(C2_fixed.subs(subs)); ca = float(C2_adapt.subs(subs))
        print(f"\n{label}: kappa_I = dlnI/dlnmu = {a1v:.3f}, a2 = I''/(2I) = {a2v:.3f}; h = eta r0 = {hv:.4f}, w = {wv:.5f} (h/w = {hv/wv:.2f})")
        print("    A    | fixed: exact ODE  1+C2 A^2 | adaptive: exact ODE  1+C2 A^2")
        for Av in (0.05, 0.1, 0.2, 0.35):
            rf = exact_ratio(a1v, a2v, Av, hv, wv, 'fixed'); ra = exact_ratio(a1v, a2v, Av, hv, wv, 'adaptive')
            print(f"   {Av:4.2f}  |   {rf:8.4f}     {1 + cf*Av*Av:8.4f}   |   {ra:8.4f}      {1 + ca*Av*Av:8.4f}", flush=True)
        print(f"    C2_fixed = {cf:.3f}, C2_adaptive = {ca:.3f}, C2_adaptive - C2_fixed = {ca - cf:.3f}")
