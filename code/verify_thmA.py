"""Independent verification of the Step 1-3 claims of Theorem A (problem10_rigorous.md), beta = 1 WLOG (all claims are scale-covariant).
 alpha = 1-eps, h(x)=e^-x-1+x, c(s)=eps*s+h(alpha s), G=s/c, D_eps = (alpha^2 s^2/2 - h)/((eps s + h)(eps + alpha^2 s/2))."""
import numpy as np, math
from scipy.integrate import quad
from p10_fisher_k1 import exact_m

h = lambda x: np.expm1(-x) + x
def D(s, eps):
    a = 1 - eps; x = a*s
    return (x*x/2 - h(x))/((eps*s + h(x))*(eps + a*a*s/2))
def c(s, eps): return eps*s + h((1-eps)*s)

print("(1) Step 1 identity  s/c(s) = 1/(eps + alpha^2 s/2) + D_eps(s):  max relative residual over a grid")
rng = np.random.default_rng(0); worst = 0
for eps in (0.9, 0.5, 0.25, 0.1, 0.01, 1e-3):
    s = np.logspace(-3, 3, 400); a = 1 - eps
    lhs = s/c(s, eps); rhs = 1/(eps + a*a*s/2) + D(s, eps); worst = max(worst, np.max(np.abs(lhs - rhs)/lhs))
print(f"    {worst:.2e}")

print("(2) h bounds: (x^2/2)(1-x/3) <= h <= x^2/2 and x^2/2 - h <= x^3/6 for x in [0,50]:")
x = np.linspace(0, 50, 200001)
print("    ", bool(np.all(h(x) >= x**2/2*(1 - x/3) - 1e-12)), bool(np.all(h(x) <= x**2/2 + 1e-12)), bool(np.all(x**2/2 - h(x) <= x**3/6 + 1e-12)))

print("(3) 'for x<=1.5: D_eps <= 4/(3 alpha)'  -> max over grid of D*alpha*3/4 (must be <= 1):")
worst = 0
for eps in (0.9, 0.5, 0.25, 0.1, 0.01, 1e-3, 1e-5):
    a = 1 - eps; s = np.linspace(1e-6, 1.5/a, 20000); worst = max(worst, np.max(D(s, eps)*a*3/4))
print(f"    {worst:.4f}")

print("(4) 'for x>=1.5, D_eps bounded by a constant' -- as written the denominator bound gives x/(h(1.5) alpha) (unbounded). With h(x) >= x h(1.5)/1.5 the bound is 3/(2 h(1.5) alpha):")
worst = 0; hb = 3/(2*h(1.5))
for eps in (0.9, 0.5, 0.25, 0.1, 0.01, 1e-3, 1e-5):
    a = 1 - eps; s = np.logspace(math.log10(1.5/a), 7, 4000); worst = max(worst, np.max(D(s, eps)*a)/hb)
print(f"    max D*alpha / (3/(2h(1.5))) = {worst:.4f}  (<=1 means the repaired bound holds; sup D*alpha over x>=1.5 = {worst*hb:.4f})")

print("(5) 'for x>=1: x^2/2 - h(x) >= 0.132 x^2':", bool(np.all((x[x>=1]**2/2 - h(x[x>=1])) >= 0.132*x[x>=1]**2 - 1e-12)))
print("    'for eps<=1/4, x>=1: denominator <= alpha x^2' (denominator/(alpha x^2), sup over x>=1 and eps in (0,1/4]):")
worst = 0; arg = None
for eps in (0.25, 0.2, 0.1, 0.01, 1e-3):
    a = 1 - eps; xx = np.logspace(0, 7, 20000); s = xx/a
    den = (eps*s + h(xx))*(eps + a*a*s/2); r = np.max(den/(a*xx*xx)); 
    if r > worst: worst, arg = r, (eps, xx[np.argmax(den/(a*xx*xx))])
print(f"    sup = {worst:.4f} at (eps, x) = {arg}  (claim needs <= 1)")
print("    crude bound used in the text's reasoning would give 20/18 =", round(20/18, 3), " (>1), so the written justification does not prove it, but the statement is numerically true.")

print("(6) 'D_eps(s) >= 0.13/beta for s >= 1/alpha, eps<=1/4': min of D*beta (beta=1):")
worst = 1e9
for eps in (0.25, 0.1, 0.01, 1e-3, 1e-5):
    a = 1 - eps; s = np.logspace(math.log10(1/a), 7, 4000); worst = min(worst, np.min(D(s, eps)))
print(f"    {worst:.4f}  (claim >= 0.13)")

print("(7) F_eps(u) <= exp(-0.13 mu (u - 1/alpha)) uniformly, using exact phi = exp(-beta mu Gamma) and P_eps = (1+au/eps)^{-k}: check F = phi/P on a grid")
worst = -1e9
for mu in (0.3, 1.0, 2.0):
    for eps in (0.25, 0.1, 0.01):
        a = 1 - eps; k = 2*mu/a**2; aa = a*a/2
        Gam = lambda u: quad(lambda s: s/c(s, eps), 0, u, limit=200)[0]
        for u in (0.5, 1, 2, 4, 8, 16, 30):
            F = math.exp(-mu*Gam(u) + k*math.log(1 + aa*u/eps))
            bound = math.exp(-0.13*mu*(u - 1/a)) if u > 1/a else 1.0
            worst = max(worst, F - bound)
print(f"    max (F - bound) = {worst:.2e}  (<= 0 means the bound holds on the grid)")

print("(8) Theorem A limits (beta=1, exact m via the ODE):")
for k0 in (2.0, 0.6, 1.0):
    mu = k0/2
    if k0 > 1: pred = 2/(2*mu - 1)
    if k0 < 1:
        F0 = lambda u: math.exp(-mu*quad(lambda s: ((s*s/2 - h(s))/(h(s)*s/2)) if s > 1e-8 else (1/3)*s*0 + (2/3 if False else 0), 0, u)[0])
        # D0(s) = (s^2/2 - h(s))/(h(s) s/2); small-s limit: (s^3/6)/((s^2/2)(s/2)) = 2/3
        D0 = lambda s: ((s*s/2 - h(s))/(h(s)*s/2)) if s > 1e-4 else 2/3
        F0 = lambda u: math.exp(-mu*quad(D0, 0, u, limit=200)[0])
        Cc = (2)**k0*quad(lambda u: u**(-k0)*F0(u), 0, 60, limit=400)[0]
    for eps in (1e-2, 1e-3, 1e-4):
        m = exact_m(mu, 1 - eps, 1.0)
        if k0 > 1: print(f"    k0={k0}: eps={eps:g}: m/eps = {m/eps:.4f}  (limit {pred:.4f})")
        elif k0 == 1.0: print(f"    k0=1: eps={eps:g}: m/(eps ln(1/eps)) = {m/(eps*math.log(1/eps)):.4f} (limit 2); local d(m/eps)/d ln(1/eps) approx -> 2")
        else: print(f"    k0={k0}: eps={eps:g}: m/eps^k0 = {m/eps**k0:.4f}  (limit C = {Cc:.4f})")
