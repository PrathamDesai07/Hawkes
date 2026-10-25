"""Numerical checks of the repaired Step 2 of Theorem A and of the imaginary-axis bounds used for TV(pi_theta, pi_theta').  beta = 1."""
import numpy as np, math
h = lambda x: np.expm1(-x) + x
c_ = 1 - math.exp(-1)

print("(R1) h(x)/x^2 decreasing: q(x) = 2 - x - (x+2)e^-x < 0 for x>0 (max over grid):")
x = np.linspace(1e-6, 60, 400001); q = 2 - x - (x + 2)*np.exp(-x); print("    max q =", q.max(), " (must be < 0);  q(0)=q'(0)=0, q''=-x e^-x<0")
print("(R2) (x^2/2 - h)/x^2 at x=1 = 1/2 - 1/e =", 0.5 - math.exp(-1), " and nondecreasing:", bool(np.all(np.diff(0.5 - h(x)/x**2) >= -1e-15)))

print("(R3) repaired denominator bound: Den/(alpha x^2) = (r x + h(x))(r + x/2)/x^2 <= 2/3 + (4/9 - c/2)/x  <= 0.795 for x>=1, r in (0,1/3]:")
xx = np.logspace(0, 7, 20000); worst = 0; w2 = 0
for r in (1e-6, 0.01, 0.1, 0.2, 1/3):
    ratio = (r*xx + h(xx))*(r + xx/2)/xx**2; bound = 2/3 + (4/9 - c_/2)/xx
    worst = max(worst, np.max(ratio - bound)); w2 = max(w2, ratio.max())
print(f"    max(ratio - bound) = {worst:.2e} (<=0 means bound holds);  max ratio = {w2:.4f} <= {2/3 + 4/9 - c_/2:.4f}")

print("(R4) global bound D_eps(s) <= 2.08/alpha for all s>0 (pieces: x<=1.5: 4/(3 alpha); x>=1.5: 3/(2 h(1.5) alpha) = %.4f/alpha):" % (3/(2*h(1.5))))
worst = 0
for eps in (0.9, 0.5, 0.25, 0.1, 0.01, 1e-3, 1e-5):
    a = 1 - eps; s = np.logspace(-6, 7, 40000); xs = a*s
    D = (xs**2/2 - h(xs))/((eps*s + h(xs))*(eps + a*a*s/2)); worst = max(worst, np.max(D*a))
print(f"    sup D*alpha = {worst:.4f}  (<= 2.08)")

print("(R5) D_eps(s) >= 0.166/alpha for s >= 1/alpha, eps <= 1/4 (from 0.132/0.795):")
worst = 1e9
for eps in (0.25, 0.1, 0.01, 1e-3, 1e-5):
    a = 1 - eps; s = np.logspace(math.log10(1/a), 7, 20000); xs = a*s
    D = (xs**2/2 - h(xs))/((eps*s + h(xs))*(eps + a*a*s/2)); worst = min(worst, np.min(D*a))
print(f"    min D*alpha = {worst:.4f} (>= 0.166)")

print("(R6) imaginary axis: |D_eps(-i tau)| <= C/alpha, claimed C = 31.4 (x>=1) and 0.873 (x<=1); actual sup (accurate series for h at small |y|; a naive exp(-y)-1+y cancels catastrophically for |tau|<1e-4 and gave a spurious 52.7):")
def hseries(y, nterms=14):
    out = np.zeros_like(y); term = np.ones_like(y)
    for n in range(1, nterms + 1):
        term = term*(-y)/n
        if n >= 2: out = out + term
    return out
def hcomplex(y):
    small = np.abs(y) < 0.5; res = np.empty_like(y); res[small] = hseries(y[small]); res[~small] = np.exp(-y[~small]) - 1 + y[~small]; return res
worst = {"x<=1": 0.0, "x>=1": 0.0}
for eps in (0.9, 0.5, 0.25, 0.1, 0.01, 1e-3, 1e-5):
    a = 1 - eps; tau = np.logspace(-6, 7, 80000); s = -1j*tau; y = a*s; hh = hcomplex(y); N = y*y/2 - hh
    D = N/((eps*s + hh)*(eps + a*a*s/2)); x = a*tau
    worst["x<=1"] = max(worst["x<=1"], np.max(np.abs(D[x <= 1]))*a); worst["x>=1"] = max(worst["x>=1"], np.max(np.abs(D[x >= 1]))*a)
print("    sup |D|*alpha:", {k: round(v, 4) for k, v in worst.items()}, " (claims 0.873 and 31.4)")
