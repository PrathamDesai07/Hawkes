"""Randomised cross-check of the INEQUALITIES used in the lemmas (fresh seeds; counts of violations over many random instances)."""
import numpy as np, math
from math import comb
rng = np.random.default_rng(2026)
rows = []
def report(name, ok, detail=""):
    rows.append((name, ok)); print(f"{'OK  ' if ok else 'FAIL'} {name} {detail}", flush=True)

def rand_params(K, eps=None):
    beta = np.sort(rng.uniform(0.3, 6.0, K)); eps = rng.uniform(0.01, 0.6) if eps is None else eps
    sh = rng.dirichlet(np.ones(K)); alpha = (1 - eps)*sh*beta; mu = rng.uniform(0.1, 4.0)
    return mu, beta, alpha, eps

# 1. Lemma 1' pathwise inequality  A Z^p <= -p eps b Z^p + sum_i a_{p,i} Z^{p-i}
viol = 0; tot = 0
for _ in range(4000):
    K = int(rng.integers(1, 5)); mu, beta, alpha, eps = rand_params(K); b, B = beta.min(), beta.max()
    c = alpha/beta; n = c.sum(); S = rng.exponential(rng.uniform(0.1, 30.0), K)*rng.integers(0, 2, K).clip(0, 1) if rng.random() < .3 else rng.exponential(rng.uniform(0.1, 30.0), K)
    Z = c@S; lam = mu + alpha@S
    for p in range(1, 7):
        A = p*Z**(p-1)*(mu - eps*lam) + lam*sum(comb(p, j)*n**j*Z**(p-j) for j in range(2, p+1))
        a = {1: p*mu + comb(p, 2)*B}
        for i in range(2, p+1): a[i] = mu*comb(p, i) + B*(comb(p, i+1) if i+1 <= p else 0)
        rhs = -p*eps*b*Z**p + sum(a[i]*Z**(p-i) for i in range(1, p+1))
        tot += 1; viol += A > rhs*(1 + 1e-9) + 1e-9
report("Lemma 1' pathwise bound A Z^p <= -p eps b Z^p + sum a_{p,i} Z^{p-i}", viol == 0, f"({tot} cases, {viol} violations)")

# 2. Lemma 4: A e^{theta Z} <= theta e^{theta Z}(mu - eps b Z/2) for theta <= eps/e
viol = 0; tot = 0
for _ in range(20000):
    K = int(rng.integers(1, 5)); mu, beta, alpha, eps = rand_params(K); b = beta.min(); c = alpha/beta; n = c.sum()
    S = rng.exponential(rng.uniform(0.1, 20.0), K); Z = c@S; lam = mu + alpha@S; theta = rng.uniform(0.05, 1.0)*eps/math.e
    A = math.exp(theta*Z)*(theta*mu + lam*(math.exp(theta*n) - 1 - theta)); rhs = theta*math.exp(theta*Z)*(mu - eps*b*Z/2)
    tot += 1; viol += A > rhs*(1 + 1e-9) + 1e-12
report("Lemma 4 pathwise inequality (exponential moment)", viol == 0, f"({tot} cases, {viol} violations)")

# 3. Lemma G sandwich (b-gamma) Z <= lam - mu <= (B-gamma) Z, gamma root
viol = 0; tot = 0
from scipy.optimize import brentq
for _ in range(3000):
    K = int(rng.integers(1, 5)); mu, beta, alpha, eps = rand_params(K)
    f = lambda g: (alpha/(beta - g)).sum() - 1.0
    g = brentq(f, 0.0, beta.min()*(1 - 1e-12))
    om = alpha/(beta - g); S = rng.exponential(rng.uniform(0.1, 20), K); Z = om@S; d = alpha@S
    tot += 1; viol += not ((beta.min() - g)*Z <= d*(1 + 1e-9) + 1e-12 and d <= (beta.max() - g)*Z*(1 + 1e-9) + 1e-12)
    assert abs(om.sum() - 1) < 1e-9 and abs(g - eps/((om/beta).sum())) < 1e-7*max(1, g), "gamma = eps/nu~"
report("Lemma G: (b-gamma)Z <= lam-mu <= (B-gamma)Z, gamma = eps/nu~", viol == 0, f"({tot} cases, {viol} violations)")

# 4. Thm T elementary bounds on d(y), E1
def h_c(y): return np.exp(1j*y) - 1 - 1j*y
ys = np.concatenate([np.geomspace(1e-3, 1, 300), np.linspace(1, 300, 60000)]); worst_d = 0; viol_E1 = 0
for r in (1e-6, 1e-3, 0.05, 0.3, 1.0, 5.0, 50.0):
    N = (1 - np.cos(ys) - ys**2/2) + 1j*(ys - np.sin(ys)); E1 = -(1 - np.cos(ys)) - 1j*((1 + r)*ys - np.sin(ys)); d = N/(E1*(r - 0.5j*ys))
    sm = ys < 0.2
    # small-y series for stability
    yy = ys[sm]; N[sm] = (-yy**4/24 + yy**6/720) + 1j*(yy**3/6 - yy**5/120); E1[sm] = -(yy**2/2 - yy**4/24 + yy**6/720) - 1j*(r*yy + yy**3/6 - yy**5/120); d = N/(E1*(r - 0.5j*ys))
    worst_d = max(worst_d, np.abs(d).max())
    m1 = ys <= 1; viol_E1 += np.any(np.abs(E1[m1]) < 0.4583*ys[m1]**2*(1 - 1e-9)) + np.any(np.abs(E1[~m1]) < 0.158*ys[~m1]*(1 - 1e-9))
    m5 = ys >= 5; viol_E1 += np.any(np.abs(E1[m5])**2 > (1 + r)**2*ys[m5]**2*(1 + 3/ys[m5])*(1 + 1e-9))
report("Thm T elementary: sup|d| <= 19, |E1| >= 0.4583 y^2 (y<=1), >= 0.158 y (y>=1), |E1|^2 <= (1+r)^2 y^2(1+3/y) (y>=5)", worst_d <= 19 and viol_E1 == 0, f"(sup|d| = {worst_d:.4f}, violations {viol_E1})")

# 5. Theorem A step 2 bounds on D_eps
def hh(x): return np.expm1(-x) + x
viol = 0; worstmax = 0; worstmin = 1e9
xs = np.concatenate([np.geomspace(1e-3, 1, 500), np.linspace(1, 400, 40000)])
for eps in (1e-4, 1e-2, 0.1, 0.25, 0.3, 0.5, 0.9):
    r = eps/(1 - eps); alpha = 1 - eps; s = xs/alpha                       # beta = 1
    # D_eps * alpha as function of x = alpha*s  (stable for small x via series of h)
    hx = np.where(xs < 0.1, xs**2/2 - xs**3/6 + xs**4/24 - xs**5/120, hh(xs))
    num = xs**2/2 - hx; den = (r*xs + hx)*(r + xs/2)                       # alpha*D = num/den after dividing alpha factors
    D_alpha = np.where(xs < 0.1, (xs**3/6 - xs**4/24 + xs**5/120)/den, num/den)
    worstmax = max(worstmax, D_alpha.max())
    if eps <= 0.25: m = xs >= 1; worstmin = min(worstmin, D_alpha[m].min())
report("Lemma 4.3: alpha D_eps <= 2.08 for all s; >= 0.166 for s >= 1/alpha (eps <= 1/4)", worstmax <= 2.08 and worstmin >= 0.166, f"(max alpha D = {worstmax:.4f}, min over s>=1/alpha = {worstmin:.4f})")

# 6. Lemma F2 positive-part trick and det Q
viol = 0
for _ in range(100000):
    x, y = rng.normal(size=2)*rng.uniform(0.1, 10, 2); viol += (x + y)**2 < 0.5*x**2 - y**2 - 1e-9
report("elementary (x+y)^2 >= x^2/2 - y^2", viol == 0, f"({viol} violations)")
viol = 0
for _ in range(20000):
    c1, c2 = np.sort(rng.uniform(0, 5, 2)); d = rng.uniform(0.01, 2);
    if (c2 - c1) > 2*math.sqrt(2)*d:
        Q = np.array([[1, (c1 + c2)/2], [(c1 + c2)/2, (c1**2 + c2**2)/2 - 2*d**2]]); viol += np.linalg.eigvalsh(Q)[0] <= 0
        v = rng.normal(size=2); Lval = 0
        t = sum(max(0.0, 0.5*(v[0] + v[1]*c)**2 - v[1]**2*d**2) for c in (c1, c2)); viol += t < v@Q@v - 1e-9
report("Lemma F2: Q positive definite iff gap > 2 sqrt2 d; sum of positive parts >= v'Qv", viol == 0, f"({viol} violations)")

# 7. Bernstein / MGF inequalities
ysd = np.linspace(-30, 1, 200001)
report("e^y <= 1 + y + (e-2) y^2 for y <= 1", np.all(np.exp(ysd) <= 1 + ysd + (math.e - 2)*ysd**2 + 1e-12))
thw = np.linspace(0, 1, 100001)
report("e^{-t} - 1 + t <= (e-2) t^2 for 0 <= t <= 1 (Lemma E4 use, with t = theta*w)", np.all(np.exp(-thw) - 1 + thw <= (math.e - 2)*thw**2 + 1e-14))
viol = 0
for _ in range(2000):
    b = rng.uniform(0.1, 5); vals = rng.uniform(0, b, 6)*rng.integers(0, 2, 6) ; pr = rng.dirichlet(np.ones(6)); th = rng.uniform(0.01, 1)/b
    Eell = pr@vals; lhs = pr@np.exp(-th*vals); rhs = math.exp(-th*(1 - th*b/2)*Eell); viol += lhs > rhs*(1 + 1e-12)
report("Lemma E5c MGF inequality E e^{-theta l} <= exp(-theta (1 - theta b/2) E l)", viol == 0, f"({viol} violations)")
# psi expansion
for dl in (0.05, 0.1, 0.25):
    r = np.concatenate([np.linspace(1 - dl, 1 - 1e-6, 4000), np.linspace(1 + 1e-6, 1 + dl, 4000)]); ratio = (r - 1 - np.log(r))/((r - 1)**2/r)
    report(f"F2: psi(r) <= (1/2 + delta)(r-1)^2/r on [1-delta,1+delta], delta={dl}", ratio.max() <= 0.5 + dl, f"(max ratio {ratio.max():.4f})")
# r(x) concavity step
viol = 0
for _ in range(20000):
    xb = rng.uniform(1e-3, 5); d = int(rng.integers(1, 5)); y = rng.uniform(0, xb, d)
    viol += y.sum() > (xb/math.log(1 + xb))*np.log1p(y).sum()*(1 + 1e-12)
report("Lemma 7c: sum y_j <= r(xbar) sum log(1+y_j) for y_j in [0,xbar]", viol == 0, f"({viol} violations)")
# tangent inequality
lam = rng.uniform(1e-3, 50, 100000); l0 = rng.uniform(1e-2, 50, 100000)
report("tangent minorant 1/lam >= 2/l0 - lam/l0^2", np.all(1/lam >= 2/l0 - lam/l0**2 - 1e-12))
# Gamma tail (burn-in part of Theorem HP)
bad = 0; N = 400000
for a in (3, 20, 200):
    g = rng.gamma(a, 1.0, N)
    for x in (0.5, 1.0, 2.0, 4.0):
        p = np.mean(g - a >= math.sqrt(2*a*x) + x); bad += p > math.exp(-x)*(1 + 0.05)
report("Gamma(a) tail P(G - a >= sqrt(2ax) + x) <= e^{-x}", bad == 0, f"(violations {bad})")
# path-length bounds
viol = 0
for _ in range(2000):
    T = 400; eta = rng.uniform(0.005, 0.5); d = rng.normal(size=T)*rng.uniform(0, 1)*(rng.random(T) < rng.uniform(0.01, 1)); e = np.zeros(T); x = 0.0
    for t in range(T): x = (1 - eta)*x - d[t]; e[t] = x
    V = np.abs(d).sum(); viol += np.mean(e**2) > np.mean(d**2)/eta**2*(1 + 1e-9) or np.mean(e**2) > V**2/(eta*T)*(1 + 1e-9)
report("Prop path-length: mean e^2 <= mean d^2/eta^2 and <= V^2/(eta T)", viol == 0, f"({viol} violations)")
# Gaussian max entropy / AM-GM on determinant (lower-bound proof)
viol = 0
for _ in range(2000):
    d = int(rng.integers(1, 5)); Aa = rng.normal(size=(d, d)); C = Aa@Aa.T + 0.01*np.eye(d); mean = rng.normal(size=d)
    hG = 0.5*np.log(((2*math.pi*math.e)**d)*np.linalg.det(C)); bound = 0.5*d*np.log(2*math.pi*math.e*(np.trace(C) + mean@mean)/d)
    viol += hG > bound + 1e-12
report("entropy step: (1/2)log((2 pi e)^d det C) <= (d/2) log(2 pi e (tr C + |m|^2)/d)", viol == 0, f"({viol} violations)")
print(f"\n{len(rows)} inequality families checked; failures: {[r[0] for r in rows if not r[1]]}")
