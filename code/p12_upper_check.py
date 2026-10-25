"""Problem 12, coding stage (after problems12_upper_bounds.md):
 (1) sympy: Lemma 0 (A y_j = -a_j y_j + mu, each event adds 1) and Lemma 1 (generator of y_j^2/W) for K = 2 with symbolic beta_1, beta_2, gamma, a_2;
 (2) random-state test of the pathwise inequality of Lemma 2 (K = 2, 3, 4);
 (3) constants: C_Y, C_mid, rigorous upper bound on lam_{K-1}(I) vs Monte Carlo, for several eps (eps-independence);
 (4) eigenvector: sin^2 bound rho_eps/g vs measured overlaps."""
import numpy as np, sympy as sp, random, math, itertools
from scipy.optimize import brentq

# ---------- (1) symbolic ----------
b1, b2, g, a2, mu, s1, s2, z = sp.symbols('b1 b2 g a2 mu s1 s2 z', positive=True)
al1, al2 = sp.symbols('al1 al2')
sol = sp.solve([al1/(b1 - g) + al2/(b2 - g) - 1, al1/(b1 - a2) + al2/(b2 - a2) - 1], [al1, al2], dict=True)[0]
A1, A2 = sol[al1], sol[al2]
lam = mu + A1*s1 + A2*s2
def Agen(f): return -b1*s1*sp.diff(f, s1) - b2*s2*sp.diff(f, s2) + lam*(f.subs({s1: s1 + 1, s2: s2 + 1}, simultaneous=True) - f)
y1 = A1*s1/(b1 - g) + A2*s2/(b2 - g); y2 = A1*s1/(b1 - a2) + A2*s2/(b2 - a2)
random.seed(5)
def rr(lo, hi): return sp.Rational(random.randint(int(lo*100), int(hi*100)), 100)
def worst(expr, trials=12):
    w = 0
    for _ in range(trials):
        b1v = rr(0.5, 1.5); b2v = b1v + rr(0.5, 3); gv = rr(0.02, 0.3)*b1v; a2v = b1v + rr(0.1, 0.9)*(b2v - b1v)
        sub = {b1: b1v, b2: b2v, g: gv, a2: a2v, mu: rr(0.3, 3), s1: rr(0.1, 9), s2: rr(0.1, 9), z: rr(0.5, 5)}
        w = max(w, abs(sp.N(expr.subs(sub), 40)))
    return float(w)
print("(1) symbolic checks (K = 2, symbolic beta, gamma, a_2; alpha solved from the two secular equations)")
print("    Lemma 0: A y_1 + gamma y_1 - mu          :", worst(Agen(y1) + g*y1 - mu))
print("    Lemma 0: A y_2 + a_2 y_2 - mu             :", worst(Agen(y2) + a2*y2 - mu))
print("    Lemma 0: y_1(s+1) - y_1(s) = 1, y_2(s+1) - y_2(s) = 1:", worst(y1.subs({s1: s1 + 1, s2: s2 + 1}, simultaneous=True) - y1 - 1), worst(y2.subs({s1: s1 + 1, s2: s2 + 1}, simultaneous=True) - y2 - 1))
W = y1 + z
f2 = y2**2/W
formula = -2*a2*y2**2/W + y2**2*((g*y1 - mu)/W**2 + lam/(W**2*(W + 1))) + lam/(W + 1) + 2*y2*(mu/W - lam/(W*(W + 1)))
print("    Lemma 1: A(y_2^2/W) - formula             :", worst(Agen(f2) - formula))

# ---------- helpers: parameters ----------
def modes(mu, beta, alpha):
    K = len(beta); Am = np.diag(beta) - np.outer(np.ones(K), alpha); ev, V = np.linalg.eig(Am); a = np.sort(ev.real)
    L = np.array([alpha/(beta - aj) for aj in a])          # rows l_j, sigma_j = 1
    assert np.allclose(L.sum(1), 1.0, atol=1e-8), L.sum(1)
    return a, L

# ---------- (2) random-state test of Lemma 2 ----------
rng = np.random.default_rng(11); viol = tot = 0
for _ in range(3000):
    K = int(rng.integers(2, 5)); beta = np.sort(rng.uniform(0.5, 6.0, K)); eps = rng.uniform(0.005, 0.05)
    alpha = (1 - eps)*rng.dirichlet(np.ones(K))*beta; mu = rng.uniform(0.2, 3.0)
    a, L = modes(mu, beta, alpha); gam, b, B = a[0], beta.min(), beta.max()
    aj = a[1:]; z = max(4*(B - gam)/aj.min(), 1.0); Bz = max(mu/z, B - gam); c1 = max(mu, Bz)
    if gam > aj.min()/4 or Bz/z > aj.min()/4: continue
    S = rng.exponential(rng.uniform(0.1, 50.0), K)*rng.integers(0, 2, K).clip(0, 1) if rng.random() < 0.3 else rng.exponential(rng.uniform(0.1, 50.0), K)
    y = L@S; Z = y[0]; Wv = Z + z; lam = mu + alpha@S
    for j in range(1, K):
        yj = y[j]
        Af = -2*aj[j-1]*yj**2/Wv + yj**2*((gam*Z - mu)/Wv**2 + lam/(Wv**2*(Wv + 1))) + lam/(Wv + 1) + 2*yj*(mu/Wv - lam/(Wv*(Wv + 1)))
        rhs = -aj[j-1]*yj**2/Wv + Bz + 2*c1**2/(aj[j-1]*z)
        tot += 1; viol += Af > rhs*(1 + 1e-9) + 1e-9
print(f"(2) Lemma 2 pathwise inequality: {tot} (state, mode) pairs, {viol} violations")

# ---------- (3) constants and comparison with Monte Carlo ----------
def Cmid_bound(mu, beta, alpha):
    K = len(beta); a, L = modes(mu, beta, alpha); gam, b, B = a[0], beta.min(), beta.max(); aj = a[1:]
    best = None
    for z in np.geomspace(max(1.0, 4*(B - gam)/aj.min()), 1e4, 80):
        Bz = max(mu/z, B - gam); c1 = max(mu, Bz); clam = min(mu/z, b - gam)
        if gam > aj.min()/4 or Bz/z > aj.min()/4: continue
        CY = max((Bz + 2*c1**2/(a_*z))/(a_*clam) for a_ in aj)
        Lm = L[1:].T                                    # K x (K-1), columns l_j (j >= 2)
        smin = np.linalg.svd(Lm, compute_uv=False).min()
        Cm = (K - 1)*CY/smin**2
        if best is None or Cm < best[0]: best = (Cm, CY, smin, z)
    return best, a

def mc_Ey2_over_lam(mu, beta, alpha, nev, seed):
    rng = random.Random(seed); K = len(beta); S = np.zeros(K); a, L = modes(mu, beta, alpha)
    acc = np.zeros(K); T = 0.0
    for _ in range(3000):
        lam0 = mu + alpha@S; w = 0.0
        while True:
            w += rng.expovariate(lam0)
            if rng.random()*lam0 <= mu + alpha@(S*np.exp(-beta*w)): break
        S = S*np.exp(-beta*w) + 1.0
    for _ in range(nev):
        lam0 = mu + alpha@S; w = 0.0
        while True:
            w += rng.expovariate(lam0)
            if rng.random()*lam0 <= mu + alpha@(S*np.exp(-beta*w)): break
        n = 17; u = np.linspace(0, w, n); wt = np.ones(n); wt[1:-1:2] = 4; wt[2:-1:2] = 2; wt *= (w/(n - 1))/3
        Su = S[None, :]*np.exp(-np.outer(u, beta)); y = Su@L.T; lm = mu + Su@alpha
        acc += wt@(y**2/lm[:, None]); T += w
        S = S*np.exp(-beta*w) + 1.0
    return acc/T

print("\n(3) K = 2, beta = (1, 2), mu = 1.5, alpha_k = (1-eps) beta_k / 2:  E[y_2^2/lam] (MC, 3e5 events) vs the rigorous constant C_Y, and the bound on lam_{K-1}(I)")
for eps in (0.1, 0.03, 0.01):
    beta = np.array([1.0, 2.0]); alpha = (1 - eps)*beta/2; mu = 1.5
    best, a = Cmid_bound(mu, beta, alpha)
    if best is None: print(f"    eps={eps}: conditions of Lemma 2 not met"); continue
    Cm, CY, smin, zb = best; mcv = mc_Ey2_over_lam(mu, beta, alpha, 300000, 1)
    print(f"    eps={eps}: a = {np.round(a, 4)}, MC E[y_2^2/lam] = {mcv[1]:.4f}, rigorous C_Y = {CY:.3f} (z = {zb:.1f}), s_min = {smin:.3f}, C_mid = {Cm:.3f} (MC middle eigenvalue of Problem 12: 0.064-0.066)", flush=True)
