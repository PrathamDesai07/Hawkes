"""Independent audit checks of the 'Done' problems (fresh random parameters)."""
import numpy as np, math
from scipy.optimize import brentq
from scipy.linalg import expm
from p10_fisher_k1_exact import exact_m
rng = np.random.default_rng(2026)

print("A. Problem 11 Prop 1: eff = (m - eps/mu)/(1-2eps+eps*mu*m) vs 1/(I^-1)_mumu from the exact K=1 matrix")
worst = 0
for _ in range(8):
    mu = rng.uniform(0.2, 3); beta = rng.uniform(0.5, 3); eps = rng.uniform(0.02, 0.5)
    m = exact_m(mu, beta, eps); a = (1-eps)*beta; El = mu/eps
    I = np.array([[m, (1-mu*m)/a], [(1-mu*m)/a, (El-2*mu+mu**2*m)/a**2]])
    eff_mat = 1/np.linalg.inv(I)[0, 0]; eff_f = (m-eps/mu)/(1-2*eps+eps*mu*m)
    worst = max(worst, abs(eff_mat-eff_f)/eff_mat)
print(f"   max relative discrepancy over 8 random (mu,beta,eps): {worst:.2e}")

print("B. Problem 10 / 11: Jensen floor m >= eps/mu, and Theorem C bound m <= 2 eps/(2mu - n^2 beta) when 2mu > n^2 beta (K=1, exact m)")
viol = 0; tight = 0.0; N = 0
for _ in range(25):
    mu = rng.uniform(0.2, 3); beta = rng.uniform(0.5, 3); eps = rng.uniform(0.02, 0.5); n = 1-eps
    m = exact_m(mu, beta, eps)
    if m < eps/mu*(1-1e-6): viol += 1
    if 2*mu > n*n*beta:
        N += 1; b = 2*eps/(2*mu-n*n*beta)
        if m > b*(1+1e-6): viol += 1
        tight = max(tight, m/b)
print(f"   violations: {viol} (of 25 floor checks + {N} bound checks); max m/bound = {tight:.3f}")

print("C. Problem 10 Prop (any K): E lambda_t <= mu(1 + beta_max t) from a cold start (mean ODE), random K=1..4 parameter sets")
bad = 0
for _ in range(40):
    K = rng.integers(1, 5); beta = np.sort(rng.uniform(0.3, 6, K)); w = rng.dirichlet(np.ones(K)); eps = rng.uniform(0.01, 0.5); mu = rng.uniform(0.2, 3)
    alpha = (1-eps)*w*beta; M = -np.diag(beta) + np.outer(np.ones(K), alpha)
    for t in (0.1, 1, 5, 20, 100, 500):
        # E S_t = int_0^t e^{M s} mu 1 ds = M^{-1}(e^{Mt}-I) mu 1
        ES = np.linalg.solve(M, (expm(M*t) - np.eye(K))@(mu*np.ones(K)))
        lam = mu + alpha@ES
        if lam > mu*(1+beta.max()*t)*(1+1e-9): bad += 1
print(f"   violations: {bad} of 240")

print("D. Problem 3 Lemma 1/1': E Z^2 <= bound and E Z <= mu(1-eps)/(eps*b); gamma = eps/nu_tilde exactly (random K=1..4)")
from p3_moments_mixing import moments
badg = 0; bad2 = 0
for _ in range(30):
    K = rng.integers(1, 5); beta = np.sort(rng.uniform(0.3, 6, K)); w = rng.dirichlet(np.ones(K)); eps = rng.uniform(0.02, 0.5); mu = rng.uniform(0.2, 3)
    alpha = (1-eps)*w*beta; n = 1-eps
    g = brentq(lambda x: (alpha/(beta-x)).sum()-1, 0, beta.min()*(1-1e-12)); om = alpha/(beta-g); nut = (om/beta).sum()
    if abs(g - eps/nut) > 1e-8*max(1, g): badg += 1
    m1, M2 = moments(mu, beta, alpha); c = alpha/beta; EZ2 = c@M2@c; b, B = beta.min(), beta.max()
    bound = mu/(2*b*eps**2)*((2*mu + n*n*B)/b + n*n*eps)
    if EZ2 > bound*(1+1e-8): bad2 += 1
print(f"   gamma=eps/nu_tilde violations {badg}/30;  E Z^2 bound violations {bad2}/30")
