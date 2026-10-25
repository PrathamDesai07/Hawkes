"""Problem 12 (Theorem 6): smallest eigenvector -> ridge direction.  K = 2, beta = (1, 2), mu = 1.5, alpha_k = (1-eps) beta_k / 2.
   rigorous pieces:  rho_eps <= (1+eps)^2 m_bound - (eps/mu)(1+2 eps)  divided by |v_r|^2, with m_bound = Theorem E explicit bound on E[1/lam];
                     g = lam_2(N(lam0)) from the tangent bound (best lam0);   sin^2 <= rho_eps / g.
   For comparison: MC eigenvector overlaps of Problem 12 (0.85, 0.99, 0.9992 at eps = 0.1, 0.03, 0.01) and the exact eigen-decomposition of the MC-free tangent matrix."""
import numpy as np
from complexity_bench import moments_fast
from tangent_fast import matrices_fast

def modes(mu, beta, alpha):
    K = len(beta); Am = np.diag(beta) - np.outer(np.ones(K), alpha); ev = np.sort(np.linalg.eig(Am)[0].real)
    L = np.array([alpha/(beta - aj) for aj in ev]); return ev, L

def thmE_bound(mu, beta, alpha, eps):
    K = len(beta); a, L = modes(mu, beta, alpha); gam = a[0]; om = L[0]; b, B = beta.min(), beta.max()
    nut = (om/beta).sum(); n = (alpha/beta).sum(); c0 = n/(2*nut); d0 = mu - c0
    if d0 <= 0: return None
    kap = (c0*om - alpha/2)/beta; kb = np.max(np.abs(kap)/om)
    if 2*kb*gam > d0/4: return None
    best = None
    for z in np.geomspace(1e-2, 1e4, 600):
        Bz = max(mu/z, B - gam); C2 = 2*kb*mu + 3*Bz*kb
        if z < 4*(mu/2 + C2)/d0: continue
        val = max(z/mu, 1/(b - gam))*2*gam/d0
        if best is None or val < best: best = val
    return best

def tangent_gap(mu, alpha, beta):
    M, M3, El = matrices_fast(mu, alpha, beta); best = None
    for t in np.linspace(0.05, 3.0, 300)/El:
        ev = np.linalg.eigvalsh(2*t*M - t*t*M3)
        if best is None or ev[1] > best[1]: best = (ev[0], ev[1], t)
    return best

mu = 1.5; beta = np.array([1.0, 2.0]); mc_overlap = {0.1: 0.85, 0.03: 0.99, 0.01: 0.9992}
print(f"{'eps':>6} {'k_K':>5} | {'m bound (Thm E)':>15} {'rho bound':>10} | {'g=lam2(N)':>10} | {'sin^2 bound':>11} | {'MC overlap -> sin^2':>20}")
for eps in (0.1, 0.03, 0.01, 0.003, 0.001, 0.0003):
    alpha = (1 - eps)*beta/2; a, L = modes(mu, beta, alpha); nut = (L[0]/beta).sum(); kK = mu*2*nut/((alpha/beta).sum())
    mb = thmE_bound(mu, beta, alpha, eps)
    vr = np.array([1.0, *(-eps*alpha/mu)]); nv2 = vr@vr
    g = tangent_gap(mu, alpha, beta)[1]
    if mb is None: print(f"{eps:6.4f} {kK:5.2f} | Theorem E conditions not met"); continue
    rho = ((1 + eps)**2*mb - (eps/mu)*(1 + 2*eps))/nv2
    s2 = rho/g
    mcs = f"{mc_overlap[eps]} -> {1 - mc_overlap[eps]**2:.4f}" if eps in mc_overlap else "-"
    print(f"{eps:6.4f} {kK:5.2f} | {mb:15.4g} {rho:10.4g} | {g:10.4g} | {min(s2, 9.99):11.4g} | {mcs:>20}")
