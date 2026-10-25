"""Lemma IM (inverse moments, any K).  f_a = W^-a + psi W^-(a+2),  W = Z + z (Z = omega.S, omega_k = alpha_k/(beta_k-gamma), sum omega = 1),
 psi = kappa.S,  kappa_k = ((a(a+1)/2) alpha_k - c_a omega_k)/beta_k,  c_a = a(a+1) c0,  c0 = n/(2 nu~),  sum kappa = 0.
 Claim (pathwise):  A f_a <= -[a(mu-(a+1)c0) - b kbar gamma - C_a/z] W^-(a+1) + a gamma W^-a,  b=a+2,  kbar = max|kappa_k|/omega_k,
 C_a = a(a+1) mu/2 + b kbar mu + B_z kbar b(b+1)/2,  B_z = max(mu/z, beta_max - gamma).   Test on random states, K=1..4, a in (0,1]."""
import numpy as np, math
from scipy.optimize import brentq
rng = np.random.default_rng(31)
def setup(mu, beta, alpha):
    n = (alpha/beta).sum(); g = brentq(lambda x: (alpha/(beta-x)).sum()-1, 0, beta.min()*(1-1e-12)); om = alpha/(beta-g); nut = (om/beta).sum(); return n, g, om, nut, n/(2*nut)
def Af(S, a, z, mu, beta, alpha, om, kap):
    lam = mu + alpha@S; Z = om@S; W = Z + z; psi = kap@S
    cont_W = -(om*beta)@S; cont_psi = -(kap*beta)@S                   # drift of W and of psi between events
    b = a + 2
    # f = W^-a + psi W^-b ; jump: W -> W+1 (Omega=1), psi unchanged (sum kappa = 0)
    dfc = -a*W**(-a-1)*cont_W + W**(-b)*cont_psi - b*psi*W**(-b-1)*cont_W
    jump = lam*(((W+1)**(-a) - W**(-a)) + psi*((W+1)**(-b) - W**(-b)))
    return dfc + jump, W, lam
viol = 0; N = 0; worst = np.inf; checked = []
for trial in range(60):
    K = int(rng.integers(1, 5)); beta = np.sort(rng.uniform(0.3, 6, K)); w = rng.dirichlet(np.ones(K)); eps = 10**rng.uniform(-2.2, -0.3); mu = rng.uniform(0.2, 4); alpha = (1-eps)*w*beta
    n, g, om, nut, c0 = setup(mu, beta, alpha); kk = mu/c0
    if kk <= 1.05: continue
    a = float(rng.uniform(0.05, 1.0))
    if a + 1 >= kk: a = 0.5*(kk-1)
    z = 10**rng.uniform(0.5, 3); ca = a*(a+1)*c0; kap = ((a*(a+1)/2)*alpha - ca*om)/beta; assert abs(kap.sum()) < 1e-9*max(1, np.abs(kap).max())
    kbar = np.max(np.abs(kap)/om); b = a + 2; Bz = max(mu/z, beta.max()-g); Ca = a*(a+1)*mu/2 + b*kbar*mu + Bz*kbar*b*(b+1)/2
    for _ in range(4000):
        S = rng.exponential(1, K)*10**rng.uniform(-3, 4, K)*rng.random(K)
        val, W, lam = Af(S, a, z, mu, beta, alpha, om, kap)
        bound = -(a*(mu-(a+1)*c0) - b*kbar*g - Ca/z)*W**(-a-1) + a*g*W**(-a)
        sl = (bound - val)*W**(a+2); worst = min(worst, sl); N += 1
        if sl < -1e-6: viol += 1
    checked.append((K, round(eps,4), round(kk,2), round(a,2)))
print(f"Lemma IM pathwise: {N} states over {len(checked)} parameter sets (K=1..4, k_K in (1.05, ..), a in (0,1], a+1<k_K): violations {viol}; min slack*W^(a+2) = {worst:.3e}")
print("examples (K, eps, k_K, a):", checked[:6])
