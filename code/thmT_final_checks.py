"""Remaining checks for Theorem T: (B) pointwise head/tail inequalities used, (D) KL(Gamma_theta||Gamma_theta') <= C_K delta^2 and the assembled TV bound."""
import numpy as np, math
from scipy import special
from thmT_constants import *
k0 = 2.0
for eps in (1e-2, 1e-3, 1e-4):
    th, x, d, z, G = true_h1(k0, eps); r = eps/(1-eps); k = k0/(1-eps)**2; alpha = 1-eps; a = alpha**2/2; b = (1-eps)/(2*eps)
    M1 = 19*k*alpha/2; ts = x/(alpha*eps); t1 = 1/(M1*eps); ph = G*np.exp(z)
    head = ts <= t1
    f1 = np.all(np.abs(np.expm1(z[head])) <= (math.e-1)*M1*eps*ts[head]*(1+1e-9))
    Gp = 1j*k*b*(1-1j*b*x)**(-k-1)
    Dlp = Gp*np.expm1(z) + G*(1j*(k/2)*d)*np.exp(z)
    f2 = np.all(alpha*eps*np.abs(Dlp[head]) <= M1*eps*(k*(math.e-1)+math.e)*np.abs(G[head])*(1+1e-9))
    Lp = (k/2)*x/E1_fun(x, r)                              # dlog phi/dx ; per dt multiply by alpha*eps
    Lt = np.abs(Lp)*alpha*eps
    m1 = (x <= 1) & (ts >= t1); m2 = x >= 1
    f3 = np.all(Lt[m1] <= 1.1*k/ts[m1]*(1+1e-9)); f4 = np.all(Lt[m2] <= 3.2*k*alpha*eps*(1+1e-9))
    print(f"eps={eps:g}: (B1) |F~-1|<=(e-1)M1 eps t on head {f1}; (B2) |Delta'|<=M1 eps(k(e-1)+e)|G| {f2}; (B3) |L'|<=1.1k/t (x<=1) {f3}; (B4) |L'|<=3.2k alpha eps (x>=1) {f4}; max|Lt|*t/k on x<=1: {np.max(Lt[m1]*ts[m1]/k):.3f}")
# (D) KL for the ridge pair, beta = 1, mu = 1
def kl(k1, s1, k2, s2):   # KL(Gamma(k1,scale s1) || Gamma(k2,scale s2))
    return (k1-k2)*special.digamma(k1) - special.gammaln(k1) + special.gammaln(k2) + k2*(math.log(s2)-math.log(s1)) + k1*(s1-s2)/s2
for eps in (1e-2, 1e-4):
    for mu in (1.0,):
        al = 1-eps; k = 2*mu/al**2; sc = al**2/2/eps
        row = []
        for dl in (0.01, 0.05, 0.1, 0.2, 0.3):
            mu2 = mu+dl; al2 = al*(1-eps*dl/mu); eps2 = 1-al2; k2 = 2*mu2/al2**2; sc2 = al2**2/2/eps2
            K = kl(k, sc, k2, sc2); row.append((dl, round(K/dl**2, 3), round(math.sqrt(K/2), 4)))
        print(f"eps={eps:g}: (delta, KL/delta^2, Pinsker TV bound):", row)
