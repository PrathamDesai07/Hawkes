import numpy as np, math
from p11_corrector_check import constants
rng = np.random.default_rng(777)

def Af_stable(S, c):
    """Exact generator of f=log W + psi/W^2 written with log1p and rational differences (no cancellation of O(lam) terms)."""
    mu, beta, alpha, om, kap, z0 = c['mu'], c['beta'], c['alpha'], c['om'], c['kap'], c['z0']
    lam = mu + alpha@S; Z = om@S; W = Z+z0; psi = kap@S; Ks = kap.sum()      # Ks = 0 by construction
    dZ = -(om*beta)@S; dpsi = -(kap*beta)@S
    cont = dZ/W + dpsi/W**2 - 2*psi*dZ/W**3
    Om = om.sum()
    jump = lam*(math.log1p(Om/W) + (psi+Ks)/(W+Om)**2 - psi/W**2)
    return cont + jump, W

print("E'. Lemma H, stable evaluation (slack scaled by W^2), 40 random parameter sets, K=1..4")
nbad = 0; worst = [np.inf, np.inf]; badinfo = []
for trial in range(40):
    K = int(rng.integers(1, 5)); beta = np.sort(rng.uniform(0.3, 8, K)); w = rng.dirichlet(np.ones(K)); eps = 10**rng.uniform(-2.3, -0.3); mu = rng.uniform(0.2, 4); z0 = 10**rng.uniform(0, 3)
    alpha = (1-eps)*w*beta; c = constants(mu, beta, alpha, z0)
    for _ in range(20000):
        S = rng.exponential(1, K)*10**rng.uniform(-3, 4, K)*rng.random(K)
        a, W = Af_stable(S, c)
        LB = (mu-c['c0']-2*c['kb']*c['gam'])/W - c['gam'] - (mu/2+c['C2'])/W**2
        UB = (mu-c['c0']+c['gam']*z0+2*c['kb']*c['gam'])/W - c['gam'] + (c['c0']*z0+c['C4'])/W**2
        sl, sh = (a-LB)*W**2, (UB-a)*W**2; worst[0] = min(worst[0], sl); worst[1] = min(worst[1], sh)
        if sl < -1e-6 or sh < -1e-6:
            nbad += 1
            if len(badinfo) < 3: badinfo.append((trial, K, round(eps,4), round(mu,3), round(z0,1), round(W,1), sl, sh))
print(f"   violations: {nbad}; min slack*W^2: lower {worst[0]:.3e}, upper {worst[1]:.3e}")
print("   first violating cases (trial,K,eps,mu,z0,W,slack_lo,slack_hi):", badinfo)

print("G'. Problem 3 Lemma 4 pathwise (e^{theta Z} divided out): theta*mu + lam*(e^{theta n}-1-theta) <= theta*(mu - eps*b*Z/2), theta=eps/e")
bad = 0; N = 0; worst = np.inf
for _ in range(30):
    K = int(rng.integers(1, 5)); beta = np.sort(rng.uniform(0.3, 6, K)); w = rng.dirichlet(np.ones(K)); eps = rng.uniform(0.01, 0.6); mu = rng.uniform(0.2, 3)
    alpha = (1-eps)*w*beta; cc = alpha/beta; n = cc.sum(); b = beta.min(); th = eps/math.e
    for _ in range(5000):
        S = rng.exponential(1, K)*10**rng.uniform(-3, 4, K)*rng.random(K); Z = cc@S; lam = mu + alpha@S
        lhs = th*mu + lam*math.expm1(th*n) - lam*th; rhs = th*(mu - eps*b*Z/2); N += 1; worst = min(worst, rhs-lhs)
        if lhs > rhs + 1e-12*max(1, abs(lhs)): bad += 1
print(f"   violations {bad} of {N}; min (rhs-lhs) = {worst:.3e}")
