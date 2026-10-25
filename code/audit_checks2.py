import numpy as np, math
from p11_corrector_check import constants, Af
from p10_fisher_k1_exact import h
rng = np.random.default_rng(777)

print("E. Problem 11 Lemma H pathwise (generator between the two bounds) on 40 random parameter sets, K=1..4, eps in [0.005,0.5], 20000 states each")
worst_lo = np.inf; worst_hi = np.inf; nbad = 0
for _ in range(40):
    K = int(rng.integers(1, 5)); beta = np.sort(rng.uniform(0.3, 8, K)); w = rng.dirichlet(np.ones(K)); eps = 10**rng.uniform(-2.3, -0.3); mu = rng.uniform(0.2, 4); z0 = 10**rng.uniform(0, 3)
    alpha = (1-eps)*w*beta; c = constants(mu, beta, alpha, z0)
    for _ in range(20000):
        S = rng.exponential(1, K)*10**rng.uniform(-3, 5, K)*rng.random(K)
        a, W = Af(S, c)
        LB = (mu-c['c0']-2*c['kb']*c['gam'])/W - c['gam'] - (mu/2+c['C2'])/W**2
        UB = (mu-c['c0']+c['gam']*z0+2*c['kb']*c['gam'])/W - c['gam'] + (c['c0']*z0+c['C4'])/W**2
        sl, sh = (a-LB)*W**2, (UB-a)*W**2          # scaled slack
        worst_lo = min(worst_lo, sl); worst_hi = min(worst_hi, sh)
        if sl < -1e-9 or sh < -1e-9: nbad += 1
print(f"   violations: {nbad}; min slack*W^2 lower bound {worst_lo:.2e}, upper bound {worst_hi:.2e}")

print("F. Problem 10 Theorem A constant: Den <= alpha*x^2 for x>=1 and D_eps >= 0.132/alpha, for eps<=1/5 and at eps=1/4")
for eps in (0.2, 0.25):
    worst_den = 0; minD = np.inf
    for beta in (0.5, 1.0, 3.0):
        al = (1-eps)*beta; x = np.linspace(1, 400, 200000); s = x/al
        den = (eps*beta*s + h(x))*(eps*beta + al*x/2); worst_den = max(worst_den, (den/(al*x*x)).max())
        D = (x*x/2 - h(x))/den; minD = min(minD, (D*al).min())
    print(f"   eps={eps}: max Den/(alpha x^2) = {worst_den:.4f} (claim <= 1);  min alpha*D_eps over x>=1 = {minD:.4f} (claim >= 0.132)")

print("G. Problem 3 Lemma 4 pathwise: A e^{theta Z} <= theta e^{theta Z}(mu - eps*b*Z/2), theta=eps/e")
bad = 0; N = 0
for _ in range(30):
    K = int(rng.integers(1, 5)); beta = np.sort(rng.uniform(0.3, 6, K)); w = rng.dirichlet(np.ones(K)); eps = rng.uniform(0.01, 0.6); mu = rng.uniform(0.2, 3)
    alpha = (1-eps)*w*beta; cc = alpha/beta; n = cc.sum(); b = beta.min(); th = eps/math.e
    for _ in range(5000):
        S = rng.exponential(1, K)*10**rng.uniform(-3, 4, K)*rng.random(K); Z = cc@S; lam = mu + alpha@S
        lhs = math.exp(th*Z)*(th*mu + lam*(math.exp(th*n)-1-th))        # exact generator (drift of Z is -(lam-mu), jump n at rate lam)
        rhs = th*math.exp(th*Z)*(mu - eps*b*Z/2); N += 1
        if lhs > rhs*(1+1e-9) + 1e-12: bad += 1
print(f"   violations {bad} of {N}")
