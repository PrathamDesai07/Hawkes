"""Lemma E0: sup_t E exp(theta S_t) < inf for theta < theta_*, where (e^theta-1)/theta = beta/alpha.  Generator algebra:
 L e^{theta s} = e^{theta s}[ mu(e^theta-1) - kappa_theta s ],  kappa_theta = beta theta - alpha (e^theta-1) > 0.
Check: (i) theta_* by bisection; (ii) Monte Carlo E exp(theta S_t) over time and its comparison with the Foster-Lyapunov bound."""
import math, random, sys
def theta_star(beta, al):
    lo, hi = 1e-9, 5.0
    f = lambda t: beta*t - al*(math.exp(t)-1)
    # f>0 on (0,theta_*), f(theta_*)=0 ; find the positive root
    t = 1e-6
    while f(t) > 0: t *= 1.01
    return t
def sim(seed, mu, beta, al, T, theta, dts):
    rng = random.Random(seed); s = 0.0; t = 0.0; out = []; k = 0
    while k < len(dts):
        lam0 = mu + al*s; w = 0.0
        while True:
            w += rng.expovariate(lam0)
            if rng.random()*lam0 <= mu + al*s*math.exp(-beta*w): break
        # record the state at checkpoints that fall inside the gap
        while k < len(dts) and t + w >= dts[k]:
            out.append(s*math.exp(-beta*(dts[k]-t))); k += 1
        s = s*math.exp(-beta*w) + 1.0; t += w
    return out
mu, beta = 1.0, 1.0
for eps in (0.3, 0.1):
    al = (1-eps)*beta; ts = theta_star(beta, al); print(f"eps={eps}: theta_* = {ts:.4f} (2 eps/(1-eps)~ {2*eps/(1-eps):.3f})")
    for frac in (0.5, 0.8):
        th = frac*ts; kap = beta*th - al*(math.exp(th)-1); s0 = 2*mu*(math.exp(th)-1)/kap
        cps = [x/eps for x in (0.5, 1, 2, 4, 8, 16, 32)]
        reps = 3000 if eps == 0.3 else 1500; acc = [0.0]*len(cps)
        for r in range(reps):
            v = sim(7000+r, mu, beta, al, cps[-1], th, cps)
            for j, x in enumerate(v): acc[j] += math.exp(th*x)
        print(f"  theta={frac}*theta_*={th:.3f}: E exp(theta S_t) at t*eps = .5,1,2,4,8,16,32 :", [round(a/reps, 2) for a in acc], f"  (Lyapunov threshold s0 = {s0:.1f})", flush=True)
