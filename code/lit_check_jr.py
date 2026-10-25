"""Tests against Jaisson-Rosenbaum (arXiv 1310.2033), Theorem 2.2, specialised to phi(s)=beta*exp(-beta s): m = 1/beta, a_T = 1-eps, lam_J/T = eps.
 Limit CIR (rescaled time):  dX = (lam_J/m)(mu - X)dt + (sqrt(lam_J)/m) sqrt(X) dB   [my reading of the extracted formula; the alternative would be sqrt(lam_J/m)].
 Real-time:  kappa = eps*beta,  sigma^2 = eps*beta^2  for X = eps*lambda.   Stationary law of CIR: Gamma(shape = 2 kappa mu/sigma^2 = 2 mu/beta, scale sigma^2/(2 kappa) = beta/2).
 Alternative reading sqrt(lam_J/m) would give shape 2*mu instead.  Choose beta=2, mu=0.5  ->  reading A: 0.5 ; reading B: 1.0.
 (1) Var/mean^2 of lambda at stationarity (grid-sampled) vs 1/shape;   (2) autocorrelation of lambda(t) vs exp(-eps*beta*t);  (3) E N_T (empty start) vs the CIR integral."""
import math, random, sys

def gap(rng, mu, al, beta, s):
    lam0 = mu + al*s; w = 0.0
    while True:
        w += rng.expovariate(lam0)
        if rng.random()*lam0 <= mu + al*s*math.exp(-beta*w): return w

def grid_lambda(seed, mu, beta, eps, T, dt, warm_time):
    """lambda(t) on a time grid, stationary start (burn-in warm_time)."""
    al = (1-eps)*beta; rng = random.Random(seed); s = 0.0; t = 0.0
    while t < warm_time:
        w = gap(rng, mu, al, beta, s); s = s*math.exp(-beta*w) + 1.0; t += w
    out = []; t = 0.0; nxt = 0.0
    while t < T:
        w = gap(rng, mu, al, beta, s)
        while nxt < t + w and nxt < T:
            out.append(mu + al*s*math.exp(-beta*(nxt - t))); nxt += dt
        s = s*math.exp(-beta*w) + 1.0; t += w
    return out

if __name__ == "__main__":
    reps = int(sys.argv[1]); mu, beta = 0.5, 2.0
    for eps in (0.2, 0.05):
        al = (1-eps)*beta; gam = eps*beta; dt = 0.5/gam; T = 3000/gam
        k_exact = 2*beta*mu/al**2; lam_bar = mu/eps
        means = []; vars_ = []; ac = {1: [], 2: [], 4: []}
        for r in range(reps):
            x = grid_lambda(10+r, mu, beta, eps, T, dt, 50/gam)
            n = len(x); m = sum(x)/n; v = sum((y - m)**2 for y in x)/(n-1); means.append(m); vars_.append(v)
            for lag in ac: ac[lag].append(sum((x[i]-m)*(x[i+lag]-m) for i in range(n-lag))/(n-lag)/v)
        M = sum(means)/reps; V = sum(vars_)/reps
        print(f"eps={eps}: mean lambda = {M:.3f} (theory mu/eps = {lam_bar:.3f});  shape estimate mean^2/var = {M*M/V:.3f}  [exact 2 beta mu/alpha^2 = {k_exact:.3f};  JR reading A 2mu/beta = {2*mu/beta:.3f};  reading B = {2*mu:.3f}]")
        for lag in ac:
            rho = sum(ac[lag])/reps
            print(f"      autocorr at lag {lag}*dt = {lag*0.5:.1f}/gamma:  empirical {rho:.3f}   exp(-gamma t) = {math.exp(-gam*lag*dt):.3f}")
    print("(3) E N_T for an empty start vs the CIR integral  mu - mu*(m/lam_J)*(1-exp(-lam_J/m)),  m=1/beta, lam_J = T*eps")
    for eps in (0.05, 0.01):
        for lamJ in (0.5, 2.0, 8.0):
            T = lamJ/eps; gam = eps*beta; ENT = mu*T/eps - mu*(1-eps)/eps*(1-math.exp(-gam*T))/gam
            V1 = eps*ENT/T; cir = mu - mu*(1/beta/lamJ)*(1 - math.exp(-lamJ*beta))
            print(f"      eps={eps} lam_J={lamJ}: E[V_1^T] exact = {V1:.4f}   CIR limit = {cir:.4f}")
