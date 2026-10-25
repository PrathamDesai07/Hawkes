"""Independent numerical tests of claims in Herrera-Marin (arXiv 2609.38268) against my K=1 model  lam = mu + alpha*S.
 (T1) Theorem 2.2:  KL(P0^T || P1^T) <= delta^2 (1+c) I_T / (rho0 rho1),  I_T = T mu (1-rho0), empty start (c=0), same mu and kernel shape.
      KL estimated by Monte Carlo as E_0 int_0^T phi(lam0,lam1) dt,  phi(x,y) = x log(x/y) - x + y  (compensator form), integrated exactly per gap with 8-pt Gauss.
      Also compared with the local quadratic approximation 0.5 * J_eta,eta * (d eta)^2 (Theorem 4.1 entry).
 (T2) Theorem 4.1 entries vs my I_mu,mu, I_mu,alpha, I_alpha,alpha, via the SCORE covariance (an estimator independent of E[aa'/lam]):
      score_mu = sum_i 1/lam(t_i) - T,  score_alpha = sum_i S(t_i-)/lam(t_i) - int S dt ; Cov/T -> I   (stationary start)."""
import math, random, sys
sys.path.insert(0, '.')
from p10_fisher_k1 import exact_m
xs = [-0.9602898564975363, -0.7966664774136267, -0.5255324099163290, -0.1834346424956498, 0.1834346424956498, 0.5255324099163290, 0.7966664774136267, 0.9602898564975363]
ws = [0.1012285362903763, 0.2223810344533745, 0.3137066458778873, 0.3626837833783620, 0.3626837833783620, 0.3137066458778873, 0.2223810344533745, 0.1012285362903763]

def gap(rng, mu, al, beta, s):
    lam0 = mu + al*s; w = 0.0
    while True:
        w += rng.expovariate(lam0)
        if rng.random()*lam0 <= mu + al*s*math.exp(-beta*w): return w

def kl_mc(reps, mu, beta, eps0, delta, T, seed=0):
    rho0 = 1 - eps0; eps1 = (1+delta)*eps0; rho1 = 1 - eps1; a0 = rho0*beta; a1 = rho1*beta; tot = 0.0; vals = []
    for r in range(reps):
        rng = random.Random(seed + r); s = 0.0; t = 0.0; kl = 0.0       # empty start at 0: S = 0
        while t < T:
            w = gap(rng, mu, a0, beta, s); wu = min(w, T - t)
            for x, wt in zip(xs, ws):
                u = 0.5*wu*(x+1); su = s*math.exp(-beta*u); l0 = mu + a0*su; l1 = mu + a1*su
                kl += 0.5*wu*wt*(l0*math.log(l0/l1) - l0 + l1)
            t += w; s = s*math.exp(-beta*w) + 1.0
        vals.append(kl)
    m = sum(vals)/reps; se = (sum((v - m)**2 for v in vals)/(reps-1))**.5/reps**.5
    return m, se, rho0, rho1

def fisher_exact(mu, beta, eps):
    al = (1-eps)*beta; m = exact_m(mu, al, beta); El = mu/eps
    return m, (1-mu*m)/al, (El - 2*mu + mu*mu*m)/al**2

def score_cov(reps, mu, beta, eps, T, seed=100):
    al = (1-eps)*beta; sm = []; sa = []
    for r in range(reps):
        rng = random.Random(seed + r); s = 0.0
        for _ in range(int(10/eps**2)): s = s*math.exp(-beta*gap(rng, mu, al, beta, s)) + 1.0
        t = 0.0; a = b = 0.0
        while t < T:
            w = gap(rng, mu, al, beta, s); e = math.exp(-beta*w); s_end = s*e; lam = mu + al*s_end
            a += 1.0/lam; b += s_end/lam; a -= w; b -= s*(1-e)/beta; t += w; s = s_end + 1.0
        sm.append(a); sa.append(b)
    n = len(sm); mm = sum(sm)/n; ma = sum(sa)/n
    c = lambda x, y, mx, my: sum((p - mx)*(q - my) for p, q in zip(x, y))/(n-1)
    return c(sm, sm, mm, mm)/T, c(sm, sa, mm, ma)/T, c(sa, sa, ma, ma)/T, mm/T, ma/T

if __name__ == "__main__":
    reps = int(sys.argv[1]); mu, beta = 1.0, 1.0
    print("T1: KL(P0||P1) vs Theorem 2.2 bound delta^2 I_T/(rho0 rho1)  and local approximation")
    print(f"{'eps0':>5} {'delta':>6} {'T':>7} | {'KL (MC)':>10} {'+-se':>7} | {'Thm 2.2 bound':>13} {'bound/KL':>9} | {'local approx':>12} {'approx/KL':>9}")
    for eps0 in (0.3, 0.1, 0.03):
        for delta in (0.5, 0.1):
            for Tfac in (3, 30):
                T = Tfac/(mu*eps0)           # mu*eps*T = Tfac  (I_T = Tfac)
                kl, se, rho0, rho1 = kl_mc(reps, mu, beta, eps0, delta, T)
                IT = T*mu*eps0; bound = delta**2*IT/(rho0*rho1)
                m, jma, jaa = fisher_exact(mu, beta, eps0)
                # local approx via eta = log(eps): d eta ~ log(1+delta); J_eta,eta per unit time = (a/rho)^2 (E lam - 2mu + mu^2 m)... use paper's entry
                a_ = eps0; J_eta = (a_**2/rho0**2)*((mu/eps0) - 2*mu + mu*mu*m); deta = math.log(1+delta)
                # full local approx needs the efficient information over mu (mu shared!) -> mu is NOT varied, so use the plain eta,eta entry
                approx = 0.5*J_eta*deta**2*T
                print(f"{eps0:5.2f} {delta:6.2f} {T:7.1f} | {kl:10.4f} {se:7.4f} | {bound:13.4f} {bound/kl:9.2f} | {approx:12.4f} {approx/kl:9.2f}", flush=True)
    print("\nT2: Fisher entries: exact (my ODE) vs score covariance (independent estimator), mu=beta=1")
    for eps in (0.3, 0.1):
        m, jma, jaa = fisher_exact(mu, beta, eps); T = 200/eps
        cmm, cma, caa, mm, ma = score_cov(reps, mu, beta, eps, T)
        print(f"eps={eps}: I_mm exact {m:.4f}  score-cov {cmm:.4f} | I_ma exact {jma:.4f}  score-cov {cma:.4f} | I_aa exact {jaa:.3f}  score-cov {caa:.3f} | mean scores/T: {mm:+.4f}, {ma:+.4f} (should be 0)", flush=True)
