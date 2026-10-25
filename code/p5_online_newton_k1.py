"""Problem 5, K=1 experiment: Fisher-preconditioned online Newton (recursive MLE) for lam = mu + alpha*S, dS=-beta S dt + dN.
Per event i (gap length w, state s just after event i-1):
   a=(1, s e^{-beta w}), b=(w, s(1-e^{-beta w})/beta), lam_i(theta)=a'theta, g_i = b - a/lam_i(theta_i)      (grad of l_i)
   F += int_0^w phi phi'/lam_theta du   (Fisher increment, closed form),  theta <- clip( theta - F^{-1} g_i )
Pseudo-regret R_n = sum_i int_gap [lam_th - lam* - lam* log(lam_th/lam*)] du  (compensated Bregman excess; E R_n = E sum of regret vs truth).
Normalised error  e_n = (theta_n-theta*)' J_n (theta_n-theta*),  J_n = sum of Fisher increments evaluated at theta*  (~ chi^2_2 ideal)."""
import math, random, sys
xq = [-0.9602898564975363, -0.7966664774136267, -0.5255324099163290, -0.1834346424956498,
       0.1834346424956498, 0.5255324099163290, 0.7966664774136267, 0.9602898564975363]
wq = [0.1012285362903763, 0.2223810344533745, 0.3137066458778873, 0.3626837833783620]*2
wq = [0.1012285362903763, 0.2223810344533745, 0.3137066458778873, 0.3626837833783620,
      0.3626837833783620, 0.3137066458778873, 0.2223810344533745, 0.1012285362903763]

def fisher_inc(mu, al, beta, s, w):
    c = al*s; e = math.exp(-beta*w)
    Iinv = (w + math.log((mu + c*e)/(mu + c))/beta)/mu          # int 1/lam
    Lam = mu*w + c*(1-e)/beta                                    # int lam
    return Iinv, (w - mu*Iinv)/al, (Lam - 2*mu*w + mu*mu*Iinv)/(al*al)

def run(mu_s, beta, eps, n, checkpoints, seed, th0, rho0=5.0, eps_c=None, mu_min=0.1, al_min=0.02, warm=0):
    rng = random.Random(seed); al_s = (1-eps)*beta; eps_c = eps/2 if eps_c is None else eps_c
    mu, al = th0; lam_bar = mu_s/eps
    # prior Fisher: rho0 time units at the prior scale (diag; S ~ lam_bar/beta)
    F = [[rho0/lam_bar, 0.0], [0.0, rho0*(lam_bar/beta)**2/lam_bar]]
    J = [0.0, 0.0, 0.0]; s = 0.0; R = 0.0; Tt = 0.0; out = {}; al_max = beta*(1-eps_c)
    for _ in range(warm):                                        # run the true process silently to reach stationarity
        lam0 = mu_s + al_s*s; w = 0.0
        while True:
            w += rng.expovariate(lam0)
            if rng.random()*lam0 <= mu_s + al_s*s*math.exp(-beta*w): break
        s = s*math.exp(-beta*w) + 1.0
    for i in range(1, n+1):
        lam0 = mu_s + al_s*s; w = 0.0                            # true gap by thinning
        while True:
            w += rng.expovariate(lam0)
            if rng.random()*lam0 <= mu_s + al_s*s*math.exp(-beta*w): break
        e = math.exp(-beta*w); a1 = s*e; b1 = s*(1-e)/beta
        # pseudo-regret increment (theta_i fixed over the gap)
        Lt = mu*w + al*b1; Ls = mu_s*w + al_s*b1; lg = 0.0
        for x, wt in zip(xq, wq):
            u = 0.5*w*(x+1); su = s*math.exp(-beta*u); ls = mu_s + al_s*su; lt = mu + al*su
            lg += 0.5*w*wt*ls*math.log(lt/ls)
        R += Lt - Ls - lg
        # oracle information at the truth, and Fisher increment at current theta
        i0, i1, i2 = fisher_inc(mu_s, al_s, beta, s, w); J[0] += i0; J[1] += i1; J[2] += i2
        f0, f1, f2 = fisher_inc(mu, al, beta, s, w); F[0][0] += f0; F[0][1] += f1; F[1][0] += f1; F[1][1] += f2
        # gradient and Newton step
        g0 = w - 1.0/(mu + al*a1); g1 = b1 - a1/(mu + al*a1)
        det = F[0][0]*F[1][1] - F[0][1]**2
        d0 = (F[1][1]*g0 - F[0][1]*g1)/det; d1 = (-F[0][1]*g0 + F[0][0]*g1)/det
        mu = max(mu - d0, mu_min); al = min(max(al - d1, al_min), al_max)
        s = a1 + 1.0; Tt += w
        if i in checkpoints:
            x0, x1 = mu - mu_s, al - al_s
            out[i] = (R, J[0]*x0*x0 + 2*J[1]*x0*x1 + J[2]*x1*x1, Tt)
    return out

if __name__ == "__main__":
    reps = int(sys.argv[1]) if len(sys.argv) > 1 else 16; n = int(sys.argv[2]) if len(sys.argv) > 2 else 20000
    cps = [50, 100, 200, 500, 1000, 2000, 5000, 10000, 20000]; cps = [c for c in cps if c <= n]
    mu_s, beta = 1.0, 1.0
    for eps in (0.3, 0.1, 0.03):
        res = [run(mu_s, beta, eps, n, set(cps), 1000+r, (2.0*mu_s, 0.5*beta)) for r in range(reps)]
        print(f"eps={eps}  (d=2, so (d/2)log n = log n)   reps={reps}")
        print("      n   E R_n   R_n - log n   E e_n (ideal ~2)   frac(e_n>6)")
        for c in cps:
            Rm = sum(r[c][0] for r in res)/reps; em = sum(r[c][1] for r in res)/reps; fr = sum(r[c][1] > 6 for r in res)/reps
            print(f"{c:7d} {Rm:8.2f} {Rm-math.log(c):12.2f} {em:14.2f} {fr:14.2f}", flush=True)
