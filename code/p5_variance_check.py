"""Check Var( int_0^T h(S_t) dt ) <= 2 T Lip(h) sqrt(m) sigma_S / gamma   (K=1, stationary start), h = trapezoid in [0,1].
 sigma_S = sqrt(E S^2 + 3 (E S)^2), gamma = beta*eps,  m = E_pi h  (estimated from the same runs)."""
import math, random, sys
from p5_fisher_conc import stationary_law_moments
xs = [-0.9602898564975363, -0.7966664774136267, -0.5255324099163290, -0.1834346424956498, 0.1834346424956498, 0.5255324099163290, 0.7966664774136267, 0.9602898564975363]
ws = [0.1012285362903763, 0.2223810344533745, 0.3137066458778873, 0.3626837833783620, 0.3626837833783620, 0.3137066458778873, 0.2223810344533745, 0.1012285362903763]

def trap(s, lo, hi, ramp):          # 0 outside [lo,hi]; ramps of width `ramp`; 1 on plateau
    if s <= lo or s >= hi: return 0.0
    return min(1.0, (s - lo)/ramp, (hi - s)/ramp)

def occupation(rng, mu, beta, eps, T, warm, h):
    al = (1-eps)*beta; s = 0.0
    for _ in range(warm):
        lam0 = mu + al*s; w = 0.0
        while True:
            w += rng.expovariate(lam0)
            if rng.random()*lam0 <= mu + al*s*math.exp(-beta*w): break
        s = s*math.exp(-beta*w) + 1.0
    t = 0.0; I = 0.0
    while t < T:
        lam0 = mu + al*s; w = 0.0
        while True:
            w += rng.expovariate(lam0)
            if rng.random()*lam0 <= mu + al*s*math.exp(-beta*w): break
        wu = min(w, T - t); sub = 4; acc = 0.0                   # split gap in 4 subintervals, 8-pt Gauss each
        for k in range(sub):
            a = wu*k/sub; b = wu*(k+1)/sub
            for x, wt in zip(xs, ws): acc += 0.5*(b-a)*wt*h(s*math.exp(-beta*(0.5*(a+b) + 0.5*(b-a)*x)))
        I += acc; t += w; s = s*math.exp(-beta*w) + 1.0
    return I

if __name__ == "__main__":
    reps = int(sys.argv[1]); mu, beta = 1.0, 1.0
    for eps in (0.3, 0.1):
        ES, ES2 = stationary_law_moments(mu, beta, eps); sig = math.sqrt(ES2 + 3*ES*ES); gam = beta*eps; d = 0.15*ES
        for name, c in (("I1", 0.5*ES), ("I2", 1.8*ES)):
            lo, hi, ramp = c - d, c + d, d/2; h = lambda s, lo=lo, hi=hi, ramp=ramp: trap(s, lo, hi, ramp); Lip = 1.0/ramp
            for Tg in (5, 20):
                T = Tg/gam; rng = random.Random(77)
                v = [occupation(rng, mu, beta, eps, T, int(20/eps**2), h) for _ in range(reps)]
                mean = sum(v)/reps; var = sum((x - mean)**2 for x in v)/(reps-1); m = mean/T
                bound = 2*T*Lip*math.sqrt(m)*sig/gam
                print(f"eps={eps} {name} T*gamma={Tg:3d}: m={m:.4f}  Var={var:9.2f}  bound={bound:11.1f}  Var/bound={var/bound:.4f}  Chebyshev P(O<=Tm/2) bound={min(1, 4*bound/(T*m)**2):.3f}  empirical={sum(x <= T*m/2 for x in v)/reps:.3f}", flush=True)
