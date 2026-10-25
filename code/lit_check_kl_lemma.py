"""Independent test of the KL lemma used in problem11_dichotomy.md (other session):  KL(P_theta^T || P_theta'^T) <= rho * u' I_T(theta) u,
 u = theta' - theta, rho = max(mu/mu', alpha/alpha'), I_T = T*I(theta) for a STATIONARY start (same initial state S_0 under both laws).
 KL estimated by Monte Carlo as E_theta int_0^T [lam log(lam/lam') - lam + lam'] dt (compensator form, 8-pt Gauss per gap).
 Direction: the 'ridge' u = (delta, -w*delta), w = I_ma / I_aa  (least informative direction), and the pure-mu direction u = (delta, 0)."""
import math, random, sys
from lit_check_herrera import gap, xs, ws, fisher_exact

def kl_stat(reps, mu, beta, eps, du, dal, T, seed=0):
    al = (1 - eps)*beta; mu2 = mu + du; al2 = al + dal; vals = []
    for r in range(reps):
        rng = random.Random(seed + r); s = 0.0
        for _ in range(int(10/eps**2)): s = s*math.exp(-beta*gap(rng, mu, al, beta, s)) + 1.0     # stationary start under theta
        t = 0.0; kl = 0.0
        while t < T:
            w = gap(rng, mu, al, beta, s); wu = min(w, T - t)
            for x, wt in zip(xs, ws):
                u = 0.5*wu*(x + 1); su = s*math.exp(-beta*u); l0 = mu + al*su; l1 = mu2 + al2*su
                kl += 0.5*wu*wt*(l0*math.log(l0/l1) - l0 + l1)
            t += w; s = s*math.exp(-beta*w) + 1.0
        vals.append(kl)
    m = sum(vals)/reps; se = (sum((v - m)**2 for v in vals)/(reps - 1))**.5/reps**.5
    return m, se, mu2, al2, al

if __name__ == "__main__":
    reps = int(sys.argv[1]); beta = 1.0
    print(f"{'k0':>4} {'eps':>5} {'T':>6} {'direction':>9} {'delta':>6} | {'KL (MC)':>10} {'+-se':>7} | {'rho u'+chr(39)+'I_T u':>11} {'(KL <= ?)':>9} | {'1/2 u'+chr(39)+'I_T u':>10} {'KL/(1/2 uIu)':>12}")
    for mu in (1.0, 0.3):
        for eps in (0.1,):
            m, jma, jaa = fisher_exact(mu, beta, eps); w = jma/jaa; T = 2.0/eps; al = (1 - eps)*beta
            Imm, Ima, Iaa = m, jma, jaa
            for name, dirw in (("ridge", w), ("mu-only", 0.0)):
                for delta in (0.1, 0.3):
                    du, dal = delta*mu, -dirw*delta*mu
                    kl, se, mu2, al2, al0 = kl_stat(reps, mu, beta, eps, du, dal, T)
                    uIu = T*(Imm*du*du + 2*Ima*du*dal + Iaa*dal*dal); rho = max(mu/mu2, al0/al2)
                    print(f"{2*mu/beta:4.1f} {eps:5.2f} {T:6.1f} {name:>9} {delta:6.2f} | {kl:10.4f} {se:7.4f} | {rho*uIu:11.4f} {str(kl <= rho*uIu):>9} | {0.5*uIu:10.4f} {kl/(0.5*uIu):12.2f}", flush=True)
