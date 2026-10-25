"""Problem 14: how does constant-step tracking error scale with eta under different drift models? K=1, alpha=(1-eps)beta known,
estimate mu_t by a Fisher-preconditioned constant-gain recursion, one update per event:
    mu_hat <- mu_hat + eta * (r/m) * (1/lam_i(mu_hat) - w_i),  w_i = gap since previous event,
 r = event rate, m = E[1/lam]  (so eta is a per-event contraction rate; I_ev = m/r).  Pure python, cumulative-time thinning.
Drift models: none; smooth (sinusoid, rate Delta); jumps (Poisson switching between two levels, rate rho); random walk (variance q per unit time)."""
import numpy as np, math, random
from p10_fisher_k1_exact import exact_m

def run(model, eta, T, mu0=1.0, beta=1.0, eps=0.2, seed=0, dt=0.05):
    rng = random.Random(seed); alpha = (1-eps)*beta
    m = exact_m(mu0, beta, eps); r = mu0/eps; gain = eta*r/m
    S = 0.0; t = 0.0; mu_hat = mu0; mu_t = mu0; level = 0; err2 = 0.0; cnt = 0; last = 0.0; burn = 0.1*T
    nxt_jump = rng.expovariate(model.get('rho', 1e-9)) if model['name'] == 'jump' else 1e18
    while t < T:
        tend = t + dt                                     # parameters held constant within a step
        if model['name'] == 'smooth': mu_t = mu0*(1 + model['amp']*math.sin(2*math.pi*tend/model['P']))
        elif model['name'] == 'rw': mu_t = max(0.3, mu_t + math.sqrt(model['q']*dt)*rng.gauss(0, 1))
        elif model['name'] == 'jump':
            while nxt_jump <= tend: level = 1-level; nxt_jump += rng.expovariate(model['rho'])
            mu_t = mu0*(1 + model['amp']*level)
        pos = t
        lam0 = mu_t + alpha*S*math.exp(-beta*0.0)
        while True:
            lam = mu_t + alpha*S
            w = rng.expovariate(lam)
            if pos + w > tend:
                S *= math.exp(-beta*(tend-pos)); pos = tend; break
            S *= math.exp(-beta*w); pos += w
            if rng.random()*lam < mu_t + alpha*S:        # accept (intensity at pos <= lam since decaying)
                lam_i = mu_t + alpha*S                   # approx; event-time intensity under true params
                gap = pos - last; last = pos
                lh = mu_hat + alpha*S
                mu_hat += gain*(1/lh - gap)
                mu_hat = min(max(mu_hat, 0.05), 20.0)
                if pos > burn: err2 += (mu_hat-mu_t)**2; cnt += 1
                S += 1.0
        t = tend
    return err2/max(cnt, 1)

if __name__ == "__main__":
    import sys
    T = float(sys.argv[1]) if len(sys.argv) > 1 else 6000
    etas = [0.003, 0.01, 0.03, 0.1, 0.3]
    models = {'none': dict(name='none'),
              'smooth': dict(name='smooth', amp=0.3, P=1500.0),
              'jump': dict(name='jump', amp=0.5, rho=0.004),
              'rw': dict(name='rw', q=2e-4)}
    res = {k: [] for k in models}
    for k, mod in models.items():
        for e in etas: res[k].append(np.mean([run(mod, e, T, seed=s) for s in range(3)]))
        print(k, ' '.join(f"{x:.4f}" for x in res[k]), flush=True)
    a = np.polyfit(etas, res['none'], 1)[0]          # variance part  ~ a*eta
    print(f"\nvariance part: MSE_none ~ {a:.3f}*eta; residual (drift part) = MSE - MSE_none, log-log slope vs eta over eta in [0.003,0.03]:")
    for k in ('smooth', 'jump', 'rw'):
        resid = np.array(res[k]) - np.array(res['none'])
        sl = np.polyfit(np.log(etas[:3]), np.log(np.maximum(resid[:3], 1e-9)), 1)[0]
        print(f"  {k}: residual {np.round(resid,4).tolist()}  slope {sl:.2f}")
