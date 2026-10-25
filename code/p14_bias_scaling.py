"""Problem 14 (v2): tracking BIAS vs eta. Replicates share the same deterministic drift path; the ensemble-mean of mu_hat_t is compared
with mu_t, so estimator noise averages out (noise variance/R is subtracted).  K=1, alpha known, eps=0.2, mu0=1, beta=1.
Models: smooth (sinusoid, amplitude A, period P) and jumps (fixed switching times, alternating levels).
Prediction: smooth -> bias^2 ~ (Delta_event/eta)^2 ; jumps -> bias^2 ~ rho D^2/(2 eta r)  (per-unit-time average)."""
import numpy as np, math, random
from p10_fisher_k1_exact import exact_m

def drift_fn(model, mu0):
    if model['name'] == 'smooth':
        return lambda t: mu0*(1 + model['amp']*math.sin(2*math.pi*t/model['P']))
    times = model['times']
    def f(t):
        lvl = np.searchsorted(times, t) % 2
        return mu0*(1 + model['amp']*lvl)
    return f

def run(model, eta, T, seed, mu0=1.0, beta=1.0, eps=0.2, dt=0.05, grid=5.0):
    rng = random.Random(seed); alpha = (1-eps)*beta; m = M_CACHE; r = mu0/eps; gain = eta*r/m
    f = drift_fn(model, mu0); S = 0.0; t = 0.0; mu_hat = mu0; last = 0.0; out = []; nextg = grid
    while t < T:
        tend = t + dt; mu_t = f(tend); pos = t
        while True:
            lam = mu_t + alpha*S; w = rng.expovariate(lam)
            if pos + w > tend: S *= math.exp(-beta*(tend-pos)); pos = tend; break
            S *= math.exp(-beta*w); pos += w
            if rng.random()*lam < mu_t + alpha*S:
                gap = pos - last; last = pos
                mu_hat += gain*(1/(mu_hat + alpha*S) - gap); mu_hat = min(max(mu_hat, 0.05), 20.0)
                S += 1.0
        t = tend
        if t >= nextg: out.append(mu_hat); nextg += grid
    return np.array(out)

M_CACHE = exact_m(1.0, 1.0, 0.2)
if __name__ == "__main__":
    import sys
    T = 3000.0; R = int(sys.argv[1]) if len(sys.argv) > 1 else 30; grid = 5.0
    tgrid = np.arange(1, int(T/grid)+1)*grid; burn = tgrid > 0.15*T
    rngj = np.random.default_rng(5); times = np.cumsum(rngj.exponential(300.0, 40))
    models = {'smooth': dict(name='smooth', amp=0.5, P=300.0), 'jump': dict(name='jump', amp=0.5, times=times)}
    etas = [0.003, 0.01, 0.03, 0.1]
    r = 1.0/0.2; mu0 = 1.0
    for k, mod in models.items():
        f = drift_fn(mod, mu0); target = np.array([f(t) for t in tgrid])
        print(f"\n{k}:")
        bs = []
        for e in etas:
            reps = np.array([run(mod, e, T, seed=100*s+int(1000*e)) for s in range(R)])
            n = min(reps.shape[1], len(target)); reps = reps[:, :n]; tg = target[:n]; bu = burn[:n]
            mean = reps.mean(0); varc = reps.var(0, ddof=1).mean()/R
            bias2 = np.mean(((mean - tg)**2)[bu]) - varc
            bs.append(bias2)
            print(f"  eta={e}: bias^2 = {bias2:.5f}   (noise var/R removed = {varc:.5f}; per-replicate MSE {np.mean(((reps-tg)**2)[:, bu]):.4f})", flush=True)
        bs = np.array(bs); sl = np.polyfit(np.log(etas[:3]), np.log(np.maximum(bs[:3], 1e-9)), 1)[0]
        print(f"  log-log slope of bias^2 vs eta over eta in [0.003, 0.03]: {sl:.2f}   (smooth predicts -2, jumps predict -1)")
