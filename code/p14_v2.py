"""Problem 14 (v2 estimators). K=1, alpha known, time-varying mu_t.  Two constant-gain estimators of mu_t, one update per event:
 'naive' : mu_hat += eta*(r/m)*(1/lam_hat - w)                       (preconditioner = stationary information, fixed)
 'newton': F <- (1-eta) F + dF_i ; mu_hat -= g_i/F   with dF_i = int_gap du/lam_hat(u), g_i = w - 1/lam_hat_i   (exponentially forgotten Fisher; effective window 1/eta events)
Pure python, cumulative-time thinning, parameters held constant on steps of length dt."""
import numpy as np, math, random
from p10_fisher_k1_exact import exact_m

def drift_fn(model, mu0):
    if model['name'] == 'none': return lambda t: mu0
    if model['name'] == 'smooth': return lambda t: mu0*(1 + model['amp']*math.sin(2*math.pi*t/model['P']))
    times = model['times']; return lambda t: mu0*(1 + model['amp']*(np.searchsorted(times, t) % 2))

def run(model, est, eta, T, seed, eps, mu0=1.0, beta=1.0, dt=0.1, grid=5.0):
    rng = random.Random(seed); alpha = (1-eps)*beta; r = mu0/eps; m = exact_m(mu0, beta, eps)
    f = drift_fn(model, mu0); S = 0.0; t = 0.0; mu_hat = mu0; last = 0.0; out = []; nextg = grid
    F = 1.0/eta*(m/r)                       # initial forgotten information ~ window * per-event information
    while t < T:
        tend = t + dt; mu_t = f(tend); pos = t
        while True:
            lam = mu_t + alpha*S; w = rng.expovariate(lam)
            if pos + w > tend: S *= math.exp(-beta*(tend-pos)); pos = tend; break
            S *= math.exp(-beta*w); pos += w
            if rng.random()*lam < mu_t + alpha*S:
                gap = pos - last; last = pos; lh = mu_hat + alpha*S
                if est == 'naive':
                    mu_hat += eta*(r/m)*(1/lh - gap)
                else:
                    s_prev = S/1.0                                    # state just before the event
                    # int over gap of 1/lam_hat: state decayed from s0=(S_after_prev) ; use s_pre=S at event with exp back-projection
                    s0 = S*math.exp(beta*gap)                         # state at gap start
                    dF = (gap + math.log((mu_hat + alpha*S)/(mu_hat + alpha*s0))/beta)/mu_hat
                    F = (1-eta)*F + dF
                    mu_hat -= (gap - 1/lh)/F
                mu_hat = min(max(mu_hat, 0.05), 20.0); S += 1.0
        t = tend
        if t >= nextg: out.append(mu_hat); nextg += grid
    return np.array(out)

def stats(model, est, eta, T, R, eps, seed0=0, grid=5.0, mu0=1.0, burn=0.2):
    tgrid = np.arange(1, int(T/grid)+1)*grid; f = drift_fn(model, mu0); target = np.array([f(t) for t in tgrid])
    reps = np.array([run(model, est, eta, T, seed0+s, eps, mu0=mu0) for s in range(R)])
    n = min(reps.shape[1], len(target)); reps = reps[:, :n]; tg = target[:n]; bu = tgrid[:n] > burn*T
    mse = np.mean(((reps-tg)**2)[:, bu]); mean = reps.mean(0)
    bias2 = np.mean(((mean-tg)**2)[bu]) - reps.var(0, ddof=1).mean()/R
    return mse, bias2
