"""Problem 14 (v3): test the linear-filter predictions for tracking bias at small eta (noise-induced bias negligible).
 smooth (sinusoid, amplitude A, period P):  bias^2 = (A^2/2) * w^2/(eta^2 + w^2),  w = 2 pi/(P r)  (per-event frequency), r = event rate = mu0/eps
 jumps (switching at Poisson times, mean 1/rho, level change D):  avg bias^2 = rho * D^2/(2 eta r)   [exponential relaxation e^{-eta r t}]"""
import numpy as np, math
from p14_bias_scaling import run, drift_fn
import p14_bias_scaling as B

T = 6000.0; R = 24; grid = 5.0; mu0 = 1.0; r = 5.0
tgrid = np.arange(1, int(T/grid)+1)*grid; burn = tgrid > 0.2*T
rngj = np.random.default_rng(7); times = np.cumsum(rngj.exponential(600.0, 40))
models = {'smooth': dict(name='smooth', amp=0.5, P=1200.0), 'jump': dict(name='jump', amp=0.5, times=times)}
etas = [0.0005, 0.001, 0.002, 0.004]
for k, mod in models.items():
    f = drift_fn(mod, mu0); target = np.array([f(t) for t in tgrid])
    print(f"\n{k}:")
    for e in etas:
        reps = np.array([run(mod, e, T, seed=100*s+int(1e4*e)) for s in range(R)])
        n = min(reps.shape[1], len(target)); reps = reps[:, :n]; tg = target[:n]; bu = burn[:n]
        mean = reps.mean(0); varc = reps.var(0, ddof=1).mean()/R
        bias2 = np.mean(((mean-tg)**2)[bu]) - varc
        if k == 'smooth':
            w = 2*math.pi/(mod['P']*r); pred = (mod['amp']**2/2)*w*w/(e*e + w*w)
        else:
            jt = times[times < T]; rho = len(jt)/T; D = mod['amp']*mu0
            pred = rho*D*D/(2*e*r)                       # long-run average, ignores the finite horizon
        print(f"  eta={e}: measured bias^2 = {bias2:.5f}   linear-filter prediction = {pred:.5f}", flush=True)
