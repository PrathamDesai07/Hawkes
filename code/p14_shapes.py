"""Other drift shapes (K=1, eps=0.4, mu0=1, beta=1; estimators 'known' and 'adaptive'); predictions from the linear filter, e_{t+1}=(1-eta)e_t + noise - d_t:
 ramp     mu_t = mu0 + s t        : steady lag  mu_hat - mu_t = -(s/r)/eta   (s = slope per time unit, r = event rate)
 triangle period P, amplitude A  : slope +-4 A mu0/P; mean-square bias ~ (4 A mu0/(P r eta))^2 (corners ignored)
 random walk, variance q per time: MSE = a eta + (q/r)/(2 eta),   a = noise coefficient (6.7-7 for adaptive at eps=0.4)"""
import numpy as np, math, sys
from p14_gen import run
eps = 0.4; mu0 = 1.0; beta = [1.0]; r = mu0/eps; R = int(sys.argv[1]); which = sys.argv[2]; ests = sys.argv[3].split(',')
etas = (0.002, 0.004)
for eta in etas:
    for est in ests:
        if which == 'ramp':
            T = 3000.0; s = 0.8/T; fn = lambda t: mu0*0.8 + s*t; mu_max = mu0*0.8 + s*T
            tg = np.arange(1, int(T/5)+1)*5.0; sel = (tg > 700) & (tg < T-50)
            lag = np.mean([ (run(est, beta, eps, fn, mu_max, eta, T, 100+k, mu0*0.8)[1] - fn(tg))[sel].mean() for k in range(R)])
            print(f"ramp  eta={eta} {est:8s}: measured steady lag {lag:+.4f}   predicted -(s/r)/eta = {-(s/r)/eta:+.4f}", flush=True)
        elif which == 'tri':
            P = 1200.0; A = 0.3; T = 6000.0
            tri = lambda t: (2/np.pi)*np.arcsin(np.sin(2*np.pi*np.asarray(t)/P))
            fn = lambda t: mu0*(1 + A*tri(t)); mu_max = mu0*(1+A)
            tg = np.arange(1, int(T/5)+1)*5.0; sel = tg > 0.2*T
            reps = np.array([run(est, beta, eps, fn, mu_max, eta, T, 200+k, mu0)[1][sel] for k in range(R)])
            mean = reps.mean(0); varc = reps.var(0, ddof=1).mean()/R; b2 = np.mean((mean - fn(tg[sel]))**2) - varc
            pred = (4*A*mu0/(P*r*eta))**2
            print(f"tri   eta={eta} {est:8s}: bias^2 {b2:.5f}   predicted (4A mu0/(P r eta))^2 = {pred:.5f} (corners smooth it: expect slightly below)", flush=True)
        else:
            T = 4000.0; q = 4e-5; dt = 0.1; a_pred = {'known': 2.8, 'adaptive': 7.0}[est]
            mses = []
            for k in range(R):
                rng = np.random.default_rng(300+k); path = np.clip(mu0 + np.cumsum(np.sqrt(q*dt)*rng.standard_normal(int(T/dt)+20)), 0.5, 1.8)
                fn = lambda t, p=path: p[np.minimum((np.asarray(t)/dt).astype(int), len(p)-1)]
                tg, x = run(est, beta, eps, lambda t, p=path: float(p[min(int(t/dt), len(p)-1)]), float(path.max()), eta, T, 400+k, mu0)
                sel = tg > 0.2*T; mses.append(np.mean((x[sel] - fn(tg[sel]))**2))
            pred = a_pred*eta + (q/r)/(2*eta)
            print(f"rw    eta={eta} {est:8s}: MSE {np.mean(mses):.4f} +- {np.std(mses,ddof=1)/math.sqrt(R):.4f}   predicted a*eta + (q/r)/(2 eta) = {pred:.4f}  (noise part {a_pred*eta:.4f}, drift part {(q/r)/(2*eta):.4f})", flush=True)
