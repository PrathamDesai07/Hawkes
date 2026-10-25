"""Problem 5's estimator (verbatim, eta=0) and its forgetting version under smooth drift of mu, alpha UNKNOWN (jointly estimated).
 Predictions (preconditioned Newton, contraction eta in every direction):
   bias^2 of mu_hat = (A^2/2) w^2/(eta^2+w^2),  w = 2 pi/(P r)      [same as alpha known]
   noise: MSE(mu_hat) ~ a eta,  a = r/(2 eff),  eff = (m - eps/mu)/(1-2eps+eps mu m)   [alpha unknown]  vs r/(2m) [alpha known]"""
import numpy as np, math, sys
from p14_p5_drift import run_drift
from p10_fisher_k1_exact import exact_m

eps = float(sys.argv[1]); R = int(sys.argv[2]); T = 6000.0; P = 1200.0; A = 0.5; mu_s = 1.0; beta = 1.0; grid = 5.0
r = mu_s/eps; m = exact_m(mu_s, beta, eps); eff = (m - eps/mu_s)/(1 - 2*eps + eps*mu_s*m)
tgrid = np.arange(1, int(T/grid)+1)*grid; target = mu_s*(1 + A*np.sin(2*np.pi*tgrid/P))
print(f"eps={eps}: r={r:.1f}  m={m:.4f}  eff={eff:.4f}   noise predictions a: alpha unknown r/(2 eff)={r/(2*eff):.2f}, alpha known r/(2m)={r/(2*m):.2f}")
for eta in (0.0, 0.001, 0.003):
    reps = []; reps_al = []
    for s in range(R):
        ev, _ = run_drift(mu_s, beta, eps, 500+s, (mu_s, (1-eps)*beta), T=T, eta=eta, A=A, P=P)
        idx = np.searchsorted(ev[:, 0], tgrid, side='right') - 1; idx = np.clip(idx, 0, len(ev)-1)
        reps.append(ev[idx, 1]); reps_al.append(ev[idx, 2])
    reps = np.array(reps); mean = reps.mean(0); varc = reps.var(0, ddof=1).mean()/R; bu = tgrid > 0.2*T
    bias2 = np.mean(((mean - target)**2)[bu]) - varc; mse = np.mean(((reps - target)**2)[:, bu])
    if eta > 0:
        w = 2*math.pi/(P*r); pred = (A*A/2)*w*w/(eta*eta + w*w)
        print(f" eta={eta}: bias^2 = {bias2:.5f} (filter prediction {pred:.5f});  per-replicate MSE = {mse:.4f}  ->  MSE-bias^2 = {mse-bias2:.4f}, /eta = {(mse-bias2)/eta:.1f}", flush=True)
    else:
        # Problem 5 verbatim: bias in successive time windows
        wins = [(0.2*T, 0.4*T), (0.4*T, 0.6*T), (0.6*T, 0.8*T), (0.8*T, T)]
        out = []
        for lo, hi in wins:
            sel = (tgrid > lo) & (tgrid <= hi); out.append(np.mean(((mean-target)**2)[sel]) - varc)
        print(f" eta=0 (Problem 5 verbatim): bias^2 in windows {[f'{lo/T:.1f}-{hi/T:.1f}T' for lo,hi in wins]} = {np.round(out,4).tolist()};  mean(mu_hat) range over t>0.2T: [{mean[bu].min():.3f}, {mean[bu].max():.3f}] vs truth range [{target[bu].min():.3f}, {target[bu].max():.3f}]", flush=True)
