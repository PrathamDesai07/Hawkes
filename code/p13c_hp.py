"""High-probability version: predictable leverage lp_i (a function of S_{i-1}, F_{i-1}, theta*_i only) dominates the actual leverage lev_i pathwise, so
tau = min{i : lp_i > lbar} is a predictable stopping time ({tau >= i} in F_{t_{i-1}}).  Measure: (i) max(lev - lp) (<=0), (ii) fraction of runs with tau <= n for several lbar, inside / outside the window."""
import numpy as np, math, sys
from p13c_engine import run
from p14_mech import Iev_matrix
mu_s = 1.0; A = 0.3; P = 1200.0; fn = lambda t: mu_s*(1 + A*np.sin(2*np.pi*t/P)); R = int(sys.argv[1])
print("setting: eta/eps^2 | max(lev-lp) | lp quantiles (50%,99%,max) | P(tau<=n) for lbar = 0.1, 0.2, 0.5, 1.0")
for eps, eta, n in ((0.4, 0.004, 15000), (0.2, 0.001, 30000), (0.2, 0.004, 30000), (0.4, 0.016, 15000)):
    Iev = Iev_matrix(mu_s, 1.0, eps); out = [run([1.0], eps, fn, mu_s*(1+A), eta, n, 1100+s, Iev, mu_s, project=False, cminus=0.6) for s in range(R)]
    lev = np.concatenate([o['lev'] for o in out]); lp = np.concatenate([o['lp'] for o in out])
    mx = [o['lp'].max() for o in out]; fr = [np.mean([m > L for m in mx]) for L in (0.1, 0.2, 0.5, 1.0)]
    print(f" eps={eps} eta={eta}: {eta/eps**2:.3f} | {np.max(lev-lp):+.2e} | {np.median(lp):.4f}, {np.quantile(lp,0.99):.4f}, {lp.max():.4f} | " + ", ".join(f"{f:.2f}" for f in fr), flush=True)
