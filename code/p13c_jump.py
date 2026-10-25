"""Sharp jump rate WITH exact projection.  Baseline matched to the noise level: average of two drift-free runs at the two plateau levels (same seeds, same reference metric)."""
import numpy as np, math, sys
from p13c_engine import run
from p14_mech import Iev_matrix
eps = 0.4; mu_s = 1.0; D = 0.3; Lt = 800; R = int(sys.argv[1]); Iev = Iev_matrix(mu_s, 1.0, eps); r = mu_s/eps; Le = Lt*r; n = int(8*Le)
fj = lambda t: mu_s*(1 + D*(np.floor(np.asarray(t)/Lt).astype(int) % 2))
print(f"exact projection, K=1, eps={eps}, jump D={D}, regime {Lt} time units")
res = []
for eta in (0.002, 0.004, 0.008):
    mj = np.concatenate([run([1.0], eps, fj, mu_s*(1+D), eta, n, 400+s, Iev, mu_s, burn=0.2)['vI'] for s in range(R)]).mean()
    base = []
    for lvl in (mu_s, mu_s*(1+D)):
        # events per regime differ with the level, so use the per-event mean error (event-weighted); the jump run is also event-weighted
        base.append(np.concatenate([run([1.0], eps, (lambda t, L=lvl: L), L if False else lvl, eta, int(4*Le), 400+s, Iev, mu_s, burn=0.2)['vI'] for s in range(R)]).mean())
    B = 0.5*(base[0] + base[1]); cost = (mj - B)*Le; res.append((eta, cost))
    print(f" eta={eta}: jumps {mj:.5f}; baselines (level 1.0, 1.3) {base[0]:.5f}, {base[1]:.5f} -> matched {B:.5f}; per-jump excess cost {cost:.3f}", flush=True)
res = np.array(res); print(f" log-log slope of the per-jump excess vs eta: {np.polyfit(np.log(res[:,0]), np.log(np.maximum(res[:,1],1e-9)), 1)[0]:.2f}  (sharp rate -1)")
