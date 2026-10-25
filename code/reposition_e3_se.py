"""E3 with spreads: K=2, beta=(1,4), w=(.5,.5), m=.625.  A(eps)=E[1/lam]E[lam] per replicate, mean +- se, and bootstrap se of the log-log slope."""
import math, random, sys
import numpy as np
from reposition_tests import E3_run
reps = int(sys.argv[1]); beta = np.array([1.0, 4.0]); w = np.array([0.5, 0.5]); m = 0.625; rng = random.Random(5)
for mu in (0.4, 1.6):
    k = 2*m*mu; per = {}
    for eps in (0.3, 0.1, 0.03):
        gam = eps/m; T = 150/gam; per[eps] = [E3_run(100 + r, mu, beta, w, eps, T)*mu/eps for r in range(reps)]
    eps_l = sorted(per, reverse=True)
    print(f"mu={mu} (k={k:.2f}): " + ", ".join(f"eps={e:g}: A={np.mean(per[e]):.2f}+-{np.std(per[e], ddof=1)/math.sqrt(reps):.2f}" for e in eps_l))
    sl = []
    for _ in range(400):
        mm = {e: np.mean([per[e][rng.randrange(reps)] for _ in range(reps)]) for e in eps_l}
        sl.append(math.log(mm[0.03]/mm[0.3])/math.log(0.3/0.03))
    print(f"      d log A / d log(1/eps) over eps 0.3->0.03: {np.mean(sl):.2f} +- {np.std(sl):.2f} (bootstrap);  prediction: " + (f"{1-k:.2f} (asymptotic)" if k < 1 else f"-> 0, A -> {k/(k-1):.2f}"), flush=True)
