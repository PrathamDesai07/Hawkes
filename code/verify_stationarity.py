"""E_pi[A f_z] = 0 with an honest error bar: per-replicate means from independent replicates, SE across replicates."""
import numpy as np, math, sys
from verify_lemmaH_EF import setup, Af
from verify_EF_mc import sample_stationary
reps = int(sys.argv[1])
beta = np.array([1.0, 2.0]); alpha = 0.95*beta/2; P = setup(beta, alpha, 1.5); g = P['gam']
for z in (1.0, 10.0):
    means = []
    for r in range(reps):
        S = sample_stationary(P, 300/g, 0.5, 500 + r, 30/g); A, W = Af(P, S, z); means.append(A.mean())
    m = float(np.mean(means)); se = float(np.std(means, ddof=1)/math.sqrt(reps))
    print(f"z={z}: E[A f_z] = {m:+.2e} +- {se:.1e} (across {reps} independent replicates), z-score {m/se:+.2f}")
