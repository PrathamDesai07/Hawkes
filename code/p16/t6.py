import numpy as np
from p16lib import *
np.set_printoptions(precision=4, suppress=True)
r, w = powerlaw_truth(0.8, 0.5, 1.0); betas=np.array([0.3,1.5])
ref=[]
for s in range(3):
    P=simulate(100+s, 2_000_000, 20000, 1.0, r, w, betas); ref.append(solve(P)[0])
ref=np.array(ref); thr=ref.mean(0); print('ref', thr, 'spread of 3 ref runs', ref.std(0))
for nev in [10_000, 40_000, 160_000, 640_000]:
    errs=[]
    for s in range(16):
        P=simulate(1000+s, nev, 5000, 1.0, r, w, betas); errs.append(solve(P)[0]-thr)
    errs=np.array(errs); rms=np.sqrt((errs**2).mean(0)); print(nev, 'rms err', rms, 'mean err (bias)', errs.mean(0))
