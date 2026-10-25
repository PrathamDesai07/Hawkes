import numpy as np
from p16lib import *
np.set_printoptions(precision=4, suppress=True)
for nstar in [0.8, 0.95]:
  r, w = powerlaw_truth(nstar, 0.5, 1.0); betas=np.array([0.1,0.5,2.5])
  P=simulate(5, 3_000_000, 50000, 1.0, r, w, betas); th=solve(P)[0]; tq=quad_surrogate(P); Lam=P.lam_mean()
  print('n*',nstar,'Lam %.3f'%Lam,'KL',th,'n=%.4f'%nvec(th,betas),'| quad',tq,'n=%.4f'%nvec(tq,betas),'| rel diff %.3f'%(np.linalg.norm(th-tq)/np.linalg.norm(th)))
