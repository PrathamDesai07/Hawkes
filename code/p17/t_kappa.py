import numpy as np
from lib17 import *
np.set_printoptions(precision=4, suppress=True, linewidth=150)
r,w=powerlaw_truth(0.8,0.5,1.0); cps=[300,1000,3000,10000,30000,100000]
for label,betas,thd,th0 in [('K=1 (mu_dag/mu*=2.7)',np.array([1.0]),np.array([2.7032,0.466]),[0.5,0.9]),('K=2',np.array([0.3,1.5]),np.array([1.8859,0.139,0.2462]),[0.5,0.2,0.6])]:
    print(label,' start',th0)
    for pc,name in [(0,'EF'),(1,'OI')]:
        for kappa in (0.5,1.0,2.0,4.0):
            A=run_online(40,300,100000,5000,1.0,r,w,betas,0.0,0.3,pc,1.0,th0,cps,kappa=kappa)
            d=len(betas)+1; err=np.sqrt(((A[:,:,:d]-thd)**2).sum(2))      # reps x cps
            print(f'  {name} kappa={kappa}: median |err| by n={cps}: {np.median(err,0)}   max over reps at last cp {err[:,-1].max():.3f}',flush=True)
