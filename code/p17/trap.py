import numpy as np
from lib17 import *
betas=np.array([0.3,1.5]); r=np.array([0.3,1.5]); w=np.array([0.12,0.45]); mu=1.0
for mm in [1e-3,0.1,0.3]:
    for pc,name in [(0,'EF'),(1,'OI')]:
        A=run_online(300,100,30000,3000,mu,r,w,betas,0.0,mm,pc,1.0,[1.0,0.1,0.3],[30000])
        th=A[:,0,:3]; trapped=np.sum(np.abs(th[:,0]-mm)<1e-6); print(f'mu_min={mm} {name}: trapped at mu=mu_min in {trapped}/300 runs; max err mu {np.abs(th[:,0]-1).max():.3f}')
