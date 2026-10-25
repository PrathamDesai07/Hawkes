import numpy as np
from fd18 import *
np.set_printoptions(precision=4, linewidth=170, suppress=True)
ns=0.8; Lam=5.0; Gf,_=truth_gamma(ns,0.5)
for r in [2,3,5]:
    print('r=',r,'  K: bias(nonneg) | bias(signed) | #zero alphas | min signed alpha')
    for K in [1,2,3,4,5,6,8,10,12,16]:
        betas=3.0*float(r)**(-np.arange(K)[::-1].astype(float))
        b1,al1,_=bias(Gf,ns,Lam,betas); b2,al2,_=bias(Gf,ns,Lam,betas,nonneg=False)
        print(f'  K={K:2d}: {b1:9.5f} | {b2:9.5f} | {int(np.sum(al1==0))} | {al2.min():.2e}')
