import numpy as np
from fd18 import *
ns=0.8; Lam=5.0; Gf,_=truth_gamma(ns,0.5)
print('bias n*-n_dagger (signed LS) vs ratio r for fixed K (top beta=3); optimum r marked')
rs=[2,3,4,5,6,8,10,13,17,22,30,40]
print('   r:  '+' '.join(f'{r:7d}' for r in rs))
for K in [2,3,4,5,6,8]:
    row=[]
    for r in rs:
        betas=3.0*float(r)**(-np.arange(K)[::-1].astype(float)); b,_,_=bias(Gf,ns,Lam,betas,nonneg=False); row.append(b)
    i=int(np.argmin(np.abs(row)))
    print(f'K={K}: '+' '.join(f'{x:7.4f}' for x in row)+f'   min|b| at r={rs[i]}')
