import numpy as np
from fd18 import *
np.set_printoptions(precision=4, suppress=True, linewidth=170)
nstar=0.8; Lam=1/(1-nstar)
for eta in [0.5]:
    Gf,ns=truth_gamma(nstar,eta); print('eta',eta,'n*',ns)
    print('F1: fixed top beta=3, ratio r, K extends downward; bias n*-n_dagger')
    for r in [2,3,5,10]:
        row=[]
        for K in [1,2,3,4,6,8,10,14,18]:
            betas=3.0*r**(-np.arange(K)[::-1].astype(float)); b,al,q=bias(Gf,ns,Lam,betas); row.append(b)
        print(f' r={r:2d}',' '.join(f'{x:.4f}' for x in row))
    print('   K=      1      2      3      4      6      8     10     14     18')
    print('F2: fixed range [1e-3,10], K geometric points')
    row=[]
    for K in [2,3,4,6,8,10,14,18,24]:
        betas=np.geomspace(1e-3,10,K); b,al,q=bias(Gf,ns,Lam,betas); row.append((K,b,np.exp(np.log(1e4)/(K-1))))
    print(' '.join(f'K={K}:b={b:.4f}(r={r:.2f})' for K,b,r in row))
