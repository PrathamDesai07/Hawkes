import numpy as np
from math import gamma
from scipy.integrate import quad
from fd18 import *
def kap(eta,r,b1,ns=0.8,signed=True):
    Lam=1/(1-ns); Gf,_=truth_gamma(ns,eta); K=int(round(np.log(3/b1)/np.log(r)))+1; betas=b1*float(r)**np.arange(K)
    b,al,q=bias(Gf,ns,Lam,betas,nonneg=not signed); return b/(ns*b1**eta/gamma(1+eta))
print('kappa(eta,r) = bias / n_low, n_low = n* (tau beta1)^eta/Gamma(1+eta); signed LS; columns beta1=1e-4,1e-5,1e-6')
for eta in [0.15,0.25,0.35,0.5]:
    k1=1-2*eta*quad(lambda x:x**eta/(1+x),0,1)[0]
    for r in [2,3,4]:
        print(f' eta={eta} r={r}:',[round(kap(eta,r,b1),3) for b1 in (1e-4,1e-5,1e-6)],f'  (single-pole heuristic {k1:.3f})')
