import numpy as np
from math import gamma
from fd18 import *
np.set_printoptions(precision=4, linewidth=170, suppress=True)
eta=0.5; ns=0.8; Lam=5.0; Gf,_=truth_gamma(ns,eta); C0=ns/gamma(eta)
r=1.3; b1=1e-5; K=int(np.log(3/b1)/np.log(r))+1; betas=b1*r**np.arange(K)
b,al,q=bias(Gf,ns,Lam,betas)
wq=C0*betas**eta*np.exp(-betas)*betas*np.log(r)
print('K',K,'bias',b)
print('k, beta_k, alpha_k/quadrature weight (first 25):')
print(np.round(al[:25]/wq[:25],3))
print('later:',np.round(al[25:60:3]/wq[25:60:3],3))
print('extra mass sum alpha_k/beta_k - sum wq/beta_k over k (in units of n_low):', (np.sum(al/betas)-np.sum(wq/betas)), ' n_low=',C0*b1**eta/eta)
# sum of the excess mass by index
ex=(al-wq)/betas; print('cumulative excess/ n_low:',np.round(np.cumsum(ex)[:25]/(C0*b1**eta/eta),3))
