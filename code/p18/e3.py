import numpy as np
from math import gamma
from fd18 import *
def kappa(eta,nstar,r,beta1,top=3.0):
    Lam=1/(1-nstar); Gf,ns=truth_gamma(nstar,eta)
    K=int(round(np.log(top/beta1)/np.log(r)))+1; betas=beta1*r**np.arange(K)
    b,al,q=bias(Gf,ns,Lam,betas)
    c=ns/gamma(eta)                                 # nu(beta)~c beta^eta ; int_0^b1 nu/beta = c b1^eta/eta
    return b/(c*beta1**eta/eta), b, K
np.set_printoptions(precision=4)
print('exponent check eta=0.5,n*=0.8,r=2: b vs beta1')
for b1 in [1e-3,1e-4,1e-5,1e-6]: k,b,K=kappa(0.5,0.8,2.0,b1); print(f' beta1={b1:g} K={K} b={b:.3e} b/beta1^eta={b/b1**0.5:.4f} kappa={k:.4f}')
print('kappa(eta) at n*=0.8, r=1.5, beta1=1e-5')
for eta in [0.2,0.35,0.5,0.75,1.0]:
    print(' eta',eta,[round(kappa(eta,0.8,r,1e-5)[0],4) for r in (1.3,1.6,2.0,3.0)])
print('kappa vs n* (eta=.5, r=1.5)')
for ns in [0.3,0.5,0.8,0.9,0.95]: print(' n*',ns,round(kappa(0.5,ns,1.5,1e-5)[0],4))
