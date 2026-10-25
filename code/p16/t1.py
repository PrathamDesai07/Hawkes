import numpy as np, time
from p16lib import *
# V1/V2: well-specified K=1: mu*=1, alpha*=0.7, beta=1 (n=0.7)
mu, al, be = 1.0, 0.7, 1.0
t0=time.time(); P = simulate(1, 400000, 5000, mu, [be], [al], [be]); print('sim', time.time()-t0, P.n, P.T)
th, k, tau = solve(P); print('theta', th, 'truth', mu, al, 'n', nvec(th,P.betas))
print('rate identity: mu_hat/(1-n_hat)=', th[0]/(1-nvec(th,P.betas)), ' N/T=', P.lam_mean(), 'sum lam dt / N =', (P.b@th).sum()/P.n)
