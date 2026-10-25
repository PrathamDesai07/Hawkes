import numpy as np
from lib17 import *
np.set_printoptions(precision=4, suppress=True)
for nstar in [0.8,0.95]:
    r,w=powerlaw_truth(nstar,0.5,1.0)
    for betas in [np.array([0.05]),np.array([1.0]),np.array([10.0]),np.array([0.05,10.0])]:
        P=simulate(3,800_000,30000,1.0,r,w,betas); th=solve(P,mu_min=0.01)[0]; g,H=P.grad_hess(th); IF=IF_matrix(P,th)
        A=np.linalg.inv(IF)@H; print(f'n*={nstar} betas={betas} theta={th} n={nvec(th,betas):.3f} eig(I_F^-1 H)={np.sort(np.linalg.eigvals(A).real)}')
