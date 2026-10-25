import numpy as np
from lib17 import *
np.set_printoptions(precision=5, suppress=True, linewidth=150)
betas=np.array([0.3,1.5]); r=np.array([0.3,1.5]); w=np.array([0.12,0.45]); mu=1.0; K=2
cps=[100,300,1000,3000,10000,30000]
# truth = model: theta*=(1,0.12,0.45); Lambda*=mu/(1-n)
th_true=np.array([mu,0.12,0.45]); Lam=mu/(1-0.7)
P=simulate(1,400000,5000,mu,r,w,betas); IF=IF_matrix(P,th_true); g,H=P.grad_hess(th_true)
print('Lambda',Lam,P.lam_mean(),' I_F vs H max rel diff',np.abs(IF-H).max()/np.abs(H).max())
for pc,name in [(0,'EF'),(1,'OI')]:
    A=run_online(150,100,30000,3000,mu,r,w,betas,0.0,0.3,pc,1.0,[1.0,0.1,0.3],cps)
    print(name)
    for i,n in enumerate(cps):
        e=A[:,i,:3]-th_true; C=np.cov(A[:,i,:3].T)*n; print(' n=%6d  n*var diag'%n,np.diag(C),' mean err',e.mean(0))
    print(' predicted Lam*I^-1 diag',Lam*np.diag(np.linalg.inv(H)))
