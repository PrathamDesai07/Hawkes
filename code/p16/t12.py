import numpy as np
from p16lib import *
np.set_printoptions(precision=5, suppress=True, linewidth=150)
r,w=powerlaw_truth(0.8,0.5,1.0); betas=np.array([0.3,1.5]); eps=0.5; K=2
c=np.concatenate([[0],1/betas]); nev=40000
est=[];G=[];Hs=[];Ts=[]
# reference theta: constrained solve on a long run
ref=np.mean([solve(simulate(100+s,1_000_000,5000,1.0,r,w,betas),eps=eps)[0] for s in range(2)],0); print('ref',ref)
for s in range(300):
    P=simulate(9000+s,nev,5000,1.0,r,w,betas); est.append(solve(P,eps=eps)[0])
    g,H=P.grad_hess(ref); G.append(g); Hs.append(H); Ts.append(P.T)
est=np.array(est); G=np.array(G); H=np.mean(Hs,0); T=np.mean(Ts)
J=np.cov(G.T)*T                     # long-run variance of per-time gradient
u,s_,vt=np.linalg.svd(c[None,:]); V=vt[1:].T
M=np.linalg.inv(V.T@H@V); C=V@M@V.T@J@V@M@V.T/T
emp=np.cov(est.T); print('empirical var diag',np.diag(emp)); print('formula (replicate J) ',np.diag(C))
print('corr(mu,a1) emp %.3f formula %.3f'%(emp[0,1]/np.sqrt(emp[0,0]*emp[1,1]),C[0,1]/np.sqrt(C[0,0]*C[1,1])))
print('E grad at ref (should be ~ -kappa*c):',G.mean(0), ' kappa from this:', -G.mean(0)[1:]*betas)
