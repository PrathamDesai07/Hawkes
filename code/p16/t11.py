import numpy as np, sys
from p16lib import *
np.set_printoptions(precision=5, suppress=True, linewidth=150)
r,w=powerlaw_truth(0.8,0.5,1.0); betas=np.array([0.3,1.5]); eps=0.5; K=2
def blocks(P,th,B=20):
    sc=(P.a/(P.a@th)[:,None]-P.b)
    idx=np.array_split(np.arange(P.n),B); Tb=np.array([P.w[i].sum() for i in idx]); Gb=np.array([sc[i].sum(0) for i in idx])
    Gb=Gb-Tb[:,None]*(Gb.sum(0)/Tb.sum()); return (Gb.T@Gb)/(B-1)/Tb.mean()
def cons_cov(P,th,act):
    g,H=P.grad_hess(th); J=blocks(P,th)
    C=np.array(act)                       # rows = active constraint normals
    if len(C)==0: V=np.eye(K+1)
    else:
        u,s,vt=np.linalg.svd(C); V=vt[len(C):].T
    M=np.linalg.inv(V.T@H@V); return V@M@V.T@J@V@M@V.T/P.T, np.linalg.inv(H)@J@np.linalg.inv(H)/P.T
nev=40000; R=300; est=[];cc=[];sw=[];kap=[]
c=np.concatenate([[0],1/betas])
for s in range(R):
    P=simulate(9000+s,nev,5000,1.0,r,w,betas); th,k,_=solve(P,eps=eps); est.append(th); kap.append(k)
    a,b=cons_cov(P,th,[c]); cc.append(a); sw.append(b)
est=np.array(est); emp=np.cov(est.T); cc=np.mean(cc,0); sw=np.mean(sw,0)
print('n_hat range', min(nvec(t,betas) for t in est), max(nvec(t,betas) for t in est))
print('empirical var diag   ',np.diag(emp)); print('active-set formula   ',np.diag(cc)); print('unconstrained sandwich',np.diag(sw))
print('empirical corr(mu,a1) %.3f formula %.3f'%(emp[0,1]/np.sqrt(emp[0,0]*emp[1,1]), cc[0,1]/np.sqrt(cc[0,0]*cc[1,1])))
print('mean kappa',np.mean(kap),'sd',np.std(kap))
print('rank check: empirical eigenvalues',np.linalg.eigvalsh(emp)); print('formula eigenvalues',np.linalg.eigvalsh(cc))
