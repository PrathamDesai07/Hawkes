import numpy as np
from lib17 import *
np.set_printoptions(precision=4, suppress=True)
gx,gw=np.polynomial.legendre.leggauss(10)
def Gcum(P):
    """cumulative G_t = int a a'/(1'a) ds over gaps (rows: gap index) -> array (n, d, d)"""
    d=P.K+1; bt=P.betas; out=np.zeros((P.n,d,d))
    for x,wt in zip(gx,gw):
        u=0.5*P.w*(x+1); f=0.5*P.w*wt
        S=P.s*np.exp(-np.outer(u,bt)); a=np.hstack([np.ones((P.n,1)),S]); den=a.sum(1)
        out+= (a[:,:,None]*a[:,None,:])*(f/den)[:,None,None]
    return np.cumsum(out,axis=0), np.cumsum(P.w)
for label,mu,r,w in [('K=2 n=0.97 (eps=0.03)',1.0,[0.3,1.5],[0.12,0.855])]:
    n=sum(a/b for a,b in zip(w,r)); betas=np.array(r)
    P=simulate(3,1_500_000,50000,mu,np.array(r),np.array(w),betas); G,T=Gcum(P)
    print(label,'true n=%.3f'%n,' events',P.n,' time',T[-1])
    for idx in [int(x) for x in (1e3,1e4,1e5,1e6,P.n-1)]:
        Gt=G[idx]; ev=np.linalg.eigvalsh(Gt); lam_min,lam_max=ev[0],ev[-1]; t=T[idx]
        f=(np.log(lam_max))**1.5            # f(x)=x^{1+a}, a=1/2
        print(f'  t={t:10.1f}  lam_min/t={lam_min/t:.4f} lam_max/t={lam_max/t:.3f}  f(log lam_max)/lam_min={f/lam_min:.2e}')
