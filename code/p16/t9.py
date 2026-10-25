import numpy as np
from p16lib import *
np.set_printoptions(precision=4, suppress=True)
# (a) a-posteriori mu lower bound using Z=alpha'S at theta-dagger
for nstar,K,ratio,b1 in [(0.8,2,5,0.3),(0.8,3,5,0.1),(0.95,3,5,0.1)]:
    r,w=powerlaw_truth(nstar,0.5,1.0); betas=b1*ratio**np.arange(K); P=simulate(7,1_500_000,30000,1.0,r,w,betas)
    th=solve(P)[0]; Zi=P.a[:,1:]@th[1:]; Zs=np.sort(Zi); xs=np.quantile(Zi,np.linspace(0.001,0.9,400))
    lb=max(np.searchsorted(Zs,x,side='right')/P.T-x for x in xs)
    print(f'n*={nstar} K={K}: mu_dag={th[0]:.4f} a-posteriori bound sup_x[R_Z(x)-x]={lb:.4f}')
# (b) sandwich variance, K=2, n*=0.8
r,w=powerlaw_truth(0.8,0.5,1.0); betas=np.array([0.3,1.5])
def sandwich(P,th,B=20):
    g,H=P.grad_hess(th)
    sc=(P.a/(P.a@th)[:,None]-P.b)            # per-gap score contributions (sum = T * (-grad))
    idx=np.array_split(np.arange(P.n),B); Tb=np.array([P.w[i].sum() for i in idx]); Gb=np.array([sc[i].sum(0) for i in idx])
    Gb=Gb-Tb[:,None]*(Gb.sum(0)/Tb.sum())     # center per unit time
    J=(Gb.T@Gb)/(B-1)/Tb.mean()
    Hi=np.linalg.inv(H); return Hi@J@Hi/P.T, Hi/P.T, H, J
def study(truth_r,truth_w,mu,betas,nev,R,label):
    ests=[];sw=[];nv=[]
    for s in range(R):
        P=simulate(5000+s,nev,5000,mu,truth_r,truth_w,betas); th=solve(P)[0]; ests.append(th)
        cs,cn,_,_=sandwich(P,th); sw.append(cs); nv.append(cn)
    ests=np.array(ests); emp=np.cov(ests.T); sw=np.mean(sw,0); nv=np.mean(nv,0)
    print(label,'\n  empirical var diag',np.diag(emp),'\n  sandwich  var diag',np.diag(sw),'\n  naive H^-1/T diag',np.diag(nv))
    print('  empirical corr(mu,a1)=%.3f sandwich %.3f'%(emp[0,1]/np.sqrt(emp[0,0]*emp[1,1]),sw[0,1]/np.sqrt(sw[0,0]*sw[1,1])))
study(r,w,1.0,betas,40000,300,'misspecified power-law truth, K=2')
# control: well-specified (K=2 truth exactly the model)
study(np.array([0.3,1.5]),np.array([0.12,0.45]),1.0,betas,40000,300,'well-specified K=2 (alpha=(0.12,0.45), n=%.2f)'%((0.12/0.3)+(0.45/1.5)))
