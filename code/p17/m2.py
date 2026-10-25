import numpy as np, pickle
from lib17 import *
np.set_printoptions(precision=4, suppress=True, linewidth=160)
r,w=powerlaw_truth(0.8,0.5,1.0); betas=np.array([0.3,1.5]); mu=1.0; d=3
ref=pickle.load(open('m1_cache.pkl','rb')); cps=[1000,3000,10000,30000,100000]; R=600
for key,eps in [('free',0.0),('act',0.5)]:
    th=ref[key]['th']; H=ref[key]['H']; J=ref[key]['J']; Lam=ref[key]['Lam']
    c=np.concatenate([[0],1/betas]); V=np.linalg.svd(c[None,:])[2][1:].T if key=='act' else np.eye(d)
    Hi=np.linalg.inv(V.T@H@V); S=V@(Lam*Hi@(V.T@J@V)@Hi)@V.T
    out=run_online(R,2000,100000,5000,mu,r,w,betas,eps,0.3,1,1.0,[2.0,0.05,0.2],cps)
    print('---',key,'OI, R=600; predicted n*Var diag',np.diag(S),' trace',np.trace(S))
    for i,n in enumerate(cps):
        e=out[:,i,:d]-th; C=n*np.cov(out[:,i,:d].T)
        print(f' n={n:6d} n*Var diag {np.diag(C)}  n*MSE {n*np.mean(np.sum(e**2,1)):.2f}  se(var mu)~{np.sqrt(2/(R-1))*C[0,0]:.0f}')
    C=n*np.cov(out[:,-1,:d].T); print(' corr emp %.3f pred %.3f'%(C[0,1]/np.sqrt(C[0,0]*C[1,1]),S[0,1]/np.sqrt(S[0,0]*S[1,1])))
