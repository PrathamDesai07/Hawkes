import numpy as np
from lib17 import *
np.set_printoptions(precision=4, suppress=True, linewidth=150)
betas=np.array([0.3,1.5]); r=np.array([0.3,1.5]); w=np.array([0.12,0.45]); mu=1.0
cps=[1000,10000,30000]; th_true=np.array([mu,0.12,0.45])
A=run_online(150,100,30000,3000,mu,r,w,betas,0.0,1e-3,1,1.0,[1.0,0.1,0.3],cps)
for i,n in enumerate(cps):
    th=A[:,i,:3]; e=th-th_true
    print(n,'sd per coord',e.std(0),'MAD-sd',1.4826*np.median(np.abs(e-np.median(e,0)),0),'max|e|',np.abs(e).max(0),'n_active mean',A[:,i,3].mean())
    idx=np.argsort(-np.abs(e[:,0]))[:5]; print('  worst mu errors',e[idx,0], ' their active counts',A[idx,i,3])
