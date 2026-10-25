import numpy as np, pickle
from lib17 import *
np.set_printoptions(precision=4, suppress=True, linewidth=160)
r,w=powerlaw_truth(0.8,0.5,1.0); betas=np.array([0.3,1.5]); mu=1.0; d=3
ref=pickle.load(open('m1_cache.pkl','rb'))
for key in ['free','act']:
    th=ref[key]['th']; H=ref[key]['H']; Lam=ref[key]['Lam']
    J2,H2=J_population(r,w,mu,betas,th,nev=100000,R=600,seed0=70000)
    c=np.concatenate([[0],1/betas]); V=np.linalg.svd(c[None,:])[2][1:].T if key=='act' else np.eye(d)
    Hi=np.linalg.inv(V.T@H@V); S=V@(Lam*Hi@(V.T@J2@V)@Hi)@V.T
    print(key,'J old diag',np.diag(ref[key]['J']),'J new diag',np.diag(J2),'\n  pred n*Var diag (new J)',np.diag(S),' corr %.3f'%(S[0,1]/np.sqrt(S[0,0]*S[1,1])),flush=True)
    ref[key]['J2']=J2
pickle.dump(ref,open('m1_cache.pkl','wb'))
