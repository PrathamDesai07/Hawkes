import numpy as np, pickle, os
from lib17 import *
np.set_printoptions(precision=5, suppress=True, linewidth=160)
r,w=powerlaw_truth(0.8,0.5,1.0); betas=np.array([0.3,1.5]); mu=1.0; K=2; d=3
cache='m1_cache.pkl'
if os.path.exists(cache): ref=pickle.load(open(cache,'rb'))
else:
    ref={}
    Ps=[simulate(300+s,1_500_000,20000,mu,r,w,betas) for s in range(2)]
    for key,eps in [('free',0.0),('act',0.5)]:
        ths=[solve(P,eps=(eps if eps>0 else None),mu_min=0.3)[0] for P in Ps]; th=np.mean(ths,0); ref[key]=dict(th=th)
        P=Ps[0]; ref[key]['IF']=IF_matrix(P,th); g,H=P.grad_hess(th); ref[key]['H']=H; ref[key]['Lam']=P.lam_mean()
        ref[key]['J'],_=J_population(r,w,mu,betas,th,R=200)
        print(key,'theta',th,'Lam',P.lam_mean(),'n',nvec(th,betas),flush=True)
    pickle.dump(ref,open(cache,'wb'))
def lyap(A,Ibar,J,Lam):
    Ii=np.linalg.inv(Ibar); Q=Lam*Ii@J@Ii; n=A.shape[0]
    return solve_continuous_lyapunov(A-0.5*np.eye(n),Q)  # (A-1/2)S+S(A-1/2)'=Q
cps=[100,300,1000,3000,10000,30000,100000]
R=200
for key,eps in [('free',0.0),('act',0.5)]:
    th=ref[key]['th']; H=ref[key]['H']; IF=ref[key]['IF']; J=ref[key]['J']; Lam=ref[key]['Lam']
    c=np.concatenate([[0],1/betas])
    if key=='act':
        u,s_,vt=np.linalg.svd(c[None,:]); V=vt[1:].T
    else: V=np.eye(d)
    for pc,name,Ibar in [(1,'OI',H),(0,'EF',IF)]:
        A=V.T@Ibar@V; A=np.linalg.inv(A)@V.T@H@V          # restricted drift matrix
        ev=np.linalg.eigvals(A).real
        IV=V.T@Ibar@V; Ji=np.linalg.inv(IV); JV=V.T@J@V
        S=solve_continuous_lyapunov(A-0.5*np.eye(len(A)),Lam*Ji@JV@Ji)
        Spred=V@S@V.T
        out=run_online(R,500,100000,5000,mu,r,w,betas,eps,0.3,pc,1.0,[2.0,0.05,0.2],cps)
        print(f'--- {key} {name}: eig(A)={ev}  (need > 1/2)')
        for i,n in enumerate(cps):
            e=out[:,i,:d]-th; C=n*np.cov(out[:,i,:d].T)
            print(f' n={n:6d} n*Var diag {np.diag(C)}  bias*sqrt(n) {e.mean(0)*np.sqrt(n)}  active frac {out[:,i,d].mean()/n:.4f}  n*dist {np.mean(n*out[:,i,d+1]):.3f}')
        print('  predicted n*Var diag',np.diag(Spred))
        C=n*np.cov(out[:,-1,:d].T); print('  corr(mu,a1) emp %.3f pred %.3f'%(C[0,1]/np.sqrt(C[0,0]*C[1,1]),Spred[0,1]/np.sqrt(Spred[0,0]*Spred[1,1])),flush=True)
