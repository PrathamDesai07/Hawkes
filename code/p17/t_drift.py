import numpy as np
from lib17 import *
np.set_printoptions(precision=4, suppress=True, linewidth=150)
gx,gw=np.polynomial.legendre.leggauss(10)
def prep(P):
    P._gl=[]
    for x,wt in zip(gx,gw):
        u=0.5*P.w*(x+1); f=0.5*P.w*wt
        S=P.s*np.exp(-np.outer(u,P.betas)); P._gl.append((f,np.hstack([np.ones((P.n,1)),S])))
def pieces(P,th,thd,pc):
    v=th-thd; T=P.T
    int_vphi=(P.b@v).sum()
    lam_i=P.a@th; ev_term=((P.a@v)/lam_i).sum()
    if pc==0:
        q=sum((f*((a@v)**2/(a@th))).sum() for f,a in P._gl)
    else:
        q=(((P.a@v)/lam_i)**2).sum()
    return q/T, (-2*int_vphi+2*ev_term)/T
def sample_theta(rng,betas,nmax,mu_min,n):
    K=len(betas); out=[]
    for _ in range(n):
        mu=rng.uniform(mu_min,6.0); al=rng.dirichlet(np.ones(K)*rng.choice([0.2,1.0]))*rng.uniform(0,nmax)*betas   # sum alpha/beta ~ n*...
        # rescale so that n = sum alpha/beta in [0,nmax]
        nn=np.sum(al/betas); 
        if nn>nmax: al*=nmax/nn
        out.append(np.concatenate([[mu],al]))
    # vertices / extremes
    for mu in (mu_min,6.0):
        out.append(np.concatenate([[mu],np.zeros(K)]))
        for k in range(K):
            al=np.zeros(K); al[k]=nmax*betas[k]; out.append(np.concatenate([[mu],al]))
    return np.array(out)
r,w=powerlaw_truth(0.8,0.5,1.0); rng=np.random.default_rng(1)
for label,betas,eps in [('K=1 free (strong misspec, mu_dag/mu*~2.7)',np.array([1.0]),0.0),('K=2 free',np.array([0.3,1.5]),0.0),('K=2 active n<=0.5',np.array([0.3,1.5]),0.5)]:
    P=simulate(21,400_000,30000,1.0,r,w,betas); prep(P)
    thd=solve(P,eps=(eps if eps>0 else None),mu_min=0.3)[0]
    nmax=1.0-eps if eps>0 else 0.999
    ratio=(P.a@thd)/P.lstar
    print(f'== {label}: theta_dag={thd}  lam_dag/lam* at events: max {ratio.max():.2f} q99 {np.quantile(ratio,.99):.2f} median {np.median(ratio):.2f}; quiet limit mu_dag/mu*={thd[0]:.2f}')
    ths=sample_theta(rng,betas,nmax,0.3,500)
    for pc,name in [(0,'EF'),(1,'OI')]:
        qs=[];rest=[];nvs=[]
        for th in ths:
            v=th-thd; nv=v@v
            if nv<1e-6: continue
            q,rr=pieces(P,th,thd,pc); qs.append(q); rest.append(rr); nvs.append(nv)
        qs=np.array(qs); rest=np.array(rest); nvs=np.array(nvs); row=[]
        for kappa in (0.25,0.5,1.0,2.0):
            vals=(kappa*qs+rest)/nvs; row.append(f'k={kappa}: max {vals.max():+.4f} (frac>0 {np.mean(vals>0):.3f})')
        print('  ',name,' | '.join(row),flush=True)
