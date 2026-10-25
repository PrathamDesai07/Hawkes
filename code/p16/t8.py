import numpy as np
from p16lib import *
np.set_printoptions(precision=4, suppress=True)
gl_x, gl_w = np.polynomial.legendre.leggauss(10)
def G_Y(P, Lam, mu_max):
    """G = E[a a'/(mu_max+Y)^2], Y=sum beta_k S_k, time average via Gauss-Legendre on each gap."""
    K=P.K; bt=P.betas; w=P.w; s=P.s; G=np.zeros((K+1,K+1))
    for x,wt in zip(gl_x,gl_w):
        u=0.5*w*(x+1); f=0.5*w*wt
        S=s*np.exp(-np.outer(u,bt)); a=np.hstack([np.ones((len(w),1)),S]); Y=S@bt
        G+= (a.T*(f/(mu_max+Y)**2))@a
    return G/P.T
for nstar,Ks in [(0.8,[(1,1,1.0),(2,5,0.3),(3,5,0.1),(4,4,0.05)]),(0.95,[(3,5,0.1)])]:
    r,w=powerlaw_truth(nstar,0.5,1.0)
    for K,ratio,b1 in Ks:
        betas=b1*ratio**np.arange(K); P=simulate(7,1_500_000,30000,1.0,r,w,betas); Lam=P.lam_mean()
        th=solve(P)[0]; g,H=P.grad_hess(th)
        # mu lower bound: sup_x [ rate of events with Y(t-)<=x  - x ],  Y=sum beta_k S_k
        Yi=P.a[:,1:]@betas; Ys=np.sort(Yi); xs=np.quantile(Yi,np.linspace(0.001,0.9,400))
        vals=[(np.searchsorted(Ys,x,side='right')/P.T - x) for x in xs]; lb=max(vals)
        # strong convexity constant m = lam_lo * lambda_min(G_Y), lam_lo = mu* = 1
        Gm=G_Y(P,Lam,Lam); m=1.0*np.linalg.eigvalsh(Gm)[0]
        print(f'n*={nstar} K={K}: mu_dag={th[0]:.4f}  a-priori lower bound={lb:.4f}   lambda_min(H(th))={np.linalg.eigvalsh(H)[0]:.3e}  certified m={m:.3e}  ratio={m/np.linalg.eigvalsh(H)[0]:.3f}')
