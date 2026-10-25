import numpy as np
from scipy.optimize import minimize
from p16lib import *
np.set_printoptions(precision=4, suppress=True)
def Phi(P, mu, eps):
    """min over alpha>=0, sum alpha/beta <= 1-eps of E[lam*/(mu+alpha'S)] = (1/T) sum_i 1/(mu+alpha.a_i[1:])"""
    A=P.a[:,1:]; K=P.K; c=1/P.betas
    f=lambda al: (1/(mu+A@al)).sum()/P.T
    g=lambda al: -(A.T@(1/(mu+A@al)**2))/P.T
    best=None
    for al0 in [np.full(K,(1-eps)/K)*P.betas]+[np.eye(K)[k]*(1-eps)*P.betas[k] for k in range(K)]:
        r=minimize(f,al0,jac=g,bounds=[(0,None)]*K,constraints=[{'type':'ineq','fun':lambda al:(1-eps)-c@al,'jac':lambda al:-c}],method='SLSQP',options={'ftol':1e-12,'maxiter':300})
        if best is None or r.fun<best.fun: best=r
    return best.fun, best.x
for nstar,K,ratio,b1 in [(0.8,1,1,1.0),(0.8,2,5,0.3),(0.8,3,5,0.1),(0.95,3,5,0.1)]:
    r,w=powerlaw_truth(nstar,0.5,1.0); betas=b1*ratio**np.arange(K); P=simulate(7,600_000,30000,1.0,r,w,betas)
    th=solve(P)[0]
    print(f'n*={nstar} K={K}: mu_dag={th[0]:.4f}, n_dag={nvec(th,betas):.3f}')
    for eps in [0.3, 0.05]:
        out=[]
        for mu in [0,1e-3,0.05,0.2,0.5,1.0]:
            if mu==0: mu=1e-9
            out.append((mu,Phi(P,mu,eps)[0]))
        print('   eps=%.2f Phi(mu):'%eps, ' '.join(f'{m:.0e}:{v:.3f}' for m,v in out))
