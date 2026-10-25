"""Check E exp(eps Z/e) <= 2 exp(4 mu/(e b)) and the tail P(Z>t) <= 2 exp(4mu/(eb)) exp(-eps t/e) by time-weighted simulation."""
import numpy as np, sys
xq, wq = np.polynomial.legendre.leggauss(24); xq = 0.5*(xq+1); wq = 0.5*wq
def run(mu, beta, eps, n_ev, seed=7):
    K=len(beta); alpha=(1-eps)*beta/K; c=alpha/beta; rng=np.random.default_rng(seed)
    th=eps/np.e; ts=np.array([0.5,1,2,4,8])/eps          # tail thresholds in units of 1/eps
    S=np.zeros(K); T=0.0; E=0.0; tail=np.zeros(len(ts)); burn=int(0.1*n_ev)
    for i in range(n_ev):
        lam=mu+alpha@S
        w=0.0
        while True:                                      # thinning with cumulative candidate times
            w+=rng.exponential(1/lam); Sw=S*np.exp(-beta*w)
            if rng.random()*lam < mu+alpha@Sw: break
        if i>=burn:
            s=w*xq; Z=(S[None,:]*np.exp(-np.outer(s,beta)))@c
            E+=w*(wq@np.exp(th*Z)); T+=w
            tail+=np.array([w*(wq@(Z>t)) for t in ts])
        S=Sw+1.0
    return E/T, tail/T, ts, th
mu=1.0
for beta in (np.array([1.]), np.array([1.,2.])):
    for eps in (0.2, 0.1):
        b=beta.min(); E,tail,ts,th=run(mu,beta,eps,150000)
        bd=2*np.exp(4*mu/(np.e*b))
        print(f"beta={beta.tolist()} eps={eps}: E exp(eps Z/e) sim={E:.3f}  bound={bd:.3f}")
        print("   t*eps  P_sim(Z>t)   bound")
        for t,p in zip(ts,tail): print(f"   {t*eps:5.1f}  {p:10.4f}  {min(1,bd*np.exp(-th*t)):8.4f}")
