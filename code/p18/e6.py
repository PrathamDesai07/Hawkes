import sys, os, numpy as np
sys.path.insert(0,'../p16')
from p16lib import *
from fd18 import *
np.set_printoptions(precision=4, linewidth=170, suppress=True)
r_,w_=powerlaw_truth(0.8,0.5,1.0); Gf,ns=truth_discrete(r_,w_); Lam=5.0
print('truncated 24-exp truth, n*=%.3f. Grid beta_k = 3*3^-(K-k); KL by batch QMLE on 1.2e6 events, paired with sample-LS surrogate'%ns)
for K in [2,3,4,5,6,7]:
    betas=3.0*3.0**(-np.arange(K)[::-1].astype(float))
    bQ,alQ,_=bias(Gf,ns,Lam,betas,nonneg=False)
    out=[]
    for seed in (11,12):
        P=simulate(seed,1_200_000,30000,1.0,r_,w_,betas)
        th=solve(P,mu_min=0.3)[0]; tq=quad_surrogate(P)
        out.append((nvec(th,betas)-nvec(tq,betas), nvec(th,betas), nvec(tq,betas)))
    d=np.mean([o[0] for o in out])
    print(f'K={K}: Q exact n={ns-bQ:.4f} (bias {bQ:.4f}); paired KL-Q shift {out[0][0]:+.4f},{out[1][0]:+.4f}; KL estimate bias {bQ-d:.4f}   [raw KL n: {out[0][1]:.4f},{out[1][1]:.4f}; raw sampleLS n: {out[0][2]:.4f},{out[1][2]:.4f}]',flush=True)
