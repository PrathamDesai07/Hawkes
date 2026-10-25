import numpy as np, sys
from p16lib import *
np.set_printoptions(precision=4, suppress=True)
nstar, eta, tau = 0.8, 0.5, 1.0
r, w = powerlaw_truth(nstar, eta, tau)
mu = 1.0
print('check n*', (w/r).sum())
NEV = 1_500_000
for K, ratio, b1 in [(1,1,1.0),(2,5,0.3),(3,5,0.1),(4,4,0.05)]:
    betas = b1*ratio**np.arange(K)
    P = simulate(7, NEV, 20000, mu, r, w, betas)
    Lam = P.lam_mean()
    th, _, _ = solve(P)
    th2, _, _ = solve(P, th0=np.concatenate([[3*th[0]], np.full(K, 0.01)]))  # different start
    tq = quad_surrogate(P)
    nk = nvec(th, betas)
    print(f'K={K} betas={betas} Lam={Lam:.4f} (true {mu/(1-nstar):.3f})')
    print('  KL theta  ', th, ' n=%.4f mu/(1-n)=%.4f' % (nk, th[0]/(1-nk)), ' restart diff %.1e' % np.abs(th-th2).max())
    print('  quad theta', tq, ' n=%.4f mu/(1-n)=%.4f' % (nvec(tq,betas), tq[0]/(1-nvec(tq,betas))))
    g,H = P.grad_hess(th); print('  |grad|=%.1e  eig(H)=' % np.abs(g).max(), np.linalg.eigvalsh(H))
