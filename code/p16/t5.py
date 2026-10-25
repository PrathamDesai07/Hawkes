import numpy as np
from p16lib import *
np.set_printoptions(precision=5, suppress=True)
nstar, eta, tau = 0.8, 0.5, 1.0
r, w = powerlaw_truth(nstar, eta, tau); mu=1.0
betas = np.array([0.1,0.5,2.5]); K=3
P = simulate(11, 2_000_000, 20000, mu, r, w, betas); Lam=P.lam_mean()
thf,_,_ = solve(P); print('free   ', thf, 'n=%.4f'%nvec(thf,betas), 'mu/Lam=%.4f'%(thf[0]/Lam))
print('--- V5 active constraint')
for eps in [0.4, 0.5, 0.7]:
    th,kap,_ = solve(P, eps=eps)
    g,_ = P.grad_hess(th)
    kap_g = -g[1:]*betas
    print(f'eps={eps} n={nvec(th,betas):.4f} (1-eps={1-eps}) active={nvec(thf,betas)>1-eps}  kappa_barrier={kap:.5f} kappa_from_grad={kap_g} g_mu={g[0]:.1e}')
    print('   mu=%.5f  eps*Lam-kappa(1-eps)=%.5f   E lam_theta=%.4f Lam=%.4f' % (th[0], eps*Lam-kap*(1-eps), th[0]+nvec(th,betas)*Lam, Lam))
print('--- V7 Newton decrement bound (mu*=%g)'%mu)
rng=np.random.default_rng(0)
g0,H0 = P.grad_hess(thf)
for sc in [0.02,0.05,0.1,0.2]:
    th0 = thf*(1+sc*rng.standard_normal(K+1)); th0=np.abs(th0)
    g,H = P.grad_hess(th0); d=np.sqrt(g@np.linalg.solve(H,g)); dist=thf-th0; hn=np.sqrt(dist@H@dist)
    bound = d/(1-d/np.sqrt(mu)) if d<np.sqrt(mu) else np.inf
    print(f'sc={sc} decrement={d:.4f} actual H-dist={hn:.4f} bound={bound:.4f} ok={hn<=bound}')
