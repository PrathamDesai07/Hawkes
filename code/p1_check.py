import sys; sys.path.insert(0,'.')
import numpy as np
from p10_fisher_k1 import exact_m, simulate
# (1) exp-concavity counterexample: H = a a^T/lam^2 ; need H >= eta * g g^T, g=b-a/lam
mu,al=1.0,0.5; s=0.7; dt=0.8; beta=1.0
a=np.array([1,s]); th=np.array([mu,al]); lam=a@th
Sint=s*(1-np.exp(-beta*dt))/beta; b=np.array([dt,Sint])
H=np.outer(a,a)/lam**2; g=b-a/lam
print("b parallel to a?", a[0]*b[1]-a[1]*b[0], " rank H:", np.linalg.matrix_rank(H))
for eta in [1e-1,1e-3,1e-6]:
    print(eta, "min eig of H - eta g g^T:", np.linalg.eigvalsh(H-eta*np.outer(g,g)).min())
# (2) K=1 stationary Fisher lambda_min vs eps
mu,b_=1.0,1.0
print("eps  lam_min  lam_min/eps  Imm  Iaa  Ima")
for eps in [.5,.2,.1,.05,.02,.01]:
    a_=(1-eps)*b_; El=mu/eps; m=exact_m(mu,a_,b_)
    I=np.array([[m,(1-mu*m)/a_],[(1-mu*m)/a_,(El-2*mu+mu**2*m)/a_**2]])
    print(eps, np.linalg.eigvalsh(I).min(), np.linalg.eigvalsh(I).min()/eps, I[0,0], I[1,1], I[0,1])
