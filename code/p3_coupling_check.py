"""Ordered coupling of two Hawkes (K=2) copies on one Poisson random measure; check D=S^max-S^min stays >=0
and E[Delta Z_t] = c^T expm(M t) D0 (exact linear ODE), where Delta Z = sum (alpha_k/beta_k) D_k."""
import numpy as np
from scipy.linalg import expm
mu=1.0; beta=np.array([1.,2.]); eps=0.2; alpha=(1-eps)*beta/2; c=alpha/beta
M=-np.diag(beta)+np.outer(np.ones(2),alpha)
rng=np.random.default_rng(5); R=3000; T=12.0; ts=np.array([2,4,6,8,10,12.])
S0max=np.array([6.,6.]); S0min=np.array([1.,1.])
acc=np.zeros(len(ts)); neg=0
for _ in range(R):
    Sx=S0max.copy(); Sn=S0min.copy(); t=0.0; vals=[]; ti=0
    while ti<len(ts):
        lx=mu+alpha@Sx; w=rng.exponential(1/lx)
        tn=min(t+w, ts[ti])
        d=tn-t; Sx=Sx*np.exp(-beta*d); Sn=Sn*np.exp(-beta*d)
        if tn==ts[ti] and t+w>ts[ti]:
            acc[ti]+=c@(Sx-Sn); neg+=(Sx<Sn-1e-12).any(); ti+=1; t=tn; continue   # memoryless: redraw
        t=tn
        u=rng.random()*lx
        # intensities at event time (state already decayed)
        lx2=mu+alpha@Sx; ln2=mu+alpha@Sn
        if u<=lx2: Sx=Sx+1
        if u<=ln2: Sn=Sn+1
acc/=R
exact=[c@expm(M*t)@(S0max-S0min) for t in ts]
print("t   sim E[dZ]   ODE exact   (order violations:",neg,")")
for t,a,e in zip(ts,acc,exact): print(f"{t:4.0f} {a:10.4f} {e:10.4f}")
