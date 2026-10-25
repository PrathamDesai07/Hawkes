import sys, os, numpy as np
from scipy.optimize import nnls
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'p16'))
def truth_discrete(r, w):
    """G(s)=sum w_j/(r_j+s); returns function of s (complex array), n*, phi(0)"""
    r=np.asarray(r); w=np.asarray(w)
    return (lambda s: (w[None,:]/(r[None,:]+s[:,None])).sum(1)), float((w/r).sum())
def truth_gamma(nstar, eta, tau=1.0, lo=-30.0, hi=8.0, h=0.02):
    """nu(d beta) ~ beta^eta exp(-beta*tau) dbeta rescaled so that int nu/beta = nstar. log-grid trapezoid."""
    x=np.arange(lo,hi,h); b=np.exp(x); wt=b**(eta+1)*np.exp(-b*tau)*h        # nu(b)*b dx ,  nu=b^eta e^{-b tau}
    wt*=nstar/np.sum(wt/b)
    return truth_discrete(b,wt)
def omega_grid(lo=-24.0,hi=14.0,h=0.05):
    y=np.arange(lo,hi,h); om=np.exp(y); return om, om*h      # d omega = omega dy
def fit(Gfun, nstar, Lam, betas, om=None, dom=None, nonneg=True):
    if om is None: om,dom=omega_grid()
    s=1j*om; G=Gfun(s); wgt=Lam/np.abs(1-G)**2*dom/np.pi            # (1/pi) int_0^inf
    Gk=1.0/(betas[None,:]+s[:,None])
    sw=np.sqrt(wgt)
    A=np.vstack([(sw[:,None]*Gk).real,(sw[:,None]*Gk).imag]); y=np.concatenate([(sw*G).real,(sw*G).imag])
    if nonneg: al,res=nnls(A,y,maxiter=20000)
    else: al=np.linalg.lstsq(A,y,rcond=None)[0]; res=np.linalg.norm(A@al-y)
    return al, res**2
def bias(Gfun,nstar,Lam,betas,**kw):
    al,q=fit(Gfun,nstar,Lam,betas,**kw); return nstar-np.sum(al/betas), al, q
if __name__=='__main__':
    # E1: validation against simulation (truncated 24-exponential truth, n*=0.8, Lambda=5)
    from p16lib import powerlaw_truth
    r,w=powerlaw_truth(0.8,0.5,1.0); Gf,ns=truth_discrete(r,w); Lam=1.0/(1-ns)
    print('n*',ns,'Lambda',Lam)
    for betas,sim in [([1.0],0.463),([0.3,1.5],0.621),([0.1,0.5,2.5],0.697),([0.05,0.2,0.8,3.2],0.733)]:
        b,al,q=bias(Gf,ns,Lam,np.array(betas)); n_dag=ns-b
        print(f'betas={betas}: n_Q (frequency domain)={n_dag:.4f}  simulated KL n_dagger={sim}   mu={1+Lam*b:.4f}  alpha={al}')
