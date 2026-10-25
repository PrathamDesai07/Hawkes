"""Problem 10 caveat 2: total Fisher information over [0,T] from a cold start (S_0=0, lam_0=mu), K=1.
 I_T = E int_0^T a a^T/lam dt.  With M_T=int E[1/lam_t]dt, L_T=int E lam_t dt:
   I_mm=M_T, I_ma=(T-mu M_T)/alpha, I_aa=(L_T-2 mu T+mu^2 M_T)/alpha^2.
 Exact transient (characteristics of d_t phi=-u beta mu phi-c(u) d_u phi, c(u)=eps beta u+h(alpha u)):
   M_T = int_0^inf du0 e^{-mu u0}/(c(u0) phi(u0)) * int_{u0}^{u_T(u0)} phi(u)du,  Theta(u_T)-Theta(u0)=T, Theta'=1/c."""
import numpy as np
from p10_fisher_k1_exact import h

def grids(mu, beta, eps, N=60000, lo=1e-12, hi=None):
    alpha = (1-eps)*beta; hi = 40/mu if hi is None else hi
    lx = np.linspace(np.log(lo), np.log(hi), N); x = np.exp(lx); dl = lx[1]-lx[0]
    c = eps*beta*x + h(alpha*x); G = x/c
    cum = lambda f: np.concatenate([[0], np.cumsum(0.5*(f[1:]*x[1:]+f[:-1]*x[:-1])*dl)])
    Gam = cum(G) + lo*G[0]                    # int_0^x G (tiny head added)
    phi = np.exp(-beta*mu*Gam)
    Th = cum(1/c)                             # Theta(x) up to constant
    Cphi = cum(phi) + lo
    return dict(x=x, lx=lx, c=c, phi=phi, Th=Th, Cphi=Cphi, alpha=alpha, dl=dl)

def M_T(mu, beta, eps, T, g=None, lam0=None):
    g = g or grids(mu, beta, eps); lam0 = mu if lam0 is None else lam0
    x, lx, c, phi, Th, Cphi = g['x'], g['lx'], g['c'], g['phi'], g['Th'], g['Cphi']
    tgt = Th + np.atleast_1d(T)[:, None]
    uT_l = np.interp(tgt, Th, lx)            # clamps at grid max
    inner = np.interp(uT_l, lx, Cphi) - Cphi
    f = np.exp(-lam0*x)/(c*phi)*inner
    integ = 0.5*(f[:,1:]*x[1:]+f[:,:-1]*x[:-1])*g['dl']
    return integ.sum(axis=1)

def L_T(mu, beta, eps, T):
    return mu/eps*T - (mu/eps-mu)*(1-np.exp(-beta*eps*T))/(beta*eps)

def info_matrix(mu, beta, eps, T, g=None):
    a = (1-eps)*beta; M = M_T(mu, beta, eps, T, g)[0]; L = L_T(mu, beta, eps, T)
    return np.array([[M, (T-mu*M)/a], [(T-mu*M)/a, (L-2*mu*T+mu**2*M)/a**2]]), M, L

def mc_M(mu, beta, eps, T, R=4000, seed=3):
    rng = np.random.default_rng(seed); a = (1-eps)*beta; tot = 0.0
    for _ in range(R):
        S = 0.0; t = 0.0; I = 0.0
        while t < T:
            lam = mu + a*S
            w = 0.0
            while True:                                   # thinning with cumulative candidate times
                w += rng.exponential(1/lam); Sw = S*np.exp(-beta*w)
                if rng.random()*lam < mu + a*Sw: break
            w2 = min(w, T-t); Se = S*np.exp(-beta*w2)
            # int_0^w2 dt/(mu+a S e^{-beta t}) = (w2 + ln((mu+a Se)/(mu+a S))/beta)/mu
            I += (w2 + np.log((mu+a*Se)/(mu+a*S))/beta)/mu
            t += w2; S = Se + 1.0 if w2 == w else Se
        tot += I
    return tot/R

if __name__ == "__main__":
    mu, beta = 1.0, 1.0
    print("validation: exact transient formula vs Monte Carlo (mu=beta=1)")
    for eps, T in ((0.2, 5.0), (0.2, 20.0), (0.5, 10.0)):
        g = grids(mu, beta, eps)
        print(f" eps={eps} T={T}: formula {M_T(mu,beta,eps,T,g)[0]:.4f}   MC {mc_M(mu,beta,eps,T):.4f}   (T/mu={T/mu})")
    g = grids(mu, beta, 0.1); print(" T large: M_T/T =", M_T(mu,beta,0.1,5000.0,g)[0]/5000, " vs stationary m =", 0.15506)
