"""K=1: exact stationary law of X = eps*lambda via its characteristic function, phi(t) = E exp(i t eps lam) = exp(L(t)),
   L(t) = beta*mu*eps^2 * int_0^t eta d eta / c(-i eps eta),  c(u) = eps*beta*u + h(alpha u),  h(y) = e^-y - 1 + y   (complex series at small |y|).
 Compares: TV(pi_theta, pi_theta') for the ridge pair theta' = (mu+delta, alpha(1 - eps delta/mu)); TV(pi_theta, Gamma comparator); Gamma-Gamma TV;
 and the H^1 distance between characteristic functions  ||phi - G||_{H^1}  (G = (1 - i a t)^{-k}) that enters  TV <= (1/(2 sqrt 2)) ||.||_{H^1}."""
import numpy as np, math, sys
from scipy import stats
def hc(y, nterms=16):
    small = np.abs(y) < 0.5; res = np.empty_like(y); ys = y[small]
    out = np.zeros_like(ys); term = np.ones_like(ys)
    for n in range(1, nterms + 1):
        term = term*(-ys)/n
        if n >= 2: out = out + term
    res[small] = out; res[~small] = np.exp(-y[~small]) - 1 + y[~small]; return res
def phi_grid(mu, beta, eps_theta, eps_scale, t):
    alpha = (1 - eps_theta)*beta; dt = t[1] - t[0]; u = -1j*eps_scale*t
    c = eps_theta*beta*u + hc(alpha*u); g = np.empty_like(t, dtype=complex); g[1:] = t[1:]/c[1:]; g[0] = 1j/(eps_scale**2*beta)*(eps_scale**2/(eps_scale**2))  # limit eta/c at 0
    g[0] = -1.0/(eps_scale**2*(1j*beta))*(-1) if False else (1j/(eps_scale**2*beta))
    cum = np.concatenate([[0], np.cumsum(0.5*(g[1:] + g[:-1])*dt)]); L = beta*mu*eps_scale**2*cum
    Lp = beta*mu*eps_scale**2*g; return np.exp(L), L, Lp           # phi, log phi, d/dt log phi
def density(phi_pos, dt, N, dx):
    # f(x_m) = (dt/2pi) sum_{j=-N/2}^{N/2-1} phi(t_j) e^{-i t_j x_m},  phi(-t)=conj phi(t)
    full = np.empty(N, dtype=complex); half = N//2; full[0:half] = phi_pos[:half]; full[half] = 0.0; full[half+1:] = np.conj(phi_pos[1:half][::-1])
    f = np.real(np.fft.fft(full))*dt/(2*np.pi); return f
if __name__ == "__main__":
    beta, mu = 1.0, float(sys.argv[1]) if len(sys.argv) > 1 else 1.0
    print(f"beta=1, mu={mu} (k0={2*mu}); columns: eps, delta | TV(pi,pi') | TV(pi,Gamma) | TV(pi',Gamma') | TV(Gamma,Gamma') | ||phi-G||_H1 | L1-mass check")
    for eps in (0.1, 0.03, 0.01, 0.003):
        N = 2**21; dx = eps/10; dt = 2*np.pi/(N*dx); t = np.arange(N//2 + 1)*dt; x = np.arange(N)*dx
        alpha = (1 - eps)*beta; a = alpha**2/(2*beta); k = 2*beta*mu/alpha**2
        ph, L, Lp = phi_grid(mu, beta, eps, eps, t); f = density(ph, dt, N, dx)
        G = (1 - 1j*a*t)**(-k); fg = stats.gamma.pdf(x, k, scale=a)
        Gp = 1j*a*k*(1 - 1j*a*t)**(-k - 1)
        h1 = math.sqrt(2*(np.sum(np.abs(ph - G)**2)*dt + np.sum(np.abs(ph*Lp - Gp)**2)*dt))
        tv_g = 0.5*np.sum(np.abs(f - fg))*dx
        for delta in (0.1, 0.3):
            mu2 = mu + delta; al2 = alpha*(1 - eps*delta/mu); eps2 = 1 - al2/beta
            ph2, _, _ = phi_grid(mu2, beta, eps2, eps, t); f2 = density(ph2, dt, N, dx)
            a2 = al2**2/(2*beta); k2 = 2*beta*mu2/al2**2
            # Gamma comparator for theta' lives in the variable eps*lambda' = eps*(X'/eps2): scale a2*eps/eps2
            fg2 = stats.gamma.pdf(x, k2, scale=a2*eps/eps2)
            # NOTE: phi2 is the char. function of eps*lambda' (scale eps), so the Gamma comparator for eps*lambda' has scale a2*eps/eps2
            print(f"  eps={eps:6.3f} delta={delta:.1f} | {0.5*np.sum(np.abs(f - f2))*dx:8.4f} | {tv_g:8.4f} | {0.5*np.sum(np.abs(f2 - fg2))*dx:8.4f} | {0.5*np.sum(np.abs(fg - fg2))*dx:8.4f} | {h1:8.4f} | {np.sum(f)*dx:.4f}", flush=True)
