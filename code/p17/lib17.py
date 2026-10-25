import os, sys, subprocess, tempfile, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'p16'))
os.environ['PATH'] = r'C:\msys64\ucrt64\bin;' + os.environ['PATH']
from p16lib import *
from scipy.linalg import solve_continuous_lyapunov
EXE17 = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'online17.exe')
gl_x, gl_w = np.polynomial.legendre.leggauss(10)

def run_online(nrep, seed0, nev, warm, mu, r, w, betas, eps, mu_min, precond, rho0, th0, cps, kappa=1.0):
    K = len(betas); cfg = tempfile.mktemp(suffix='.txt'); out = tempfile.mktemp(suffix='.bin')
    toks = [nrep, seed0, nev, warm, float(mu), len(r)] + [x for p in zip(r, w) for x in p] + [K] + list(betas) + [eps, mu_min, precond, rho0, kappa] + list(th0) + [len(cps)] + list(cps)
    open(cfg, 'w').write(' '.join(repr(float(t)) if isinstance(t, (float, np.floating)) else str(int(t)) for t in toks))
    subprocess.run([EXE17, cfg, out], check=True)
    A = np.fromfile(out, dtype=np.float64).reshape(nrep, len(cps), K + 3); os.remove(cfg); os.remove(out)
    return A  # [..., :K+1]=theta, K+1 = nactive, K+2 = dist

def IF_matrix(P, th):
    """I_F = E[a a'/lam_theta] time average, Gauss-Legendre per gap."""
    K, bt = P.K, P.betas; d = K + 1; IF = np.zeros((d, d))
    for x, wt in zip(gl_x, gl_w):
        u = 0.5 * P.w * (x + 1); f = 0.5 * P.w * wt
        S = P.s * np.exp(-np.outer(u, bt)); a = np.hstack([np.ones((P.n, 1)), S]); lam = a @ th
        IF += (a.T * (f / lam)) @ a
    return IF / P.T

def J_population(truth_r, truth_w, mu, betas, th, nev=40000, R=300, seed0=9000):
    G = []; Hs = []; Ts = []
    for s in range(R):
        P = simulate(seed0 + s, nev, 5000, mu, truth_r, truth_w, betas)
        g, H = P.grad_hess(th); G.append(g); Hs.append(H); Ts.append(P.T)
    return np.cov(np.array(G).T) * np.mean(Ts), np.mean(Hs, 0)
