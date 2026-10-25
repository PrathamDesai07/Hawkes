import os, subprocess, numpy as np, tempfile
HERE = os.path.dirname(os.path.abspath(__file__))
EXE = os.path.join(HERE, 'sim16.exe')
os.environ['PATH'] = r'C:\msys64\ucrt64\bin;' + os.environ['PATH']

def powerlaw_truth(nstar, eta, tau, M=24, umin=1e-3, umax=60.0):
    """phi(t)=sum w_j exp(-r_j t) approximating nstar*eta/tau*(1+t/tau)^-(1+eta); rescaled so sum w/r = nstar."""
    u = np.geomspace(umin, umax, M); lu = np.log(u)
    wq = u**(eta + 1) * np.exp(-u) * np.gradient(lu)          # int u^eta e^-u du in log-grid
    r = u / tau; w = wq.copy()
    w *= nstar / np.sum(w / r)
    return r, w

def simulate(seed, nev, warm, mu, r, w, betas, tag=''):
    out = os.path.join(tempfile.gettempdir(), f'p16_{os.getpid()}_{seed}{tag}.bin')
    args = [EXE, str(seed), str(nev), str(warm), repr(mu), str(len(r))]
    for a, b in zip(r, w): args += [repr(float(a)), repr(float(b))]
    args += [str(len(betas))] + [repr(float(b)) for b in betas] + [out]
    subprocess.run(args, check=True)
    K = len(betas)
    A = np.fromfile(out, dtype=np.float64).reshape(-1, 2 + K); os.remove(out)
    return Path(A, np.asarray(betas, float))

class Path:
    """Per-gap features. w_i gap length, lam_star at t_i-, s_i = S(t_{i-1}+). a_i=(1,s e^{-b w}), b_i=(w, s(1-e^{-b w})/b)."""
    def __init__(self, A, betas):
        self.betas = betas; self.w = A[:, 0]; self.lstar = A[:, 1]; s = A[:, 2:]
        e = np.exp(-np.outer(self.w, betas))
        n = len(self.w)
        self.a = np.hstack([np.ones((n, 1)), s * e])
        self.b = np.hstack([self.w[:, None], s * (1 - e) / betas])
        self.T = self.w.sum(); self.n = n; self.K = len(betas)
        self.s = s
    def lam_mean(self): return self.n / self.T
    def obj(self, th):
        lam = self.a @ th
        return (self.b @ th).sum() / self.T - np.log(lam).sum() / self.T
    def grad_hess(self, th):
        lam = self.a @ th
        g = (self.b.sum(0) - (self.a / lam[:, None]).sum(0)) / self.T
        H = (self.a.T * (1 / lam**2)) @ self.a / self.T
        return g, H
    def oracle_grad_hess(self, th):
        """population-form gradient/Hessian using lambda* at the event times weighted by 1/lambda*: E[a(1-l*/l)] ~ sum over time. Uses the Palm identity only, so equals empirical (check only)."""
        raise NotImplementedError

def nvec(th, betas): return float(np.sum(th[1:] / betas))

def solve(P, eps=None, mu_min=1e-6, th0=None, verbose=False):
    """Minimise empirical KL objective on mu>=mu_min, alpha>=0, (n<=1-eps if eps) by log-barrier Newton. Returns th, kappa estimate, info."""
    K, betas = P.K, P.betas
    th = np.array(th0, float) if th0 is not None else np.concatenate([[0.5 * P.lam_mean()], np.full(K, 0.3 * P.lam_mean() * betas.min() / K)])
    if eps is not None:
        # strictly feasible start
        c = 1 / betas; sc = (1 - eps) * 0.5 / (c @ th[1:]) if c @ th[1:] > 0 else 1
        th[1:] = th[1:] * min(1, sc)
    def feas(x):
        if x[0] <= mu_min or np.any(x[1:] <= 0): return False
        return eps is None or (1 - eps - nvec(x, betas)) > 0
    def phi(x, tau):
        v = P.obj(x) - tau * (np.log(x[0] - mu_min) + np.log(x[1:]).sum())
        if eps is not None: v -= tau * np.log(1 - eps - nvec(x, betas))
        return v
    tau = 1e-2
    while True:
        for it in range(100):
            g, H = P.grad_hess(th)
            g = g.copy(); H = H.copy()
            g[0] -= tau / (th[0] - mu_min); H[0, 0] += tau / (th[0] - mu_min)**2
            for k in range(K): g[1 + k] -= tau / th[1 + k]; H[1 + k, 1 + k] += tau / th[1 + k]**2
            if eps is not None:
                sl = 1 - eps - nvec(th, betas); cv = np.concatenate([[0], 1 / betas])
                g += tau * cv / sl; H += tau * np.outer(cv, cv) / sl**2
            d = np.linalg.solve(H, g); dec = g @ d
            if dec < 1e-14: break
            t = 1.0
            while not feas(th - t * d): t *= 0.5
            f0 = phi(th, tau)
            while phi(th - t * d, tau) > f0 - 0.25 * t * dec and t > 1e-12: t *= 0.5
            th = th - t * d
        if tau < 1e-11: break
        tau *= 0.1
    kappa = None
    if eps is not None:
        kappa = tau / (1 - eps - nvec(th, betas))
    return th, kappa, tau

def quad_surrogate(P):
    """Minimise E(lam*-theta'a)^2 : G theta = (1/T) sum a_i  (E[lam* a]=E[a dN]/dt)."""
    K, bt = P.K, P.betas
    G = np.zeros((K + 1, K + 1)); w = P.w; s = P.s
    G[0, 0] = w.sum()
    for k in range(K):
        G[0, 1 + k] = G[1 + k, 0] = (s[:, k] * (1 - np.exp(-bt[k] * w)) / bt[k]).sum()
        for j in range(K):
            G[1 + j, 1 + k] = (s[:, j] * s[:, k] * (1 - np.exp(-(bt[j] + bt[k]) * w)) / (bt[j] + bt[k])).sum()
    G /= P.T; rhs = P.a.sum(0) / P.T
    return np.linalg.solve(G, rhs)
