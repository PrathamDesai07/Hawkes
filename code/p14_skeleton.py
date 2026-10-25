"""Problem 14, amplitude effect: DETERMINISTIC skeleton of the mean error dynamics (no noise).
 Truth: mu_t = 1 + A sin(2 pi t/P), alpha fixed (eps fixed), beta = 1.  Per-event information at the truth: f_mat(mu) = I(mu)/r(mu), r = mu/eps, I from the exact Laplace-ODE m(mu).
 Mean error e (2-vector, alpha unknown) obeys   e_{i+1} = (Id - F_i^{-1} f_mat(mu_t)) e_i - (d_i, 0),
   'fixed'    : F_i = f_mat(1)/eta                        (constant preconditioner at the reference parameter)
   'adaptive' : F_i = (1-eta) F_{i-1} + f_mat(mu_t)       (forgetting window)
   'known'    : scalar, alpha known, I_mumu only.
 Output: fundamental error^2 of e_mu(t), as a ratio to the constant-gain filter value (A^2/2) w^2/(eta^2+w^2), w = 2 pi/(P r0)."""
import numpy as np, math, sys
from p10_fisher_k1 import exact_m

def build(eps, beta=1.0, grid=np.linspace(0.3, 1.9, 33)):
    ms = np.array([exact_m(mu, (1 - eps)*beta, beta) for mu in grid]); return grid, ms

def make_f(eps, beta=1.0):
    grid, ms = build(eps, beta); alpha = (1 - eps)*beta
    def fmat(mu):
        m = np.interp(mu, grid, ms); El = mu/eps
        I = np.array([[m, (1 - mu*m)/alpha], [(1 - mu*m)/alpha, (El - 2*mu + mu*mu*m)/alpha**2]]); return I/(mu/eps)
    return fmat

def run(kind, A, eps, eta, P=1200.0, ncyc=14):
    fmat = make_f(eps); mu0 = 1.0; r0 = mu0/eps
    f_ref = fmat(mu0); F = f_ref/eta; e = np.zeros(2); t = 0.0; ts = []; es = []
    T = ncyc*P
    while t < T:
        mu = 1 + A*math.sin(2*math.pi*t/P); r = mu/eps; fm = fmat(mu); tn = t + 1.0/r; d = A*(math.sin(2*math.pi*tn/P) - math.sin(2*math.pi*t/P))
        if kind == 'fixed': G = np.linalg.solve(F, fm)
        elif kind == 'adaptive': F = (1 - eta)*F + fm; G = np.linalg.solve(F, fm)
        elif kind == 'known': F = f_ref[0, 0]/eta; G = np.array([[fm[0, 0]/F, 0], [0, 0]])
        elif kind == 'known_adaptive':
            F = (1 - eta)*F + fm if np.ndim(F) == 0 else (1 - eta)*F[0, 0] + fm[0, 0]; G = np.array([[fm[0, 0]/F, 0], [0, 0]]); F = F
        e = (np.eye(2) - G)@e - np.array([d, 0.0]); t = tn; ts.append(t); es.append(e[0])
    ts = np.array(ts); es = np.array(es); sel = ts > 0.5*T; t = ts[sel]; y = es[sel]
    X = np.column_stack([np.ones_like(t)] + [f(2*math.pi*k*t/P) for k in (1, 2, 3) for f in (np.sin, np.cos)])
    c, *_ = np.linalg.lstsq(X, y, rcond=None); fund = (c[1]**2 + c[2]**2)/2; w = 2*math.pi/(P*r0); pred = (A*A/2)*w*w/(eta*eta + w*w)
    return fund/pred, c[0], math.hypot(c[3], c[4])

if __name__ == "__main__":
    eps = float(sys.argv[1]) if len(sys.argv) > 1 else 0.4; eta = float(sys.argv[2]) if len(sys.argv) > 2 else 0.004
    print(f"deterministic skeleton, eps={eps}, eta={eta}: ratio of fundamental error^2 to the constant-gain filter (offset, 2nd harmonic amplitude in brackets)")
    print("  A      | fixed               | adaptive (alpha unknown) | adaptive/fixed | known (alpha known, fixed)")
    for A in (0.05, 0.1, 0.2, 0.35, 0.5):
        rf = run('fixed', A, eps, eta); ra = run('adaptive', A, eps, eta); rk = run('known', A, eps, eta)
        print(f"  {A:5.2f}  | {rf[0]:6.3f} ({rf[1]:+.3f},{rf[2]:.3f}) | {ra[0]:6.3f} ({ra[1]:+.3f},{ra[2]:.3f})     | {ra[0]/rf[0]:6.3f}         | {rk[0]:6.3f}", flush=True)
    print("  measured (Monte Carlo, eps=0.4, eta=0.004): A=0.2: fixed 1.07, adaptive 1.41, ratio 1.32; A=0.35: 0.98, 1.48, 1.48; A=0.5: 0.92, 1.61, 1.75;  known: 1.07, 1.15, 1.24")
