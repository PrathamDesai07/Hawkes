"""Problem 14: analytic (second-order) explanation of the amplitude effect, scalar case (alpha known), checked against the deterministic skeleton.
 Mean error:  e' = -eta Phi(s) e - d(s)  (s = event index), Phi = f(mu_t)/f_ref  ('fixed')  or  f(mu_t)/(eta Ftilde)  ('adaptive', Ftilde = exponentially lagged f).
 Linear response:  Phi = 1 + A (g_+ e^{i w s} + conj(g_+) e^{-i w s}),  g_+ = (kappa/2i) * c(w),  c = 1 (fixed),  c = 1 - H(w) = i w/(eta + i w) (adaptive),  kappa = dln f/dln mu.
 Harmonic balance:  E_1 = e1 [1 + (eta A^2/(eta + i w)) ( g_+^2 conj(e1)/e1 + |g_+|^2 (1 + eta/(eta + 2 i w)) )] + O(A^4),  e1 = -D_1/(eta + i w) (constant-gain filter).
 Predicted ratio of fundamental error^2 to the constant-gain filter:  |E_1/e1|^2.  Forcing modulation (time dilation) is ignored here: it enters at the same order and explains the remaining gap."""
import numpy as np, math, sys
from p10_fisher_k1 import exact_m

eps = float(sys.argv[1]) if len(sys.argv) > 1 else 0.4; eta = float(sys.argv[2]) if len(sys.argv) > 2 else 0.004; P = 1200.0
alpha = 1 - eps
def f(mu): return exact_m(mu, alpha, 1.0)/(mu/eps)                    # information about mu per event (alpha known)
h = 1e-3; kappa = (math.log(f(1 + h)) - math.log(f(1 - h)))/(math.log(1 + h) - math.log(1 - h))
w = 2*math.pi/(P/eps)                                                  # per-event frequency at the reference rate r0 = 1/eps
print(f"eps={eps}, eta={eta}: kappa = dln f/dln mu = {kappa:.3f} (alpha known), per-event frequency w = {w:.5f}, eta/w = {eta/w:.2f}")

def predicted(A, kind):
    """continuous-time harmonic balance: de/dt = -Gamma(t) e - dmu/dt, forcing exactly first-harmonic; Gamma = eta r_t Phi_t,
    Gamma/(eta r0) = I(mu_t)/I_ref ('fixed')  or  I(mu_t)/I_lagged ('adaptive'),  I = information per unit time (kappa_I = kappa + 1)."""
    eta_t = eta/eps; wt = 2*math.pi/P; kI = kappa + 1.0
    c = 1.0 if kind == 'fixed' else (1j*wt/(eta_t + 1j*wt))
    gp = kI/(2j)*c; e1 = -1.0/(eta_t + 1j*wt)        # complex amplitude of the constant-gain response (forcing amplitude A w/2 real): fixes the phase in conj(e1)/e1
    corr = (eta_t*A*A/(eta_t + 1j*wt))*(gp**2*np.conj(e1)/e1 + abs(gp)**2*(1 + eta_t/(eta_t + 2j*wt)))
    return abs(1 + corr)**2

# deterministic scalar skeleton
def skeleton(kind, A, ncyc=14):
    mu0 = 1.0; r0 = mu0/eps; fref = f(mu0); F = fref/eta; e = 0.0; t = 0.0; T = ncyc*P; ts = []; es = []
    grid = np.linspace(0.3, 1.9, 33); fv = np.array([f(m) for m in grid])
    while t < T:
        mu = 1 + A*math.sin(2*math.pi*t/P); r = mu/eps; fm = np.interp(mu, grid, fv); tn = t + 1.0/r; d = A*(math.sin(2*math.pi*tn/P) - math.sin(2*math.pi*t/P))
        if kind == 'adaptive': F = (1 - eta)*F + fm
        e = (1 - fm/F)*e - d; t = tn; ts.append(t); es.append(e)
    ts = np.array(ts); es = np.array(es); sel = ts > 0.5*T; tt = ts[sel]; y = es[sel]
    X = np.column_stack([np.ones_like(tt)] + [fn(2*math.pi*k*tt/P) for k in (1, 2, 3) for fn in (np.sin, np.cos)])
    cc, *_ = np.linalg.lstsq(X, y, rcond=None); fund = (cc[1]**2 + cc[2]**2)/2; ww = 2*math.pi/(P*r0); pred = (A*A/2)*ww*ww/(eta*eta + ww*ww)
    return fund/pred

print("   A    | fixed: skeleton  analytic | adaptive: skeleton  analytic")
for A in (0.1, 0.2, 0.35, 0.5):
    print(f"  {A:4.2f}  |   {skeleton('fixed', A):7.3f}   {predicted(A, 'fixed'):7.3f}  |   {skeleton('adaptive', A):7.3f}    {predicted(A, 'adaptive'):7.3f}", flush=True)
