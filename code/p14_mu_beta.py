"""Problem 14: dependence of the noise constant a and the linear-window onset on k = 2 mu/beta (beta=1; only k matters by time rescaling).
(a) a(mu) = per-replicate MSE/eta (no drift, small eta) vs prediction mu/(2 eps m) with exact m.
(b) smooth drift with fixed per-event drift frequency omega=0.001, A=0.3: ratio bias^2/filter vs eta; onset where ratio > 1.5."""
import numpy as np, math, sys
from p14_v2 import stats
from p10_fisher_k1_exact import exact_m
eps = 0.2; mu = float(sys.argv[1]); R = int(sys.argv[2]); part = sys.argv[3]
r = mu/eps; m = exact_m(mu, 1.0, eps); k = 2*mu
if part == 'a':
    T = 30000.0/r
    for eta in (0.003, 0.006):
        mse, _ = stats(dict(name='none'), 'newton', eta, T, R, eps, seed0=11, mu0=mu)
        print(f"mu={mu} (k={k:.2f}) eta={eta}: a_meas = MSE/eta = {mse/eta:.2f}   predicted mu/(2 eps m) = {r/(2*m):.2f}   a*eps^2 = {mse/eta*eps**2:.3f} (pred {r/(2*m)*eps**2:.3f})", flush=True)
else:
    om = 0.001; P = 2*math.pi/(om*r); A = 0.3; T = 4*P; a_pred = r/(2*m)
    print(f"mu={mu} (k={k:.2f}): P={P:.0f} time units, T={T:.0f}, noise a~{a_pred:.1f}")
    for eta in (0.001, 0.002, 0.004, 0.008, 0.016):
        pred = (A*A*mu**2/2)*om**2/(eta**2 + om**2)
        mse, b2 = stats(dict(name='smooth', amp=A, P=P), 'newton', eta, T, R, eps, seed0=21, mu0=mu)
        print(f"   eta={eta}: bias^2 {b2:.5f} vs filter {pred:.5f}  ratio {b2/pred:.2f}   a*eta/mu^2 = {a_pred*eta/mu**2:.4f}  (eta/eps^2 = {eta/eps**2:.3f})", flush=True)
