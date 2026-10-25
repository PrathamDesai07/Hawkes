"""Problem 15: optimal eta when alpha is also estimated (K=1, adaptive recursion, p14_gen), smooth drift.
Prediction MSE(eta) = 1.15*r/(2 eff) * eta + bias2_tv(eta) ;  eff = (m - eps/mu)/(1-2 eps + eps mu m)  (information about mu with alpha unknown)."""
import numpy as np, math, sys
from p14_gen import run
from p14_tv_filter import tv_bias2
from p10_fisher_k1_exact import exact_m
eps = float(sys.argv[1]); R = int(sys.argv[2]); A = float(sys.argv[4]) if len(sys.argv) > 4 else 0.2; P = 1200.0; T = 6000.0; mu0 = 1.0; etas = [float(x) for x in sys.argv[3].split(',')]
r = mu0/eps; m = exact_m(mu0, 1.0, eps); eff = (m - eps/mu0)/(1 - 2*eps + eps*mu0*m); a_u = 1.15*r/(2*eff); a_k = 1.14*r/(2*m)
fn = lambda t: mu0*(1 + A*math.sin(2*math.pi*t/P)); fn_arr = lambda t: mu0*(1 + A*np.sin(2*np.pi*t/P)); tg = np.arange(1, int(T/5.0)+1)*5.0; sel = tg > 0.2*T
rows = []
for eta in etas:
    mses = [np.mean((run('adaptive', [1.0], eps, fn, mu0*(1+A), eta, T, 90+k, mu0)[1][sel] - fn_arr(tg[sel]))**2) for k in range(R)]
    pred = a_u*eta + tv_bias2('sin', A, eta, eps, P=P, mu0=mu0, cycles=6, n_per=12000); rows.append((eta, np.mean(mses), pred))
    print(f"  eps={eps} eta={eta}: MSE_sim {np.mean(mses):.5f} +- {np.std(mses,ddof=1)/math.sqrt(R):.5f}   MSE_pred {pred:.5f}  (noise {a_u*eta:.5f})", flush=True)
rows = np.array(rows); le = np.log(rows[:, 0]); i = int(np.argmin(rows[:, 1]))
c = np.polyfit(le[max(i-1,0):i+2], rows[max(i-1,0):i+2, 1], 2) if 0 < i < len(etas)-1 else None
eta_sim = math.exp(-c[1]/(2*c[0])) if c is not None and c[0] > 0 else rows[i, 0]
fine = np.logspace(np.log10(etas[0]/2), np.log10(etas[-1]*2), 50)
pf_u = [a_u*e + tv_bias2('sin', A, e, eps, P=P, mu0=mu0, cycles=6, n_per=12000) for e in fine]; pf_k = [a_k*e + tv_bias2('sin', A, e, eps, P=P, mu0=mu0, cycles=6, n_per=12000) for e in fine]
print(f"eps={eps}: noise coefficients a(alpha unknown)={a_u:.2f}, a(alpha known)={a_k:.2f}; eta_opt simulated (alpha unknown) {eta_sim:.5f}; predicted argmin alpha unknown {fine[int(np.argmin(pf_u))]:.5f}, alpha known {fine[int(np.argmin(pf_k))]:.5f}")
