"""Problem 15, jump drift: alternating levels mu0 <-> (1+D)mu0 with fixed regime length L_t (time units), alpha known, forgetting-Newton.
Prediction MSE_pred(eta) = c_a a0 eta + bias2_tv(eta) (deterministic time-varying-gain filter on the jump path), a0 = r/(2m), c_a = 1.14.
Asymptotic law (relaxation fast vs regime):  eta*_asym = D mu0 / sqrt(2 a L_e),  L_e = L_t r events per regime."""
import numpy as np, math, sys
from p14_v2 import stats
from p10_fisher_k1_exact import exact_m
eps = float(sys.argv[1]); Lt = float(sys.argv[2]); R = int(sys.argv[3]); etas = [float(x) for x in sys.argv[4].split(',')]; D = float(sys.argv[5]) if len(sys.argv) > 5 else 0.3
mu0 = 1.0; r = mu0/eps; m = exact_m(mu0, 1.0, eps); a0 = 1.14*r/(2*m); Le = Lt*r; T = 8*Lt
times = Lt*np.arange(1, 12)
def pred_bias2(eta):
    tt = np.linspace(0, T, 300001); dt = tt[1]-tt[0]; mu = mu0*(1 + D*(np.floor(tt/Lt).astype(int) % 2)); y = np.empty_like(mu); y[0] = mu0
    for i in range(len(tt)-1): y[i+1] = y[i] + dt*eta*(mu[i]/eps)*(mu[i]-y[i])
    s = tt > 2*Lt; return np.mean((y[s]-mu[s])**2)
rows = []
for eta in etas:
    mse, b2 = stats(dict(name='jump', amp=D, times=times), 'newton', eta, T, R, eps, seed0=71, mu0=mu0, burn=0.25)
    pb = pred_bias2(eta); rows.append((eta, mse, a0*eta + pb)); print(f"  eps={eps} L_t={Lt:.0f} (L_e={Le:.0f}) eta={eta}: MSE_sim {mse:.5f}  MSE_pred {a0*eta+pb:.5f} (noise {a0*eta:.5f}, bias {pb:.5f}); eta*L_e = {eta*Le:.1f}", flush=True)
rows = np.array(rows); le = np.log(rows[:, 0])
def am(y):
    i = int(np.argmin(y))
    if 0 < i < len(y)-1:
        c = np.polyfit(le[i-1:i+2], y[i-1:i+2], 2)
        return math.exp(-c[1]/(2*c[0])) if c[0] > 0 else rows[i, 0]
    return rows[i, 0]
fine = np.logspace(np.log10(etas[0]/2), np.log10(etas[-1]*2), 40); pf = [a0*e + pred_bias2(e) for e in fine]
print(f"eps={eps}, L_t={Lt}: eta_opt simulated {am(rows[:,1]):.5f}; predicted argmin {fine[int(np.argmin(pf))]:.5f}; asymptotic D mu/sqrt(2 a L_e) = {D*mu0/math.sqrt(2*a0*Le):.5f}; a0={a0:.2f}; min MSE sim {rows[:,1].min():.5f}, pred {min(pf):.5f}")
