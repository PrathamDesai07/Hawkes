"""Problem 15: simulated MSE(eta) vs prediction  MSE_pred(eta) = c_a * a0 * eta + bias2_tv(eta),  a0 = r/(2m) (alpha known), c_a=1.14 (measured),
bias2_tv = deterministic time-varying-gain filter response (p14_tv_filter).  alpha known, forgetting-Newton (p14_v2).
Prints the simulated and predicted argmin and the asymptotic law eta*_asym = (A^2 mu^2 w_e^2 / a)^{1/3}  (eta >> w_e limit of  aη + (A^2 mu^2/2) w_e^2/η^2),  w_e = 2pi/(P r)."""
import numpy as np, math, sys
from p14_v2 import stats
from p14_tv_filter import tv_bias2
from p10_fisher_k1_exact import exact_m
eps = float(sys.argv[1]); R = int(sys.argv[2]); A = float(sys.argv[3]); P = float(sys.argv[4]); T = float(sys.argv[5]); etas = [float(x) for x in sys.argv[6].split(',')]
mu0 = 1.0; r = mu0/eps; m = exact_m(mu0, 1.0, eps); a0 = 1.14*r/(2*m); we = 2*math.pi/(P*r)
rows = []
for eta in etas:
    mse, b2 = stats(dict(name='smooth', amp=A, P=P), 'newton', eta, T, R, eps, seed0=51, mu0=mu0)
    pred = a0*eta + tv_bias2('sin', A, eta, eps, P=P, mu0=mu0, cycles=6, n_per=12000)
    rows.append((eta, mse, pred)); print(f"  eps={eps} A={A} eta={eta}: MSE_sim {mse:.5f}   MSE_pred {pred:.5f}   (noise {a0*eta:.5f}, tv-bias {pred-a0*eta:.5f})", flush=True)
rows = np.array(rows); le = np.log(rows[:, 0])
def argmin_log(y):
    i = int(np.argmin(y)); 
    if 0 < i < len(y)-1:
        c = np.polyfit(le[i-1:i+2], y[i-1:i+2], 2); return math.exp(-c[1]/(2*c[0])) if c[0] > 0 else rows[i, 0]
    return rows[i, 0]
fine = np.logspace(np.log10(etas[0]/2), np.log10(etas[-1]*2), 60)
pred_fine = [a0*e + tv_bias2('sin', A, e, eps, P=P, mu0=mu0, cycles=6, n_per=12000) for e in fine]
eta_pred = fine[int(np.argmin(pred_fine))]; eta_asym = (A*A*mu0*mu0*we*we/a0)**(1/3)
print(f"eps={eps}: eta_opt simulated (parabola in log eta) = {argmin_log(rows[:,1]):.5f};  predicted argmin = {eta_pred:.5f};  asymptotic (A^2 mu^2 w^2/a)^(1/3) = {eta_asym:.5f};  w_e={we:.5f}, a0={a0:.2f}")
print(f"        MSE at simulated optimum {rows[:,1].min():.5f}; MSE_pred min {min(pred_fine):.5f}")
