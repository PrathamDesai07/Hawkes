"""Jump class: alternating levels mu0 <-> 1.3 mu0 with fixed spacing L = Lev/r events->time; ensemble bias^2 vs deterministic time-varying-gain filter
dy/dt = eta r_t (mu_t - y), r_t = mu_t/eps  (and the constant-gain closed form D^2 tanh(kL/2)/(2kL)).  est in {known, adaptive}."""
import numpy as np, math, sys
from p14_gen import run
cfgs = {'K1_e04': ([1.0], 1.0, 0.4), 'K1_e02': ([1.0], 1.0, 0.2), 'K1_k06': ([1.0], 0.3, 0.4), 'K1_k4': ([1.0], 2.0, 0.4), 'K2a': ([1.0, 2.0], 1.5, 0.4)}
name = sys.argv[1]; R = int(sys.argv[2]); est = sys.argv[3]; eta = float(sys.argv[4]); Lev = float(sys.argv[5]) if len(sys.argv) > 5 else 1250.0
beta, mu0, eps = cfgs[name]; beta = np.array(beta); r0 = mu0/eps; L = Lev/r0; T = 8*L; D = 0.3*mu0
level = lambda t: np.floor(np.asarray(t)/L).astype(int) % 2
fn_arr = lambda t: mu0 + D*level(t); fn = lambda t: float(mu0 + D*(int(t//L) % 2)); mu_max = mu0 + D
tg = np.arange(1, int(T/5.0)+1)*5.0; sel = tg > 2*L
reps = np.array([run(est, beta, eps, fn, mu_max, eta, T, 800+k, mu0)[1][sel] for k in range(R)])
mean = reps.mean(0); varc = reps.var(0, ddof=1).mean()/R; tgt = fn_arr(tg[sel]); bias2 = np.mean((mean - tgt)**2) - varc
# deterministic prediction
tt = np.linspace(0, T, 400001); dtt = tt[1]-tt[0]; mu = fn_arr(tt); y = np.empty_like(mu); y[0] = mu0
for i in range(len(tt)-1): y[i+1] = y[i] + dtt*eta*(mu[i]/eps)*(mu[i]-y[i])
s2 = tt > 2*L; pred_tv = np.mean((y[s2]-mu[s2])**2)
kappa = eta*r0*(1+0.15); pred_cg = D**2*math.tanh(kappa*L/2)/(2*kappa*L)
print(f"{name} est={est} eta={eta} L={L:.0f} (eta*Lev={eta*Lev:.1f}): measured bias^2 {bias2:.5f} | time-varying-gain filter {pred_tv:.5f} (ratio {bias2/pred_tv:.2f}) | const-gain tanh formula {pred_cg:.5f} | noise var/R removed {varc:.5f}", flush=True)
