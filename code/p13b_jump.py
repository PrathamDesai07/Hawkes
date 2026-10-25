"""Sharp jump rate: alternating levels mu_s <-> (1+D) mu_s with regime length L_t.  Drift part  v^d:  per-jump cumulative
  S_U = sum over a regime of U^d_i = ||v^d_i||^2_{F_{i-1}}   vs the Minkowski bound  (1/eta) * delta^2,  delta^2 = (1+lev/c-)^2 * sum_i ||d_i||^2_{F_i}  (one jump gap per regime)
  S_I = sum over a regime of ||v^d_i||^2_Iev  (Fisher-metric drift cost per jump)  -> should scale like 1/eta (old bound: 1/eta^2)."""
import numpy as np, math, sys
from p13b_engine import run
from p14_mech import Iev_matrix
eps = 0.4; mu_s = 1.0; D = float(sys.argv[1]); Lt = float(sys.argv[2]); R = int(sys.argv[3]); etas = [float(x) for x in sys.argv[4].split(',')]
Iev = Iev_matrix(mu_s, 1.0, eps); r = mu_s/eps; Le = Lt*r; nreg = 8; n = int(nreg*Le)
fn = lambda t: mu_s*(1 + D*(np.floor(np.asarray(t)/Lt).astype(int) % 2)); mu_max = mu_s*(1+D)
print(f"jump drift D={D} (|d|={D*mu_s}), regime {Lt} time units = {Le:.0f} events, eps={eps}")
res = []
for eta in etas:
    out = [run([1.0], eps, fn, mu_max, eta, n, 500+s, Iev, mu_s, burn=0.2) for s in range(R)]
    cat = lambda k: np.concatenate([o[k] for o in out]); N = len(cat('Ud')); per_reg = Le   # events per regime (post burn average)
    SU = cat('Ud').mean()*per_reg; SI = cat('vdI').mean()*per_reg; sd2 = cat('dF2').mean()*per_reg; lev = cat('lev').max()
    bound = (1+lev/0.6)**2*sd2/eta
    res.append((eta, SI)); print(f" eta={eta}: per-jump sum U^d = {SU:.3f}  <=  bound (1/eta)(1+lev/c-)^2 sum||d||^2_F = {bound:.3f} (ratio {SU/bound:.3f});   per-jump Fisher-metric drift cost S_I = {SI:.5f}", flush=True)
res = np.array(res)
if len(res) > 1:
    sl = np.polyfit(np.log(res[:, 0]), np.log(res[:, 1]), 1)[0]; print(f" log-log slope of S_I vs eta over the range: {sl:.2f}  (sharp rate: -1; the old Young-step bound has -2)")
