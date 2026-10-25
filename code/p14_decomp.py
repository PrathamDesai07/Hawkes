"""Decompose the ensemble-mean tracking error into: offset, fundamental-mismatch, 2nd harmonic, residual (low-frequency noise etc.).
Error series d(t)=mean_reps mu_hat(t) - mu_t on the window t>0.2T (integer number of periods). Contributions to mean(d^2)-noise_var/R."""
import numpy as np, math, sys
from p14_transfer import one
def run(est, eps, eta, R, A, P=1200.0, T=6000.0, mu_s=1.0, seed0=4000):
    grid = 5.0; tg = np.arange(1, int(T/grid)+1)*grid; sel = tg > 0.2*T; t = tg[sel]; tgt = mu_s*(1+A*np.sin(2*np.pi*t/P))
    reps = np.array([one(est, eps, eta, T, seed0+s, A, P)[:len(tg)][sel] for s in range(R)])
    d = reps.mean(0) - tgt; varc = reps.var(0, ddof=1).mean()/R
    X = np.column_stack([np.ones_like(t)] + [f(2*np.pi*k*t/P) for k in (1, 2) for f in (np.sin, np.cos)])
    c, *_ = np.linalg.lstsq(X, d, rcond=None); fit = X@c; parts = {}
    parts['offset^2'] = c[0]**2; parts['fund error^2'] = (c[1]**2 + c[2]**2)/2; parts['2nd harm^2'] = (c[3]**2 + c[4]**2)/2
    res = d - fit; parts['residual'] = np.mean(res**2) - varc
    total = np.mean(d**2) - varc
    return total, parts
if __name__ == "__main__":
    eps = float(sys.argv[1]); R = int(sys.argv[2]); eta = float(sys.argv[3]); A = float(sys.argv[4]); r = 1/eps; w = 2*math.pi/(1200.0*r)
    pred = (A*A/2)*w*w/(eta*eta + w*w); print(f"eps={eps} eta={eta} A={A}: filter bias^2 {pred:.5f}")
    for est in sys.argv[5].split(','):
        tot, p = run(est, eps, eta, R, A)
        print(f"  {est:8s}: total {tot:.5f} (ratio {tot/pred:.2f}) = " + ", ".join(f"{k} {v:.5f}" for k, v in p.items()), flush=True)
