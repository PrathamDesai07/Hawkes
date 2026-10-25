"""Ramp lag with the offset removed: lag = (mean error of up-ramp - mean error of down-ramp)/2 in a window around the midpoint (same mu range, same offset)."""
import numpy as np, math, sys
from p14_gen import run
eps = 0.4; mu0 = 1.0; beta = [1.0]; r = mu0/eps; R = int(sys.argv[1]); T = 3000.0; s = 0.8/T
tg = np.arange(1, int(T/5)+1)*5.0; win = (tg > T/2-400) & (tg < T/2+400)
up = lambda t: 0.8 + s*t; dn = lambda t: 1.6 - s*t
for eta in (0.002, 0.004):
    for est in sys.argv[2].split(','):
        eu = np.mean([(run(est, beta, eps, up, 1.6, eta, T, 100+k, 1.2)[1] - up(tg))[win].mean() for k in range(R)])
        ed = np.mean([(run(est, beta, eps, dn, 1.6, eta, T, 100+k, 1.2)[1] - dn(tg))[win].mean() for k in range(R)])
        print(f"eta={eta} {est:8s}: up-ramp err {eu:+.4f}, down-ramp err {ed:+.4f} -> lag {(eu-ed)/2:+.4f}, offset {(eu+ed)/2:+.4f};  predicted lag -(s/r)/eta = {-(s/r)/eta:+.4f}", flush=True)
