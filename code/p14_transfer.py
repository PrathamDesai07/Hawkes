"""Transfer-function estimate of the tracking filter: demodulate mu_hat(t) at the drift frequency.
H_hat = mean over replicates of  (2/N) sum_t (mu_hat_t - mean) e^{-i w t} / (A mu_s),  compared with the first-order filter H = eta/(eta + i omega_event).
Standard error from the spread over replicates.  Estimators: 'known' (alpha known, forgetting-Newton, p14_v2), 'adaptive', 'fixed' (alpha unknown, p14_mech)."""
import numpy as np, math, sys
import p14_mech as M
from p14_v2 import run as run_known

def one(est, eps, eta, T, seed, A, P, mu_s=1.0, beta=1.0):
    if est == 'known':
        out = run_known(dict(name='smooth', amp=A, P=P), 'newton', eta, T, seed, eps, mu0=mu_s, beta=beta)
        return out
    return M.run(est, mu_s, beta, eps, eta, T, seed, A, P)

def transfer(est, eps, eta, R, A=0.3, P=1200.0, T=6000.0, mu_s=1.0, seed0=4000):
    grid = 5.0; tg = np.arange(1, int(T/grid)+1)*grid; sel = tg > 0.2*T; t = tg[sel]; e = np.exp(-2j*np.pi*t/P)
    cs = []; offs = []
    for s in range(R):
        x = one(est, eps, eta, T, seed0+s, A, P)[:len(tg)][sel]
        cs.append(1j*2*np.mean((x - x.mean())*e)/(A*mu_s)); offs.append(x.mean() - mu_s)
    cs = np.array(cs); H = cs.mean(); se = cs.std(ddof=1)/math.sqrt(R); r = mu_s/eps; w = 2*math.pi/(P*r)
    Hf = eta/(eta + 1j*w)
    return H, se, Hf, np.mean(offs), np.std(offs, ddof=1)/math.sqrt(R)

if __name__ == "__main__":
    eps = float(sys.argv[1]); R = int(sys.argv[2]); etas = [float(x) for x in sys.argv[3].split(',')]; ests = sys.argv[4].split(',')
    print(f"eps={eps}, A=0.3, P=1200: complex response H (re, im) +- se (modulus se) vs first-order filter; mean offset of mu_hat")
    for eta in etas:
        for est in ests:
            H, se, Hf, off, offse = transfer(est, eps, eta, R)
            print(f" eta={eta} {est:8s}: H = {H.real:.3f} {H.imag:+.3f}i (|H|={abs(H):.3f}, arg={np.angle(H):+.3f})  se {se:.3f} | filter: {Hf.real:.3f} {Hf.imag:+.3f}i (|H|={abs(Hf):.3f}, arg={np.angle(Hf):+.3f}) | offset {off:+.3f} +- {offse:.3f}", flush=True)
