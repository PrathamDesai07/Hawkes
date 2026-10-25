"""Problem 6: summarise code/p6_out/*.txt (output of p6_scaling.exe) against the predictions P1-P7 of problem06_epsilon_explicit.md."""
import glob, math, os, re
import numpy as np

def load(fn):
    hdr = open(fn).readline()
    a = np.loadtxt(fn, comments='#')
    return hdr, a   # cols: T ER se R-(d/2)lnT lndetF/2 lmin lmax N R-(d/2)lnN pred

def slope(T, R, lo, hi):
    m = (T >= lo) & (T <= hi)
    return np.polyfit(np.log(T[m]), R[m], 1)[0] if m.sum() >= 2 else float('nan')

def sat_time(T, R, d, frac=0.8, win=3):
    """first T where the local slope of R against ln T (over `win` grid steps) reaches frac*d/2 and stays there"""
    ok = []
    for i in range(len(T) - win):
        ok.append((R[i + win] - R[i]) / math.log(T[i + win] / T[i]) >= frac * d / 2)
    for i in range(len(ok)):
        if all(ok[i:]): return T[i]
    return float('nan')

rows = []
for fn in sorted(glob.glob('p6_out/*.txt')):
    if os.path.basename(fn) == 'log.txt': continue
    hdr, a = load(fn)
    K = int(re.search(r'K=(\d+)', hdr).group(1)); d = K + 1
    eps = float(re.search(r'eps=([\d.]+)', hdr).group(1)); mu = float(re.search(r' mu=([\d.]+)', hdr).group(1))
    beta0 = float(re.search(r'beta=([\d.]+)', hdr).group(1))
    T, R, se, c_t, hd, lmin, lmax, N, c_n, pred = a.T
    tail = T >= T[-1] / 10
    rows.append(dict(f=os.path.basename(fn)[:-4], K=K, eps=eps, mu=mu, d=d, k0=2 * mu / beta0,
                     slope_last=slope(T, R, T[-1] / 30, T[-1]), c_time=c_t[-1], c_time_se=se[-1], c_ev=c_n[-1],
                     Tsat=sat_time(T, R, d), lmin=lmin[tail].mean(), lmax=lmax[tail].mean(),
                     detI=(lmin * lmax)[tail].mean() if K == 1 else float('nan'),
                     ratio_pred=(R[-4:] / pred[-4:]).mean(), R_t_star=np.interp((2 * (2 * mu / beta0 - 1)) / (mu * eps) if K == 1 else 0, T, R)))
print(f"{'run':34s} {'slope':>6s} {'c_time':>8s} {'c_events':>8s} {'T_sat':>8s} {'lmin':>9s} {'pred lmin':>9s} {'lmax*eps':>8s} {'det I':>7s} {'R/pred':>6s}")
for r in rows:
    pl = r['eps'] / (r['mu'] * (r['k0'] - 1)) if r['K'] == 1 and r['k0'] > 1 else float('nan')
    print(f"{r['f']:34s} {r['slope_last']:6.2f} {r['c_time']:8.2f} {r['c_ev']:8.2f} {r['Tsat']:8.1f} {r['lmin']:9.4f} {pl:9.4f} {r['lmax']*r['eps']:8.3f} {r['detI']:7.3f} {r['ratio_pred']:6.2f}")
