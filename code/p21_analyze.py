"""Problem 21: summarise code/p21_out (raw rows from p21_baselines.exe) against B1-B5 of problem21_baselines.md; figures -> figures/p21_*.png; numbers -> code/p21_summary.txt"""
import glob, os, re, sys, types
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.modules.setdefault('mpmath', types.ModuleType('mpmath')); sys.modules['mpmath'].mp = types.SimpleNamespace(dps=60)
from p2_tangent_bound import matrices

HERE = os.path.dirname(os.path.abspath(__file__)); FIG = os.path.join(HERE, '..', 'figures'); os.makedirs(FIG, exist_ok=True)
OUT = open(os.path.join(HERE, 'p21_summary.txt'), 'w')
def say(s=''):
    print(s); OUT.write(s + '\n')

def load(fn):
    head = open(fn).readline(); g = lambda k: re.search(k + r'=([\d.,]+)', head).group(1)
    K = int(g('K')); d = K + 1; eps = float(g('eps')); mu = float(g(' mu')); beta = np.array([float(x) for x in g('beta').split(',')])
    data = {}; Jacc = {}
    for line in open(fn):
        if line.startswith('# '): continue
        t = line.split()
        if t[0] == '#I':
            rep, n = int(t[1]), int(t[2]); Jacc.setdefault(n, []).append([float(x) for x in t[3:]])
        else:
            rep, name, n = int(t[0]), t[1], int(t[2]); T = float(t[3]); R = float(t[4]); e = float(t[5]); v = [float(x) for x in t[6:6 + d]]
            data.setdefault(name, {}).setdefault(n, []).append((R, e, *v))
    arr = {m: {n: np.array(rows) for n, rows in dd.items()} for m, dd in data.items()}
    ns = sorted(next(iter(arr.values())).keys())
    return dict(K=K, d=d, eps=eps, mu=mu, beta=beta, arr=arr, ns=ns, J={n: np.array(v) for n, v in Jacc.items()})

def mean_ci(x): return x.mean(), 1.96 * x.std(ddof=1) / np.sqrt(len(x))
def Imat(r, n):
    d = r['d']; J = r['J'][n]; T = J[:, 0]; Jm = J[:, 1:].mean(0); Tm = T.mean(); M = np.zeros((d, d)); c = 0
    for p in range(d):
        for q in range(p + 1): M[p, q] = M[q, p] = Jm[c]; c += 1
    return M / Tm, Tm

runs = {os.path.basename(f)[4:-4]: load(f) for f in sorted(glob.glob(os.path.join(HERE, 'p21_out', 'acc_*.csv')))}
families = {'ours (online Newton)': ['newton'], 'batch MLE (refit at n)': ['batch'], 'Solo-type EWLS, rho=1': ['ls_r1'], 'Solo-type EWLS, rho=1-1/2000': ['ls_r2000'], 'Solo-type EWLS, rho=1-1/200': ['ls_r200']}
def hw_names(r, kind): return sorted([m for m in r['arr'] if m.startswith('hw_' + kind)], key=lambda s: float(s[4:]))
def best_hw(r, kind, n):
    best = None
    for m in hw_names(r, kind):
        e = r['arr'][m][n][:, 1].mean()
        if best is None or e < best[1]: best = (m, e)
    return best[0]

say("== Problem 21 summary: mean Fisher-normalised error e_n = (theta-theta*)' J_n (theta-theta*) (target d for an efficient estimator), identical paths, 95 % CI over replicates")
for key, r in runs.items():
    d = r['d']; ns = r['ns']; nlast = ns[-1]; nsel = [n for n in (100, 316, 1000, 3162, 10000, 100000) if n in ns or any(abs(n - m) / n < 0.03 for m in ns)]
    nsel = [min(ns, key=lambda m: abs(m - n)) for n in (100, 1000, 10000, 100000)]
    say(f"\n{key} (K={r['K']}, eps={r['eps']}, mu={r['mu']}, reps={len(r['arr']['newton'][nlast])})   checkpoints n = {nsel}")
    rows = []
    for name, ms in families.items(): rows.append((name, ms[0]))
    for kind, lab in (('c', 'HW-type, constant step'), ('s', 'HW-type, eta0/sqrt(i)'), ('i', 'HW-type, eta0/i')):
        rows.append((lab + ' [tuned on the truth: ' + best_hw(r, kind, nlast) + ']', best_hw(r, kind, nlast)))
        rows.append((lab + ' [default eta0 = 0.1]', 'hw_%s0.1' % kind))
    say(f"  {'method':62s} " + " ".join(f"{'e @ n='+str(n):>16s}" for n in nsel) + f" {'E R @ n=max':>12s} {'1e3 R/n':>9s}")
    for lab, m in rows:
        cells = []
        for n in nsel:
            mm, ci = mean_ci(r['arr'][m][n][:, 1]); cells.append(f"{mm:8.2f}+-{ci:5.2f}")
        R = r['arr'][m][nlast][:, 0]; Rs = (f"{np.nanmean(R):12.2f} {1e3*np.nanmean(R)/nlast:9.3f}") if not np.isnan(R).all() else f"{'-':>12s} {'-':>9s}"
        say(f"  {lab:62s} " + " ".join(f"{c:>16s}" for c in cells) + Rs)
    # LS efficiency prediction tr(I B^-1 C B^-1)
    I, Tm = Imat(r, nlast); beta = r['beta']; eps = r['eps']; mu = r['mu']; alpha = (1 - eps) * beta / r['K']
    M, M3, El = matrices(mu, alpha, beta); B = M; C = M3; Bi = np.linalg.inv(B); eff_ls = np.trace(I @ Bi @ C @ Bi)
    say(f"  B2 least-squares efficiency: predicted limit E e_n = tr(I B^-1 C B^-1) = {eff_ls:.2f}  (MLE: d = {d});  measured Solo-type rho=1 at n={nlast}: {r['arr']['ls_r1'][nlast][:,1].mean():.2f}")
    # regret slopes: slope of E R_n against ln n on the last decade
    for m in ['newton', 'ls_r1', 'ls_r2000', 'ls_r200'] + [best_hw(r, 's', nlast), best_hw(r, 'c', nlast)]:
        n1, n0 = nlast, min(ns, key=lambda x: abs(x - nlast / 10)); R1 = r['arr'][m][n1][:, 0].mean(); R0 = r['arr'][m][n0][:, 0].mean()
        say(f"     regret slope vs ln n, last decade, {m:12s}: {(R1 - R0) / np.log(n1 / n0):6.2f}   (ours target d/2 = {d/2})")

# ---- figures
cols = {'newton': 'k', 'batch': 'tab:blue', 'ls_r1': 'tab:green', 'ls_r2000': 'tab:olive', 'ls_r200': 'tab:cyan'}
for key, r in runs.items():
    if r['K'] != 1 and key != 'k2_eps0.1': continue
    d = r['d']; ns = np.array(r['ns']); nlast = ns[-1]
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.4))
    for name, m in (('ours (online Newton)', 'newton'), ('batch MLE', 'batch'), ('EWLS rho=1', 'ls_r1'), ('EWLS rho=1-1/2000', 'ls_r2000'), ('EWLS rho=1-1/200', 'ls_r200'), ('HW-type const step (tuned)', best_hw(r, 'c', nlast)), ('HW-type eta0/sqrt(i) (tuned)', best_hw(r, 's', nlast)), ('HW-type eta0/i (tuned)', best_hw(r, 'i', nlast))):
        e = np.array([r['arr'][m][n][:, 1].mean() for n in ns]); ls_ = '-' if not m.startswith('hw') else '--'
        ax[0].loglog(ns, e, ls_, label=name, color=cols.get(m, None))
    ax[0].axhline(d, color='gray', ls=':'); ax[0].set_xlabel('events n'); ax[0].set_ylabel("E e_n  (Fisher-normalised squared error)"); ax[0].set_title(f"K={r['K']}, eps={r['eps']}: accuracy on identical paths (dotted: d)", fontsize=9); ax[0].legend(fontsize=7)
    for name, m in (('ours', 'newton'), ('EWLS rho=1', 'ls_r1'), ('EWLS rho=1-1/2000', 'ls_r2000'), ('EWLS rho=1-1/200', 'ls_r200'), ('HW-type eta0/sqrt(i) (tuned)', best_hw(r, 's', nlast)), ('HW-type const step (tuned)', best_hw(r, 'c', nlast))):
        R = np.array([r['arr'][m][n][:, 0].mean() for n in ns]); ax[1].semilogx(ns, R, '-' if not m.startswith('hw') else '--', label=name, color=cols.get(m, None))
    ax[1].set_xlabel('events n'); ax[1].set_ylabel('E R_n (regret)'); ax[1].set_title('sequential regret (batch MLE is not a sequential predictor)', fontsize=9); ax[1].legend(fontsize=7)
    plt.tight_layout(); plt.savefig(os.path.join(FIG, f'p21_{key}.png'), dpi=140); plt.close()

# ---- timing table
say("\n== Timing (microseconds per event unless stated; single thread, minimum of 3 passes, data generation excluded)")
for fn in sorted(glob.glob(os.path.join(HERE, 'p21_out', 'time_*.txt'))):
    say(os.path.basename(fn)); [say('  ' + l.rstrip()) for l in open(fn)]
OUT.close()
