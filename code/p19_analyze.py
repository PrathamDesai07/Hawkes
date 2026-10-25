"""Problem 19: fit slopes and draw the figures from code/p19_out/*.csv (raw output of p19_sim.exe).
Predictions tested (problem19_simulation_study.md):
  P1 regret:   E R_T = (d/2) ln T + c(eps) + o(1)      -> slope of E R against ln T tends to d/2 (late window), c flat in eps (time units)
  P2 raw error: E (theta_hat_j - theta_j)^2 = (I^{-1})_jj / T (1+o(1))  -> log-log slope -1 and MSE*T -> (I^{-1})_jj, I measured from J_T/T
  P3 normalised error: E D_n -> d and E e_n -> d (chi^2_d-like), flat in T after burn-in
  P4 burn-in:  T_exit(delta) ~ c_exit/(beta eps delta^2), c_exit ~ 2 (Problem 6); start cost of a far start = P_0/2 = (rho_0/2)[eps (Dmu)^2/mu + (mu/eps) sum (Dalpha_k/beta_k)^2]
Slope CIs: bootstrap over replicates (B = 400). Figures -> figures/p19_*.png ; numbers -> code/p19_summary.txt"""
import glob, os, re, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

rng = np.random.default_rng(19)
HERE = os.path.dirname(os.path.abspath(__file__)); FIG = os.path.join(HERE, '..', 'figures'); os.makedirs(FIG, exist_ok=True)
OUT = open(os.path.join(HERE, 'p19_summary.txt'), 'w')
def say(s=''):
    print(s); OUT.write(s + '\n')

def load(fn):
    head = open(fn).readline(); g = lambda k: re.search(k + r'=([\d.,]+)', head).group(1)
    K = int(g('K')); d = K + 1; eps = float(g('eps')); mu = float(g(' mu')); beta = [float(x) for x in g('beta').split(',')]
    rows = []; texit = []
    for line in open(fn):
        if line.startswith('#TEXIT'): texit.append(float(line.split()[2]))
        elif not line.startswith('#'): rows.append(line.split())
    a = np.array(rows, float); rep = a[:, 0].astype(int); T = a[:, 1]; ut = np.unique(T.round(6))
    reps = np.unique(rep); nr = len(reps); nT = len(a) // nr
    A = a.reshape(nr, nT, -1)           # replicate x checkpoint x columns  (rows are written replicate-major)
    return dict(K=K, d=d, eps=eps, mu=mu, beta=np.array(beta), start=re.search(r'start=(\w+)', head).group(1), A=A, T=A[:, :, 1].mean(0), texit=np.array(texit), nr=nr)

def slope_fit(T, y, lo, hi):
    m = (T >= lo) & (T <= hi)
    return np.polyfit(np.log(T[m]), y[m], 1)[0]
def boot_slope(T, Y, lo, hi, B=400):
    """Y: replicate x checkpoint; slope of the replicate-mean against ln T on [lo, hi]; bootstrap over replicates."""
    est = slope_fit(T, Y.mean(0), lo, hi); n = Y.shape[0]; bs = []
    for _ in range(B):
        idx = rng.integers(0, n, n); bs.append(slope_fit(T, Y[idx].mean(0), lo, hi))
    return est, np.percentile(bs, [2.5, 97.5])
def loglog_slope(T, Y, lo, hi, B=400):
    m = (T >= lo) & (T <= hi); lt = np.log(T[m]); n = Y.shape[0]
    f = lambda Z: np.polyfit(lt, np.log(Z.mean(0)[m]), 1)[0]
    est = f(Y); bs = [f(Y[rng.integers(0, n, n)]) for _ in range(B)]
    return est, np.percentile(bs, [2.5, 97.5])
def ci(Y):   # mean and 95% CI half-width over replicates
    return Y.mean(0), 1.96 * Y.std(0, ddof=1) / np.sqrt(Y.shape[0])

runs = {}
for fn in sorted(glob.glob(os.path.join(HERE, 'p19_out', '*.csv'))):
    r = load(fn); runs[os.path.basename(fn)[:-4]] = r

def Imat(r, idx_T):
    """information matrix estimate J_T/T averaged over replicates at checkpoint idx_T"""
    d = r['d']; J = r['A'][:, idx_T, 6 + d:].mean(0); I = np.zeros((d, d)); c = 0
    for p in range(d):
        for q in range(p + 1):
            I[p, q] = I[q, p] = J[c]; c += 1
    return I / r['T'][idx_T]

cols = {0.3: 'tab:blue', 0.1: 'tab:green', 0.03: 'tab:orange', 0.01: 'tab:red'}
say("== Problem 19 summary (all CIs: 95 %, bootstrap over replicates)")
res = {}
for key, r in runs.items():
    d = r['d']; T = r['T']; A = r['A']; eps = r['eps']; Tmax = T[-1]; R = A[:, :, 3]
    late = (Tmax / 10, Tmax); mid = (Tmax / 100, Tmax / 10); i_end = len(T) - 1
    s_late, ci_late = boot_slope(T, R, *late); s_mid, ci_mid = boot_slope(T, R, *mid)
    c = (R - d / 2 * np.log(T)[None, :]).mean(0)[T >= Tmax / 10].mean()
    I = Imat(r, i_end); Iinv = np.linalg.inv(I)
    mse = np.stack([A[:, :, 6 + j] ** 2 for j in range(d)])          # d x reps x T
    ll = [loglog_slope(T, mse[j], *late) for j in range(d)]
    mseT = [(mse[j].mean(0) * T)[T >= Tmax / 10].mean() for j in range(d)]      # pointwise MSE*T averaged over the late window (not <MSE><T>, which is biased by Jensen)
    D = A[:, :, 4]; E = A[:, :, 5]
    Dm = D.mean(0)[T >= Tmax / 10].mean(); Em = E.mean(0)[T >= Tmax / 10].mean()
    tex = r['texit']; res[key] = dict(s_late=s_late, ci_late=ci_late, s_mid=s_mid, c=c, ll=ll, mseT=mseT, Iinv=np.diag(Iinv), Dm=Dm, Em=Em, texit_med=np.median(tex), texit_p90=np.percentile(tex, 90))
    say(f"\n{key}  (K={r['K']}, d={d}, eps={eps}, start={r['start']}, reps={r['nr']}, T up to {Tmax:.0f})")
    say(f"  P1 regret slope vs ln T: late [{late[0]:.0f},{late[1]:.0f}] {s_late:.2f}  CI [{ci_late[0]:.2f},{ci_late[1]:.2f}]   mid [{mid[0]:.0f},{mid[1]:.0f}] {s_mid:.2f} [{ci_mid[0]:.2f},{ci_mid[1]:.2f}]   target d/2 = {d/2}")
    say(f"     constant c = E R - (d/2) ln T (late window mean) = {c:.2f}")
    say("  P2 raw error: log-log slope (late) " + ", ".join(f"v{j}: {ll[j][0]:.2f} [{ll[j][1][0]:.2f},{ll[j][1][1]:.2f}]" for j in range(d)) + "   target -1")
    say("     MSE*T (late) / (I^-1)_jj : " + ", ".join(f"v{j}: {mseT[j]:.3g}/{Iinv[j,j]:.3g} = {mseT[j]/Iinv[j,j]:.2f}" for j in range(d)))
    say(f"  P3 E D_n = {Dm:.2f}, E e_n = {Em:.2f}   (target d = {d})")
    thstar = np.concatenate([[r['mu']], (1 - eps) * r['beta'] / r['K']]); ratio_j = np.diag(Iinv) / (0.25 * thstar) ** 2; jb = int(np.argmax(ratio_j))
    tpred = 1.4 * ratio_j[jb]; res[key]['tpred'] = tpred; res[key]['jb'] = jb
    say(f"  P4 T_exit(0.25): median {np.median(tex):.0f}, p90 {np.percentile(tex,90):.0f};  eps*median = {eps*np.median(tex):.1f}   (mixing form c_exit/(beta delta^2) = {2/(r['beta'][0]*0.25**2):.0f})")
    say(f"     information form  z^2 max_j (I^-1)_jj/(delta theta*_j)^2 with z^2 = 1.4: binding coordinate {jb}, prediction {tpred:.0f}, observed/predicted = {np.median(tex)/tpred:.2f}")

# ---- Figure 1: regret
k1 = [k for k in runs if k.startswith('k1_') and k.endswith('_near')]
k1.sort(key=lambda k: -runs[k]['eps'])
fig, ax = plt.subplots(1, 2, figsize=(11.5, 4.4))
for k in k1:
    r = runs[k]; T = r['T']; R = r['A'][:, :, 3]; m, h = ci(R); e = r['eps']
    ax[0].plot(T, m, color=cols[e], label=f"eps = {e}"); ax[0].fill_between(T, m - h, m + h, color=cols[e], alpha=0.2)
    s, c95 = res[k]['s_late'], res[k]['ci_late']; lo, hi = T[-1] / 10, T[-1]
    ax[0].plot([lo, hi], [m[T >= lo][0], m[T >= lo][0] + s * np.log(hi / lo)], color=cols[e], ls='--', lw=1)
    ax[0].annotate(f"{s:.2f}", (hi, m[-1]), textcoords='offset points', xytext=(4, -2), color=cols[e], fontsize=8)
    ax[1].plot(T, (R - np.log(T)[None, :]).mean(0), color=cols[e], label=f"eps = {e}")
    ax[1].fill_between(T, *(lambda mm, hh: (mm - hh, mm + hh))(*ci(R - np.log(T)[None, :])), color=cols[e], alpha=0.2)
Tg = np.array([10, 3e4]); ax[0].plot(Tg, np.log(Tg) - np.log(10) + 0.3, 'k:', lw=1, label='slope d/2 = 1')
ax[0].set_xscale('log'); ax[0].set_xlabel('T (time units)'); ax[0].set_ylabel('E R_T  (pseudo-regret = expected regret)'); ax[0].set_title('K=1: regret vs T (dashed: fit on last decade; number = fitted slope)', fontsize=9); ax[0].legend(fontsize=8)
ax[1].set_xscale('log'); ax[1].set_xlabel('T'); ax[1].set_ylabel('E R_T - ln T'); ax[1].set_title('constant c(eps) in time units: flat in eps', fontsize=9); ax[1].legend(fontsize=8); ax[1].axhline(0, color='k', lw=0.5)
plt.tight_layout(); plt.savefig(os.path.join(FIG, 'p19_k1_regret.png'), dpi=140); plt.close()

# ---- Figure 2: raw error
fig, ax = plt.subplots(1, 2, figsize=(11.5, 4.4))
for j, name in ((0, 'mu'), (1, 'alpha')):
    for k in k1:
        r = runs[k]; T = r['T']; e = r['eps']; mse = (r['A'][:, :, 6 + j] ** 2); m = mse.mean(0)
        ax[j].loglog(T, m, color=cols[e], label=f"eps = {e}  (slope {res[k]['ll'][j][0]:.2f})")
        ax[j].loglog(T[T > T[-1] / 30], res[k]['Iinv'][j] / T[T > T[-1] / 30], color=cols[e], ls=':', lw=1)
    Tg = np.array([1, 3e4]); ax[j].loglog(Tg, 3 * Tg ** -1.0 * (0.3 if j == 0 else 0.03), 'k--', lw=0.8, label='slope -1')
    ax[j].set_xlabel('T'); ax[j].set_ylabel(f'E (hat {name} - {name}*)^2'); ax[j].legend(fontsize=7); ax[j].set_title(f'{name}: MSE (dotted: (I^-1)_jj / T from measured J_T/T)', fontsize=9)
plt.tight_layout(); plt.savefig(os.path.join(FIG, 'p19_k1_error.png'), dpi=140); plt.close()

# ---- Figure 3: normalised error
fig, ax = plt.subplots(1, 2, figsize=(11.5, 4.2))
for k in k1:
    r = runs[k]; T = r['T']; e = r['eps']
    for a_, col, nm in ((ax[0], 4, 'D_n = ||.||_F^2'), (ax[1], 5, 'e_n = ||.||_J^2')):
        m, h = ci(r['A'][:, :, col]); a_.semilogx(T, m, color=cols[e], label=f"eps = {e}"); a_.fill_between(T, m - h, m + h, color=cols[e], alpha=0.2)
        a_.set_title(f'E {nm}  (target d = 2)', fontsize=9)
for a_ in ax: a_.axhline(2, color='k', ls=':'); a_.set_xlabel('T'); a_.legend(fontsize=8); a_.set_ylim(0, 8)
plt.tight_layout(); plt.savefig(os.path.join(FIG, 'p19_k1_normalised.png'), dpi=140); plt.close()

# ---- Figure 4: burn-in and start cost
fig, ax = plt.subplots(1, 2, figsize=(11.5, 4.3))
ee = np.array([runs[k]['eps'] for k in k1]); med = np.array([np.median(runs[k]['texit']) for k in k1]); p90 = np.array([np.percentile(runs[k]['texit'], 90) for k in k1])
ax[0].loglog(1 / ee, med, 'o-', label='median T_exit (delta=0.25)'); ax[0].loglog(1 / ee, p90, 's--', label='p90')
ax[0].loglog(1 / ee, 32 / ee, 'k:', label='c_exit/(beta eps delta^2), c_exit = 2'); ax[0].loglog(1 / ee, [res[k]['tpred'] for k in k1], 'g^-.', label='information form: 1.4 max_j (I^-1)_jj/(delta theta_j)^2'); ax[0].set_xlabel('1/eps'); ax[0].set_ylabel('time'); ax[0].legend(fontsize=8); ax[0].set_title('burn-in: last exit from R(0.25)', fontsize=9)
far = [k for k in runs if k.startswith('k1_') and k.endswith('_far')]; far.sort(key=lambda k: -runs[k]['eps'])
for k in far:
    e = runs[k]['eps']; kn = f"k1_eps{e}_near"
    if kn not in runs: continue
    T = runs[k]['T']; Rf = runs[k]['A'][:, :, 3]; Rn = runs[kn]['A'][:, :, 3]
    ax[1].semilogx(T, Rf.mean(0) - Rn.mean(0), color=cols[e], label=f"eps = {e}  ((P0 far - P0 near)/2 = {0.5 * (e * 1.0 + (1 / e) * ((0.5 - (1 - e)) ** 2) - (e * 0.09 + (1 / e) * (0.2 * (1 - e)) ** 2)):.2f})")
ax[1].set_xlabel('T'); ax[1].set_ylabel('E R(far start) - E R(near start)'); ax[1].set_title('start cost of the far start (legend: predicted (P0_far - P0_near)/2)', fontsize=9); ax[1].legend(fontsize=7)
plt.tight_layout(); plt.savefig(os.path.join(FIG, 'p19_burnin.png'), dpi=140); plt.close()
say("\n== far-start excess over near start at the end of the run:")
for k in far:
    e = runs[k]['eps']; kn = f"k1_eps{e}_near"
    if kn not in runs: continue
    Rf = runs[k]['A'][:, :, 3][:, -1]; Rn = runs[kn]['A'][:, :, 3][:, -1]
    p0f = e * 1.0 + (1 / e) * (0.5 - (1 - e)) ** 2; p0n = e * 0.09 + (1 / e) * (0.2 * (1 - e)) ** 2
    se = np.sqrt(Rf.var(ddof=1) / len(Rf) + Rn.var(ddof=1) / len(Rn))
    say(f"  eps={e}: excess {Rf.mean() - Rn.mean():.2f} +- {1.96*se:.2f};  (P_0^far - P_0^near)/2 = {(p0f - p0n)/2:.2f}")

# ---- Figure 5: K=2
k2 = [k for k in runs if k.startswith('k2_')]; k2.sort(key=lambda k: -runs[k]['eps'])
if k2:
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.3))
    for k in k2:
        r = runs[k]; T = r['T']; e = r['eps']; R = r['A'][:, :, 3]; m, h = ci(R)
        ax[0].plot(T, m, color=cols[e], label=f"eps = {e}  (late slope {res[k]['s_late']:.2f})"); ax[0].fill_between(T, m - h, m + h, color=cols[e], alpha=0.2)
        for j, mk in ((0, '-'), (1, '--')):
            ax[1].loglog(T, (r['A'][:, :, 6 + j] ** 2).mean(0), color=cols[e], ls=mk, label=f"eps = {e}" if j == 0 else None)
        m2, h2 = ci(r['A'][:, :, 4]); ax[2].semilogx(T, m2, color=cols[e], label=f"eps = {e}")
    Tg = np.array([10, 3e4]); ax[0].plot(Tg, 1.5 * np.log(Tg) - 1.5 * np.log(10) + 0.5, 'k:', label='slope d/2 = 1.5')
    ax[0].set_xscale('log'); ax[0].set_xlabel('T'); ax[0].set_ylabel('E R_T'); ax[0].set_title('K=2 (beta = 1,2; mu = 1.5): regret', fontsize=9); ax[0].legend(fontsize=7)
    ax[1].loglog(Tg, 0.02 / Tg * 10, 'k:', label='slope -1'); ax[1].set_xlabel('T'); ax[1].set_title('K=2: MSE of mu (solid) and of alpha_1 (dashed)', fontsize=9); ax[1].legend(fontsize=7)
    ax[2].axhline(3, color='k', ls=':'); ax[2].set_ylim(0, 9); ax[2].set_xlabel('T'); ax[2].set_title('K=2: E D_n (target d = 3)', fontsize=9); ax[2].legend(fontsize=7)
    plt.tight_layout(); plt.savefig(os.path.join(FIG, 'p19_k2.png'), dpi=140); plt.close()
OUT.close()
