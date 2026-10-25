"""Problem 21, real-data half: analyse code/p21_real_out/run.txt, time.txt (output of p21_real.exe) -> tables in code/p21_real_summary.txt, figure figures/p21_real.png.
Score = prequential point-process log-loss per event (lower is better), regret per event relative to the best FIXED parameter in hindsight over the evaluation period (pooled_test, in-sample)."""
import os, re, collections
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__)); OUT = open(os.path.join(HERE, 'p21_real_summary.txt'), 'w')
def say(s=''):
    print(s); OUT.write(s + '\n')
loss = collections.defaultdict(dict); nev = {}; dist = collections.defaultdict(dict); theta = {}; tim = []
def lines():
    for fn in ('run.txt', 'run_fast.txt', 'run_fast2.txt'):
        for l in open(os.path.join(HERE, 'p21_real_out', fn)): yield l
for l in lines():
    t = l.split()
    if not t: continue
    if t[0] == 'LOSS': d, n, m, v = int(t[1]), int(t[2]), t[3], float(t[4]); loss[m][d] = v; nev[d] = n
    elif t[0] == 'DIST': d, m = int(t[1]), t[2]; dist[m][d] = (float(t[3]), float(t[4]))
    elif t[0] == 'THETA': theta[t[1]] = l.strip()
for l in open(os.path.join(HERE, 'p21_real_out', 'time.txt')):
    if l.startswith('TIME'): tim.append(l.strip())
days = np.array(sorted(nev)); N = np.array([nev[d] for d in days]); TEST0 = 60
tmask = days >= TEST0
def tot(m, mask): return sum(loss[m][d] for d in days[mask])
Ntest = N[tmask].sum(); Nall = N.sum()
ref = tot('pooled_test', tmask)
say(f"== Problem 21 on the Deribit BTC option-trade stream (K=3, beta=(0.1,1,10)/s, windows 13-20 UTC, {len(days)} days, {Nall} events; evaluation period days >= {TEST0}: {Ntest} events)")
say(f"best fixed parameter in hindsight over the evaluation period (in-sample): {theta.get('pooled_test','')}")
say(f"best fixed parameter over the whole year: {theta.get('pooled_all','')}")
say(f"best per-day parameter (oracle, in-sample): loss/event {tot('oracle_day', tmask)/Ntest:.4f}; best fixed (pooled_test): {ref/Ntest:.4f}; the gap {1e3*(ref - tot('oracle_day', tmask))/Ntest:.1f} x 1e-3 nats/event is the room that tracking can exploit")
fam = {'ours (online Newton, eta = 1/(m N_prevday))': sorted([m for m in loss if m.startswith('newton_')]),
       'Solo-type EWLS (rho)': sorted([m for m in loss if m.startswith('ls_')]),
       'HW-type constant step': sorted([m for m in loss if m.startswith('hw_c')]), 'HW-type eta0/sqrt(i)': sorted([m for m in loss if m.startswith('hw_s')]), 'HW-type eta0/i': sorted([m for m in loss if m.startswith('hw_i')])}
sel = {}; best_test = {}
for f, ms in fam.items():
    dev = {m: sum(loss[m][d] for d in days[days < TEST0]) for m in ms}; sel[f] = min(dev, key=dev.get); best_test[f] = min(ms, key=lambda m: tot(m, tmask))
rows = []
for f in fam:
    rows.append((f'{f}: chosen on days 0-59 -> {sel[f]}', sel[f])); rows.append((f'{f}: best on the test period [oracle] -> {best_test[f]}', best_test[f]))
rows += [('batch MLE refit daily on all previous days', 'batch_all'), ('batch MLE refit daily on the previous 30 days', 'batch_30d'), ('ours growing window (eta=0, no hyperparameter)', 'newton_g')]
rows += [('best per-day parameter (oracle, in-sample)', 'oracle_day'), ('best fixed parameter, test period (in-sample)', 'pooled_test')]
say(f"\n{'method':70s} {'loss/event':>11s} {'excess vs best fixed (1e-3 nats/event)':>40s} {'excess, days 60-364 total (nats)':>34s}")
res = {}
for lab, m in rows:
    L = tot(m, tmask); res[m] = L; say(f"{lab:70s} {L/Ntest:11.4f} {1e3*(L - ref)/Ntest:40.2f} {L - ref:34.1f}")
say("\nrank by total test-period log-loss among the methods chosen on days 0-59 (plus batch and growing): " + " < ".join(m for m, _ in sorted([(m, res[m]) for lab, m in rows if ('chosen' in lab or lab.startswith('batch') or 'growing' in lab)], key=lambda x: x[1])))
say("\n-- distance to the pooled batch MLE of the same days: e = (theta - theta_batch)' H (theta - theta_batch) (H observed information of the pooled fit; ~ chi^2_4 = 4 would be one standard error), and n_hat - n_hat(batch)")
cols = [5, 30, 120, 364]
say(f"{'method':32s} " + " ".join(f"{'day '+str(c):>22s}" for c in cols))
for m in ['newton_m5', 'newton_m1', 'newton_g', 'hw_c0.1', 'hw_s1', 'hw_i1', 'ls_r1', 'ls_r20000', 'ls_r2000', 'ls_r200']:
    say(f"{m:32s} " + " ".join(f"{dist[m][c][0]:12.1f} ({dist[m][c][1]:+.3f})" for c in cols))
say("   pooled batch n_hat at those days: " + ", ".join(f"day {c}: {dist['batch_pooled'][c][1]:.3f}" for c in cols))
say("\n-- CPU on the stream"); [say('   ' + t) for t in tim]

# figure
fig, ax = plt.subplots(1, 2, figsize=(13, 4.6))
dt = days[tmask]; Nt = N[tmask]
def cum(m): return np.cumsum([loss[m][d] - loss['pooled_test'][d] for d in dt])
for lab, m, st in (('ours, m=0.05 (chosen on days 0-59)', 'newton_m0.05', 'r-'), ('ours, eta=1/N', 'newton_m1', 'm-'), ('ours, growing', 'newton_g', 'b-'), ('batch all history', 'batch_all', 'c-'), ('batch 30 days', 'batch_30d', 'g-'),
                   ('EWLS rho=1-1/100 (chosen)', 'ls_r100', 'y-'), ('EWLS rho=1', 'ls_r1', 'y--'), ('HW-type eta0/sqrt(i), eta0=3 (chosen)', 'hw_s3', 'k--'), ('oracle per-day', 'oracle_day', 'k:')):
    ax[0].plot(dt, cum(m), st, lw=1.2, label=lab)
ax[0].axhline(0, color='0.5', lw=0.5); ax[0].set_xlabel('day'); ax[0].set_ylabel('cumulative log-loss minus best fixed theta (nats)'); ax[0].set_title('prequential regret vs the best fixed parameter (days >= 60); below 0 = beats a fixed parameter', fontsize=8); ax[0].legend(fontsize=6, ncol=2)
for m, lab, st in (('newton_m5', 'ours eta=1/(5N)', 'r-o'), ('newton_m1', 'ours eta=1/N', 'm-o'), ('newton_g', 'ours growing', 'b-o'), ('ls_r1', 'EWLS rho=1', 'y-o'), ('ls_r2000', 'EWLS 1-1/2000', 'y--o'), ('hw_s1', 'HW-type eta0/sqrt(i)', 'k--o')):
    cs = sorted(dist[m]); ax[1].semilogy(cs, [max(dist[m][c][0], 1e-3) for c in cs], st, ms=3, lw=1, label=lab)
ax[1].axhline(4, color='0.5', ls=':'); ax[1].set_xlabel('day'); ax[1].set_ylabel("e = (theta - theta_batch)' H (theta - theta_batch)"); ax[1].set_title('distance to the pooled batch MLE of the same days (dotted: chi^2_4 mean)', fontsize=8); ax[1].legend(fontsize=7)
plt.tight_layout(); plt.savefig(os.path.join(HERE, '..', 'figures', 'p21_real.png'), dpi=140); plt.close(); OUT.close()
