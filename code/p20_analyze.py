"""Problem 20: figures and numbers from code/p20_out (fits of the exponential-dictionary Hawkes model to the Deribit BTC option-trade stream).
Figures -> figures/p20_*.png ; numbers -> code/p20_summary.txt"""
import os, re, datetime as dt
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__)); OUTD = os.path.join(HERE, 'p20_out'); FIG = os.path.join(HERE, '..', 'figures'); os.makedirs(FIG, exist_ok=True)
SUM = open(os.path.join(HERE, 'p20_summary.txt'), 'w')
def say(s=''):
    print(s); SUM.write(s + '\n')
D0 = dt.date(2024, 11, 1)
def date(i): return D0 + dt.timedelta(days=int(i))

def read_daily(fn):
    res = {}
    for l in open(fn):
        if l.startswith('#'): continue
        a, b, c = l.split('|')[0], l.split('|')[1], l.split('|')[2]; t = a.split(); day = int(t[0]); dic = t[1]; nev = int(t[2])
        nh, se = [float(x) for x in b.split()]; th = [float(x) for x in c.split()]; q = [float(x) for x in l.split('|')[3].split()]
        res.setdefault(dic, {})[day] = dict(nev=nev, n=nh, se=se, th=th, nll=q[0], ks=q[1], rho=q[2], me=q[3])
    return res
def read_online(fn):
    r = {}
    for l in open(fn):
        if l.startswith('#'): continue
        a, b, c = l.split('|'); t = a.split(); r[int(t[0])] = (int(t[1]), float(b), [float(x) for x in c.split()])
    return r
def read_rolling(fn):
    rows = []
    for l in open(fn):
        if l.startswith('#'): continue
        a, b, c = l.split('|'); t = a.split(); n, se = [float(x) for x in b.split()]; th = [float(x) for x in c.split()]; rows.append((float(t[0]), int(t[1]), n, se, th))
    return rows

dly = read_daily(os.path.join(OUTD, 'daily_13_20.txt')); d24 = read_daily(os.path.join(OUTD, 'daily_0_24.txt'))
onl5 = read_online(os.path.join(OUTD, 'online_eta5.txt')); onl1 = read_online(os.path.join(OUTD, 'online_eta1.txt')); onlg = read_online(os.path.join(OUTD, 'online_growing.txt'))
K3 = '0.1,1,10'; days = np.array(sorted(dly[K3].keys())); n3 = np.array([dly[K3][d]['n'] for d in days]); se3 = np.array([dly[K3][d]['se'] for d in days]); nev3 = np.array([dly[K3][d]['nev'] for d in days])
mu3 = np.array([dly[K3][d]['th'][0] for d in days])

say("== Problem 20: exponential-dictionary Hawkes fits to the Deribit BTC option-trade stream (events = distinct milliseconds, window 13:00-20:00 UTC, 2024-11-01 .. 2025-10-31)")
# --- dictionary selection (from select.txt)
say("\n-- dictionary selection (train = days 0-29, validation = days 30-59; pooled fit; log-lik per event, higher is better)")
for l in open(os.path.join(OUTD, 'select.txt')):
    if not l.startswith('#'): say('  ' + l.rstrip())

say(f"\n-- main dictionary beta = (0.1, 1, 10) /s (K = 3), {len(days)} daily windows")
say(f"   events per window: median {np.median(nev3):.0f}, range {nev3.min()}-{nev3.max()}")
say(f"   n_hat per day: median {np.median(n3):.3f}, 10-90 % {np.percentile(n3,10):.3f}-{np.percentile(n3,90):.3f}, min {n3.min():.3f} ({date(days[n3.argmin()])}), max {n3.max():.3f} ({date(days[n3.argmax()])});  median s.e. {np.median(se3):.3f}")
z = (n3 - np.median(n3)) / se3
say(f"   day-to-day dispersion of n_hat relative to its own s.e.: sd(z) = {z.std():.1f} (1 would mean no regime variation beyond estimation noise)")
cc = np.corrcoef(np.log(nev3), n3)[0, 1]; cm = np.corrcoef(np.log(nev3), np.log(mu3))[0, 1]
say(f"   corr(log events/day, n_hat) = {cc:+.2f};  corr(log events/day, log mu_hat) = {cm:+.2f}")
wk = np.array([date(d).weekday() >= 5 for d in days]); say(f"   weekend median n_hat {np.median(n3[wk]):.3f} vs weekday {np.median(n3[~wk]):.3f}; events/day weekend {np.median(nev3[wk]):.0f} vs weekday {np.median(nev3[~wk]):.0f}")
top = np.argsort(-nev3)[:6]; say("   six busiest days: " + "; ".join(f"{date(days[i])} ({nev3[i]} ev, n_hat {n3[i]:.3f}, mu {mu3[i]:.3f})" for i in top))
say(f"   monthly median n_hat: " + ", ".join(f"{date(days[m]).strftime('%Y-%m')} {np.median(n3[[date(d).month == date(days[m]).month and date(d).year == date(days[m]).year for d in days]]):.2f}" for m in sorted({np.where([(date(d).year, date(d).month) == (date(days[i]).year, date(days[i]).month) for d in days])[0][0] for i in range(len(days))})))

# --- misspecification: n_hat vs dictionary
say("\n-- misspecification: n_hat and fit quality across dictionaries (same days, same windows)")
say(f"   {'dictionary':14s} {'median n_hat':>12s} {'10-90 %':>14s} {'nll/ev (median)':>16s} {'median KS':>10s} {'frac KS>crit':>13s} {'median rho1':>12s}")
sumd = {}
for dic in ['1', '0.3,3', '0.1,1,10', '0.1,1,10,100']:
    ds = sorted(dly[dic].keys()); nn = np.array([dly[dic][d]['n'] for d in ds]); ks = np.array([dly[dic][d]['ks'] for d in ds]); ne = np.array([dly[dic][d]['nev'] for d in ds]); crit = 1.358 / np.sqrt(ne)
    nl = np.array([dly[dic][d]['nll'] for d in ds]); rho = np.array([dly[dic][d]['rho'] for d in ds]); sumd[dic] = (nn, nl, ks)
    say(f"   {dic:14s} {np.median(nn):12.3f} {np.percentile(nn,10):6.3f}-{np.percentile(nn,90):6.3f} {np.median(nl):16.3f} {np.median(ks):10.3f} {np.mean(ks > crit):13.2f} {np.median(rho):12.3f}")
dd = np.array([sumd['0.1,1,10'][0][i] - sumd['1'][0][i] for i in range(min(len(sumd['1'][0]), len(sumd['0.1,1,10'][0])))])
say(f"   median daily difference n_hat(K=3) - n_hat(K=1, beta=1/s) = {np.median(dd):+.3f} (IQR {np.percentile(dd,25):+.3f} .. {np.percentile(dd,75):+.3f})")

# --- seasonality: 13-20 vs 24 h
d24k = d24[K3]; common = [d for d in days if d in d24k]; n24 = np.array([d24k[d]['n'] for d in common]); n7 = np.array([dly[K3][d]['n'] for d in common])
ks24 = np.array([d24k[d]['ks'] for d in common]); nl24 = np.array([d24k[d]['nll'] for d in common])
say(f"\n-- seasonality: K=3 fitted to the full 24 h day vs the 13-20 UTC window: median n_hat {np.median(n24):.3f} vs {np.median(n7):.3f}; median KS {np.median(ks24):.3f} vs {np.median([dly[K3][d]['ks'] for d in common]):.3f}; mu_hat(24h) median {np.median([d24k[d]['th'][0] for d in common]):.3f} vs {np.median(mu3):.3f}")

# --- online vs batch
say("\n-- online forgetting Newton vs batch daily MLE (K=3), end-of-window n_hat")
for name, o in (('forgetting eta = 1/(5 N_prevday)', onl5), ('forgetting eta = 1/N_prevday', onl1), ('growing window (eta = 0)', onlg)):
    ds = [d for d in days if d in o]; no = np.array([o[d][1] for d in ds]); nb = np.array([dly[K3][d]['n'] for d in ds]); df = no - nb
    cs = np.cumsum(df) / (np.arange(len(df)) + 1)
    say(f"   {name:36s}: mean(online - batch) = {df.mean():+.3f}, sd = {df.std():.3f}, mean |.| = {np.abs(df).mean():.3f}; corr = {np.corrcoef(no, nb)[0,1]:+.2f}; first 30 days mean diff {df[:30].mean():+.3f}, last 120 days {df[-120:].mean():+.3f}; sd(online n_hat) = {no.std():.3f} vs sd(batch) = {nb.std():.3f}")

say("\n-- growing-window online Newton vs the POOLED batch MLE over the same days (K=3): same estimator target, so the difference is estimator error")
for l in open(os.path.join(OUTD, 'pooled.txt')):
    m = re.match(r'pooled days \[0,(\d+)\):.*\| ([\d.]+)\s*$', l.strip()); 
    if m:
        last = int(m.group(1)) - 1; ob = onlg.get(last)
        if ob: say(f"   days 0..{last}: pooled batch n_hat = {float(m.group(2)):.4f}, growing online n_hat = {ob[1]:.4f}, difference {ob[1] - float(m.group(2)):+.4f}")
say(f"   daily-window median n_hat = {np.median(n3):.3f}; event-weighted mean of daily n_hat = {np.average(n3, weights=nev3):.3f}; pooled MLE over the year 0.759: heterogeneity (busy days have larger n_hat) explains the offset between 'pooled' and 'median daily'")

# ---------- figures
cm_ = plt.get_cmap('tab10')
x = np.array([date(d) for d in days])
fig, ax = plt.subplots(2, 1, figsize=(12, 7.4), sharex=True, gridspec_kw={'height_ratios': [3, 1.2]})
ax[0].plot(x, n3, '.', ms=3, color='0.45', label='batch MLE per day (13-20 UTC), K=3, beta=(0.1,1,10)/s'); ax[0].fill_between(x, n3 - 1.96 * se3, n3 + 1.96 * se3, color='0.7', alpha=0.5, lw=0)
ax[0].plot([date(d) for d in sorted(onl5)], [onl5[d][1] for d in sorted(onl5)], color='tab:red', lw=1.4, label='online forgetting Newton, eta = 1/(5 N)')
ax[0].plot([date(d) for d in sorted(onlg)], [onlg[d][1] for d in sorted(onlg)], color='tab:blue', lw=1.2, ls='--', label='online growing window (eta = 0)')
ax[0].axvline(date(343), color='k', ls=':', lw=0.8); ax[0].text(date(343), ax[0].get_ylim()[0] if False else 0.38, ' 2025-10-10', fontsize=7, rotation=90, va='bottom')
ax[0].set_ylabel('estimated branching ratio  n_hat(t)'); ax[0].set_ylim(0.35, 0.95); ax[0].legend(fontsize=8, loc='lower left'); ax[0].set_title('Deribit BTC option trades: n_hat over one year (95 % bands: batch delta-method)', fontsize=10)
ax[1].semilogy(x, nev3, color='tab:green', lw=0.9); ax[1].set_ylabel('events / window'); ax[1].set_xlabel('date')
plt.tight_layout(); plt.savefig(os.path.join(FIG, 'p20_nhat_year.png'), dpi=140); plt.close()

fig, ax = plt.subplots(1, 3, figsize=(15, 4.4))
labs = {'1': 'K=1 (1/s)', '0.3,3': 'K=2 (0.3,3)', '0.1,1,10': 'K=3 (0.1,1,10)', '0.1,1,10,100': 'K=4 (0.1..100)'}
ax[0].boxplot([sumd[k][0] for k in labs], tick_labels=[labs[k] for k in labs], showfliers=False); ax[0].set_ylabel('n_hat per day'); ax[0].set_title('misspecification: n_hat grows with the dictionary', fontsize=9); ax[0].tick_params(axis='x', labelsize=7)
ax[1].boxplot([sumd[k][1] for k in labs], tick_labels=[labs[k] for k in labs], showfliers=False); ax[1].set_ylabel('in-sample -loglik / event'); ax[1].set_title('fit improves with K (log-lik per event)', fontsize=9); ax[1].tick_params(axis='x', labelsize=7)
ks_ = [sumd[k][2] for k in labs]; ax[2].boxplot(ks_, tick_labels=[labs[k] for k in labs], showfliers=False); ax[2].set_ylabel('KS distance of rescaled gaps to Exp(1)'); ax[2].set_title('time-rescaling residuals (all dictionaries reject)', fontsize=9); ax[2].tick_params(axis='x', labelsize=7)
plt.tight_layout(); plt.savefig(os.path.join(FIG, 'p20_misspecification.png'), dpi=140); plt.close()

fig, ax = plt.subplots(1, 2, figsize=(12, 4.4))
for nm, fn, c in (('2025-10-06 .. 2025-10-16 (includes 2025-10-10)', 'rolling_oct2025.txt', 'tab:red'), ('2025-02-28 .. 2025-03-07 (includes 2025-03-02)', 'rolling_mar2025.txt', 'tab:orange'), ('2025-09-15 .. 2025-09-25 (quiet reference)', 'rolling_calm.txt', 'tab:blue')):
    r = read_rolling(os.path.join(OUTD, fn)); t = np.array([q[0] for q in r]); n = np.array([q[2] for q in r]); se = np.array([q[3] for q in r]); ev_ = np.array([q[1] for q in r])
    ax[0].plot((t - t[0]) / 24, n, color=c, lw=1, label=nm); ax[0].fill_between((t - t[0]) / 24, n - 1.96 * se, n + 1.96 * se, color=c, alpha=0.2, lw=0)
    ax[1].semilogy((t - t[0]) / 24, ev_, color=c, lw=1, label=nm)
    say(f"-- rolling 1-hour fits, {nm}: median n_hat {np.median(n):.3f} (10-90 % {np.percentile(n,10):.2f}-{np.percentile(n,90):.2f}), median s.e. {np.median(se):.3f}; events/hour median {np.median(ev_):.0f}, max {ev_.max()}; sd(n_hat)/median s.e. = {n.std()/np.median(se):.1f}")
ax[0].set_xlabel('days from window start'); ax[0].set_ylabel('hourly n_hat (K=3)'); ax[0].legend(fontsize=7); ax[0].set_title('intra-day tracking: hourly batch fits', fontsize=9); ax[1].set_xlabel('days from window start'); ax[1].set_ylabel('events per hour'); ax[1].legend(fontsize=7)
plt.tight_layout(); plt.savefig(os.path.join(FIG, 'p20_rolling.png'), dpi=140); plt.close()

fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
ax[0].plot(x, np.array([d24k[d]['n'] if d in d24k else np.nan for d in days]), color='tab:purple', lw=0.8, label='full 24 h day (seasonality ignored)'); ax[0].plot(x, n3, color='0.4', lw=0.8, label='13-20 UTC window'); ax[0].set_ylabel('n_hat (K=3)'); ax[0].legend(fontsize=8); ax[0].set_title('seasonality as misspecification', fontsize=9)
ax[1].scatter(nev3, n3, s=6, c=[date(d).month for d in days], cmap='viridis'); ax[1].set_xscale('log'); ax[1].set_xlabel('events in window'); ax[1].set_ylabel('n_hat'); ax[1].set_title('n_hat vs activity level (colour = month)', fontsize=9)
plt.tight_layout(); plt.savefig(os.path.join(FIG, 'p20_seasonality_activity.png'), dpi=140); plt.close()
SUM.close()
