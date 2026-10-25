"""Problem 20: build the event stream from the Deribit BTC option trades (code/p20_data/events.bin, float64 seconds since 2024-11-01 00:00:00 UTC).
Event = distinct millisecond timestamp over all instruments (trades sharing a millisecond, typically fills of one order, are one event).
Also writes p20_data/prep_summary.txt (counts, ties, daily counts) used by the note."""
import glob, os, numpy as np, pandas as pd
base = r"C:\Users\PrathamDesai\Desktop\market-data\btc-1y\options_trades"
EPOCH_MS = int(pd.Timestamp('2024-11-01', tz='UTC').timestamp() * 1000)
out = []; stats = []
for f in sorted(glob.glob(os.path.join(base, '*.csv.gz'))):
    d = pd.read_csv(f, usecols=['timestamp', 'block_trade_id', 'liquidation'])
    ts = d.timestamp.values.astype(np.int64); u = np.unique(ts)
    stats.append((os.path.basename(f)[-14:-7], len(ts), len(u), int(d.block_trade_id.notna().sum()), int(d.liquidation.notna().sum()))); out.append(u)
t = np.unique(np.concatenate(out)); sec = (t - EPOCH_MS) / 1000.0
os.makedirs(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'p20_data'), exist_ok=True)
here = os.path.dirname(os.path.abspath(__file__)); sec.astype('<f8').tofile(os.path.join(here, 'p20_data', 'events.bin'))
day = (sec // 86400).astype(int); cnt = np.bincount(day)
with open(os.path.join(here, 'p20_data', 'prep_summary.txt'), 'w') as g:
    g.write("month  trades  distinct_ms  block  liquidation\n")
    for s in stats: g.write("%s %d %d %d %d\n" % s)
    g.write("total trades %d, distinct-ms events %d, tie fraction %.3f, days %d\n" % (sum(s[1] for s in stats), len(sec), 1 - len(sec) / sum(s[1] for s in stats), len(cnt)))
    g.write("daily events: min %d median %d max %d (day index of max %d)\n" % (cnt[cnt > 0].min(), np.median(cnt[cnt > 0]), cnt.max(), cnt.argmax()))
    g.write("day,events\n" + "\n".join("%d,%d" % (i, c) for i, c in enumerate(cnt)))
print(open(os.path.join(here, 'p20_data', 'prep_summary.txt')).read()[:1500])
