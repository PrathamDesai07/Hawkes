"""Problem 8: summarise code/p8_out/*.txt against T1-T6 of problem08_regret_notion.md (last row = n = 20000 events)."""
import glob, os, re, numpy as np
print(f"{'run':20s} {'(d/2)ln n':>9s} {'E r':>7s} {'E R':>7s} {'r-R (se)':>14s} {'R/(Q/2)':>7s} {'Dn':>5s} {'avg-err':>7s} {'avg/(2cp R)':>12s} {'hind':>6s} {'hind-r':>6s} {'wilks':>6s} {'viol':>4s} {'gap R/(Q/2) range':>18s}")
for fn in sorted(glob.glob('p8_out/*.txt')):
    h = open(fn).readline(); d = int(re.search(r'd=(\d+)', h).group(1))
    rows = [l for l in open(fn) if not l.startswith('#')]
    t = rows[-1].replace('[', ' ').replace(']', ' ').replace(',', ' ').replace('minres=', ' ').replace('|', ' ').split()
    n = int(t[0]); v = [float(x) for x in t[1:]]
    Er, ser, ER, seR, dm, sd, Qh, ratio, Dn, en, av, avfrac, hind, wil, glo, ghi, viol = v[:17]
    print(f"{os.path.basename(fn)[:-4]:20s} {d/2*np.log(n):9.2f} {Er:7.2f} {ER:7.2f} {dm:+7.2f}({sd:.2f}) {ratio:7.3f} {Dn:5.2f} {av:7.2f} {avfrac:12.3f} {hind:6.2f} {hind-Er:6.2f} {wil:6.2f} {int(viol):4d} [{glo:.2f},{ghi:.2f}]")
