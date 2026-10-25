"""Problem 6, second analysis: potential defect, count of x>1 gaps, exit time T_exit, constants for k0<=1 and K=2,3 (reads p6_out2/*.txt)."""
import glob, os, re, numpy as np
def sl(T,R,lo,hi):
    m=(T>=lo)&(T<=hi); return np.polyfit(np.log(T[m]),R[m],1)[0]
print(f"{'run':26s} {'k0':>4s} {'d/2':>4s} | {'R-(d/2)lnT':>10s} {'se':>5s} {'slope':>5s} | {'defect':>7s} {'ln det':>7s} {'#x>1':>6s} {'events':>8s} | {'Texit med':>9s} {'p90':>7s} {'(k0-1)/(mu eps d^2)':>19s} {'1/(b eps d^2)':>13s}")
for fn in sorted(glob.glob('p6_out2/*.txt')):
    h=open(fn).readline(); txt=open(fn).read()
    K=int(re.search(r'K=(\d+)',h).group(1)); d=K+1
    eps=float(re.search(r'eps=([\d.]+)',h).group(1)); mu=float(re.search(r' mu=([\d.]+)',h).group(1)); b0=float(re.search(r'beta=([\d.]+)',h).group(1))
    a=np.loadtxt(fn,comments='#'); T=a[:,0]
    te=re.search(r'median=([\d.]+) p90=([\d.]+)',txt); med,p90=float(te.group(1)),float(te.group(2))
    k0=2*mu/b0; delta=0.25
    lev=(k0-1)/(mu*eps*delta**2) if k0>1 else float('nan'); mix=1/(b0*eps*delta**2)
    print(f"{os.path.basename(fn)[:-4]:26s} {k0:4.1f} {d/2:4.1f} | {a[-1,3]:10.2f} {a[-1,2]:5.2f} {sl(T,a[:,1],T[-1]/30,T[-1]):5.2f} | {a[-1,10]:7.2f} {a[-1,4]*2:7.1f} {a[-1,12]:6.1f} {a[-1,7]:8.0f} | {med:9.1f} {p90:7.1f} {lev:19.1f} {mix:13.1f}")
