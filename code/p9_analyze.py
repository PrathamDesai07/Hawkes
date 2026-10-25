"""Problem 9: summarise code/p9_out/ (output of p9_bayes.exe) against T1-T7 of problem09_lower_bound.md."""
import glob, os, re, sys, types, numpy as np
sys.modules.setdefault('mpmath', types.ModuleType('mpmath')); sys.modules['mpmath'].mp = types.SimpleNamespace(dps=60)
from p2_tangent_bound import matrices

def load(fn):
    h = open(fn).readline(); a = np.loadtxt(fn, comments='#')
    g = lambda k: float(re.search(k + r'=([\d.]+)', h).group(1))
    return dict(eps=g('eps'), mu=g(' mu'), beta=g('beta'), w=g(' w'), lnvol=float(re.search(r'ln\|Theta0\|=([-\d.]+)', h).group(1)), a=a)

def tang_N(mu, beta, eps):
    alpha = np.array([(1-eps)*beta]); M, M3, El = matrices(mu, alpha, np.array([beta])); best = None
    for f in np.concatenate([np.linspace(0.8, 3.0, 45), np.linspace(3.2, 8.0, 17)]):
        l0 = f*El; N = 2*M/l0 - M3/l0**2; N = (N+N.T)/2; ev = np.linalg.eigvalsh(N)
        if ev.min() > 0 and (best is None or np.linalg.det(N) > np.linalg.det(best)): best = N
    return best

runs = {}
for fn in sorted(glob.glob('p9_out/b_*.txt')):
    r = load(fn); runs[os.path.basename(fn)[2:-4]] = r

print("== (A) Bayes regret of Newton vs benchmark I_G and measured-covariance bound L_emp, selected T*eps   (mu=1, beta=1; d/2=1)")
print(f"{'eps':>6s} {'T':>8s} {'T*eps':>6s} | {'E R':>7s} {'se':>5s} {'L_emp':>7s} {'gap':>6s} {'I_G':>7s} {'E R-I_G':>8s} | {'kappa':>6s} {'L_emp<=ER':>9s}")
for key in sorted([k for k in runs if k.startswith('mu1_')], key=lambda k: -runs[k]['eps']):
    r = runs[key]; a = r['a']; T = a[:, 0]
    for te in (0.3, 1, 3, 10, 30, 100, 300):
        i = np.argmin(abs(T*r['eps'] - te))
        if abs(np.log(T[i]*r['eps']/te)) > 0.3: continue
        print(f"{r['eps']:6.3f} {T[i]:8.1f} {T[i]*r['eps']:6.2f} | {a[i,1]:7.3f} {a[i,2]:5.3f} {a[i,3]:7.3f} {a[i,4]:6.3f} {a[i,5]:7.3f} {a[i,6]:+8.3f} | {a[i,10]:6.3f} {str(a[i,3] <= a[i,1] + 2*a[i,2]):>9s}")
print("\n== (B) additive constant  c = E R - (d/2) ln T - (1/2) ln det Ibar - ln|Theta0|   at T*eps ~ 30 and ~ 100 (should be flat in eps, T)")
print(f"{'mu':>4s} {'eps':>6s} {'k0':>4s} | {'T*eps=30: c':>12s} {'T*eps=100: c':>13s} | {'ln|Th0|':>8s} {'1/2 ln detI':>11s} {'(1-k0)/2 ln(1/eps)':>19s}")
for key in sorted(runs, key=lambda k: (runs[k]['mu'], -runs[k]['eps'])):
    r = runs[key]; a = r['a']; T = a[:, 0]; out = []
    for te in (30, 100):
        i = np.argmin(abs(T*r['eps'] - te))
        if abs(np.log(T[i]*r['eps']/te)) > 0.35: out.append(float('nan')); continue
        out.append(a[i, 1] - np.log(T[i]) - 0.5*a[i, 11] - r['lnvol'])
    k0 = 2*r['mu']/r['beta']; ldI = a[np.argmin(abs(T*r['eps']-30)), 11]
    print(f"{r['mu']:4.1f} {r['eps']:6.3f} {k0:4.1f} | {out[0]:12.3f} {out[1]:13.3f} | {r['lnvol']:8.3f} {0.5*ldI:11.3f} {max(0,(1-k0)/2)*np.log(1/r['eps']):19.3f}")
print("\n== (C) rigorous-form lower bound  L_N = ln|Theta0| + (d/2) ln(T/(2 pi e kappa_N)) + (1/2) ln det N   (Prop 9.1 with A = tangent minorant N), at the largest T")
print(f"{'eps':>6s} {'T':>8s} | {'E R':>7s} {'L_emp':>7s} {'L_N':>7s} | {'kappa_N':>8s} {'ln det N':>9s} {'ln det I(MC)':>12s}")
for key in sorted([k for k in runs if k.startswith('mu1_')], key=lambda k: -runs[k]['eps']):
    r = runs[key]; a = r['a']; i = len(a)-1; T = a[i, 0]; N = tang_N(r['mu'], r['beta'], r['eps'])
    C = np.array([[a[i, 7], a[i, 9]], [a[i, 9], a[i, 8]]]); Nh = np.linalg.cholesky(N)          # N = Lh Lh'
    kN = np.linalg.eigvalsh(Nh.T @ (T*C) @ Nh).max()
    LN = r['lnvol'] + np.log(T/(2*np.pi*np.e*kN)) + 0.5*np.log(np.linalg.det(N))
    print(f"{r['eps']:6.3f} {T:8.1f} | {a[i,1]:7.3f} {a[i,3]:7.3f} {LN:7.3f} | {kN:8.3f} {np.log(np.linalg.det(N)):9.3f} {a[i,11]:12.3f}")
print("\n== (D) local slope of E R against ln T (crossover): slope over one decade ending at T*eps = x")
print(f"{'eps':>6s} " + " ".join(f"{'T*eps='+str(x):>11s}" for x in (1, 3, 10, 30, 100)))
for key in sorted([k for k in runs if k.startswith('mu1_')], key=lambda k: -runs[k]['eps']):
    r = runs[key]; a = r['a']; T = a[:, 0]; row = []
    for te in (1, 3, 10, 30, 100):
        j = np.argmin(abs(T*r['eps'] - te)); lo = np.argmin(abs(T - T[j]/10))
        row.append(f"{(a[j,1]-a[lo,1])/np.log(T[j]/T[lo]):11.2f}" if T[lo] < T[j]*0.5 and abs(np.log(T[j]*r['eps']/te)) < 0.3 and T[lo] >= T[0]*0.99 else f"{'-':>11s}")
    print(f"{r['eps']:6.3f} " + " ".join(row))
print("\n== (E) two-point bound, Prop 9.0:  E R_T/T  vs  u'Iu delta^2/2  (ratio must lie in [1/c~+, 1/c~-^2], -> 1 as delta -> 0)")
for fn in sorted(glob.glob('p9_out/tp_*.txt')):
    print(os.path.basename(fn)[:-4]); [print('   ' + l.rstrip()) for l in open(fn) if l.startswith('u=')]
print("\n== (F) scaling collapse: Bayes regret E R over the stability band as a function of x = T*eps alone (mu = beta = 1)")
xs = (0.3, 1, 3, 10, 30, 100, 300)
print(f"{'eps':>6s} " + " ".join(f"{'x='+str(x):>8s}" for x in xs) + "   (E R; asymptote ln x + const)")
for key in sorted([k for k in runs if k.startswith('mu1_')], key=lambda k: -runs[k]['eps']):
    r = runs[key]; a = r['a']; T = a[:, 0]; row = []
    for x in xs:
        j = np.argmin(abs(T*r['eps'] - x)); row.append(f"{a[j,1]:8.3f}" if abs(np.log(T[j]*r['eps']/x)) < 0.2 else f"{'-':>8s}")
    print(f"{r['eps']:6.3f} " + " ".join(row))
