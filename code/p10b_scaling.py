"""Cold-start total information over [0,T], T=c/eps, K=1. M_T = info about mu (alpha known);
eff = det/I_aa = efficient info about mu (alpha unknown). Jensen floor: ln(1+beta T)/(mu beta)."""
import numpy as np
from p10b_transient import grids, info_matrix
beta = 1.0; c = 2.0
for mu in (1.0, 0.5, 0.3):
    k = 2*mu/beta; A = 2/(beta*(2*mu-beta)) if k > 1 else np.nan
    print(f"\n mu={mu} (k={k:.2f})  T=c/eps, c={c}   predicted log-slope A={A:.3f} (k>1);  Jensen coefficient 1/(mu beta)={1/(mu*beta):.3f}")
    print(f"{'eps':>8} {'T':>9} {'M_T':>9} {'eff info':>9} {'Jensen':>8} {'M_T/lnT':>8} {'dM/dlnT':>8} {'M_T/T^(1-k)':>12}")
    prev = None
    for eps in (3e-2, 1e-2, 3e-3, 1e-3, 3e-4):
        T = c/eps; g = grids(mu, beta, eps, N=80000, lo=1e-13)
        I, M, L = info_matrix(mu, beta, eps, T, g)
        eff = np.linalg.det(I)/I[1, 1]
        sl = "" if prev is None else f"{(M-prev[1])/(np.log(T)-np.log(prev[0])):8.3f}"
        pw = M/T**(1-k) if k < 1 else np.nan
        print(f"{eps:8.0e} {T:9.0f} {M:9.4f} {eff:9.4f} {np.log(1+beta*T)/(mu*beta):8.3f} {M/np.log(T):8.3f} {sl:>8} {pw:12.4f}")
        prev = (T, M)
