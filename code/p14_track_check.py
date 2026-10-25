"""Problem 14, tracking theorem: check of the proof's recursion and of the final bound (K = 1, both parameters unknown, forgetting Newton, exact F-projection).
 comparator of gap i: th*_i = th*(t_{i-1}) (value at the gap start);  v_i = th_i - th*_i;  D_{i-1} = ||v_i||^2_{F_{i-1}};  D_i = ||th_{i+1} - th*_{i+1}||^2_{F_i};  d_i = th*_{i+1} - th*_i.
 (P) pathwise:  D_i <= (1+rho) [ (1-eta) D_{i-1} + X_i - 2<g_i,v_i> + G_i ] + (1+1/rho) ||d_i||^2_{F_i},  X_i = v' dF v, G_i = g' F_i^-1 g  (any rho>0; here rho = eta/2).
 (E) expectation: E[X_i - 2<g_i,v_i>] <= E Xi_i,  Xi_i = int (phi'(th*_i - th*_{t_{i-1}+u}))^2/lam_th du.
 (S) steady state:  mean D <= (2/eta)[(1+eta/2)(mean G + mean Xi) + (1+2/eta) mean ||d||^2_F]  (measured means on the RHS);
 (F) final:  mean ||v||^2_{I_ev} <= 2 eta d (1+..) + 4 (||d||_{I_ev}/eta)^2  with I_ev = I(th*)/r  (exact Fisher per event)."""
import math, random, sys
import numpy as np
from p5_online_newton_k1 import fisher_inc, xq, wq
from p5_projected_exact import project_F, fnorm2
from p10_fisher_k1 import exact_m

def run(seed, eps, eta, A, P, T, beta=1.0, mu_s=1.0, rho0=1.0):
    rng = random.Random(seed); al_s = (1 - eps)*beta; mu_t = lambda t: mu_s*(1 + A*math.sin(2*math.pi*t/P))
    lam_bar = mu_s/eps; F = [[rho0/lam_bar, 0.0], [0.0, rho0*(lam_bar/beta)**2/lam_bar]]
    lo = (0.2*mu_s, 0.3*al_s); hi = (3*mu_s, beta*(1 - eps/2)); th = (mu_s, al_s); s = 0.0; t = 0.0
    D_prev = fnorm2(F, th[0] - mu_t(0.0), th[1] - al_s)
    acc = dict(viol=0, n=0, D=0.0, G=0.0, Xi=0.0, d2F=0.0, XmG=0.0, vI=[], Dsum=0.0)
    ms = {}  # cache exact info matrices
    def Iev(mu):
        key = round(mu, 3)
        if key not in ms:
            m = exact_m(mu, al_s, beta); El = mu/eps
            ms[key] = np.array([[m, (1 - mu*m)/al_s], [(1 - mu*m)/al_s, (El - 2*mu + mu*mu*m)/al_s**2]])/(mu/eps)
        return ms[key]
    n = 0
    while t < T:
        # thinning with the TRUE time-varying intensity: upper bound mu_max + al_s*s
        lam0 = mu_s*(1 + abs(A)) + al_s*s; w = 0.0
        while True:
            w += rng.expovariate(lam0)
            if rng.random()*lam0 <= mu_t(t + w) + al_s*s*math.exp(-beta*w): break
        e = math.exp(-beta*w); a1 = s*e; b1 = s*(1 - e)/beta; ths = (mu_t(t), al_s)       # comparator at gap start
        v0, v1 = th[0] - ths[0], th[1] - ths[1]
        # gradient at th
        lam_th = th[0] + th[1]*a1; g0 = w - 1.0/lam_th; g1 = b1 - a1/lam_th
        f0, f1, f2 = fisher_inc(th[0], th[1], beta, s, w); Fp = [[F[0][0], F[0][1]], [F[1][0], F[1][1]]]
        D_im1 = fnorm2(Fp, v0, v1)                                              # D_{i-1} = ||th_i - th*_i||^2_{F_{i-1}}
        Xi_ = f0*v0*v0 + 2*f1*v0*v1 + f2*v1*v1
        gv = g0*v0 + g1*v1
        F = [[(1 - eta)*F[0][0] + f0, (1 - eta)*F[0][1] + f1], [(1 - eta)*F[1][0] + f1, (1 - eta)*F[1][1] + f2]]
        det = F[0][0]*F[1][1] - F[0][1]**2
        d0 = (F[1][1]*g0 - F[0][1]*g1)/det; d1 = (-F[0][1]*g0 + F[0][0]*g1)/det; G = g0*d0 + g1*d1
        th_new, _ = project_F(F, (th[0] - d0, th[1] - d1), lo, hi)
        t_next = t + w; ths_next = (mu_t(t_next), al_s); dvec = (ths_next[0] - ths[0], 0.0)
        D_i = fnorm2(F, th_new[0] - ths_next[0], th_new[1] - ths_next[1]); d2F = fnorm2(F, dvec[0], dvec[1])
        rho = eta/2
        rhs = (1 + rho)*((1 - eta)*D_im1 + Xi_ - 2*gv + G) + (1 + 1/rho)*d2F
        # in-gap drift energy  Xi_gap = int (phi'(th*_i - th*(t+u)))^2/lam_th du  (only mu drifts)
        Xg = 0.0
        for x, wt in zip(xq, wq):
            u = 0.5*w*(x + 1); lt = th[0] + th[1]*s*math.exp(-beta*u); dm = ths[0] - mu_t(t + u); Xg += 0.5*w*wt*dm*dm/lt
        if t > 0.3*T:
            acc['viol'] += D_i > rhs*(1 + 1e-9) + 1e-12; acc['n'] += 1; acc['D'] += D_i; acc['G'] += G; acc['Xi'] += Xg; acc['d2F'] += d2F; acc['XmG'] += Xi_ - 2*gv
            Ie = Iev(ths_next[0]); vv = np.array([th_new[0] - ths_next[0], th_new[1] - ths_next[1]]); acc['vI'].append((vv@Ie@vv, dvec[0]**2*Ie[0, 0]))
        th = th_new; s = a1 + 1.0; t = t_next
    return acc

if __name__ == "__main__":
    eps = float(sys.argv[1]); eta = float(sys.argv[2]); A = float(sys.argv[3]); R = int(sys.argv[4]); P = 1200.0; T = 8000.0
    out = [run(100 + r, eps, eta, A, P, T) for r in range(R)]
    n = sum(o['n'] for o in out)
    viol = sum(o['viol'] for o in out); mD = sum(o['D'] for o in out)/n; mG = sum(o['G'] for o in out)/n; mXi = sum(o['Xi'] for o in out)/n
    md2 = sum(o['d2F'] for o in out)/n; mXmG = sum(o['XmG'] for o in out)/n
    vI = np.array([x for o in out for x in o['vI']]); mv = vI[:, 0].mean(); mdI = vI[:, 1].mean()
    d = 2
    rhsS = (2/eta)*((1 + eta/2)*(mG + mXi) + (1 + 2/eta)*md2)
    print(f"eps={eps} eta={eta} A={A}  ({n} gaps after burn-in, {R} replicates)")
    print(f"  (P) pathwise recursion violations: {viol} of {n}")
    print(f"  (E) mean(X - 2<g,v>) = {mXmG:.4g}  <=  mean Xi_gap = {mXi:.4g} ?   {'yes' if mXmG <= mXi + 3*abs(mXmG)*0 + 1e-12 + 0.02*abs(mG) else 'check'}  (mean G = {mG:.4g}; eta*d = {eta*d:.4g})")
    print(f"  (S) mean D = {mD:.4g}  <=  steady-state RHS = {rhsS:.4g}  (ratio {mD/rhsS:.3f})   [mean ||d||^2_F = {md2:.3g}]")
    bound = 2*eta*d + 4*mdI/eta**2
    print(f"  (F) mean ||v||^2_{{I_ev}} = {mv:.4g}  vs closed-form bound 2 eta d + 4 ||d||^2_I/eta^2 = {bound:.4g} (ratio {mv/bound:.3f});  noise-only prediction eta d/2 = {eta*d/2:.4g}, filter drift prediction ||d||^2_I/eta^2-ish = {mdI/eta**2:.4g}")
