"""Problem 11, K=1: two-point (Le Cam) lower bounds from Fisher information.
KL(P_th^T || P_th'^T) = E_th int [lam log(lam/lam') - (lam-lam')]dt <= E_th int (lam-lam')^2/lam' dt <= rho * u^T I_T(th) u,
u = th - th', rho = sup lam/lam' <= max(mu/mu', alpha/alpha').  Direction u=(delta, -w*delta), w = I_ma/I_aa (efficient).
Here: u^T I_T u for the cold-start exact I_T (p10b_transient) and the stationary T*I(theta*) (Theorem A)."""
import numpy as np
from p10b_transient import grids, info_matrix
from p10_fisher_k1_exact import exact_m

def quad_form(I, delta, w): u = np.array([delta, -w*delta]); return u@I@u

def table(mu, beta, c, Ts):
    k = 2*mu/beta
    print(f"\n mu={mu}, beta={beta} (k={k:.2f}), eps=c/T with c={c}")
    print(f"{'T':>7} {'eps':>8} | stationary: T*eff  | cold start: M_T  eff_T | KL/delta^2 (cold, w=eps*beta/mu) | delta for KL<=0.5/rho")
    for T in Ts:
        eps = c/T; a = (1-eps)*beta
        # stationary start: I(theta*) from Theorem A quantities
        m = exact_m(mu, beta, eps); El = mu/eps
        Ist = np.array([[m, (1-mu*m)/a], [(1-mu*m)/a, (El-2*mu+mu**2*m)/a**2]])
        eff_st = 1/np.linalg.inv(Ist)[0, 0]
        g = grids(mu, beta, eps, N=80000, lo=1e-14)
        Ic, M, L = info_matrix(mu, beta, eps, T, g); eff_c = 1/np.linalg.inv(Ic)[0, 0]
        w = eps*beta/mu                                # = I_ma/I_aa to leading order, cheap, explicit
        qf = quad_form(Ic, 1.0, w)
        dmax = np.sqrt(0.5/qf)                         # rho ~ 1 for small delta
        print(f"{T:7.0f} {eps:8.1e} | {T*eff_st:18.4f} | {M:9.3f} {eff_c:8.3f} | {qf:30.3f} | {dmax:8.3f}")

if __name__ == "__main__":
    for mu in (1.0, 0.5, 0.3): table(mu, 1.0, 2.0, (100, 1000, 10000))
