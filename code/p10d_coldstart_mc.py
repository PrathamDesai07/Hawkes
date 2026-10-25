"""Cold-start (S_0=0) total Fisher information matrix I_T=E int_0^T a a^T/lam dt by path Monte Carlo, any K.
Reports I_mumu (alpha known), efficient info for mu (alpha unknown) = 1/(I^{-1})_mumu, vs Jensen floor ln(1+B T)/(mu B)."""
import numpy as np, sys
xq, wq = np.polynomial.legendre.leggauss(16); xq = 0.5*(xq+1); wq = 0.5*wq

def path(mu, beta, alpha, cps, rng):
    K = len(beta); S = np.zeros(K); t = 0.0; A = np.zeros((K+1, K+1)); out = []; ci = 0
    def acc(S0, w0, w1):
        s = w0 + (w1-w0)*xq; Sn = S0[None, :]*np.exp(-np.outer(s, beta))
        a = np.hstack([np.ones((len(s), 1)), Sn]); l = mu + Sn@alpha
        return ((w1-w0)*wq/l*a.T)@a
    while ci < len(cps):
        lam = mu + alpha@S; w = 0.0
        while True:                                       # thinning, cumulative candidate times
            w += rng.exponential(1/lam)
            if rng.random()*lam < mu + alpha@(S*np.exp(-beta*w)): break
        tn = t + w; w0 = 0.0
        while ci < len(cps) and cps[ci] <= tn:            # checkpoints inside this gap
            A += acc(S, w0, cps[ci]-t); w0 = cps[ci]-t; out.append(A.copy()); ci += 1
        if ci >= len(cps): break
        A += acc(S, w0, w); S = S*np.exp(-beta*w) + 1.0; t = tn
    return out

def run(mu, beta, eps, cps, R, seed=11, w=None):
    K = len(beta); w = np.ones(K)/K if w is None else w
    alpha = (1-eps)*w*beta; rng = np.random.default_rng(seed)
    acc = [np.zeros((K+1, K+1)) for _ in cps]
    for _ in range(R):
        for i, A in enumerate(path(mu, beta, alpha, cps, rng)): acc[i] += A
    return [A/R for A in acc]

if __name__ == "__main__":
    from p10b_transient import grids, info_matrix
    mu = 1.0
    # validation, K=1
    eps = 0.2; cps = [5.0, 20.0]
    est = run(mu, np.array([1.0]), eps, cps, 400)
    g = grids(mu, 1.0, eps)
    print("validation K=1 eps=0.2: (I_mm, I_ma, I_aa)  MC vs exact")
    for T, I in zip(cps, est):
        ex = info_matrix(mu, 1.0, eps, T, g)[0]
        print(f"  T={T}: MC ({I[0,0]:.3f},{I[0,1]:.3f},{I[1,1]:.2f})  exact ({ex[0,0]:.3f},{ex[0,1]:.3f},{ex[1,1]:.2f})")
    # K=2 and K=1 comparison at eps=0.05
    eps = 0.05; cps = [5.0, 20.0, 80.0, 320.0]; R = int(sys.argv[1]) if len(sys.argv) > 1 else 60
    print(f"\ncold start, eps={eps}, mu={mu}: Imu (alpha known), eff (alpha unknown), Jensen floor, T/ln-slope of eff")
    for beta in (np.array([1.0]), np.array([1.0, 2.0])):
        B = beta.max(); est = run(mu, beta, eps, cps, R if len(beta) > 1 else 60)
        print(f" beta={beta.tolist()}")
        prev = None
        for T, I in zip(cps, est):
            eff = 1/np.linalg.inv(I)[0, 0]
            sl = "" if prev is None else f"{(eff-prev[1])/(np.log(T)-np.log(prev[0])):7.3f}"
            print(f"   T={T:6.0f}  I_mm={I[0,0]:8.3f}  eff={eff:8.3f}  floor={np.log(1+B*T)/(mu*B):6.3f}  d eff/d lnT={sl}", flush=True)
            prev = (T, eff)
