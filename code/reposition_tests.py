"""Evidence for repositioning C3/C4.
 E1  exact K=1 stationary:  1/2 log det I(theta*; eps)  vs log(1/eps), for k0 = 2 mu/beta in {0.2,0.6,1,2,4}.   det I = (A-1)/alpha^2, A = E[1/lam] E[lam].
 E3  K=2 sum of exponentials, weights w, rates beta: Jaisson-Rosenbaum predicts a Gamma/CIR limit with shape 2*m*mu, m = sum w_k/beta_k (mean lifetime).
     Measure A_T = E[1/lam]E[lam] (time-average by quadrature) for 2 m mu < 1 and > 1; prediction A -> k/(k-1) (k>1) and A ~ eps^{-(1-k)} (k<1)."""
import math, random, sys
import numpy as np
from p10_fisher_k1 import exact_m

def E1():
    print("E1: exact K=1, beta=1, mu = k0/2. Columns: eps, A=E[1/lam]E[lam], 1/2 log det I; then local slope d(1/2 log det I)/d log(1/eps) between successive eps")
    for k0 in (0.2, 0.6, 1.0, 2.0, 4.0):
        mu = k0/2; rows = []
        for eps in (0.3, 0.1, 0.03, 0.01, 0.003, 0.001):
            al = 1 - eps; m = exact_m(mu, al, 1.0); A = m*mu/eps; det = (A - 1)/al**2
            rows.append((eps, A, 0.5*math.log(det)))
        sl = [(rows[i+1][2] - rows[i][2])/math.log(rows[i][0]/rows[i+1][0]) for i in range(len(rows)-1)]
        pred = (1 - k0)/2 if k0 < 1 else 0.0
        print(f"  k0={k0}: " + "  ".join(f"[{e:g}: A={a:.3g}, {h:+.2f}]" for e, a, h in rows))
        print(f"          slopes {[round(s, 3) for s in sl]}   asymptotic prediction (1-k0)/2 for k0<1, 0 for k0>1: {pred:.2f}" + ("  (k0=1: expected ~ (1/2) d log log(1/eps), i.e. -> 0 slowly)" if k0 == 1.0 else ""))

xq, wq = np.polynomial.legendre.leggauss(24); xq = 0.5*(xq + 1); wq = 0.5*wq
def E3_run(seed, mu, beta, w, eps, T):
    K = len(beta); alpha = (1 - eps)*beta*w; rng = np.random.default_rng(seed); S = np.zeros(K); t = 0.0; I = 0.0
    for _ in range(int(20*50/eps**2*0.02)):                           # burn-in
        lam = mu + alpha@S; ww = 0.0
        while True:
            ww += rng.exponential(1/lam)
            if rng.random()*lam < mu + alpha@(S*np.exp(-beta*ww)): break
        S = S*np.exp(-beta*ww) + 1.0
    while t < T:
        lam = mu + alpha@S; ww = 0.0
        while True:
            ww += rng.exponential(1/lam)
            if rng.random()*lam < mu + alpha@(S*np.exp(-beta*ww)): break
        s = ww*xq; l = mu + (S[None, :]*np.exp(-np.outer(s, beta)))@alpha
        I += ww*np.sum(wq/l); t += ww; S = S*np.exp(-beta*ww) + 1.0
    return I/t

def E3(reps):
    beta = np.array([1.0, 4.0]); w = np.array([0.5, 0.5]); m = float(np.sum(w/beta))
    print(f"E3: K=2, beta={beta.tolist()}, w={w.tolist()}, mean lifetime m = sum w/beta = {m}; shape k = 2 m mu; threshold mu* = 1/(2m) = {1/(2*m):.3f}")
    for mu in (0.4, 0.8, 1.6):
        k = 2*m*mu; res = []
        for eps in (0.3, 0.1, 0.03):
            gam = eps/m; T = 150/gam; vals = [E3_run(7 + r, mu, beta, w, eps, T) for r in range(reps)]
            Einv = sum(vals)/reps; A = Einv*mu/eps; res.append((eps, A))
        sl = [math.log(res[i+1][1]/res[i][1])/math.log(res[i][0]/res[i+1][0]) for i in range(len(res)-1)]
        pred = ("A -> k/(k-1) = %.2f" % (k/(k-1))) if k > 1 else ("A ~ eps^-(1-k): slope %.2f" % (1 - k))
        print(f"  mu={mu} (k={k:.2f}): A(eps) = " + ", ".join(f"{e:g}:{a:.2f}" for e, a in res) + f";  d log A/d log(1/eps) = {[round(s, 2) for s in sl]};  prediction: {pred}", flush=True)

if __name__ == "__main__":
    E1()
    E3(int(sys.argv[1]) if len(sys.argv) > 1 else 6)
