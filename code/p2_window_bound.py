"""Problem 2: window lower bound  I >= J_h  (Loewner order), where
 J_h = (1/h) E_pi int_0^h phi_s(u) phi_s(u)^T / lam_s(u) * exp(-Lambda_s(u)) du,
 phi_s(u)=(1, s_k e^{-beta_k u}), lam_s(u)=mu+sum alpha_k s_k e^{-beta_k u}, Lambda_s=int_0^u lam_s,
 s ~ stationary law of the state. Proof: I=(1/h)E int_0^h a a^T/lam dt over a stationary window; the integrand is >=0, so keep only the
 event {no arrival in [0,t)}, which has conditional probability exp(-Lambda_s(t)) given F_0 and on which a,lam are the deterministic phi_s, lam_s.
 Also reports cond of the exponential dictionary (Cauchy matrix 1/(beta_i+beta_j)) vs K, r."""
import sys, numpy as np
xq, wq = np.polynomial.legendre.leggauss(24); xq = .5*(xq+1); wq = .5*wq

def run(mu, beta, eps, n_ev, hs, seed=1, return_samples=False):
    K = len(beta); alpha = (1-eps)*beta/K; rng = np.random.default_rng(seed)
    S = np.zeros(K); lam = mu; A = np.zeros((K+1, K+1)); T = 0.0
    tau = 0.3*eps/mu; samples = []; nxt = 0.0
    for it in range(n_ev + 500):
        lam = mu + alpha@S; w = 0.0
        while True:
            w += rng.exponential(1/lam)
            if rng.random()*lam < mu + alpha@(S*np.exp(-beta*w)): break
        if it >= 500:                                   # burn-in 500 events
            s = w*xq; Sn = S[None, :]*np.exp(-np.outer(s, beta)); a = np.hstack([np.ones((24, 1)), Sn]); l = mu + Sn@alpha
            A += (w*wq/l*a.T)@a
            while nxt < T + w:                          # time-grid samples of the state
                samples.append(S*np.exp(-beta*(nxt-T))); nxt += tau
            T += w
        S = S*np.exp(-beta*w) + 1.0
    A /= T; Ss = np.array(samples)
    Js = {}
    for h in hs:
        u = h*xq; E = np.exp(-np.outer(u, beta))        # (Q,K)
        ph = np.concatenate([np.ones((len(Ss), 24, 1)), Ss[:, None, :]*E[None]], axis=2)   # (N,Q,K+1)
        lm = mu + (Ss[:, None, :]*E[None])@alpha
        Lam = mu*u[None] + (Ss*alpha)@((1-E)/beta[None]).T
        wt = wq*h/lm*np.exp(-Lam)/h
        Js[h] = np.einsum('nq,nqi,nqj->ij', wt, ph, ph)/len(Ss)
    return (A, Js, Ss) if return_samples else (A, Js)

def mineig(M): return np.linalg.eigvalsh(M)[0]

if __name__ == "__main__":
    n_ev = int(sys.argv[1]) if len(sys.argv) > 1 else 100000
    mu = 1.0
    print("Cauchy-matrix conditioning of the dictionary, C_ij=1/(beta_i+beta_j), geometric grid beta_k=r^(k-1)")
    for r in (1.5, 2.0, 4.0):
        print(f"  r={r}:", {K: f"{np.linalg.cond(1/(np.add.outer(r**np.arange(K), r**np.arange(K)))):.2e}" for K in (1, 2, 3, 4, 5, 6)})
    print(f"\nwindow bound: lam_min(I) (MC) vs best lam_min(J_h) over h in {{0.03..0.5}}*eps/mu ; Loewner check lam_min(I-J_h)")
    print(f"{'K':>2} {'r':>3} {'eps':>5} {'lmin(I)':>9} {'lmin(J)':>9} {'ratio':>7} {'best h/(eps/mu)':>15} {'min eig(I-J)':>12}")
    for K, r in ((1, 2.0), (2, 2.0), (2, 4.0), (3, 2.0)):
        beta = r**np.arange(K)
        for eps in (0.3, 0.1, 0.05):
            hs = [c*eps/mu for c in (0.03, 0.06, 0.12, 0.25, 0.5)]
            A, Js = run(mu, beta, eps, n_ev, hs)
            best = max(hs, key=lambda h: mineig(Js[h]))
            dc = lambda M: np.linalg.cond(M/np.outer(np.sqrt(np.diag(M)), np.sqrt(np.diag(M))))
            print(f"cond raw {np.linalg.cond(A):9.1f} diag-scaled {dc(A):7.2f} |", end=" ")
            print(f"{K:2d} {r:3.0f} {eps:5.2f} {mineig(A):9.4f} {mineig(Js[best]):9.4f} {mineig(Js[best])/mineig(A):7.3f} {best/(eps/mu):15.2f} {mineig(A-Js[best]):12.4f}", flush=True)
