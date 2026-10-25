"""Problem 14, K=2, alpha UNKNOWN: forgetting-Newton recursion for theta=(mu, alpha_1, alpha_2) under smooth drift of mu.
 Truth: beta=(1,2), alpha_k=(1-eps)beta_k/2 (n=1-eps), mu_t = mu_s(1+A sin(2 pi t/P)).
 Update per event: F <- (1-eta)F + int_gap phi phi^T/lam_hat du  (8-pt Gauss-Legendre), g = b - a/lam_hat(a),  theta -= F^{-1} g, with projection
 mu>=0.1, alpha_k>=0.01 and n_hat <= 1-eps/2 (scale alphas).   Prediction: bias^2 = (A^2/2) w^2/(eta^2+w^2) (x ~1.2-2 observed for K=1),  MSE_noise ~ a eta, a = r (I^{-1})_mumu / 2."""
import numpy as np, math, random, sys
from p2_identifiability import mc
xq, wq = np.polynomial.legendre.leggauss(8); xq = 0.5*(xq+1); wq = 0.5*wq

def run(eps, eta, T, seed, A, P, mu_s=1.5, beta=np.array([1., 2.]), grid=5.0, rho0=3.0):
    rng = random.Random(seed); K = 2; al_s = (1-eps)*beta/K; th = np.concatenate([[mu_s], al_s]); lam_bar = mu_s/eps
    F = np.diag([rho0/lam_bar] + [rho0*(lam_bar/beta[k])**2/lam_bar for k in range(K)]); S = np.zeros(K); Tt = 0.0
    mu_max = mu_s*(1+abs(A)); out = []; nextg = grid; last = (mu_s, al_s[0], al_s[1])
    while Tt < T:
        lam0 = mu_max + al_s@S; w = 0.0
        while True:
            w += rng.expovariate(lam0); mt = mu_s*(1 + A*math.sin(2*math.pi*(Tt+w)/P))
            if rng.random()*lam0 <= mt + al_s@(S*np.exp(-beta*w)): break
        e = np.exp(-beta*w); a = np.concatenate([[1.0], S*e]); b = np.concatenate([[w], S*(1-e)/beta])
        u = w*xq; Su = S[None, :]*np.exp(-np.outer(u, beta)); phi = np.hstack([np.ones((8, 1)), Su]); lam_u = phi@th
        dF = (w*wq/lam_u*phi.T)@phi
        F = (1-eta)*F + dF
        g = b - a/(a@th)
        th = th - np.linalg.solve(F, g)
        th[0] = max(th[0], 0.1); th[1:] = np.maximum(th[1:], 0.01)
        nh = (th[1:]/beta).sum()
        if nh > 1-eps/2: th[1:] *= (1-eps/2)/nh
        S = S*e + 1.0; Tt += w
        while Tt >= nextg: out.append(th[0]); nextg += grid
    return np.array(out)

if __name__ == "__main__":
    eps = float(sys.argv[1]); R = int(sys.argv[2]); T = 6000.0; P = 1200.0; A = float(sys.argv[3]); mu_s = 1.5; beta = np.array([1., 2.])
    r = mu_s/eps; tgrid = np.arange(1, int(T/grid_) + 1)*5.0 if False else np.arange(1, int(T/5.0)+1)*5.0
    target = mu_s*(1 + A*np.sin(2*np.pi*tgrid/P)); w = 2*math.pi/(P*r); bu = tgrid > 0.2*T
    I, _ = mc(mu_s, beta, eps, 250000, seed=4); Iinv_mm = np.linalg.inv(I)[0, 0]
    print(f"K=2 beta=(1,2) mu_s={mu_s} eps={eps}: r={r:.2f}, omega={w:.5f}, I_mumu={I[0,0]:.4f}, 1/eff=(I^-1)_mumu={Iinv_mm:.2f}  -> predicted noise a = r/(2 eff) = {r*Iinv_mm/2:.1f}", flush=True)
    for eta in (0.002, 0.004):
        reps = np.array([run(eps, eta, T, 40+s, A, P)[:len(tgrid)] for s in range(R)]); n = reps.shape[1]
        mean = reps.mean(0); varc = reps.var(0, ddof=1).mean()/R; tg = target[:n]; b = bu[:n]
        bias2 = np.mean(((mean-tg)**2)[b]) - varc; mse = np.mean(((reps-tg)**2)[:, b]); pred = (A*A*mu_s**2/2)*w*w/(eta*eta+w*w)
        print(f" eta={eta} (eta/eps^2={eta/eps**2:.3f}): bias^2 {bias2:.5f} vs filter {pred:.5f} (ratio {bias2/pred:.2f}); noise part {(mse-bias2):.4f} = {(mse-bias2)/eta:.1f}*eta", flush=True)
