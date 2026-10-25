"""Stationary MC (K=2, beta=(1,2), alpha_k=(1-eps)beta_k/2) for the integral steps of Theorems E and F (problem11_dichotomy.md Sec 4b).
 (0) stationarity:  E_pi[A f_z] = 0  (own generator, exact A f at sampled states)
 (E) E[1/W] <= 2 gamma/delta0 at z0 = 4(mu/2 + C2)/delta0   (k_K > 1, needs 2 kbar gamma <= delta0/4)  and  E[1/lam] <= max(z0/mu, 1/(b-gamma)) 2 eps/(nu delta0)
 (F) (c0 z + C4) E[1/W^2] >= gamma + (c0 - mu - gamma z - 2 kbar gamma) E[1/W]  for z in {1,5,10}   (k_K < 1)"""
import numpy as np, math, sys
from verify_lemmaH_EF import setup, consts, Af

def sample_stationary(P, T, dt, seed, burn):
    be, al, mu = P['beta'], P['alpha'], P['mu']; K = len(be); rng = np.random.default_rng(seed); S = np.zeros(K); t = 0.0; out = []; nxt = burn
    while t < T + burn:
        lam = mu + al@S; w = 0.0
        while True:
            w += rng.exponential(1/lam)
            if rng.random()*lam < mu + al@(S*np.exp(-be*w)): break
        while nxt < t + w and nxt < T + burn:
            out.append(S*np.exp(-be*(nxt - t))); nxt += dt
        S = S*np.exp(-be*w) + 1.0; t += w
    return np.array(out)

def run(name, beta, n, mu, K, reps, Tmult=300.0, dt=0.5):
    beta = np.array(beta, float); alpha = n*beta/K; P = setup(beta, alpha, mu); g = P['gam']; T = Tmult/g
    samples = np.vstack([sample_stationary(P, T, dt, 11 + r, 30/g) for r in range(reps)]); N = len(samples)
    print(f"{name}: eps={P['eps']:.3f} gamma={g:.4f} nu~={P['nu']:.3f} c0={P['c0']:.3f} k_K={P['kK']:.3f} kbar={P['kbar']:.3f}; {N} grid samples (correlated; SE below is naive-iid x sqrt(10) as a rough allowance)")
    Z = samples@P['om']
    print(f"   E Z = {Z.mean():.2f} (theory mu/gamma = {mu/g:.2f})")
    for z in (1.0, 10.0):
        A, W = Af(P, samples, z); se = A.std()/math.sqrt(N)*math.sqrt(10)
        print(f"   (0) z={z}: E[A f_z] = {A.mean():+.2e} +- {se:.1e}   (typical |A f| = {np.mean(np.abs(A)):.2e})")
    c0, kb, mu = P['c0'], P['kbar'], P['mu']
    if P['kK'] > 1:
        d0 = mu - c0; Bz, C2, C4 = consts(P, 1e9); z0 = 4*(mu/2 + C2)/d0; Bz, C2, C4 = consts(P, z0); z0 = max(z0, 4*(mu/2 + C2)/d0)
        W = Z + z0; lam = mu + samples@P['alpha']
        print(f"   (E) condition 2 kbar gamma <= delta0/4: {2*kb*g:.3f} <= {d0/4:.3f}: {2*kb*g <= d0/4};  z0 = {z0:.2f}")
        print(f"       E[1/W] = {np.mean(1/W):.4f}  <=? 2 gamma/delta0 = {2*g/d0:.4f};   E[1/lam] = {np.mean(1/lam):.4f}  <=? bound {max(z0/mu, 1/(P['b']-g))*2*P['eps']/(P['nu']*d0):.4f}   (E[1/lam]/eps = {np.mean(1/lam)/P['eps']:.3f})")
    else:
        for z in (1.0, 5.0, 10.0):
            Bz, C2, C4 = consts(P, z); W = Z + z; lhs = (c0*z + C4)*np.mean(1/W**2); rhs = g + (c0 - mu - g*z - 2*kb*g)*np.mean(1/W)
            print(f"   (F) z={z}: (c0 z + C4) E[1/W^2] = {lhs:.4f}  >=? gamma + (c0-mu-gamma z-2 kbar gamma) E[1/W] = {rhs:.4f}   [{lhs >= rhs}]")

if __name__ == "__main__":
    reps = int(sys.argv[1])
    run("K=2 beta=(1,2) mu=1.5 (k_K>1)", [1.0, 2.0], 0.95, 1.5, 2, reps)
    run("K=2 beta=(1,2) mu=0.4 (k_K<1)", [1.0, 2.0], 0.95, 0.4, 2, reps)
    run("K=2 beta=(1,2) mu=1.5, small eps (k_K>1)", [1.0, 2.0], 0.98, 1.5, 2, max(2, reps//2), Tmult=150.0)
