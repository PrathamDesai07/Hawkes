"""Mechanism of the extra lag for alpha unknown (K=1). Two-parameter estimator, one update per event, three preconditioners:
 'adaptive': F <- (1-eta)F + dF_i (data-driven, realised information)  [what was tested before]
 'fixed'   : step = eta * I_ev^{-1} g with I_ev = I(theta*_ref)/r constant (deterministic matrix, exact I from m)
 'mean_dF' : step = eta * Fbar_i^{-1} g with Fbar_i = (1-eta)Fbar + E[dF_i | state] ... not used
Compare ensemble-mean bias^2 with the filter formula, and the noise coefficient."""
import math, random, numpy as np, sys
from p5_online_newton_k1 import fisher_inc
from p10_fisher_k1_exact import exact_m

def Iev_matrix(mu, beta, eps):
    a = (1-eps)*beta; m = exact_m(mu, beta, eps); El = mu/eps; r = mu/eps
    I = np.array([[m, (1-mu*m)/a], [(1-mu*m)/a, (El-2*mu+mu**2*m)/a**2]]); return I/r

def run(pre, mu_s, beta, eps, eta, T, seed, A, P, grid=5.0):
    rng = random.Random(seed); al_s = (1-eps)*beta; mu, al = mu_s, al_s; r = mu_s/eps; mu_max = mu_s*(1+abs(A))
    Iev = Iev_matrix(mu_s, beta, eps); Pinv = np.linalg.inv(Iev)
    F = Iev/eta; s = 0.0; Tt = 0.0; ev = []; al_max = beta*(1-eps/2)
    while Tt < T:
        lam0 = mu_max + al_s*s; w = 0.0
        while True:
            w += rng.expovariate(lam0); mt = mu_s*(1 + A*math.sin(2*math.pi*(Tt+w)/P))
            if rng.random()*lam0 <= mt + al_s*s*math.exp(-beta*w): break
        e = math.exp(-beta*w); a1 = s*e; b1 = s*(1-e)/beta
        g = np.array([w - 1.0/(mu + al*a1), b1 - a1/(mu + al*a1)])
        if pre == 'adaptive':
            f0, f1, f2 = fisher_inc(mu, al, beta, s, w); F = (1-eta)*F + np.array([[f0, f1], [f1, f2]]); step = np.linalg.solve(F, g)
        else:
            step = eta*(Pinv@g)
        mu = max(mu - step[0], 0.1); al = min(max(al - step[1], 0.02), al_max)
        s = a1 + 1.0; Tt += w
        ev.append((Tt, mu))
    ev = np.array(ev); tg = np.arange(1, int(T/grid)+1)*grid
    idx = np.clip(np.searchsorted(ev[:, 0], tg, side='right') - 1, 0, len(ev)-1)       # last event at or before the grid time (causal)
    return ev[idx, 1]

if __name__ == "__main__":
    eps = float(sys.argv[1]); R = int(sys.argv[2]); A = float(sys.argv[3]); mu_s = 1.0; beta = 1.0; P = 1200.0; T = 6000.0
    r = mu_s/eps; w = 2*math.pi/(P*r); tgrid = np.arange(1, int(T/5.0)+1)*5.0; target = mu_s*(1 + A*np.sin(2*np.pi*tgrid/P)); bu = tgrid > 0.2*T
    print(f"eps={eps} A={A} omega={w:.5f}: alpha unknown; ratio = ensemble bias^2 / filter formula")
    for eta in (0.002, 0.004):
        pred = (A*A*mu_s**2/2)*w*w/(eta*eta + w*w)
        for pre in ('adaptive', 'fixed'):
            reps = np.array([run(pre, mu_s, beta, eps, eta, T, 600+s, A, P)[:len(tgrid)] for s in range(R)])
            mean = reps.mean(0); varc = reps.var(0, ddof=1).mean()/R; n = reps.shape[1]
            b2 = np.mean(((mean-target[:n])**2)[bu[:n]]) - varc; mse = np.mean(((reps-target[:n])**2)[:, bu[:n]])
            print(f"  eta={eta} {pre:8s}: bias^2 {b2:.5f} vs filter {pred:.5f}  ratio {b2/pred:.2f};  noise part/eta = {(mse-b2)/eta:.1f}", flush=True)
