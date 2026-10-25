"""Problem 13: step-by-step numerical check of the tracking chain (K=1, alpha unknown, step-hold truth).
 truth in gap i: theta*_i = (mu_i, alpha*) with mu_i = mu_s (1+A sin(2 pi t_{i-1}/P))  (constant within the gap, predictable)
 recursion: F_i = (1-eta)F_{i-1} + dF_i(theta_i);  theta_{i+1} = theta_i - F_i^{-1} g_i  (clipped as in Problem 5; clip fraction reported)
 recorded per event: U_i=v'F_{i-1}v, <g,v>, v'dF v, |s|^2_{F_i^-1}, q_i=tr(F_{i-1}^-1 dJ_i), |h|^2_{F_i^-1}, xbar_i = lam_max(((1-eta)F_{i-1})^-1 dF_i), ||d_i||^2_{F_i}, c_F sample."""
import math, random, numpy as np, sys
from p5_online_newton_k1 import fisher_inc
from p14_mech import Iev_matrix

def run(eps, eta, n, seed, A, P, mu_s=1.0, beta=1.0, burn=0.25):
    rng = random.Random(seed); al_s = (1-eps)*beta; r = mu_s/eps; Iev = Iev_matrix(mu_s, beta, eps); Ih = np.linalg.cholesky(np.linalg.inv(Iev))   # Iev^{-1} = Ih Ih^T
    th = np.array([mu_s, al_s]); F = Iev/eta; s = 0.0; Tt = 0.0
    rec = {k: [] for k in ('U','gv','vdFv','s2','q','h2','x','dF2','cF','dd','clip','Ud','trbound','mseI')}
    mu_prev = mu_s*(1 + A*math.sin(2*math.pi*Tt/P)); lam_logdet0 = np.linalg.slogdet(F)[1]; tr_sum = 0.0
    for i in range(n):
        mu_i = mu_s*(1 + A*math.sin(2*math.pi*Tt/P)); th_s = np.array([mu_i, al_s])
        lam0 = mu_i + al_s*s; w = 0.0
        while True:
            w += rng.expovariate(lam0)
            if rng.random()*lam0 <= mu_i + al_s*s*math.exp(-beta*w): break
        e = math.exp(-beta*w); a1 = s*e; b1 = s*(1-e)/beta
        f0, f1, f2 = fisher_inc(th[0], th[1], beta, s, w); dF = np.array([[f0, f1], [f1, f2]])
        j0, j1, j2 = fisher_inc(mu_i, al_s, beta, s, w); dJ = np.array([[j0, j1], [j1, j2]])
        lam_t = th[0] + th[1]*a1; lam_s = mu_i + al_s*a1
        g = np.array([w - 1/lam_t, b1 - a1/lam_t]); sv = np.array([w - 1/lam_s, b1 - a1/lam_s]); h = g - sv
        v = th - th_s; Fprev = F; F = (1-eta)*Fprev + dF; Fi = np.linalg.inv(F); Fpi = np.linalg.inv(Fprev)
        U = v@Fprev@v; gv = g@v; vdFv = v@dF@v
        mu_next = mu_s*(1 + A*math.sin(2*math.pi*(Tt+w)/P)); d = np.array([mu_next - mu_i, 0.0])
        xb = np.max(np.linalg.eigvals(np.linalg.inv((1-eta)*Fprev)@dF).real)
        step = Fi@g; th_new = th - step; clipped = (th_new[0] < 0.1) or (th_new[1] < 0.02) or (th_new[1] > beta*(1-eps/2))
        th_new[0] = max(th_new[0], 0.1); th_new[1] = min(max(th_new[1], 0.02), beta*(1-eps/2))
        if i > burn*n:
            for k, val in (('U', U), ('gv', gv), ('vdFv', vdFv), ('s2', sv@Fi@sv), ('q', np.trace(Fpi@dJ)), ('h2', h@Fi@h), ('x', xb), ('dF2', d@F@d),
                           ('cF', np.min(np.linalg.eigvalsh(Ih.T@(eta*Fprev)@Ih))), ('dd', d@Iev@d), ('clip', float(clipped)), ('trbound', np.trace(Fpi@dF)), ('mseI', v@Iev@v)):
                rec[k].append(val)
        th = th_new; s = a1 + 1.0; Tt += w
    return {k: np.array(v) for k, v in rec.items()}, np.linalg.slogdet(F)[1] - lam_logdet0

if __name__ == "__main__":
    eps = float(sys.argv[1]); eta = float(sys.argv[2]); A = float(sys.argv[3]); P = float(sys.argv[4]); n = int(sys.argv[5]); R = int(sys.argv[6])
    out = []
    for s in range(R): out.append(run(eps, eta, n, 100+s, A, P))
    cat = lambda k: np.concatenate([o[0][k] for o in out])
    d = 2; r = 1/eps
    gv, vdFv = cat('gv'), cat('vdFv'); se = np.std(gv - vdFv)/math.sqrt(len(gv))
    print(f"eps={eps} eta={eta} A={A} P={P}, n={n} x {R} runs (post-burn events {len(gv)}), clip fraction {cat('clip').mean():.4f}")
    print(f" C2 expectation identity: mean<g,v> = {gv.mean():.5f}, mean v'dF v = {vdFv.mean():.5f}, difference {np.mean(gv-vdFv):+.5f} +- {se:.5f}")
    U = cat('U'); q = cat('q'); dF2 = cat('dF2'); s2 = cat('s2'); h2 = cat('h2'); xb = cat('x'); cF = cat('cF'); dd = cat('dd'); trb = cat('trbound')
    zeta = eta/2; rho = 1.0; C_q = (1+zeta)*(1+rho)/(1-eta); C_d = 1+1/zeta
    bound_U = (2/eta)*(C_q*q.mean() + C_d*dF2.mean())
    print(f" steps 1-5: mean U = {U.mean():.4f}  <=  (2/eta)[C_q mean q + C_d mean ||d||^2_F] = {bound_U:.4f}  (ratio {U.mean()/bound_U:.3f}); noise part {(2/eta)*C_q*q.mean():.4f}, drift part {(2/eta)*C_d*dF2.mean():.4f}")
    print(f" E|s|^2_{{F_i^-1}} = {s2.mean():.5f} <= E q_i/(1-eta) = {q.mean()/(1-eta):.5f};  E|h|^2_{{F_i^-1}} / E v'dF v = {h2.mean()/vdFv.mean():.4f}  (assumption A1: B_h)")
    xbar = np.quantile(xb, 0.999); print(f" gap leverage: xbar (99.9% quantile) = {xbar:.4f}, max {xb.max():.4f}, mean q_i = {q.mean():.5f} vs d*(-ln(1-eta)) = {d*(-math.log(1-eta)):.5f}")
    lg = np.mean([o[1] for o in out]); kx = (1-eta)*xb.max()/math.log(1+xb.max()); print(f" discount-trace lemma: mean tr(F_(i-1)^-1 dF_i) = {trb.mean():.5f} <= (1-eta) xbar/log(1+xbar) [d(-ln(1-eta)) + (1/n) dlogdet] ~ {kx*d*(-math.log(1-eta)):.5f} (+ transient {kx*lg/ (n*(1-0.25)):.5f})")
    print(f" A4: c_F = lambda_min(I^-1/2 (eta F) I^-1/2): median {np.median(cF):.3f}, 5% {np.quantile(cF,0.05):.3f}, 1% {np.quantile(cF,0.01):.3f}, min {cF.min():.3f}")
    print(f" drift: mean eta*||d||^2_F / ||d||^2_I-ratio = {np.mean(eta*dF2)/max(np.mean(dd),1e-30):.3f} (the theorem uses c_I/c^- ), mean ||d||^2_I = {dd.mean():.3e}")
    mse = cat('mseI'); bound_I = (eta/np.maximum(cF.min(),1e-9))*bound_U
    print(f" Fisher-metric MSE mean v'I_ev v = {mse.mean():.5f};  bound (eta/c_F,min) * bound_U = {bound_I:.4f} (ratio {mse.mean()/bound_I:.4f});  using the median c_F: {(eta/np.median(cF))*bound_U:.4f}")
    print(f" noise-only prediction eta*d/2 = {eta*d/2:.5f} (efficient estimator); measured/that = {mse.mean()/(eta*d/2):.2f}")
