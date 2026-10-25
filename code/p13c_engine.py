"""Problem 13, third pass: engine with the EXACT F-norm projection onto Theta_c = {lo<=x<=hi, c'x<=1-eps_c} (Problem 4's project_dense), continuous drift in mu.
Records, per event: U_i=||v_i||^2_{F_{i-1}}, the one-step quantities of Lemma N evaluated at the ACTUAL projected iterate
 N_i = ||(I-G_i)v_i - F_i^-1 s_i||^2_{F_i}   (should satisfy E N_i <= (1-eta) E U_i + E s'F_i^-1 s_i/(1-lev) ),
non-expansiveness gap  ne_i = ||Pi(y_i)-theta*(t_i)||_{F_i} - ||y_i - theta*(t_i)||_{F_i}  (<= 0 expected),  q_i, dF2_i=||d_i||^2_{F_i}, leverage, and PREDICTABLE leverage
 lp_i = sup_u phi_i(u)'F_{i-1}^-1 phi_i(u) / ((1-eta) c- lam*_i(u)^2)   (depends only on the past: S_{i-1}, F_{i-1})."""
import numpy as np, math, random
from p13b_engine import quad_info
from p4_dense_projection import project_dense

def run(beta, eps, mu_fn, mu_max, eta, n, seed, Iev, mu_ref, mu_hi=4.0, w_share=None, burn=0.25, cminus=0.6, project=True):
    rng = random.Random(seed); beta = np.asarray(beta, float); K = len(beta); d = K+1
    ws = np.ones(K)/K if w_share is None else np.asarray(w_share); al = (1-eps)*ws*beta
    th = np.concatenate([[mu_ref], al]); F = Iev/eta; S = np.zeros(K); Tt = 0.0
    lo = np.concatenate([[0.1], np.full(K, 0.01)]); hi = np.concatenate([[mu_hi], beta*0.98]); cvec = np.concatenate([[0.0], 1/beta]); rbud = 1 - eps/2
    rec = {k: [] for k in ('U','N','sF','lev','dF2','q','ne','proj','lp','lpgap','cF')}
    Ih = np.linalg.cholesky(np.linalg.inv(Iev)); ugrid = np.concatenate([[0.0], np.logspace(-2, 1.5, 40)])
    for i in range(n):
        mu0 = mu_fn(Tt); thstar0 = np.concatenate([[mu0], al])
        lam0 = mu_max + al@S; w = 0.0
        while True:
            w += rng.expovariate(lam0)
            if rng.random()*lam0 <= mu_fn(Tt+w) + al@(S*np.exp(-beta*w)): break
        e = np.exp(-beta*w); a = np.concatenate([[1.0], S*e]); b = np.concatenate([[w], S*(1-e)/beta])
        thstar1 = np.concatenate([[mu_fn(Tt+w)], al]); dvec = thstar1 - thstar0
        lam_s = mu_fn(Tt+w) + al@(S*e); lam_t = a@th
        g = b - a/lam_t; sv = b - a/lam_s; A = np.outer(a, a)/(lam_s*lam_t)
        dF = quad_info(None, None, beta, S, w, Tt, est_theta=th); dJ = quad_info(mu_fn, al, beta, S, w, Tt)
        v = th - thstar0; Fp = F; F = (1-eta)*Fp + dF; Fi = np.linalg.inv(F); Fpi = np.linalg.inv(Fp); G = Fi@A; IG = np.eye(d) - G
        y = th - Fi@g
        thn = project_dense(F, y, lo, hi, cvec, rbud) if project else y
        if i > burn*n:
            ahat = a/math.sqrt(lam_s*lam_t); lev = ahat@Fi@ahat
            vN = IG@v - Fi@sv
            # predictable leverage: sup over the gap vectors phi(u) (depends on S_{i-1}, F_{i-1}, theta* at gap start only)
            Sg = S[None, :]*np.exp(-np.outer(ugrid, beta)); phis = np.hstack([np.ones((len(ugrid), 1)), Sg]); lam_u = mu0 + Sg@al
            lp = np.max(np.einsum('ij,jk,ik->i', phis, Fpi, phis)/((1-eta)*cminus*lam_u**2))
            lpg = lev
            for k, val in (('U', v@Fp@v), ('N', vN@F@vN), ('sF', sv@Fi@sv), ('lev', lev), ('dF2', dvec@F@dvec), ('q', np.trace(Fpi@dJ)),
                           ('ne', math.sqrt(max((thn-thstar1)@F@(thn-thstar1), 0)) - math.sqrt(max((y-thstar1)@F@(y-thstar1), 0))),
                           ('proj', float(np.max(np.abs(thn-y)) > 1e-12)), ('lp', lp), ('lpgap', lpg),
                           ('cF', np.min(np.linalg.eigvalsh(Ih.T@(eta*Fp)@Ih)))):
                rec[k].append(val)
        th = thn; S = S*e + 1.0; Tt += w
    return {k: np.array(v) for k, v in rec.items()}
