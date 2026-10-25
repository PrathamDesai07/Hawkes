"""Problem 13, second pass: general-K engine with CONTINUOUS drift in mu (inside each gap), recording everything the new lemmas need.
 exact recursion:  v_{i+1} = (I-G_i) v_i - (I-G_i) d_i - F_i^{-1} s_i,   G_i = F_i^{-1} A_i,  A_i = a a^T/(lam*(w) lam_theta(w)),
 v_i = theta_i - theta*(t_{i-1}),  d_i = theta*(t_i) - theta*(t_{i-1}),  s_i = b_i - a_i/lam*(w_i)  (true time-varying lam*).
 split with the SAME realised matrices:  v^xi_{i+1} = (I-G_i)v^xi_i - F_i^{-1}s_i  (v^xi_1 = v_1),   v^d_{i+1} = (I-G_i)v^d_i - (I-G_i)d_i  (v^d_1 = 0),  v = v^xi + v^d exactly."""
import numpy as np, math, random
xg, wg = np.polynomial.legendre.leggauss(8); xg = 0.5*(xg+1); wg = 0.5*wg

def quad_info(th_mu_fn, th_alpha, beta, S, w, t0, nsub=3, est_theta=None):
    """int_0^w phi phi^T / lam du with lam(u) = est_theta^T phi (est) or the truth mu_fn(t0+u)+alpha.S e^{-beta u} (truth)."""
    K = len(beta); F = np.zeros((K+1, K+1)); h = w/nsub
    for j in range(nsub):
        u = h*(j + xg); Su = S[None, :]*np.exp(-np.outer(u, beta)); phi = np.hstack([np.ones((len(u), 1)), Su])
        lam = phi@est_theta if est_theta is not None else th_mu_fn(t0+u) + Su@th_alpha
        F += (h*wg/lam*phi.T)@phi
    return F

def run(beta, eps, mu_fn, mu_max, eta, n, seed, Iev, mu_ref, w_share=None, burn=0.25, jump_times=None):
    rng = random.Random(seed); beta = np.asarray(beta, float); K = len(beta); d = K+1
    ws = np.ones(K)/K if w_share is None else np.asarray(w_share); al = (1-eps)*ws*beta
    th = np.concatenate([[mu_ref], al]); F = Iev/eta; S = np.zeros(K); Tt = 0.0
    Ih = np.linalg.cholesky(np.linalg.inv(Iev))       # Iev^{-1} = Ih Ih^T  (whitening: Ih^T X Ih)
    vxi = np.zeros(d); vd = np.zeros(d); first = True
    rec = {k: [] for k in ('U','Uxi','Ud','vI','vdI','vxiI','q','x','lev','cF','dF2','dI','lhs_d','rhs_d','lhs_x','rhs_x','lhs_r','rhs_r','sF','sumcheck','clip')}
    rng_w = np.random.default_rng(seed+7); nclip = [0]
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
        if first: vxi = v.copy(); vd = np.zeros(d); first = False
        U = v@Fp@v; Uxi = vxi@Fp@vxi; Ud = vd@Fp@vd
        ahat = a/math.sqrt(lam_s*lam_t); lev = ahat@Fi@ahat
        if i > burn*n:
            wr = rng_w.standard_normal(d)*np.sqrt(np.diag(np.linalg.inv(Fp)))     # a predictable random direction (scaled to F-norm ~ 1)
            for k, val in (('U', U), ('Uxi', Uxi), ('Ud', Ud), ('vI', v@Iev@v), ('vdI', vd@Iev@vd), ('vxiI', vxi@Iev@vxi), ('q', np.trace(Fpi@dJ)),
                           ('x', np.max(np.linalg.eigvals(np.linalg.inv((1-eta)*Fp)@dF).real)), ('lev', lev),
                           ('cF', np.min(np.linalg.eigvalsh(Ih.T@(eta*Fp)@Ih))), ('dF2', dvec@F@dvec), ('dI', dvec@Iev@dvec),
                           ('lhs_d', np.sum((IG@vd)*(F@(IG@vd)))), ('rhs_d', (1-eta)*Ud - (1-lev)*(vd@dF@vd)),
                           ('lhs_x', np.sum((IG@vxi)*(F@(IG@vxi)))), ('rhs_x', (1-eta)*Uxi - (1-lev)*(vxi@dF@vxi)),
                           ('lhs_r', np.sum((IG@wr)*(F@(IG@wr)))), ('rhs_r', (1-eta)*(wr@Fp@wr) - (1-lev)*(wr@dF@wr)),
                           ('sF', sv@Fi@sv), ('sumcheck', np.max(np.abs(v - vxi - vd))), ('clip', 0.0)):
                rec[k].append(val)
        # updates
        vxi_new = IG@vxi - Fi@sv; vd_new = IG@vd - IG@dvec
        th_new = th - Fi@g; vxi, vd = vxi_new, vd_new; th_raw = th_new.copy()
        th_new[0] = max(th_new[0], 0.1); th_new[1:] = np.maximum(th_new[1:], 0.01); nh = (th_new[1:]/beta).sum()
        if nh > 1-eps/2: th_new[1:] *= (1-eps/2)/nh
        clipped = bool(np.max(np.abs(th_new - th_raw)) > 0)
        if clipped:                                    # resynchronise the diagnostic split after a projection step (absorbed into v^xi)
            vxi = (th_new - thstar1) - vd; nclip[0] += 1
        if i > burn*n: rec['clip'][-1] = float(clipped) if rec['clip'] else 0.0
        th = th_new; S = S*e + 1.0; Tt += w
    out = {k: np.array(v) for k, v in rec.items()}; out['nclip'] = nclip[0]; return out
