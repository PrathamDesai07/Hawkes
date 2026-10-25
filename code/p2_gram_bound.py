"""Problem 2, Gram-matrix route. Rigorous bound (see problem02_gram_bound.md):
  lam_min(I) >= p * exp(-L h) * min(1, s_lo^2) * lam_min(G_h) / (L h),
  p >= 1 - exp(-mu t0) - mu/(eps L),  s_lo = exp(-beta_K t0),
  G_h = Gram of {1, e^{-beta_k u}} on [0,h].   Free parameters t0, L, h are optimised."""
import numpy as np

def gram(beta, h):
    g = np.concatenate([[0.0], beta]); s = np.add.outer(g, g)
    return h*np.where(s*h < 1e-12, 1.0, -np.expm1(-s*h)/np.where(s == 0, 1, s*h))

def lmin_gram(beta, h): return np.linalg.eigvalsh(gram(beta, h))[0]

def bound(mu, beta, eps):
    best = (0, None)
    for q in (0.05, 0.1, 0.2, 0.3, 0.4):                    # q = mu/(eps L)  (Markov)
        L = mu/(eps*q)
        for t0 in np.linspace(0.1, 6, 40)/mu:
            p = 1 - np.exp(-mu*t0) - q
            if p <= 0: continue
            slo = np.exp(-beta[-1]*t0)
            for c in (0.01, 0.03, 0.1, 0.3, 1.0, 3.0):      # h = c/L
                h = c/L
                val = p*np.exp(-c)*min(1, slo**2)*lmin_gram(beta, h)/(L*h)
                if val > best[0]: best = (val, dict(q=q, t0=round(t0, 2), c=c, p=round(p, 3), slo=slo))
    return best

if __name__ == "__main__":
    print("scaling of lam_min(G_h) with h (beta = 1,2,4,..): log-log slope between h=1e-3 and 1e-2; expect 2K+1")
    for K, r in ((1, 2.0), (2, 2.0), (3, 2.0)):
        beta = r**np.arange(K); h1, h2 = 1e-3, 1e-2
        sl = np.log(lmin_gram(beta, h2)/lmin_gram(beta, h1))/np.log(h2/h1)
        print(f"  K={K}: slope {sl:.2f} (2K+1 = {2*K+1})")
    mc = {(1, 2, .3): .1389, (1, 2, .1): .0654, (1, 2, .05): .0383, (2, 2, .3): .0337, (2, 2, .1): .0455, (2, 2, .05): .0435,
          (2, 4, .3): .0522, (2, 4, .1): .0677, (2, 4, .05): .0635, (3, 2, .3): .0033, (3, 2, .1): .0045, (3, 2, .05): .0051}
    print(f"\n{'K':>2} {'r':>3} {'eps':>5} {'MC lmin(I)':>10} {'Gram bound':>11} {'bound/MC':>9}  best params")
    for (K, r, eps), m in mc.items():
        beta = float(r)**np.arange(K); v, par = bound(1.0, beta, eps)
        print(f"{K:2d} {r:3d} {eps:5.2f} {m:10.4f} {v:11.3e} {v/m:9.2e}  {par}")
