"""Check the small-ball/truncation step of Step 3 (K=1):  rho_z m := E[1/lam; lam<z] <= (eps/a)(k/(k-1)) p_z^{1-1/k},  p_z <= e (eps z/(a+eps z))^{k_eps},  a=alpha^2/(2 beta),
z = c1/eps.  Simulation of the stationary law (time-weighted) at eps=0.1."""
import numpy as np, math, random, sys
from p10_fisher_k1_exact import exact_m
eps = float(sys.argv[1]); nev = int(sys.argv[2]); mu = 1.0; beta = 1.0; al = (1-eps)*beta; a = al**2/(2*beta); k = 2*beta*mu/al**2
m = exact_m(mu, beta, eps); rng = random.Random(3); S = 0.0; Ttot = 0.0
cs = [0.3, 0.6, 1.0, 2.0]; acc_p = np.zeros(len(cs)); acc_i = np.zeros(len(cs)); xg, wg = np.polynomial.legendre.leggauss(16)
for i in range(nev):
    lam0 = mu + al*S; w = 0.0
    while True:
        w += rng.expovariate(lam0)
        if rng.random()*lam0 <= mu + al*S*math.exp(-beta*w): break
    u = w*(xg+1)/2; wt = w*wg/2; lam = mu + al*S*np.exp(-beta*u)
    for j, c1 in enumerate(cs):
        ind = lam < c1/eps; acc_p[j] += np.sum(wt*ind); acc_i[j] += np.sum(wt*ind/lam)
    Ttot += w; S = S*math.exp(-beta*w) + 1.0
print(f"eps={eps}: k_eps={k:.3f}, a={a:.3f}, m={m:.5f}; stationary check: E[1/lam] over all time not recomputed")
print(f"{'c1':>5} {'p_z (sim)':>10} {'p_z bound':>10} {'rho_z (sim)':>12} {'rho_z bound':>12}")
for j, c1 in enumerate(cs):
    pz = acc_p[j]/Ttot; rho = acc_i[j]/Ttot/m; pzb = math.e*(c1/(a+c1))**k
    rb = (eps/(a*m))*(k/(k-1))*min(pzb, 1.0)**(1-1/k)
    print(f"{c1:5.2f} {pz:10.4f} {pzb:10.4f} {rho:12.4f} {rb:12.4f}")
