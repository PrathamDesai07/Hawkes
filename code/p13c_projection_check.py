import numpy as np, math, sys
from p13c_engine import run
from p14_mech import Iev_matrix
from p2_identifiability import mc
cfg = sys.argv[1]; eta = float(sys.argv[2]); A = float(sys.argv[3]); P = float(sys.argv[4]); n = int(sys.argv[5]); R = int(sys.argv[6]); eps = float(sys.argv[7])
if cfg == 'K1': beta = [1.0]; mu_s = 1.0; Iev = Iev_matrix(mu_s, 1.0, eps)
else:
    beta = [1.0, 2.0]; mu_s = 1.5; I, _ = mc(mu_s, np.array(beta), eps, 250000, seed=4); Iev = I/(mu_s/eps)
fn = lambda t: mu_s*(1 + A*np.sin(2*np.pi*t/P)); mu_max = mu_s*(1+abs(A))
out = [run(beta, eps, fn, mu_max, eta, n, 900+s, Iev, mu_s) for s in range(R)]
cat = lambda k: np.concatenate([o[k] for o in out]); m_ = lambda k: cat(k).mean()
U, N, sF, lev, dF2, q, ne, pr, lp = (cat(k) for k in ('U','N','sF','lev','dF2','q','ne','proj','lp'))
print(f"{cfg} eps={eps} eta={eta} A={A}: {len(U)} events; projected steps {pr.mean()*100:.3f}% (exact F-norm projection)")
print(f" (a) non-expansiveness about the moving truth: max_i [ ||Pi(y)-theta*(t_i)||_F - ||y-theta*(t_i)||_F ] = {ne.max():.2e}  (must be <= 0; steps where the projection is active: mean gap {ne[pr>0].mean() if (pr>0).any() else 0:.3e})")
rhs = (1-eta)*m_('U') + m_('sF')/(1-lev.max())
print(f" (b) Lemma N at the actual iterate: mean N = {m_('N'):.5f}  <=  (1-eta) mean U + mean s'F^-1 s/(1-lev_max) = {rhs:.5f}   (diff {m_('N')-rhs:+.5f}; lev_max = {lev.max():.3f})")
nu = m_('sF')/(1-lev.max()); delta = (1 + lev.max()/0.6)*math.sqrt(m_('dF2'))
Vb = (delta + math.sqrt((1-eta)*delta**2 + eta*nu))/eta
print(f" (c) combined steady-state bound V <= (delta + sqrt((1-eta)delta^2 + eta nu))/eta = {Vb:.4f}  vs measured sqrt(mean U) = {math.sqrt(m_('U')):.4f}   (ratio {math.sqrt(m_('U'))/Vb:.3f}; noise nu={nu:.5f}, delta={delta:.5f})")
print(f" mean U = {m_('U'):.4f}; Fisher-metric error mean v'Iev v not recorded here; c_F median {np.median(cat('cF')):.3f}, min {cat('cF').min():.3f}")
