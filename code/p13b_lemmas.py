import numpy as np, math, sys
from p13b_engine import run
from p14_mech import Iev_matrix
from p2_identifiability import mc
cfg = sys.argv[1]; eta = float(sys.argv[2]); A = float(sys.argv[3]); P = float(sys.argv[4]); n = int(sys.argv[5]); R = int(sys.argv[6])
if cfg == 'K1': beta = [1.0]; mu_s = 1.0; eps = float(sys.argv[7]); Iev = Iev_matrix(mu_s, 1.0, eps)
else:
    beta = [1.0, 2.0]; mu_s = 1.5; eps = float(sys.argv[7]); I, _ = mc(mu_s, np.array(beta), eps, 250000, seed=4); Iev = I/(mu_s/eps)
fn = lambda t: mu_s*(1 + A*np.sin(2*np.pi*t/P)); mu_max = mu_s*(1+abs(A))
out = [run(beta, eps, fn, mu_max, eta, n, 300+s, Iev, mu_s) for s in range(R)]
cat = lambda k: np.concatenate([o[k] for o in out]); m_ = lambda k: cat(k).mean(); se = lambda a_, b_: np.std(cat(a_)-cat(b_))/math.sqrt(len(cat(a_)))
print(f"{cfg} eps={eps} eta={eta} A={A} P={P}: {len(cat('U'))} events; continuous drift inside gaps")
sc = cat('sumcheck'); print(f" exact split v = v^xi + v^d: clipped steps (all runs) {sum(o['nclip'] for o in out)}; max |v - v^xi - v^d| = {sc.max():.2e}, median {np.median(sc):.2e}")
for nm, l, r in (('v^d', 'lhs_d', 'rhs_d'), ('v^xi', 'lhs_x', 'rhs_x'), ('random predictable w', 'lhs_r', 'rhs_r')):
    print(f" contraction lemma, w = {nm}: mean E||(I-G)w||^2_F = {m_(l):.5f} <= mean[(1-eta)||w||^2 - (1-lev) w'dF w] = {m_(r):.5f}   (diff {m_(l)-m_(r):+.5f} +- {se(l, r):.5f})")
lev = cat('lev'); print(f" point leverage a^T F^-1 a/(lam* lam_theta): mean {lev.mean():.4f}, 99.9% {np.quantile(lev,0.999):.4f}, max {lev.max():.4f}")
U, Uxi, Ud = m_('U'), m_('Uxi'), m_('Ud')
print(f" Minkowski: sqrt(mean U) = {math.sqrt(U):.4f} <= sqrt(mean U^xi) + sqrt(mean U^d) = {math.sqrt(Uxi)+math.sqrt(Ud):.4f}   (U^xi {Uxi:.4f}, U^d {Ud:.4f})")
q = m_('q'); d_ = len(beta)+1; cF = cat('cF')
print(f" noise part: mean U^xi = {Uxi:.4f}  vs  (1/eta)*mean q/(1-eta)*(1+...) = {q/(eta*(1-eta)):.4f}   (q = {q:.5f}, d*eta = {d_*eta:.5f})")
print(f" Fisher-metric errors: mean v'Iv = {m_('vI'):.5f} (noise part {m_('vxiI'):.5f}, drift part {m_('vdI'):.5f}); c_F median {np.median(cF):.3f}, min {cF.min():.3f}; eta*d/2 = {eta*d_/2:.5f}")
print(f" drift forcing: mean ||d||^2_F = {m_('dF2'):.3e}, mean ||d||^2_I = {m_('dI'):.3e}, ratio eta*||d||_F^2/||d||_I^2 = {eta*m_('dF2')/m_('dI'):.3f}")
