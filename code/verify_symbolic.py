"""Symbolic cross-check (sympy) of the algebraic identities behind the lemmas.  The generator of the Markov process S in R^K_+ is
      A f = - sum_k beta_k s_k df/ds_k + lam(s) [ f(s + 1) - f(s) ],   lam = mu + sum_k alpha_k s_k,
 built symbolically with EXACT derivatives; each claimed identity is compared with A f and the difference is evaluated at many random rational points
 (an identity between analytic functions that vanishes at many random points to 1e-30 is, for practical purposes, verified; where sympy simplify() succeeds we say so)."""
import sympy as sp, random, itertools
from sympy import Rational as R

random.seed(1)
def rr(a=1, b=9): return R(random.randint(a*10, b*10), 10)
res = []

def check(name, expr, symbols, ntry=12, lo=1, hi=9, extra=None):
    worst = 0
    for _ in range(ntry):
        sub = {s: rr(lo, hi) for s in symbols}
        if extra: sub = extra(sub)
        try: v = abs(sp.N(expr.subs(sub), 40))
        except Exception as e: v = sp.oo
        worst = max(worst, v)
    res.append((name, float(worst))); print(f"{'OK ' if worst < 1e-25 else 'FAIL'} {name:70s} max|diff| = {float(worst):.2e}", flush=True)

mu, b1, b2, a1, a2, s1, s2, th, z = sp.symbols('mu b1 b2 a1 a2 s1 s2 th z', positive=True)
S = [s1, s2]; beta = [b1, b2]; alpha = [a1, a2]
lam = mu + a1*s1 + a2*s2
def A(f):
    return -sum(beta[k]*S[k]*sp.diff(f, S[k]) for k in range(2)) + lam*(f.subs({s1: s1+1, s2: s2+1}, simultaneous=True) - f)

n = a1/b1 + a2/b2; eps = 1 - n
Zc = (a1/b1)*s1 + (a2/b2)*s2
check("Prop keyid:  A Z_c = mu - eps*lam", A(Zc) - (mu - eps*lam), [mu, b1, b2, a1, a2, s1, s2])
for p in (2, 3, 4):
    lhs = A(Zc**p)
    rhs = p*Zc**(p-1)*(mu - eps*lam) + lam*sum(sp.binomial(p, j)*n**j*Zc**(p-j) for j in range(2, p+1))
    check(f"Lemma 1': A Z^p = pZ^(p-1)(mu-eps lam) + lam sum C(p,j) n^j Z^(p-j), p={p}", lhs - rhs, [mu, b1, b2, a1, a2, s1, s2])
# exponential moment (Lemma 4): A e^{th Zc} = e^{th Zc}[th mu + lam(e^{th n} - 1 - th)]
check("Lemma 4: A e^{th Z_c} = e^{th Z_c}[th mu + lam (e^{th n} - 1 - th)]", A(sp.exp(th*Zc)) - sp.exp(th*Zc)*(th*mu + lam*(sp.exp(th*n) - 1 - th)), [mu, b1, b2, a1, a2, s1, s2, th], hi=3)
# third moments (used by the tangent bound)
for (i, j, l) in ((0, 0, 0), (0, 0, 1), (0, 1, 1), (1, 1, 1)):
    f = S[i]*S[j]*S[l]
    rhs = -(beta[i] + beta[j] + beta[l])*f + lam*(S[i]*S[j] + S[i]*S[l] + S[j]*S[l] + S[i] + S[j] + S[l] + 1)
    check(f"third moments: A(S_i S_j S_l) formula for ({i},{j},{l})", A(f) - rhs, [mu, b1, b2, a1, a2, s1, s2])
# second moments
for (i, j) in ((0, 0), (0, 1), (1, 1)):
    f = S[i]*S[j]; rhs = -(beta[i] + beta[j])*f + lam*(S[i] + S[j] + 1)
    check(f"second moments: A(S_i S_j) = -(b_i+b_j)S_iS_j + lam(S_i+S_j+1), ({i},{j})", A(f) - rhs, [mu, b1, b2, a1, a2, s1, s2])

# Perron objects: choose gamma, omega1; omega2 = 1 - omega1; alpha_k = (beta_k - gamma) omega_k
g, w1 = sp.symbols('g w1', positive=True); w2 = 1 - w1
sub_al = {a1: (b1 - g)*w1, a2: (b2 - g)*w2}
Zw = w1*s1 + w2*s2
def A_p(f): return A(f).subs(sub_al)
lam_p = lam.subs(sub_al)
check("Lemma G:  A Z_omega = mu - gamma Z_omega  (omega_k = alpha_k/(beta_k-gamma), sum omega = 1)", sp.simplify(A_p(Zw) - (mu - g*Zw)), [mu, b1, b2, s1, s2], extra=lambda d: {**d, g: rr(1, 5)/10, w1: rr(1, 9)/10})
# nu~, n, c0 and kappa
nut = w1/b1 + w2/b2; n_p = sub_al[a1]/b1 + sub_al[a2]/b2
check("Exact relation 1 = n + gamma*nu~  (n = sum alpha_k/beta_k)", 1 - (n_p + g*nut), [b1, b2], extra=lambda d: {**d, g: rr(1, 5)/10, w1: rr(1, 9)/10})
c0 = n_p/(2*nut)
kap = [(c0*w1 - sub_al[a1]/2)/b1, (c0*w2 - sub_al[a2]/2)/b2]
check("sum kappa_k = 0   (c0 = n/(2 nu~))", kap[0] + kap[1], [b1, b2], extra=lambda d: {**d, g: rr(1, 5)/10, w1: rr(1, 9)/10})
W = Zw + z
logW = sp.log(W)
lhs = A_p(logW); rhs = (mu - g*Zw)/W + lam_p*(sp.log(1 + 1/W) - 1/W)
check("Lemma H (i):  A log W = (mu-gamma Z)/W + lam[log(1+1/W) - 1/W]", lhs - rhs, [mu, b1, b2, s1, s2, z], extra=lambda d: {**d, g: rr(1, 5)/10, w1: rr(1, 9)/10})
psi = kap[0]*s1 + kap[1]*s2
lhs = A_p(psi/W**2)
rhs = -(kap[0]*b1*s1 + kap[1]*b2*s2)/W**2 + 2*psi*(w1*b1*s1 + w2*b2*s2 - lam_p)/W**3 + lam_p*psi*(3*W + 2)/(W**3*(W + 1)**2)
check("Lemma H (ii): A(psi/W^2) = -sum kappa beta S/W^2 + 2 psi (sum beta omega S - lam)/W^3 + lam psi (3W+2)/(W^3 (W+1)^2)", lhs - rhs, [mu, b1, b2, s1, s2, z], extra=lambda d: {**d, g: rr(1, 5)/10, w1: rr(1, 9)/10})
check("Lemma H (iii): sum beta_k omega_k S_k - lam = gamma Z - mu", (w1*b1*s1 + w2*b2*s2 - lam_p) - (g*Zw - mu), [mu, b1, b2, s1, s2], extra=lambda d: {**d, g: rr(1, 5)/10, w1: rr(1, 9)/10})
check("Lemma H (iv): lam/2 + sum kappa_k beta_k S_k = mu/2 + c0 Z", (lam_p/2 + kap[0]*b1*s1 + kap[1]*b2*s2) - (mu/2 + c0*Zw), [mu, b1, b2, s1, s2], extra=lambda d: {**d, g: rr(1, 5)/10, w1: rr(1, 9)/10})
check("elementary: 1/W^2 - 1/(W+1)^2 = 2/W^3 - (3W+2)/(W^3 (W+1)^2)", 1/z**2 - 1/(z+1)**2 - (2/z**3 - (3*z+2)/(z**3*(z+1)**2)), [z])

# K = 1 objects
s, mu1, be, al, e1, uu = sp.symbols('s mu1 be al e1 uu', positive=True)
lam1 = mu1 + al*s
def A1(f): return -be*s*sp.diff(f, s) + lam1*(f.subs(s, s+1) - f)
kap_th = be*th - al*(sp.exp(th) - 1)
check("Lemma E0:  A e^{th S} = e^{th S}[mu (e^th-1) - kappa_th S]", A1(sp.exp(th*s)) - sp.exp(th*s)*(mu1*(sp.exp(th) - 1) - kap_th*s), [mu1, be, al, s, th], hi=2)
# Laplace ODE:  A e^{-u lam} = e^{-u lam}[u beta (lam - mu) + lam (e^{-u alpha} - 1)]
lam_s = sp.Symbol('lam_s', positive=True)
Alam = lambda f: -be*(lam_s - mu1)*sp.diff(f, lam_s) + lam_s*(f.subs(lam_s, lam_s + al) - f)
check("Laplace ODE generator:  A e^{-u lam} = e^{-u lam}[u beta (lam-mu) + lam(e^{-u alpha}-1)]",
      Alam(sp.exp(-uu*lam_s)) - sp.exp(-uu*lam_s)*(uu*be*(lam_s - mu1) + lam_s*(sp.exp(-uu*al) - 1)), [mu1, be, al, lam_s, uu], hi=2)
# Gamma factorisation D_eps
x = sp.Symbol('x', positive=True); eps_ = sp.Symbol('eps_', positive=True)
alpha_ = be*(1 - eps_)
hfun = lambda y: sp.exp(-y) - 1 + y
c_ = eps_*be*x + hfun(alpha_*x)
D = (alpha_**2*x**2/2 - hfun(alpha_*x))/(c_*(eps_*be + alpha_**2*x/2))
check("Theorem A step 1:  s/c(s) = 1/(eps beta + alpha^2 s/2) + D_eps(s)", x/c_ - 1/(eps_*be + alpha_**2*x/2) - D, [be, x, eps_], hi=2, extra=lambda d: {**d, eps_: rr(1, 8)/10})
# Fisher reduction and efficient information
m_, El, mu_, al_, ep = sp.symbols('m El mu_ al_ ep', positive=True)
Imm = m_; Ima = (1 - mu_*m_)/al_; Iaa = (El - 2*mu_ + mu_**2*m_)/al_**2
check("Fisher reduction: det I = (m E lam - 1)/alpha^2", (Imm*Iaa - Ima**2) - (m_*El - 1)/al_**2, [m_, El, mu_, al_])
eff = (Imm*Iaa - Ima**2)/Iaa
check("Prop eff: det/I_aa = (m - eps/mu)/(1 - 2 eps + eps mu m) when E lam = mu/eps", eff.subs(El, mu_/ep) - (m_ - ep/mu_)/(1 - 2*ep + ep*mu_*m_), [m_, mu_, al_, ep], extra=lambda d: {**d, ep: rr(1, 8)/10})
# ridge Rayleigh quotient: E[(1+eps - eps lam/mu)^2/lam] = (1+eps)^2 m - (eps/mu)(1+2 eps) when E lam = mu/eps
lamv = sp.Symbol('lamv', positive=True)
expr = sp.expand((1 + ep - ep*lamv/mu_)**2/lamv)   # linear in 1/lam, constant, lam
Ex = (1 + ep)**2*m_ + (-2*(1 + ep)*ep/mu_) + (ep/mu_)**2*(mu_/ep)
check("Prop ridge: v'Iv = (1+eps)^2 m - (eps/mu)(1+2 eps)", Ex - ((1 + ep)**2*m_ - (ep/mu_)*(1 + 2*ep)), [m_, mu_, ep], extra=lambda d: {**d, ep: rr(1, 8)/10})
# Lemma 1 (curvature): Hessian of -log(a'theta) + b'theta
t0, t1, a0, a1_, bb0, bb1 = sp.symbols('t0 t1 a0 a1_ bb0 bb1', positive=True)
lth = a0*t0 + a1_*t1; ell = -sp.log(lth) + bb0*t0 + bb1*t1
Hs = sp.hessian(ell, (t0, t1)); Ht = sp.Matrix([[a0*a0, a0*a1_], [a1_*a0, a1_*a1_]])/lth**2
check("Lemma curvature:  Hessian of -log(a'theta) + b'theta = a a'/lam^2", (Hs - Ht).norm(), [t0, t1, a0, a1_, bb0, bb1])
# Newton / martingale scalar identity (m_i simplification)
P_, av, ls_, Xint = sp.symbols('P av ls Xint', positive=True)
lt_ = ls_ + av        # lam_theta(tau) = lam*(tau) + a'v
Lint = P_ - Xint      # int (phi'v) lam*/lam_theta = int phi'v (1 - phi'v/lam_theta) = P - X
lhs = 2*(P_ - av/ls_) + (2*av**2/(ls_*lt_) - 2*Xint)
rhs = -2*av/lt_ + 2*Lint
check("Lemma C1: 2<s,v> + A = -2 a'v/lam_theta(tau) + 2 int(phi'v)lam*/lam_theta", lhs - rhs, [P_, av, ls_, Xint])
# Newton step identity
G_, Fi, v_ = sp.symbols('G_ Fi v_', positive=True)
check("Newton identity (scalar): (v - g/F)^2 F = v^2 F - 2 g v + g^2/F", (v_ - G_/Fi)**2*Fi - (v_**2*Fi - 2*G_*v_ + G_**2/Fi), [G_, Fi, v_])
# tanh jump law
q_, D_, kL, tt = sp.symbols('q_ D_ kL tt', positive=True)
e0 = D_/(1 + sp.exp(-kL))
avg = e0**2*(1 - sp.exp(-2*kL))/(2*kL)
check("Prop tanh:  e0^2 (1-e^{-2kL})/(2kL) = D^2 tanh(kL/2)/(2kL),  e0 = D/(1+e^{-kL})", avg - D_**2*sp.tanh(kL/2)/(2*kL), [D_, kL], hi=5)
check("Prop tanh (periodicity): e0 = D - e0 e^{-kL}", e0 - (D_ - e0*sp.exp(-kL)), [D_, kL], hi=5)
# sinusoid
kk, om = sp.symbols('kk om', positive=True)
H = kk/(kk + sp.I*om); errtf = 1 - H
check("Prop sin:  |1 - kappa/(kappa+i omega)|^2 = omega^2/(kappa^2+omega^2)", sp.Abs(errtf)**2 - om**2/(kk**2 + om**2), [kk, om])
# time-change: tail of Exp(1) used in Lemma E5a (trivial) - skip
# Theorem T / elementary: N(y), E1(y) real/imag parts and the algebra d(y) = N / (E1 (r - i y/2))
y, r_ = sp.symbols('y r_', positive=True)
Hh = sp.exp(sp.I*y) - 1 - sp.I*y        # h(-i y) with h(x) = e^{-x} - 1 + x
Nn = (1 - sp.cos(y) - y**2/2) + sp.I*(y - sp.sin(y))
check("Thm T: -(h(y_) - y_^2/2) with y_ = -i y  equals N(y)", -(Hh - (-y**2/2)) - Nn, [y])
E1 = -(1 - sp.cos(y)) - sp.I*((1 + r_)*y - sp.sin(y))
chk = -sp.I*r_*y + Hh        # c(s) at s=-i tau, y = alpha tau: eps beta s = r*alpha*s = -i r y
check("Thm T: c(-i tau) = -i r y + h(-i y) equals E1(y)", chk - E1, [y, r_])

bad = [r for r in res if r[1] >= 1e-25]
print(f"\n{len(res)} identities checked, {len(res) - len(bad)} verified to < 1e-25, failures: {bad}")
