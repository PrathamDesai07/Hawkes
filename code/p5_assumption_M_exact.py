"""Assumption M without a Gamma comparator (K=1, any k0, any eps):
   E_pi h = (1/pi) int_0^inf H(t) Re[ e^{i t c} conj(phi(t)) ] dt,    phi(t) = E e^{i t X} = exp(L),  X = eps*lambda,
   H(t) = (4/rho) sin(A t) sin(rho t/2)/t^2  (Fourier transform of the trapezoid h centred at c: A=(s+p)/2, rho=s-p, plateau half p, support half s).
 phi is exact (Lemma R):  dL/dx = (k/2) x/E1(x),  x = alpha*eps*t,  E1 = -(1-cos x) - i((1+r)x - sin x),  r = eps/(1-eps),  k = 2 beta mu/alpha^2.
 Quadrature: composite Gauss-Legendre panels in t, with L at the nodes from the panel-wise spectral integration matrix.
 Truncation: |phi| is non-increasing (Lemma elem(c)) and |H| <= 4/(rho t^2), so |tail beyond T| <= |phi(T)| * 4/(pi rho T)   (explicit).
 Convergence is checked by refining panels and nodes."""
import numpy as np, math, sys

def gl_matrices(n):
    x, w = np.polynomial.legendre.leggauss(n); x = 0.5*(x+1); w = 0.5*w          # on [0,1]
    V = np.polynomial.legendre.legvander(2*x-1, n-1)                             # n x n Vandermonde in Legendre basis
    # integral_0^{x_j} P_m(2s-1) ds  for Legendre basis
    C = np.zeros((n, n))
    for m in range(n):
        c = np.zeros(m+1); c[m] = 1
        ic = np.polynomial.legendre.legint(c, lbnd=-1)                           # antiderivative from -1
        C[:, m] = 0.5*np.polynomial.legendre.legval(2*x-1, ic)
    S = C @ np.linalg.inv(V)                                                     # maps f(nodes) -> integral_0^{x_j} f
    return x, w, S

def f_L(x, r, k):
    x = np.asarray(x, float); out = np.empty(x.shape, complex); sm = x < 0.2
    y = x[sm]
    c1 = y**2/2 - y**4/24 + y**6/720 - y**8/40320 + y**10/3628800
    im = y**3/6 - y**5/120 + y**7/5040 - y**9/362880
    out[sm] = (k/2)*y/(-c1 - 1j*(r*y + im))
    y = x[~sm]; out[~sm] = (k/2)*y/(-(1-np.cos(y)) - 1j*((1+r)*y - np.sin(y)))
    return out

def E_h(mu, eps, cS, panel=0.4, n=16, xmax_t=None, beta=1.0):
    alpha = (1-eps)*beta; r = eps/(1-eps); k = 2*beta*mu/alpha**2
    c = mu*(eps + (1-eps)*cS); dX = 0.15*alpha*mu/beta; p = dX/2; s = dX; A = (s+p)/2; rho = s-p
    # choose T so that the explicit tail bound is < 1e-9 (cap at Tmax)
    xg, wg, Smat = gl_matrices(n)
    Tmax = 3e5 if xmax_t is None else xmax_t
    tcur = 0.0; Lc = 0.0 + 0j; J = 0.0; step = panel
    chunk = 4000
    while tcur < Tmax:
        t0 = tcur + step*np.arange(chunk)
        # t nodes
        T = t0[:, None] + step*xg[None, :]
        x = alpha*eps*T
        f = f_L(x, r, k)*(alpha*eps)                       # dL/dt at nodes
        incr = step*(f @ Smat.T)                           # integral from panel start to each node
        panel_tot = step*(f @ wg)
        Lstart = Lc + np.concatenate([[0], np.cumsum(panel_tot)[:-1]])
        L = Lstart[:, None] + incr
        phi = np.exp(L)
        H = (4/rho)*np.sin(A*T)*np.sin(rho*T/2)/T**2
        g = H*np.real(np.exp(1j*c*T)*np.conj(phi))
        J += (step*(g @ wg)).sum()/math.pi
        Lc = Lc + panel_tot.sum(); tcur = t0[-1] + step
        absphi_T = math.exp(Lc.real)
        tail = absphi_T*4/(math.pi*rho*tcur)
        if tail < 2e-6: break
    return J, tail, tcur, absphi_T

if __name__ == "__main__":
    for mu in (1.0, 0.6, 0.3, 2.0):
        print(f"mu={mu} (k0={2*mu:.1f}, beta=1)")
        for eps in (0.3, 0.1, 0.03, 0.01, 0.001):
            row = []
            for cS in (0.5, 1.8):
                J1, tail1, T1, _ = E_h(mu, eps, cS, panel=0.4, n=16)
                J2, tail2, T2, _ = E_h(mu, eps, cS, panel=0.2, n=24)
                row.append((round(J2, 6), f"conv {abs(J1-J2):.1e}", f"tail {tail2:.1e}", f"T={T2:.0f}"))
            print(f"  eps={eps}: I1 {row[0]};  I2 {row[1]}", flush=True)
