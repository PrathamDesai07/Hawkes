"""Deterministic prediction with the time-varying gain implied by a per-event gain eta:  dy/dt = eta * r_t * (mu_t - y),  r_t = mu_t/eps  (events per unit time at the current mu).
Returns bias^2 = mean over the last periods of (y - mu_t)^2 (no noise).  Compare with the constant-gain first-order filter formula and with measurements."""
import numpy as np, math
def tv_bias2(shape, A, eta, eps, P=1200.0, mu0=1.0, cycles=8, n_per=24000):
    t = np.linspace(0, cycles*P, cycles*n_per+1); dt = t[1]-t[0]
    if shape == 'sin': mu = mu0*(1 + A*np.sin(2*np.pi*t/P))
    else: mu = mu0*(1 + A*(2/np.pi)*np.arcsin(np.sin(2*np.pi*t/P)))
    y = np.empty_like(mu); y[0] = mu[0]
    for i in range(len(t)-1):
        k = eta*mu[i]/eps; y[i+1] = y[i] + dt*k*(mu[i] - y[i])
    last = t > (cycles-3)*P; return np.mean((y[last]-mu[last])**2)
if __name__ == "__main__":
    eps = 0.4; eta = 0.004; P = 1200.0; r0 = 1/eps; w = 2*math.pi/(P*r0)
    print("sinusoid, eps=0.4, eta=0.004, P=1200:  A | const-gain filter | time-varying-gain filter | ratio | measured ratios (known, adaptive, fixed) from p14_decomp fundamental-error^2")
    meas = {0.2: (1.07, 1.41, 1.07), 0.35: (1.15, 1.48, 0.98), 0.5: (1.24, 1.61, 0.92)}
    for A in (0.2, 0.35, 0.5):
        cg = (A*A/2)*w*w/(0.004**2 + w*w); tv = tv_bias2('sin', A, eta, eps)
        print(f"  A={A}: {cg:.5f} | {tv:.5f} | ratio {tv/cg:.2f} | measured {meas[A]}")
    for eta in (0.002, 0.004):
        for A in (0.3,):
            print(f"triangle eps=0.4 A={A} eta={eta}: time-varying-gain filter bias^2 = {tv_bias2('tri', A, eta, eps):.5f}")
