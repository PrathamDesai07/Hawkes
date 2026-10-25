"""Transfer function (sinusoidal mu drift) for general K, causal sampling, vs first-order filter; noise coefficient vs r/(2 eff) with eff from stationary MC."""
import numpy as np, math, sys
from p14_gen import run
from p2_identifiability import mc
cfgs = {'K2a': ([1., 2.], 1.5, 0.4), 'K2b': ([1., 4.], 1.0, 0.4), 'K3a': ([1., 2., 4.], 1.5, 0.4), 'K3b': ([1., 3., 9.], 1.5, 0.4), 'K2c': ([1., 2.], 1.5, 0.2)}
name = sys.argv[1]; R = int(sys.argv[2]); eta = float(sys.argv[3]); est = sys.argv[4]; A = float(sys.argv[5]) if len(sys.argv) > 5 else 0.3
beta, mu_s, eps = cfgs[name]; beta = np.array(beta); K = len(beta); r = mu_s/eps; P = 1200.0; T = 6000.0
w = 2*math.pi/(P*r); fn = lambda t: mu_s*(1 + A*np.sin(2*np.pi*t/P)); mu_max = mu_s*(1+A)
tg = np.arange(1, int(T/5.0)+1)*5.0; sel = tg > 0.2*T; t = tg[sel]; e = np.exp(-2j*np.pi*t/P); tgt = fn(t)
reps = []; cs = []; offs = []
for s in range(R):
    _, x = run(est, beta, eps, fn, mu_max, eta, T, 5000+s, mu_s); x = x[sel]; reps.append(x)
    cs.append(1j*2*np.mean((x - x.mean())*e)/(A*mu_s)); offs.append(x.mean() - mu_s)
cs = np.array(cs); H = cs.mean(); se = cs.std(ddof=1)/math.sqrt(R); Hf = eta/(eta + 1j*w)
reps = np.array(reps); mse = np.mean((reps - tgt)**2); mean = reps.mean(0); varc = reps.var(0, ddof=1).mean()/R
bias2 = np.mean((mean - tgt)**2) - varc; noise = (mse - bias2)
I, _ = mc(mu_s, beta, eps, 200000, seed=9); eff = 1/np.linalg.inv(I)[0, 0] if est == 'adaptive' else I[0, 0]
print(f"{name} beta={beta.tolist()} mu={mu_s} eps={eps} est={est} eta={eta} A={A}: H={H.real:.3f}{H.imag:+.3f}i (|H|={abs(H):.3f}, arg={np.angle(H):+.3f}) se {se:.3f} | filter {Hf.real:.3f}{Hf.imag:+.3f}i (|H|={abs(Hf):.3f}, arg={np.angle(Hf):+.3f}) | offset {np.mean(offs):+.3f}")
print(f"      noise part {noise:.4f} = {noise/eta:.1f}*eta vs predicted r/(2 eff) = {r/(2*eff):.1f} (eff={eff:.4f}, cond(I)={np.linalg.cond(I):.0f})", flush=True)
