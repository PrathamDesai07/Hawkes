"""One-off: record Theorem 17A* (self-contained proof under D_kappa) in the Problem 17 notes."""
import zipfile, re, shutil
from xml.sax.saxutils import escape
src = 'hawkes_research_problem_tracker.xlsx'
zin = zipfile.ZipFile(src); x = zin.read('xl/worksheets/sheet1.xml').decode()
note = (" SELF-CONTAINED PROOF 25 Oct 2026 (problem17 section 7, code/p17/t_drift.py, t_kappa.py): Theorem 17A* [P]: a.s. convergence of the EF/OI recursion with gain constant kappa to the constrained KL pseudo-truth WITHOUT citing stochastic-approximation theory, under the drift condition (D_kappa) E[kappa (v'a)^2/lam_theta] - 2 v'grad L(theta) <= 0 on Theta. "
        "Proof: F-norm Lyapunov function + non-expansive projection; every non-drift term is o(t) by the ergodic theorem / martingale SLLN (no Poisson equation, no rates, no Markov property of the truth); W_N=o(N) and W_N>=cN|v|^2 give v->0. "
        "Closed form: Psi-bar = -2v'gradL(theta-dagger) - E[(v'a)^2 (2 lam* - kappa lam-dagger)/(lam-dagger lam_theta)], so (D_kappa) holds for kappa<=2 inf lam*/lam-dagger (EF) and, for OI, on a ratio-restricted Theta; well-specified: all kappa<=2. "
        "Numerics: D_1 holds on random samples of Theta for EF (K=1 strong misspec, K=2, active constraint; max Psi-bar/|v|^2 = -0.027,-0.016,-0.026), fails at kappa=2 and, for OI, at extreme theta; all kappa (0.5-4) converge numerically (condition is sufficient, not necessary). "
        "Still open: general result without D_kappa (cited ODE route), rates, regret.")
m = re.search(r'<c r="M18"( s="\d+")[^>]*?(?:/>|><is><t[^>]*>(.*?)</t></is></c>)', x, re.S); assert m
x = x[:m.start()] + f'<c r="M18"{m.group(1)} t="inlineStr"><is><t>{m.group(2) or ""}{escape(note)}</t></is></c>' + x[m.end():]
tmp = 'tracker_new_p17b.xlsx'
with zipfile.ZipFile(tmp, 'w', zipfile.ZIP_DEFLATED) as zo:
    for i in zin.infolist(): zo.writestr(i, x.encode() if i.filename == 'xl/worksheets/sheet1.xml' else zin.read(i.filename))
zin.close(); shutil.move(tmp, src); print('ok', zipfile.ZipFile(src).testzip())
