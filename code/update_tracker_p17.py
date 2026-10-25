"""One-off: Problem 17 (row 18) status and notes."""
import zipfile, re, shutil
from xml.sax.saxutils import escape
src = 'hawkes_research_problem_tracker.xlsx'
zin = zipfile.ZipFile(src); x = zin.read('xl/worksheets/sheet1.xml').decode()
note = (" 25 Oct 2026. problem17_online_convergence.md; code/p17/ (online17.c, lib17.py, w1,m1,m2,m3,trap,ev). Draft, not independently verified. "
 "Theorem 17A [P*]: projected online Newton (Fisher or observed-information preconditioner, exact F-norm projection) converges a.s. to the unique constrained KL pseudo-truth theta-dagger_Theta of Problem 16, active constraint included (ODE method with Markov noise + matrix gain; relies on cited Benveniste-Metivier-Priouret conditions, not re-verified). "
 "New exact facts: Lyapunov drift of ||theta-theta-dagger||^2_F is -v'(2Hbar - I_F)v with a favourable-sign multiplier term when the constraint binds; A = I_F^-1 H has eigenvalue exactly 1 at an interior theta-dagger. "
 "Theorem 17B/C [P*]: sqrt(n)(theta_n - theta-dagger) -> N(0,Sigma), Lyapunov equation (A-1/2)Sigma+Sigma(A-1/2)'=Lambda* Ibar^-1 J Ibar^-1; with observed-information preconditioner Sigma = Lambda* H^-1 J H^-1 (batch sandwich) with no side condition; active constraint: same in the face subspace, n*dist to face -> const (O(1/n)), projection active in a stable 3% of steps. "
 "Numerics (power-law truth n*=0.8, K=2, 600 replicates up to 1e5 events): n*Var matches prediction to 5-14% interior and 1-12% with active constraint (MC ~6%); n*MSE flat (rate 1/n); eigenvalues of I_F^-1 H in [0.97,1.002] in every example so EF=OI; OI corner-trapping found in 1/300 runs when mu_min=1e-3, none for EF or mu_min>=0.1. "
 "MISSING for Done: self-contained proof of the Poisson-equation/regularity step, rigorous excursion bound for the distance to the face, degenerate (no strict complementarity) case, regret statement (Prop 17.4 is a conjecture), K>2.")
row=18
m = re.search(rf'<c r="L{row}"( s="\d+")[^>]*?>.*?</c>', x, re.S); assert m
x = x[:m.start()] + f'<c r="L{row}"{m.group(1)} t="inlineStr"><is><t>In progress</t></is></c>' + x[m.end():]
m = re.search(rf'<c r="M{row}"( s="\d+")[^>]*?(?:/>|><is><t[^>]*>(.*?)</t></is></c>)', x, re.S); assert m
x = x[:m.start()] + f'<c r="M{row}"{m.group(1)} t="inlineStr"><is><t>{m.group(2) or ""}{escape(note)}</t></is></c>' + x[m.end():]
tmp='tracker_new_p17.xlsx'
with zipfile.ZipFile(tmp,'w',zipfile.ZIP_DEFLATED) as zo:
    for i in zin.infolist(): zo.writestr(i, x.encode() if i.filename=='xl/worksheets/sheet1.xml' else zin.read(i.filename))
zin.close(); shutil.move(tmp,src); print('ok'); y=zipfile.ZipFile(src).read('xl/worksheets/sheet1.xml').decode(); print(re.search(r'<c r="C18"[^>]*>.*?</c>',y,re.S).group(0)[:150])
