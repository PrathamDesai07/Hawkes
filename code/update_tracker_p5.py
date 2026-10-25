"""One-off: set Status = In progress and write Notes for problem 5 (row 6) in the tracker; other cells untouched."""
import zipfile, re, shutil
from xml.sax.saxutils import escape

src = 'hawkes_research_problem_tracker.xlsx'
zin = zipfile.ZipFile(src)
x = zin.read('xl/worksheets/sheet1.xml').decode()
assert re.search(r'<c r="A6"[^>]*><v>5</v>', x)

note = ("25 Oct 2026. problem05_online_newton.md; code/p5_*.py. First pass, K=1 only. Proposition 1 (rigorous, short): for the unprojected Fisher-scoring/online-Newton update, "
        "E<g_i,v_i | past> = E[v_i' dF_i v_i | past] exactly (first-event compensator identity), hence E R_n <= ||th_1-th*||^2_{F_0} + E sum g_i' F_i^-1 g_i, no exp-concavity needed. "
        "Checked numerically (pathwise identity to 1e-13; expectation identity within SE; bound ~2x loose). Empirical: regret slope (d/2) in log n; E[R_n - log T_n] ~ flat in eps "
        "(-1.2,-2.3,-2.4,-2.8 for eps=.3,.1,.03,.01): no visible +C log(1/eps) in time units; apparent 0.5/eps burn-in cost was the prior term rho0*|th_1-th*|^2, vanishes with a weak prior. "
        "NOT proven: bound on G (needs leverage bound a'F^-1 a/lam = O(1) and local control, i.e. Problem 7), sharp constant 1/2 (needs localisation), projected algorithm "
        "(diag projection not non-expansive in F-norm; unprojected iterates exit Theta_eps in 55-85% of runs), K>=2, high-probability version. Draft, not independently verified.")

pat = r'<c r="M6"( s="\d+")[^>]*?(?:/>|><is><t>.*?</t></is></c>)'
x, n = re.subn(pat, lambda m: f'<c r="M6"{m.group(1)} t="inlineStr"><is><t>{escape(note)}</t></is></c>', x, flags=re.S)
assert n == 1, n
x, a = re.subn(r'(<c r="L6"[^>]*><is><t>)[^<]*(</t></is></c>)', lambda m: m.group(1) + 'In progress' + m.group(2), x)
assert a == 1, a
tmp = 'tracker_new.xlsx'
with zipfile.ZipFile(tmp, 'w', zipfile.ZIP_DEFLATED) as zo:
    for i in zin.infolist():
        zo.writestr(i, x.encode() if i.filename == 'xl/worksheets/sheet1.xml' else zin.read(i.filename))
zin.close()
shutil.move(tmp, src)
z = zipfile.ZipFile(src)
print(z.testzip())
y = z.read('xl/worksheets/sheet1.xml').decode()
for r in range(2, 8):
    print(r, re.search(rf'<c r="A{r}"[^>]*><v>(\d+)</v>', y).group(1), re.search(rf'<c r="L{r}"[^>]*><is><t>(.*?)</t>', y).group(1))
