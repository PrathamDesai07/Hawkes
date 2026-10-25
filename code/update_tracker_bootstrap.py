"""One-off: append a short update to Notes of problems 5 (row 6) and 7 (row 8); statuses untouched."""
import zipfile, re, shutil
from xml.sax.saxutils import escape

src = 'hawkes_research_problem_tracker.xlsx'
zin = zipfile.ZipFile(src)
x = zin.read('xl/worksheets/sheet1.xml').decode()
for r, i in ((6, 5), (8, 7)):
    assert re.search(rf'<c r="A{r}"[^>]*><v>{i}</v>', x), (r, i)

e5 = (" UPDATE 25 Oct (problem05_bootstrap.md): TRUE BOOTSTRAP, K=1, two-stage, polynomial tails. Stage 1 with box constants: pathwise D_j <= D0' - X - M + G (M martingale, E m^2 <= (4/c-) E[v'dF v]); "
      "Doob L2 + Markov => P(sup D_j >= 3x) <= (1+16/c0-) Gamma/x, Gamma=O(d log T); with F >= lam_T I this forces the iterates into R(delta) for t>=T1 when delta^2 >= 24(1+16/c0-)c0+ L Gamma/(q m lamQ rho*^2 T1) (T1 gamma ~ d log T1/delta^2). "
      "Stage 2: with c+-=1+-delta, A->1, B->0, psi-factor 2 -> 1+O(delta): regret (d/2+o(1)) log T on the good event. Numerics: pathwise inequality margin 1e-14; E sum m^2/E X ~4; sup D_j <=14; "
      "last exit time from R(delta) ~ 1-3/(delta^2 gamma), slope of E R_n in log T ~0.9 (target d/2=1). NOT proved: expected regret with sharp constant (needs exponential tails; bad event only ~1/T1), localisation bookkeeping, K>=2, a-priori box.")
e7 = (" UPDATE 25 Oct: bootstrap done at polynomial-tail level (problem05_bootstrap.md); remaining for the leverage assumptions: exponential tails, localisation, K>=2.")

def append(x, r, extra):
    m = re.search(rf'<c r="M{r}"( s="\d+")[^>]*><is><t>(.*?)</t></is></c>', x, flags=re.S)
    assert m, r
    return x[:m.start()] + f'<c r="M{r}"{m.group(1)} t="inlineStr"><is><t>{m.group(2)}{escape(extra)}</t></is></c>' + x[m.end():]

x = append(x, 6, e5); x = append(x, 8, e7)
tmp = 'tracker_new.xlsx'
with zipfile.ZipFile(tmp, 'w', zipfile.ZIP_DEFLATED) as zo:
    for i in zin.infolist():
        zo.writestr(i, x.encode() if i.filename == 'xl/worksheets/sheet1.xml' else zin.read(i.filename))
zin.close()
shutil.move(tmp, src)
z = zipfile.ZipFile(src)
print(z.testzip())
y = z.read('xl/worksheets/sheet1.xml').decode()
for r in (6, 8):
    print(r, re.search(rf'<c r="L{r}"[^>]*><is><t>(.*?)</t>', y).group(1), len(re.search(rf'<c r="M{r}"[^>]*><is><t>(.*?)</t>', y, flags=re.S).group(1)))
