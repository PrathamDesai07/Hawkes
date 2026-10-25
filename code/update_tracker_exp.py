"""One-off: append a short update to Notes of problems 5 (row 6) and 7 (row 8); statuses untouched."""
import zipfile, re, shutil
from xml.sax.saxutils import escape

src = 'hawkes_research_problem_tracker.xlsx'
zin = zipfile.ZipFile(src)
x = zin.read('xl/worksheets/sheet1.xml').decode()
for r, i in ((6, 5), (8, 7)):
    assert re.search(rf'<c r="A{r}"[^>]*><v>{i}</v>', x), (r, i)

e5 = (" UPDATE 25 Oct (problem05_exponential_tails.md): exponential tails. E1 exact identities (eps*Lam_T = mu T + n M'_T - n dS; Poisson eq u=(s-ES)/gamma), verified 1e-14. "
      "E2 Lam_T tails exp(-eta^2 mu eps T/2n^2). E3 occupation functionals via Poisson equation + Freedman: P(int h <= Tm/2) <= 3 exp(-c T gamma), c=m^2 mu/(32 kappa^2 (1+eta) beta) ~3e-7 "
      "(empirical decay ~0.04-0.056 per unit T*gamma: bound valid but 1e5x pessimistic). E4 Newton martingale: with prob>=1-e^-x, -M_t <= X_t + x/theta for all t (theta=min(1/(2k1), c-/2.87)) so sup D_j <= D0' + x/theta + G_n. "
      "Theorem C-exp: bootstrap with T1 = O((ln(1/q)+d log T)/(gamma delta^2)). NOT done: E5 (exp bound on G_n, martingale part of score sum) in detail; expected regret with sharp constant needs episode/restart argument "
      "(union bound would pay crude linear rate on [0, T_c ~ log T_n]); constants; K>=2.")
e7 = " UPDATE 25 Oct: exponential-tail versions of Theorem F and the bootstrap in problem05_exponential_tails.md (G_n control sketched only)."

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
