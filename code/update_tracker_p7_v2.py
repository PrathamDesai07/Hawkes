"""One-off: append a sentence about the high-probability result to the Notes of problems 5 (row 6) and 7 (row 8); statuses untouched."""
import zipfile, re, shutil
from xml.sax.saxutils import escape

src = 'hawkes_research_problem_tracker.xlsx'
zin = zipfile.ZipFile(src)
x = zin.read('xl/worksheets/sheet1.xml').decode()
assert re.search(r'<c r="A6"[^>]*><v>5</v>', x) and re.search(r'<c r="A8"[^>]*><v>7</v>', x)

extra = (" UPDATE 25 Oct: problem05_fisher_concentration.md (K=1, oracle information): Lemma F1 var(int h(S_t)dt) <= 2T Lip sqrt(m) sigma_S/gamma via the Problem-3 contraction; "
         "Lemma F2 two-interval bound J_T >= (O/L) Q (exact, deterministic); Theorem F: J_T >= (mT/4L) Q w.p. >=1-2q for all T>=T0=16 sigma_S/(rho gamma m^1.5 q) = O(1/(beta eps)), "
         "so Theorem-7 leverage assumptions hold w.h.p. given (L1). Numerics: variance bound holds with 500-1500x margin; pathwise leverage bound 0 violations in ~1300 checks; "
         "rigorous T0*gamma ~ 4e5 vs empirical ~20-50 (Chebyshev, 1e4x conservative). OPEN: circularity (L1 needs consistency which needs the information bound), Assumption M (pi-mass of bumps, numerical), exponential tails, projection, K>=2.")

def append(x, r):
    m = re.search(rf'<c r="M{r}"( s="\d+")[^>]*><is><t>(.*?)</t></is></c>', x, flags=re.S)
    assert m, r
    new = f'<c r="M{r}"{m.group(1)} t="inlineStr"><is><t>{m.group(2)}{escape(extra)}</t></is></c>'
    return x[:m.start()] + new + x[m.end():]

x = append(x, 6)
x = append(x, 8)
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
