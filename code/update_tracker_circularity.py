"""One-off: append a short update to Notes of problems 4 (row 5), 5 (row 6) and 7 (row 8); statuses untouched."""
import zipfile, re, shutil
from xml.sax.saxutils import escape

src = 'hawkes_research_problem_tracker.xlsx'
zin = zipfile.ZipFile(src)
x = zin.read('xl/worksheets/sheet1.xml').decode()
for r, i in ((5, 4), (6, 5), (8, 7)):
    assert re.search(rf'<c r="A{r}"[^>]*><v>{i}</v>', x), (r, i)

e4 = (" UPDATE 25 Oct (problem05_circularity.md): for K=1 the diagonal clip is NOT adequate: it violates non-expansiveness in 0.7-4.6% of steps and inflates regret ~10x when the box is tight "
      "(E R_n 16 vs 3.5 at eps=.1; 79 vs 5.4 at eps=.03). Exact F-norm projection onto a box is O(1) for K=1 (edge minimisers + corners) and had 0 violations in 24 configs. Use it.")
e5 = (" UPDATE 25 Oct (problem05_circularity.md): circularity avoided by design, not by a consistency bootstrap: compact box Theta_c known a priori + exact F-norm projection => (L1) holds deterministically "
      "(c+=max(mu_hi/mu*,al_hi/al*), c-=min(...)); Prop 1 survives projection (non-expansiveness in F-norm); leverage w.h.p. after T0=O(1/(beta eps)) via Theorem F. "
      "Numerics (K=1): R_n - log T_n between -4.4 and -0.5 in 24 configs, 0 non-expansiveness violations. Sketched only: stopping-time localisation of the w.h.p. event; burn-in term crude. "
      "OPEN: true bootstrap (c->1), K>=2 (QP projection), a priori scale knowledge, exponential tails.")
e7 = (" UPDATE 25 Oct: (L1) now holds by construction for the compact-box exact-projection algorithm (problem05_circularity.md); F >= F0 + J/c+ verified PSD in all runs. Remaining: localisation bookkeeping.")

def append(x, r, extra):
    m = re.search(rf'<c r="M{r}"( s="\d+")[^>]*><is><t>(.*?)</t></is></c>', x, flags=re.S)
    assert m, r
    return x[:m.start()] + f'<c r="M{r}"{m.group(1)} t="inlineStr"><is><t>{m.group(2)}{escape(extra)}</t></is></c>' + x[m.end():]

x = append(x, 5, e4); x = append(x, 6, e5); x = append(x, 8, e7)
tmp = 'tracker_new.xlsx'
with zipfile.ZipFile(tmp, 'w', zipfile.ZIP_DEFLATED) as zo:
    for i in zin.infolist():
        zo.writestr(i, x.encode() if i.filename == 'xl/worksheets/sheet1.xml' else zin.read(i.filename))
zin.close()
shutil.move(tmp, src)
z = zipfile.ZipFile(src)
print(z.testzip())
y = z.read('xl/worksheets/sheet1.xml').decode()
for r in (5, 6, 8):
    print(r, re.search(rf'<c r="L{r}"[^>]*><is><t>(.*?)</t>', y).group(1), len(re.search(rf'<c r="M{r}"[^>]*><is><t>(.*?)</t>', y, flags=re.S).group(1)))
