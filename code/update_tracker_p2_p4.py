"""One-off: update tracker notes for problems 2 and 4 (inline-string cells in sheet1.xml)."""
import zipfile, re, shutil
from xml.sax.saxutils import escape

src = 'hawkes_research_problem_tracker.xlsx'
zin = zipfile.ZipFile(src)
x = zin.read('xl/worksheets/sheet1.xml').decode()
assert re.search(r'<c r="A3"[^>]*><v>2</v>', x) and re.search(r'<c r="A5"[^>]*><v>4</v>', x)

n2 = ("25 Oct 2026. problem02_identifiability.md, code/p2_window_bound.py, code/p2_lemma3_check.py. Found and fixed a thinning bug in p2_identifiability.py "
      "(waiting time redrawn after rejection) that had made lam_min ~10x too small for K>=2. Results (mu=1, beta=r^(k-1), MC 60-80k events): "
      "Lemma 3 bound (eps/mu)*inf_v(E|v'a|)^2 is within 2-4x of lam_min for K<=3. K>=2: lam_min does NOT shrink with eps (K=2,r=2: ~0.045; r=4: ~0.065), "
      "cond ~ 1/eps (vs eps^-2 for K=1); lam_min falls ~10x per added exponential at r=2 (0.069, 0.046, 0.0045 for K=1,2,3 at eps=0.1); wider grid helps. "
      "Rigorous Lemma W: I >= J_h (first-arrival window); J_h->I as h->0 so it is a route, not a number. "
      "NOT done: closed-form lower bound in K and grid ratio (route: Gram matrix of {1,e^{-beta_k u}} on [0,h]). Draft, not independently verified.")
n4 = ("24-25 Oct 2026. problem04_projection.md, code/p4_projection.py. Exact projection onto {mu>=mu_min, alpha>=0, sum alpha_k/beta_k<=1-eps} in a DIAGONAL metric: "
      "mu separable; alpha_k = max(0, y_k - tau*w_k/d_k), tau from sorted breakpoints, O(K log K); matches SLSQP on 300 random instances (4e-6); <1 ms for K up to 1e4. "
      "Regret-constant effect: non-expansive in the update's own metric, so no extra constant; the cost is cond(D^-1/2 I D^-1/2). K=1: universal ~5.83 (rho^2->1/2), "
      "independent of eps (raw ~1e6 at eps=1e-3). K>=2 (corrected 25 Oct after fixing a simulator bug): scaled cond only 1.2-3x below raw and still grows as eps falls "
      "(K=2,r=2: 118->368; K=3: 765->2477), so a diagonal metric is NOT enough for K>=2. Open: dense-metric projection (QP, ~O(K^2)), O(K) version, K>=2 preconditioner. "
      "Status In progress because 'Done when' asks for the Fisher/ONS metric.")


def setrow(x, r, note, status=None):
    pat = rf'<c r="M{r}"( s="\d+")[^>]*?(?:/>|><is><t>.*?</t></is></c>)'
    x, b = re.subn(pat, lambda m: f'<c r="M{r}"{m.group(1)} t="inlineStr"><is><t>{escape(note)}</t></is></c>', x, flags=re.S)
    assert b == 1, (r, b)
    if status:
        x, a = re.subn(rf'(<c r="L{r}"[^>]*><is><t>)[^<]*(</t></is></c>)', lambda m: m.group(1) + status + m.group(2), x)
        assert a == 1, (r, a)
    return x


x = setrow(x, 3, n2, 'In progress')
x = setrow(x, 5, n4)
tmp = 'tracker_new.xlsx'
with zipfile.ZipFile(tmp, 'w', zipfile.ZIP_DEFLATED) as zo:
    for i in zin.infolist():
        zo.writestr(i, x.encode() if i.filename == 'xl/worksheets/sheet1.xml' else zin.read(i.filename))
zin.close()
shutil.move(tmp, src)
z = zipfile.ZipFile(src)
print(z.testzip())
y = z.read('xl/worksheets/sheet1.xml').decode()
for r in (2, 3, 4, 5, 11):
    print(r, re.search(rf'<c r="L{r}"[^>]*><is><t>(.*?)</t>', y).group(1))
