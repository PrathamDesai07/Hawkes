"""One-off: Status = In progress + Notes for problem 7 (row 8) in the tracker; other cells untouched."""
import zipfile, re, shutil
from xml.sax.saxutils import escape

src = 'hawkes_research_problem_tracker.xlsx'
zin = zipfile.ZipFile(src)
x = zin.read('xl/worksheets/sheet1.xml').decode()
assert re.search(r'<c r="A8"[^>]*><v>7</v>', x)

note = ("25 Oct 2026. problem07_leverage.md; code/p5_leverage.py, p5_theorem_check.py, p5_component_check.py. K=1 only. Split gradient g = s (score at theta*) + h (drift). "
        "Lemma 7a (rigorous): E[s|past]=0, E[s s'|past]=E[dJ] (stopped martingale). Lemma 7b: E[h'F^-1 h] <= (ell/mu_min) E[v'dF v], ell = sup_u phi'F^-1 phi/lam* (point leverage). "
        "Lemma 7c: tr(F^-1 dF) <= r(x) logdet increment, r(x)=x/log(1+x) (gap leverage x). Theorem 7 (conditional on a.s. L1-L3): E R_n <= [P0 + (1+eta) c+ r(x) E logdet(F_n/F_0)]/(1-B), B=(1+1/eta) ell/mu_min<1, "
        "closing the gradient-norm term of Problem 5 (growth ~ d log T, constant ~3-6x too big). Leverage scale: worst case is the quiet gap, ell ~ [F^-1]_mumu/mu* ~ (k-1)/(eps T), so B<1 needs ~4(k-1)/(eps^2 mu_min) events: leverage burn-in ~1/eps^2 events (matches Problem 3 mixing; observed medians 1520 at eps=.1, 10086 at eps=.03). "
        "Numerics: all three lemmas verified; drift term is negligible in practice (0.02 vs bound 3); full bound valid but 15-40x loose. "
        "NOT done: L1-L3 as stopping-time/high-probability statements (needs Fisher concentration, Problems 2-3), projection, K>=2, sharp constant. Draft, not independently verified.")

pat = r'<c r="M8"( s="\d+")[^>]*?(?:/>|><is><t>.*?</t></is></c>)'
x, n = re.subn(pat, lambda m: f'<c r="M8"{m.group(1)} t="inlineStr"><is><t>{escape(note)}</t></is></c>', x, flags=re.S)
assert n == 1, n
x, a = re.subn(r'(<c r="L8"[^>]*><is><t>)[^<]*(</t></is></c>)', lambda m: m.group(1) + 'In progress' + m.group(2), x)
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
for r in range(2, 10):
    print(r, re.search(rf'<c r="A{r}"[^>]*><v>(\d+)</v>', y).group(1), re.search(rf'<c r="L{r}"[^>]*><is><t>(.*?)</t>', y).group(1))
