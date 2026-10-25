"""One-off: refresh the Notes cell of problem 2 (row 3) in the tracker; leave all other cells as found."""
import zipfile, re, shutil
from xml.sax.saxutils import escape

src = 'hawkes_research_problem_tracker.xlsx'
zin = zipfile.ZipFile(src)
x = zin.read('xl/worksheets/sheet1.xml').decode()
assert re.search(r'<c r="A3"[^>]*><v>2</v>', x)

note = ("25 Oct 2026. problem02_identifiability.md; code/p2_window_bound.py, p2_lemma3_check.py, p2_gram_bound.py, p2_truncation_check.py, p2_truncation_bound.py. "
        "Fixed a thinning bug in p2_identifiability.py (had made lam_min ~10x too small for K>=2). Findings (mu=1, beta=r^(k-1), MC 60-80k events): "
        "K>=2 lam_min does NOT shrink with eps (K=2: ~0.045 at r=2, ~0.065 at r=4), cond ~ 1/eps (K=1: eps^-2); lam_min falls ~10x per added exponential (r=2: 0.069, 0.046, 0.0045 for K=1,2,3 at eps=0.1). "
        "Lower bounds: (1) Gram-matrix route is rigorous but useless (1e-6 to 1e-16 of truth, decays like eps^(2K+1)). "
        "(2) Theorem T (truncation + Lemma-B moments): lam_min(I) >= (1/L)[lam_min(M) - sqrt(E|a|^4) sqrt(E lam^p / L^p)], M = E[a a^T] exact via a K^2 linear system (matches MC to 1-3%); "
        "bound is 2e-3 to 5e-2 of the MC lam_min for K<=3 and non-vacuous; pure truncation without tail penalty would be within 2-8x. "
        "Open: closed-form formula in K and grid ratio, sharper tail bound, eps-sharpness for K>=2. Depends on Lemma B (Problem 3 draft). Draft, not independently verified.")

pat = r'<c r="M3"( s="\d+")[^>]*?(?:/>|><is><t>.*?</t></is></c>)'
x, n = re.subn(pat, lambda m: f'<c r="M3"{m.group(1)} t="inlineStr"><is><t>{escape(note)}</t></is></c>', x, flags=re.S)
assert n == 1, n
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
