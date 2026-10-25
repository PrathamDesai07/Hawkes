import zipfile, re, shutil
from xml.sax.saxutils import escape
src = 'hawkes_research_problem_tracker.xlsx'
zin = zipfile.ZipFile(src); x = zin.read('xl/worksheets/sheet1.xml').decode()
assert re.search(r'<c r="A6"[^>]*><v>5</v>', x)
m = re.search(r'<c r="M6"( s="\d+")[^>]*?(?:/>|><is><t[^>]*>(.*?)</t></is></c>)', x, re.S); assert m
add = (" E5 + EXPECTED REGRET 25 Oct (problem05_E5_expected_regret.md, K1_main_results.md): E5 PROVED (Lemmas E5a-d, Prop E5: G_n <= (1+o(1)) logdet + o(1)x + B_X X w.p. >= 1-4e^-x; pathwise inequalities checked, 0 violations in 1.4e6 gaps). "
       "Found+fixed a flaw in Theorem 7/Lemma 7c: 'gap leverage <= xbar a.s.' cannot hold (depends on Exp(1) variable); truncation fixes it. Theorem R (K=1): E R_n <= (1/2+O(delta+eta)) d ln T + C, C independent of n, via charging bad events to their own time (no episode union bound); "
       "PROVED MODULO Thm F-exp, Assumption M (numerical), known box, and unwritten optional-stopping across random restarts. Rate constants ~1e5 worse than observed. K>=2 open.")
new = f'<c r="M6"{m.group(1)} t="inlineStr"><is><t>{m.group(2) or ""}{escape(add)}</t></is></c>'
x = x[:m.start()] + new + x[m.end():]
tmp = 'tracker_new.xlsx'
with zipfile.ZipFile(tmp, 'w', zipfile.ZIP_DEFLATED) as zo:
    for i in zin.infolist(): zo.writestr(i, x.encode() if i.filename == 'xl/worksheets/sheet1.xml' else zin.read(i.filename))
zin.close(); shutil.move(tmp, src); print(zipfile.ZipFile(src).testzip())
