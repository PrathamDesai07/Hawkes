"""One-off: append second-pass note to Problem 16 (row 17) Notes cell."""
import zipfile, re, shutil
from xml.sax.saxutils import escape
src = 'hawkes_research_problem_tracker.xlsx'
zin = zipfile.ZipFile(src); x = zin.read('xl/worksheets/sheet1.xml').decode()
note = (" SECOND PASS 25 Oct 2026 (problem16_kl_pseudotruth.md sec. 8, code/p16/t8.py, t9.py): explicit strong-convexity constant Hess L >= lam_lo*E[aa'/(Lambda*+sum beta_k S_k)^2] (certified, valid and non-vacuous in all runs, conservative by 17-300x, worse with K and n*); "
        "consistency of the QMLE written out; sandwich variance H^-1 J H^-1 with J=H+misspecification terms, reproduced numerically (naive H^-1/T underestimates var(mu) by 40%, sandwich by 17%; equal in the well-specified control); CLT tagged P* (needs mixing). "
        "Still open: an a priori proof that mu-dagger>mu_min (the natural bounds are vacuous; only the a-posteriori convexity certificate holds), active-constraint limit law. Replaces the earlier 'NOT proved' list except those two items.")
row=17
m = re.search(rf'<c r="M{row}"( s="\d+")[^>]*?(?:/>|><is><t[^>]*>(.*?)</t></is></c>)', x, re.S); assert m
x = x[:m.start()] + f'<c r="M{row}"{m.group(1)} t="inlineStr"><is><t>{m.group(2) or ""}{escape(note)}</t></is></c>' + x[m.end():]
tmp = 'tracker_new_p16.xlsx'
with zipfile.ZipFile(tmp, 'w', zipfile.ZIP_DEFLATED) as zo:
    for i in zin.infolist(): zo.writestr(i, x.encode() if i.filename == 'xl/worksheets/sheet1.xml' else zin.read(i.filename))
zin.close(); shutil.move(tmp, src); print('ok', zipfile.ZipFile(src).testzip())
