import zipfile, re, shutil
from xml.sax.saxutils import escape
src = 'hawkes_research_problem_tracker.xlsx'
zin = zipfile.ZipFile(src); x = zin.read('xl/worksheets/sheet1.xml').decode()
notes = {
 14: " LaTeX: part4_resolutions.tex subsection added 2026-10-25. Drift check: coupling contraction matches exact rate gamma under fast mu-drift (A=0.4,P=60); K=1 window c_win under mu-drift 0.56/0.67 (fixed metric, C=16/64), 0.66/0.80 (local metric). THIRD PASS 2026-10-25 (problem13_third_pass.md): projection proved for exact F-norm projection (diagonal clip NOT covered); high-prob A1' via predictable leverage/stopping (A2 overshoot step open); Step-3 constant derived line by line (valid, ~1e10 conservative); A4 for K>=2 CONDITIONAL on U2 (kurtosis-type) and H1 (basis conditioning), K=2 numerics support; drift covered for mu-drift only (alpha fixed, box constant), alpha-drift open; projected jump-rate slope numerically noisy.",
}
for row, t in notes.items():
    m = re.search(rf'<c r="M{row}"( s="\d+")[^>]*?(?:/>|><is><t[^>]*>(.*?)</t></is></c>)', x, re.S); assert m, row
    x = x[:m.start()] + f'<c r="M{row}"{m.group(1)} t="inlineStr"><is><t>{m.group(2) or ""}{escape(t)}</t></is></c>' + x[m.end():]
tmp = 'tracker_new.xlsx'
with zipfile.ZipFile(tmp, 'w', zipfile.ZIP_DEFLATED) as zo:
    for i in zin.infolist():
        zo.writestr(i, x.encode() if i.filename == 'xl/worksheets/sheet1.xml' else zin.read(i.filename))
zin.close(); shutil.move(tmp, src); print(zipfile.ZipFile(src).testzip())
