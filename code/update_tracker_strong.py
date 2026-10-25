import zipfile, re, shutil
from xml.sax.saxutils import escape
src = 'hawkes_research_problem_tracker.xlsx'
zin = zipfile.ZipFile(src); x = zin.read('xl/worksheets/sheet1.xml').decode()
notes = {
 3: " STRENGTHENED 25 Oct (tangent_bound.md): tangent-line minorant I >= (2/l0)M - (1/l0^2)M3 with exact moment systems (orders<=3) gives a rigorous lower bound within 1.3-3.5x of Monte Carlo (K=1,2,3), correct eps-scaling (~0.25 eps); was 1e3-1e8 loose.",
 13: " STRENGTHENED 25 Oct (tangent_bound.md): by Weyl, EVERY eigenvalue of I is bounded below rigorously (bound/MC = 0.23-0.85): top ~1/eps, middle O(1) independent of eps, ridge >= ~0.25 eps. Upper bounds (ridge Rayleigh) and eigenvector convergence still numerical.",
 6: " STRENGTHENED 25 Oct: (1) theoremR_proof_complete.md writes every step of the expected-regret theorem (assembly needs only predictable good-gap indicators + one per-time failure bound H3); (2) information level without Assumption M via the tangent minorant: lambda_*/2 ~ 0.013 for the K=2 example vs ~1e-6 before.",
}
for row, t in notes.items():
    m = re.search(rf'<c r="M{row}"( s="\d+")[^>]*?(?:/>|><is><t[^>]*>(.*?)</t></is></c>)', x, re.S); assert m, row
    x = x[:m.start()] + f'<c r="M{row}"{m.group(1)} t="inlineStr"><is><t>{m.group(2) or ""}{escape(t)}</t></is></c>' + x[m.end():]
tmp = 'tracker_new.xlsx'
with zipfile.ZipFile(tmp, 'w', zipfile.ZIP_DEFLATED) as zo:
    for i in zin.infolist():
        zo.writestr(i, x.encode() if i.filename == 'xl/worksheets/sheet1.xml' else zin.read(i.filename))
zin.close(); shutil.move(tmp, src); print(zipfile.ZipFile(src).testzip())
