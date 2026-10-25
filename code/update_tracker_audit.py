import zipfile, re, shutil
from xml.sax.saxutils import escape
src = 'hawkes_research_problem_tracker.xlsx'
zin = zipfile.ZipFile(src); x = zin.read('xl/worksheets/sheet1.xml').decode()
notes = {
 12: " FINAL AUDIT 25 Oct (final_audit.md): the 'no consistent estimator given S0' statements (Thm B, Cor D, Cor E') were stronger than the KL lemma proves: the correct bound is TV(Q) <= TV(pi_theta, pi_theta') + sqrt(rho T u'Iu/2). Rigorous impossibility: K=1 stationary record (Thm T), any K stationary record with unobserved past (Thm V). K>=2 with S0 observed: open.",
 6: " FINAL AUDIT 25 Oct (final_audit.md): fixed constant in Lemma E5c (xbar = lbar m/(c^- mu*), not c^+), replaced loose general-K Lemma E1 by the exact Z_c identity; K=2 pathwise Newton/martingale inequalities verified (code/p5_k2_pathwise_check.py).",
}
for row, t in notes.items():
    m = re.search(rf'<c r="M{row}"( s="\d+")[^>]*?(?:/>|><is><t[^>]*>(.*?)</t></is></c>)', x, re.S); assert m, row
    x = x[:m.start()] + f'<c r="M{row}"{m.group(1)} t="inlineStr"><is><t>{m.group(2) or ""}{escape(t)}</t></is></c>' + x[m.end():]
tmp = 'tracker_new.xlsx'
with zipfile.ZipFile(tmp, 'w', zipfile.ZIP_DEFLATED) as zo:
    for i in zin.infolist():
        zo.writestr(i, x.encode() if i.filename == 'xl/worksheets/sheet1.xml' else zin.read(i.filename))
zin.close(); shutil.move(tmp, src); print(zipfile.ZipFile(src).testzip())
