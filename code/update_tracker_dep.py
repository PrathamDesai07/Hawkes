import zipfile, re, shutil
from xml.sax.saxutils import escape
src = 'hawkes_research_problem_tracker.xlsx'
zin = zipfile.ZipFile(src); x = zin.read('xl/worksheets/sheet1.xml').decode()
m = re.search(r'<c r="M6"( s="\d+")[^>]*?(?:/>|><is><t[^>]*>(.*?)</t></is></c>)', x, re.S); assert m
add = (" DEPENDENCIES 25 Oct (K1_completion.md): optional stopping closed (predictable good-gap indicators replace restarts); integrability in Thm F-exp closed; "
       "Assumption M closed for eps<=eps_M via the Gamma comparator of Thm T (E_Gamma h = 0.164, 0.044 vs simulated 0.15, 0.05); C(eps)=O(mu d ln(1/(eps delta))/(beta eps^2 delta^2)) (order only). "
       "Still assumed: known box; K>=2 open; constants conservative; Thm T itself is a proof outline with numerically verified estimates.")
x = x[:m.start()] + f'<c r="M6"{m.group(1)} t="inlineStr"><is><t>{m.group(2) or ""}{escape(add)}</t></is></c>' + x[m.end():]
tmp = 'tracker_new.xlsx'
with zipfile.ZipFile(tmp, 'w', zipfile.ZIP_DEFLATED) as zo:
    for i in zin.infolist(): zo.writestr(i, x.encode() if i.filename == 'xl/worksheets/sheet1.xml' else zin.read(i.filename))
zin.close(); shutil.move(tmp, src)
y = zipfile.ZipFile(src).read('xl/worksheets/sheet1.xml').decode()
print(zipfile.ZipFile(src).testzip(), y.count('E5 + EXPECTED REGRET 25 Oct'), y.count('DEPENDENCIES 25 Oct'))
