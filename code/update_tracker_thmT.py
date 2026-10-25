import zipfile, re, shutil
from xml.sax.saxutils import escape
src='hawkes_research_problem_tracker.xlsx'; zin=zipfile.ZipFile(src); x=zin.read('xl/worksheets/sheet1.xml').decode()
add={12:" THM T CONSTANTS 25 Oct (thmT_explicit_constants.md): Theorem T now has an explicit chain: |d|<=19 (true sup 1.143), tail exponent q=k0/(2(1+3/X0))>1/2 (k/(1+r)^2=k0 exactly), certified TV(pi,Gamma)<=Phi/(2 sqrt2); k0=2: 2e-2 (eps=1e-4), 1e-3 (1e-5), 8.6e-5 (1e-6); bound 1e2-1e3 above truth, valid in all tests. Ridge pair: TV<=2.4e-3+0.39 delta (eps<=1e-5, k0=2, delta<=0.3). Phi by floating-point quadrature of hand-proved majorants, not interval arithmetic; KL step evaluated on a grid; informative only for eps<~1e-4.",
     6:" ASSUMPTION M certified 25 Oct (thmT_explicit_constants.md): K=1, k0=2, mu=beta=1: m>=0.0248 for eps<=1e-4 (E_Gamma h = 0.164, 0.0444 minus TV<=0.0196)."}
for row,t in add.items():
    m=re.search(rf'<c r="M{row}"( s="\d+")[^>]*?(?:/>|><is><t[^>]*>(.*?)</t></is></c>)',x,re.S); assert m
    x=x[:m.start()]+f'<c r="M{row}"{m.group(1)} t="inlineStr"><is><t>{m.group(2) or ""}{escape(t)}</t></is></c>'+x[m.end():]
tmp='tracker_new.xlsx'
with zipfile.ZipFile(tmp,'w',zipfile.ZIP_DEFLATED) as zo:
    for i in zin.infolist(): zo.writestr(i, x.encode() if i.filename=='xl/worksheets/sheet1.xml' else zin.read(i.filename))
zin.close(); shutil.move(tmp,src); print(zipfile.ZipFile(src).testzip())
