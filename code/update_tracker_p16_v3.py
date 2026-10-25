"""One-off: append third-pass note to Problem 16 (row 17)."""
import zipfile, re, shutil
from xml.sax.saxutils import escape
src = 'hawkes_research_problem_tracker.xlsx'
zin = zipfile.ZipFile(src); x = zin.read('xl/worksheets/sheet1.xml').decode()
note = (" THIRD PASS 25 Oct 2026 (sec. 9, code/p16/t10.py, t11.py, t12.py): closed the last two open items. (i) A priori lower bound mu-dagger >= mu_0, mu_0 = sup{mu: min_{alpha in A_eps} E[lam*/(mu+alpha'S)] >= 1} [P], computable from the data law alone; "
        "non-vacuous in all runs (Phi(0)=1.06-1.63) but conservative 2-7x (mu_0~0.28 vs mu-dagger 1.55 at K=3, eps=0.05). (ii) Limit law with active constraints [P*]: sqrt(T)(theta-hat - theta-dagger) -> V(V'HV)^-1 V'Z, rank d-|active|; "
        "tested K=2, n<=0.5 active: var(mu) 3.41e-3 empirical vs 3.45e-3 formula (population J), corr -0.209 vs -0.210; unconstrained sandwich overstates 5x. "
        "Remaining caveats: CLT needs mixing (P*), strict complementarity assumed, batch-means J is 17-27% low at B=20.")
row=17
m = re.search(rf'<c r="M{row}"( s="\d+")[^>]*?(?:/>|><is><t[^>]*>(.*?)</t></is></c>)', x, re.S); assert m
x = x[:m.start()] + f'<c r="M{row}"{m.group(1)} t="inlineStr"><is><t>{m.group(2) or ""}{escape(note)}</t></is></c>' + x[m.end():]
tmp = 'tracker_new_p16.xlsx'
with zipfile.ZipFile(tmp, 'w', zipfile.ZIP_DEFLATED) as zo:
    for i in zin.infolist(): zo.writestr(i, x.encode() if i.filename == 'xl/worksheets/sheet1.xml' else zin.read(i.filename))
zin.close(); shutil.move(tmp, src); print('ok', zipfile.ZipFile(src).testzip())
