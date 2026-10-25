"""One-off: update Read status on the Literature sheet (sheet4) for the papers actually read, and append prior-art notes to problems 9 and 11 (sheet1)."""
import zipfile, re, shutil
from xml.sax.saxutils import escape

src = 'hawkes_research_problem_tracker.xlsx'
zin = zipfile.ZipFile(src)
lit = zin.read('xl/worksheets/sheet4.xml').decode()
pro = zin.read('xl/worksheets/sheet1.xml').decode()

# sanity: row n+1 of the Literature sheet is read-order n
def lit_cell(x, ref):
    m = re.search(rf'<c r="{ref}"([^>]*?)(?:/>|>(.*?)</c>)', x, re.S); assert m, ref; return m
checks = {"B2": "Learning the exogenous rate", "B3": "Minimax and adaptive", "B4": "Tracking Dynamic Point", "B5": "Online Learning for Multivariate",
          "B6": "Nonparametric Hawkes", "B10": "Limit theorems for nearly unstable", "B11": "Logarithmic regret"}
for ref, frag in checks.items():
    assert frag in re.sub(r'<[^>]+>', '', lit_cell(lit, ref).group(0)), (ref, frag)

status = {
 "H2": "Read 25 Oct 2026 (full text). Thm 2.2, Cor 2.3, Thm 4.1, Rem 4.8 re-derived and tested numerically; overlap with C3 and the kappa<=1 half of C4 confirmed. See literature_crosscheck.md.",
 "H3": "Partly read 25 Oct 2026 (abstract, Sec 5, Sec 6.2 only). Thm 3.1 and its proof NOT checked.",
 "H4": "Partly read 25 Oct 2026 (Thm 1 statement, Appendix C parametrisation). Proof NOT checked.",
 "H5": "Not read (the arXiv version in the next row was read instead).",
 "H6": "Read 25 Oct 2026 (assumptions, eq. (8)-(9), Thm 2, Cor 1). Proof of Thm 2 only skimmed. See literature_crosscheck.md.",
 "H10": "Read 25 Oct 2026 (assumptions, Prop 2.1-2.2, Thm 2.2, heuristic derivation); tested against exact K=1 Hawkes. Proof NOT checked.",
 "H11": "Original returned HTTP 403. Read Hazan's monograph arXiv 1909.05207 Sec 4.2/4.4 instead (criterion, Thm 4.5 statement); proof not re-derived.",
}
for ref, txt in status.items():
    m = lit_cell(lit, ref); s = re.search(r' s="(\d+)"', m.group(1)).group(1)
    new = f'<c r="{ref}" s="{s}" t="inlineStr"><is><t>{escape(txt)}</t></is></c>'
    lit = lit[:m.start()] + new + lit[m.end():]

notes = {
 (10, 9): " LITERATURE 25 Oct: the two-point KL bound is already Herrera-Marin arXiv 2609.38268 Thm 2.2 + Cor 2.3 (empty start, shared mu and kernel shape, any normalised shape incl. sums of exponentials; unrefereed). Re-derived and tested on K=1 (bound held in 12/12 configs, 2-7x loose). Novelty here must come from online-estimator regret or other directions (individual alpha_k, stationary start).",
 (12, 11): " LITERATURE 25 Oct: divergence of mu-information for kappa<=1 is Herrera-Marin Remark 4.8 (Fatou); he states the kappa>1 direction is NOT proved. Gamma shape/Feller threshold k=2*beta*mu/alpha^2 confirmed against Jaisson-Rosenbaum Thm 2.2 (shape 2m*mu; tested, autocorrelation exp(-eps*beta*t) also confirmed).",
}
for (row, pid), extra in notes.items():
    assert re.search(rf'<c r="A{row}"[^>]*><v>{pid}</v>', pro), (row, pid)
    m = re.search(rf'<c r="M{row}"( s="\d+")[^>]*?(?:/>|><is><t>(.*?)</t></is></c>)', pro, re.S); assert m, row
    old = m.group(2) or ''
    pro = pro[:m.start()] + f'<c r="M{row}"{m.group(1)} t="inlineStr"><is><t>{old}{escape(extra)}</t></is></c>' + pro[m.end():]

tmp = 'tracker_new.xlsx'
with zipfile.ZipFile(tmp, 'w', zipfile.ZIP_DEFLATED) as zo:
    for i in zin.infolist():
        if i.filename == 'xl/worksheets/sheet4.xml': zo.writestr(i, lit.encode())
        elif i.filename == 'xl/worksheets/sheet1.xml': zo.writestr(i, pro.encode())
        else: zo.writestr(i, zin.read(i.filename))
zin.close()
shutil.move(tmp, src)
z = zipfile.ZipFile(src); print(z.testzip())
y = z.read('xl/worksheets/sheet4.xml').decode()
for r in (2, 3, 4, 5, 6, 10, 11):
    print(r, re.sub(r'<[^>]+>', '', re.search(rf'<c r="H{r}".*?</c>', y, re.S).group(0))[:90])
