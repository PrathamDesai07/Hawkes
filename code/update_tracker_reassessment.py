"""One-off: append reassessment pointers to the Notes of problems 5, 9, 11 and 12 (rows 6, 10, 12, 13); statuses untouched."""
import zipfile, re, shutil
from xml.sax.saxutils import escape

src = 'hawkes_research_problem_tracker.xlsx'
zin = zipfile.ZipFile(src)
x = zin.read('xl/worksheets/sheet1.xml').decode()
rows = {6: 5, 10: 9, 12: 11, 13: 12}
for r, i in rows.items():
    assert re.search(rf'<c r="A{r}"[^>]*><v>{i}</v>', x), (r, i)

extra = {
 6: " REASSESSMENT 25 Oct (novelty_reassessment.md): C2's 'C log(1/eps)' should be made explicit: conjecture regret = (d/2) log T + (1/2) log det I(theta*;eps) + O(1) in time units, with (1/2) log det I ~ ((1-k)^+/2) log(1/eps), k=2*mu*m (m=mean lifetime). Exact K=1 slopes confirm the information side (0.40 for k0=0.2); the regret side is only consistent with it so far (underpowered experiment).",
 10: " REASSESSMENT 25 Oct: as written this problem (two-point bound) is prior art (Herrera-Marin Thm 2.2). Suggested replacement: (a) regret lower bound for sequential estimators >= (d/2) log T + (1/2) log det I - O(1); (b) identifiability map / matching UPPER bound on lambda_min(I) for the K-dictionary. Neither is proved; literature on regret lower bounds for point processes not yet searched.",
 12: " REASSESSMENT 25 Oct: divergence half = Herrera-Marin Remark 4.8; he leaves the converse (kappa>1) unproved. Other session's Theorem A (K=1, stationary) claims it with rates (I have NOT verified its proof). Exact slopes of (1/2) log det I vs log(1/eps): 0.347,0.397,0.404,0.403,0.402 -> (1-k0)/2=0.40 for k0=0.2; -> 0 for k0>=2. K=2 dictionary: A_T matches the K=1 value at equal k=2*mu*m (1.29/1.53/1.79 vs 1.26/1.55/1.78): supports threshold 2*mu*sum(w_k/beta_k)=1. Observation model matters (stationary vs empty start).",
 13: " REASSESSMENT 25 Oct: Herrera-Marin II (2609.38267) Sec 5 reports numerically an 'almost flat direction' from K~5 when fitting sums of exponentials to a long-memory kernel (slow poles trade mass with the baseline; not proved, misspecified setting). Same mechanism as my Problem 2 for the well-specified dictionary (lambda_min(I) falls ~10x per added exponential at r=2). Only a LOWER bound on lambda_min is proved here (Theorem T); a matching UPPER bound as K or the grid ratio grows would turn the flat direction into a theorem.",
}

def append(x, r, text):
    m = re.search(rf'<c r="M{r}"( s="\d+")[^>]*?(?:/>|><is><t[^>]*>(.*?)</t></is></c>)', x, re.S); assert m, r
    old = m.group(2) or ''
    return x[:m.start()] + f'<c r="M{r}"{m.group(1)} t="inlineStr"><is><t>{old}{escape(text)}</t></is></c>' + x[m.end():]

for r, t in extra.items():
    x = append(x, r, t)
tmp = 'tracker_new.xlsx'
with zipfile.ZipFile(tmp, 'w', zipfile.ZIP_DEFLATED) as zo:
    for i in zin.infolist():
        zo.writestr(i, x.encode() if i.filename == 'xl/worksheets/sheet1.xml' else zin.read(i.filename))
zin.close()
shutil.move(tmp, src)
z = zipfile.ZipFile(src); print(z.testzip())
y = z.read('xl/worksheets/sheet1.xml').decode()
for r in rows:
    print(r, re.search(rf'<c r="L{r}"[^>]*><is><t>(.*?)</t>', y).group(1), len(re.search(rf'<c r="M{r}"[^>]*><is><t>(.*?)</t>', y, flags=re.S).group(1)))
