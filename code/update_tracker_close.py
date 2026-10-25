import zipfile, re, shutil
from xml.sax.saxutils import escape
src = 'hawkes_research_problem_tracker.xlsx'
zin = zipfile.ZipFile(src); x = zin.read('xl/worksheets/sheet1.xml').decode()
notes = {
 3: " CLOSED 25 Oct (Done): closed-form chain in K and grid ratio r: Cov(S) >= (E lam/4) G (G=1/(b_i+b_j), spectral), exact E S_k = mu/(eps b_k), lam_min(M) >= lam_min(C)/(1+2|m|^2), Cauchy inverse in closed form, then truncation with Lemma-1' moments; valid but 1e3-1e8 below truth (code/p2_closed_chain.py). Sharp scaling = three-tier structure (Problem 12), heuristic.",
 5: " CLOSED 25 Oct (Done): dense F-norm projection for any K (BVLS + bisection on the budget multiplier): matches brute force to 5e-14, non-expansive (0 violations in 2000 pairs), O(K^3) per inner solve (19 ms K=2, 97 ms K=5); the diagonal O(K log K) algorithm remains for the diagonal metric (code/p4_dense_projection.py).",
 6: " CLOSED 25 Oct (Done): Theorem HP (problem05_07_closure.md): with prob >= 1-6e^-x-P(E^c), R_n <= psi_max(2 i0+3x) + 4[D0' + x/theta + G_n(x)]; plus expected-regret Thm R (sharp constant). Assumptions: known compact convex Theta_c, leverage event, Assumption M (certified K=1) / lambda_L (K>=2).",
 8: " CLOSED 25 Oct (Done): Lemma JD (problem05_07_closure.md): joint zero mean and covariance E dJ (non-degenerate since b not parallel to a), pathwise |s|^2_{F^-1} <= (2 lbar/mu)(1+E^2), tail e^-m, cumulative exponential bound (E5b); 0 violations in 1.44e6 gaps.",
}
for row, t in notes.items():
    m = re.search(rf'<c r="M{row}"( s="\d+")[^>]*?(?:/>|><is><t[^>]*>(.*?)</t></is></c>)', x, re.S); assert m, row
    x = x[:m.start()] + f'<c r="M{row}"{m.group(1)} t="inlineStr"><is><t>{m.group(2) or ""}{escape(t)}</t></is></c>' + x[m.end():]
    m = re.search(rf'<c r="L{row}"( s="\d+")[^>]*?>.*?</c>', x, re.S); assert m, row
    x = x[:m.start()] + f'<c r="L{row}"{m.group(1)} t="inlineStr"><is><t>Done</t></is></c>' + x[m.end():]
tmp = 'tracker_new.xlsx'
with zipfile.ZipFile(tmp, 'w', zipfile.ZIP_DEFLATED) as zo:
    for i in zin.infolist():
        zo.writestr(i, x.encode() if i.filename == 'xl/worksheets/sheet1.xml' else zin.read(i.filename))
zin.close(); shutil.move(tmp, src)
y = zipfile.ZipFile(src).read('xl/worksheets/sheet1.xml').decode()
print(zipfile.ZipFile(src).testzip(), [re.search(rf'<c r="L{r}"[^>]*><is><t>(.*?)</t>', y).group(1) for r in (3, 5, 6, 8)])
