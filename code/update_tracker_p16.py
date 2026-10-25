"""One-off: mark Problem 16 (row 17) Done with notes; leave other cells as found."""
import zipfile, re, shutil
from xml.sax.saxutils import escape
src = 'hawkes_research_problem_tracker.xlsx'
zin = zipfile.ZipFile(src); x = zin.read('xl/worksheets/sheet1.xml').decode()
note = (" 25 Oct 2026. problem16_kl_pseudotruth.md; code/p16/ (sim16.c, p16lib.py, t1..t7.py). Theorem 16 [P, draft, not independently verified]: under finite rate and E[lam*(1+|a|^2)]<inf, "
 "L(theta)=E[lam* log(lam*/lam_theta)-lam*+lam_theta] is convex, strictly convex for distinct beta_k (linear independence of exponentials on a gap), coercive on the cone, so a unique theta-dagger exists on Theta_eps and on the cone; KKT: "
 "dL/dmu=E[1-lam*/lam]=rho_0, dL/dalpha_k=E[S_k(1-lam*/lam)]=rho_k-kappa/beta_k. Rate matching (new): mu-dagger=(1-n-dagger)Lambda* minus multiplier terms, so n(theta-dagger)<1 automatically and the constraint n<=1-eps is active iff the unconstrained "
 "immigrant share mu/Lambda* < eps; when active mu-dagger=eps*Lambda*-kappa(1-eps). Quasi-MLE converges a.s. to theta-dagger. Numerics (power-law truth n*=0.8,0.95, K=1..4): rate matching to 4 digits, unique minimiser from several starts, "
 "KKT/kappa identities hold under active constraint, error of QMLE ~T^-1/2, Newton-decrement bound never violated. Found: n(theta-dagger)<n*, increasing with K (0.46,0.62,0.70,0.73 at n*=0.8) - input to Problem 18. "
 "KL and quadratic surrogate differ by only 0.1-0.3% in these examples, so the advantage is structural not numerical. NOT proved: explicit strong-convexity constant, sandwich variance, a priori mu-dagger>mu_min.")
for pid, status in [(16, 'Done')]:
    row = pid + 1
    m = re.search(rf'<c r="L{row}"( s="\d+")[^>]*?>.*?</c>', x, re.S); assert m
    x = x[:m.start()] + f'<c r="L{row}"{m.group(1)} t="inlineStr"><is><t>{status}</t></is></c>' + x[m.end():]
    m = re.search(rf'<c r="M{row}"( s="\d+")[^>]*?(?:/>|><is><t[^>]*>(.*?)</t></is></c>)', x, re.S); assert m
    x = x[:m.start()] + f'<c r="M{row}"{m.group(1)} t="inlineStr"><is><t>{m.group(2) or ""}{escape(note)}</t></is></c>' + x[m.end():]
tmp = 'tracker_new_p16.xlsx'
with zipfile.ZipFile(tmp, 'w', zipfile.ZIP_DEFLATED) as zo:
    for i in zin.infolist(): zo.writestr(i, x.encode() if i.filename == 'xl/worksheets/sheet1.xml' else zin.read(i.filename))
zin.close(); shutil.move(tmp, src); print('ok', zipfile.ZipFile(src).testzip())
y = zipfile.ZipFile(src).read('xl/worksheets/sheet1.xml').decode()
print(re.search(r'<c r="B17"[^>]*>.*?</c>', y, re.S).group(0)[:200]); print(re.search(r'<c r="L17"[^>]*>.*?</c>', y, re.S).group(0))
