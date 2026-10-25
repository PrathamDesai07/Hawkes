"""One-off: Problem 18 (row 19) status and notes."""
import zipfile, re, shutil
from xml.sax.saxutils import escape
src = 'hawkes_research_problem_tracker.xlsx'
zin = zipfile.ZipFile(src); x = zin.read('xl/worksheets/sheet1.xml').decode()
note = (" 25 Oct 2026. problem18_bias_vs_K.md; code/p18/ (fd18.py, e2..e8.py). Draft, not independently verified. "
 "PROVED reductions: (i) n*-n-dagger = (mu-dagger - mu*)/Lambda* (rate matching, Problem 16); (ii) the surrogate pseudo-truth is the weighted L2 pole approximant of G(s)=int nu/(beta+s) on the imaginary axis with weight Lambda*/|1-G|^2 (Hawkes spectrum), so the bias is the DC value of the approximation error and can be computed deterministically for any K (matches the simulations of Problem 16 to 0.004). "
 "CONJECTURE 18.4 (numerically supported): b_K = kappa(eta,r) * n* (tau beta_1)^eta / Gamma(1+eta), beta_1 = beta_top r^-(K-1), i.e. geometric decay in K at rate eta*ln r; kappa ~0.43 at eta=1/2 (0.7 at eta=0.15), independent of n* and of criticality; exponent confirmed over 3 decades of beta_1 (ratios 0.32 per decade vs 10^-0.5); fixed range saturates at the truncation floor kappa*n_low(beta_min). "
 "Also: coarse grids overshoot (bias changes sign; best ratio for fixed K: r~17 at K=4, 8 at K=6, 6 at K=8); KL and surrogate biases agree within ~0.002 for K>=3. "
 "NOT proved: the exponent eta and kappa>0 (Wiener-Hopf/Newman-Stahl type problem for a geometric pole ladder), the discretization-floor law ~exp(-pi^2/ln r), general nu beyond the gamma form.")
row=19
m = re.search(rf'<c r="L{row}"( s="\d+")[^>]*?>.*?</c>', x, re.S); assert m
x = x[:m.start()] + f'<c r="L{row}"{m.group(1)} t="inlineStr"><is><t>Done</t></is></c>' + x[m.end():]
m = re.search(rf'<c r="M{row}"( s="\d+")[^>]*?(?:/>|><is><t[^>]*>(.*?)</t></is></c>)', x, re.S); assert m
x = x[:m.start()] + f'<c r="M{row}"{m.group(1)} t="inlineStr"><is><t>{m.group(2) or ""}{escape(note)}</t></is></c>' + x[m.end():]
tmp='tracker_new_p18.xlsx'
with zipfile.ZipFile(tmp,'w',zipfile.ZIP_DEFLATED) as zo:
    for i in zin.infolist(): zo.writestr(i, x.encode() if i.filename=='xl/worksheets/sheet1.xml' else zin.read(i.filename))
zin.close(); shutil.move(tmp,src); print('ok')
