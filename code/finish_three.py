import io, zipfile, re, shutil
from xml.sax.saxutils import escape

def rd(p): return io.open(p, encoding='utf8').read()
def wr(p, s): io.open(p, 'w', encoding='utf8', newline='\n').write(s)

m = rd('latex/hawkes_master.tex')
if 'part4_resolutions' not in m:
    m = m.replace("\\input{part3_tracking_review}", "\\input{part4_resolutions}\n\\input{part3_tracking_review}", 1)
wr('latex/hawkes_master.tex', m)
t = rd('latex/part3_tracking_review.tex')
t = t.replace(r"Spectrum: three tiers & \Pt{}/\Nt{} & all tiers bounded below (Thm.~\ref{thm:tangent}); upper bounds only for the ridge; eigenvector convergence numerical\\",
              r"Spectrum: three tiers & \Ps{} & all tiers bounded above and below (Thm.~\ref{thm:tangent}, Thm.~\ref{thm:tiers}); eigenvector convergence rate $O(\eps/g)$; constants loose; $k_K>1$\\")
t = t.replace(r"Tracking formulas for the linear filter & \Pt{} & relevance to the estimators is \Nt{}\\",
              r"Tracking formulas for the linear filter & \Pt{} & --\\" "\n" r"Tracking theorem for forgetting Newton (Thm.~\ref{thm:track}) & \Ps{} & window hypothesis (W); factor $4$ loose; jump class not tight\\" "\n" r"Amplitude effect: $O(A^2)$ closed form (scalar) & \Pt{} & matrix ($\alpha$ unknown) case numerical skeleton only\\")
t = t.replace(r"Divergence for $k_K<1$ (lower bound) & \Pt{} & no matching upper bound / estimator\\",
              r"Divergence of expected information, $k_K<1$ & \Pt{} & does not imply consistency\\" "\n" r"Matching upper bound $\mathcal I_T^{-1/2}$ (Thm.~\ref{thm:UB}) and non-consistency at $\eps=c/T$ for every $k$ (Thm.~\ref{thm:NC}) & \Pt{}/\Ps{} & $K=1$, $\alpha$ known; NC needs $\TV(\pi_\theta,\pi_{\theta'})$\\")
wr('latex/part3_tracking_review.tex', t)

src = 'hawkes_research_problem_tracker.xlsx'
zin = zipfile.ZipFile(src); x = zin.read('xl/worksheets/sheet1.xml').decode()
notes = {13: " THREE-TIER THEOREM 25 Oct (problem12_upper_bounds.md): fast modes y_j = sum alpha_k S_k/(beta_k - a_j) satisfy A y_j = -a_j y_j + mu exactly; Lyapunov y_j^2/W gives E[y_j^2/lam] = O(1) uniformly in eps => middle eigenvalues bounded above and below (both eps-independent), top Theta(1/eps), ridge O(eps); smallest eigenvector -> ridge direction at rate O(eps/g) (Davis-Kahan); symbolic + pathwise checks pass; constants ~20x loose.",
         12: " MATCHING RATE / NON-CONSISTENCY 25 Oct (problem11_observed_information.md): Theorem UB: |mu_hat - mu| <= C sqrt(loglog I_T)/sqrt(I_T) on {I_T >= M0} (observed information, any k); Theorem NC: at eps=c/T NO estimator is consistent for ANY k (stopped-KL bounded pathwise on an event of probability >= p > 0), so divergence of expected information for k<1 does not give consistency; simulations: P(|err|>mu/2) = 0.40 (k=0.6), 0.34 (k=2), flat in T.",
         15: " TRACKING THEOREM 25 Oct (problem14_tracking_theorem.md): pathwise recursion D_i <= (1+rho)[(1-eta)D_{i-1} + X - 2<g,v> + G] + (1+1/rho)||d||^2_F, E[X - 2<g,v>] <= E Xi; steady state E||v||^2_I <~ 2 eta d + 4(||d||_I/eta)^2 (modulo window hypothesis W); checked in 57,616 gaps (0 violations; measured 0.0077 vs bound 0.0463). AMPLITUDE EFFECT derived: O(A^2) gain-modulation with closed forms C2_fixed, C2_adaptive (sympy), verified vs exact ODE to 4 digits; alpha-unknown matrix case: deterministic skeleton reproduces adaptive/fixed 1.11, 1.33, 1.67 (MC 1.32, 1.48, 1.75)."}
for row, tx in notes.items():
    mm = re.search(rf'<c r="M{row}"( s="\d+")[^>]*?(?:/>|><is><t[^>]*>(.*?)</t></is></c>)', x, re.S); assert mm, row
    x = x[:mm.start()] + f'<c r="M{row}"{mm.group(1)} t="inlineStr"><is><t>{mm.group(2) or ""}{escape(tx)}</t></is></c>' + x[mm.end():]
tmp = 'tracker_new.xlsx'
with zipfile.ZipFile(tmp, 'w', zipfile.ZIP_DEFLATED) as zo:
    for i in zin.infolist(): zo.writestr(i, x.encode() if i.filename == 'xl/worksheets/sheet1.xml' else zin.read(i.filename))
zin.close(); shutil.move(tmp, src); print('tracker ok', zipfile.ZipFile(src).testzip())

with io.open('final_audit.md', 'a', encoding='utf8', newline='\n') as f:
    f.write("""

---
## Addendum 2 (25 Oct 2026): the three items that were "not strengthened"
- **Problem 12:** now a theorem (k_K > 1): middle tier bounded above by C_mid (eps-independent, Lyapunov argument on the exact fast modes), top Theta(1/eps), ridge O(eps), smallest eigenvector -> ridge at rate O(eps/g). Constants about 20x loose; eigenvector bound informative for eps <~ 1e-3 with the rigorous E[1/lambda] bound.
- **Problem 14:** tracking recursion and bound proved modulo the window hypothesis (W); checked in 57,616 gaps; amplitude effect derived in closed form for the scalar case and reproduced by a two-parameter deterministic skeleton.
- **Problem 11:** the matching upper bound exists in terms of the observed information (Theorem UB). New finding: at eps = c/T no estimator is consistent for ANY k (Theorem NC, K = 1): the expected-information divergence for k < 1 does not translate into consistency. This changes the headline of the information dichotomy: it concerns the typical accuracy (bounded vs growing like T^{1-k}), not consistency.
""")
print('done')
