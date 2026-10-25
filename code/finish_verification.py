import io, zipfile, re, shutil
from xml.sax.saxutils import escape

def rd(p): return io.open(p, encoding='utf8').read()
def wr(p, s): io.open(p, 'w', encoding='utf8', newline='\n').write(s)

s = rd('verification_and_complexity.md')
a = "pays only for K ≳ 100 in NumPy (measured: dense faster for d = 3–51, 1.14× slower→faster at d = 101); in compiled code the crossover is d ≈ q"
assert a in s
s = s.replace(a, "pays only for K ≳ 100 in NumPy (measured: dense faster for d = 3–51, Sherman–Morrison 1.14× faster at d = 101); in compiled code the crossover should be near d ≈ q")
wr('verification_and_complexity.md', s)

p3 = rd('latex/part3_tracking_review.tex')
marker = "\\subsection{Open problems, in order of value}"
assert marker in p3
new = r"""\subsection{Verification record and complexity}\label{sec:verif}
\textbf{Algebraic identities.} Thirty-four identities behind the lemmas (generator computations for $Z_c$, $Z^p$, $e^{\theta Z}$, $S_iS_jS_l$, $e^{-u\lambda}$, $e^{\theta S}$; the Perron objects $\omega,\gamma,\tilde\nu,\kappa$; Lemma~\ref{lem:H}(i)--(iv); Theorem~\ref{thm:A} step 1; the Fisher reduction, efficient information and ridge Rayleigh quotient; the identity behind Lemma~\ref{lem:C1}; the tanh and sinusoid laws) were verified symbolically with exact derivatives of the generator (\texttt{code/verify\_symbolic.py}). \textbf{Inequalities.} Eighteen families (Lemma~\ref{lem:moments} pathwise, Lemma~\ref{lem:expmoment}, Lemma~\ref{lem:G}, the elementary bounds of Lemma~\ref{lem:elem}, Lemma~\ref{lem:gammafact}, Lemma~\ref{lem:F2}, the Bernstein/MGF inequalities, the $\psi$ expansion, the tangent minorant, the Gamma tail, the path-length bounds, the entropy step) had no violation over $10^4$--$10^5$ random cases each (\texttt{code/verify\_inequalities.py}). The probabilistic theorems are chains of these steps with standard tools; they were re-derived step by step but cannot be machine-checked, and Theorem~\ref{thm:V} and the input (H3) of Theorem~\ref{thm:R} are the best targets for independent review.

\textbf{Complexity.} The estimator needs $O(K^2)$ memory and $O(K)$ state; per event the quadrature costs $O(qK^2)$ and the Newton direction $O(K^3)$ (or $O(qK^2)$ with rank-one updates of $F^{-1}$, which pays only for $K\gtrsim100$ in NumPy). The dense $F$-projection, needed by the regret proof, was improved from $\approx55$--$90$ box-QP solves (bisection on the budget multiplier) to $2$--$7$ (safeguarded Newton on the multiplier): $15$--$30\times$ faster, same solution to $1.7\times10^{-13}$. The moment systems of order $2$ and $3$ are Kronecker sums of $A=\mathrm{diag}(\beta)-\mathbf 1\alpha\T$ and are solved by diagonalising $A$: $O(K^4)$ time and $O(K^3)$ memory instead of $O(K^9)$ and $O(K^6)$ ($K=40$ in $21$ ms where the dense system would need $33$~GB). The tangent bound is concave in $t=1/\lambda_0$ ($N(t)=2tM-t^2M_3$ is operator-concave), so a Brent search with $13$--$26$ evaluations replaces a $75$-point grid. See \texttt{verification\_and\_complexity.md}.

"""
p3 = p3.replace(marker, new + marker, 1)
wr('latex/part3_tracking_review.tex', p3)

src = 'hawkes_research_problem_tracker.xlsx'
zin = zipfile.ZipFile(src); x = zin.read('xl/worksheets/sheet1.xml').decode()
notes = {5: " COMPLEXITY 25 Oct (verification_and_complexity.md): dense F-projection now uses safeguarded Newton on the budget multiplier (2-7 inner box-QP solves instead of 55-90; x15-30 faster, same solution to 1.7e-13). Estimator memory O(K^2), state O(K), per-event O(qK^2 + K^3).",
         3: " COMPLEXITY 25 Oct: moments of order <=3 by Kronecker-sum diagonalisation, O(K^4) time / O(K^3) memory (was O(K^9) / O(K^6)); lambda0-search by Brent (concavity in 1/lambda0), 13-26 evaluations (was 75)."}
for row, t in notes.items():
    m = re.search(rf'<c r="M{row}"( s="\d+")[^>]*?(?:/>|><is><t[^>]*>(.*?)</t></is></c>)', x, re.S); assert m, row
    x = x[:m.start()] + f'<c r="M{row}"{m.group(1)} t="inlineStr"><is><t>{m.group(2) or ""}{escape(t)}</t></is></c>' + x[m.end():]
tmp = 'tracker_new.xlsx'
with zipfile.ZipFile(tmp, 'w', zipfile.ZIP_DEFLATED) as zo:
    for i in zin.infolist(): zo.writestr(i, x.encode() if i.filename == 'xl/worksheets/sheet1.xml' else zin.read(i.filename))
zin.close(); shutil.move(tmp, src); print('ok', zipfile.ZipFile(src).testzip())
