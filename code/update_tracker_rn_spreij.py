"""One-off: record Rong-Nair overlap resolution and the Spreij-for-Hawkes corollary."""
import zipfile, re, shutil
from xml.sax.saxutils import escape
src = 'hawkes_research_problem_tracker.xlsx'
zin = zipfile.ZipFile(src)
s1 = zin.read('xl/worksheets/sheet1.xml').decode()
s4 = zin.read('xl/worksheets/sheet4.xml').decode()


def set_cell(x, ref, text, append=False):
    m = re.search(rf'<c r="{ref}"( s="\d+")[^>]*?(?:/>|><is><t[^>]*>(.*?)</t></is></c>)', x, re.S)
    assert m, ref
    old = (m.group(2) or '') if append else ''
    return x[:m.start()] + f'<c r="{ref}"{m.group(1)} t="inlineStr"><is><t>{old}{escape(text)}</t></is></c>' + x[m.end():]


def find_row(x, needle):
    for m in re.finditer(r'<row r="(\d+)"[^>]*>.*?</row>', x, re.S):
        if needle in m.group(0):
            return int(m.group(1))
    raise KeyError(needle)


s1 = set_cell(s1, 'M17', " NOVELTY CORRECTION 25 Oct 2026 (Rong-Nair PDF checked): the surrogate (least-squares) rate matching c* = Lambda(1-Gamma*) and the spectral form of the LS pseudo-truth are KNOWN (their Lemma 4(b), Lemma 13, Thm 15, eq. 5.7). Ours: the KL rate matching and n-dagger<1 (their Gamma*<1 is unproved), binding criterion, active-constraint limit law, strong-convexity constant. Their CLT-without-mixing uses a linear martingale representation of lambda via the Hawkes resolvent and does not transfer to the KL score, so the mixing caveat stays.", append=True)
s1 = set_cell(s1, 'M19', " NOVELTY CORRECTION 25 Oct 2026: Prop 18.2 (frequency-domain weighted-L2 form with the Hawkes Bartlett spectrum) is KNOWN for least squares: Rong-Nair arXiv 2602.20795 Lemma 4(b), Lemma 13(b), Thm 15. Not claimed. The bias law in K for a power-law truth (Conjecture 18.4) is not in their paper (their numerics: sum-of-exponentials truth, Laguerre P=1..5, non-monotone, no rate).", append=True)
s1 = set_cell(s1, 'M18', " SPREIJ COROLLARY 25 Oct 2026 (refpapers_conclusions.md 5.2): Spreij 1991 Thm 3 applies to the Hawkes dictionary (conditions verified: G_t/t -> positive definite limit by the ergodic theorem and Wasserstein contraction; f(x)=x^2 growth condition holds), so the well-specified AML recursion is a.s. consistent at every fixed n<1 without the Problem 5 concentration machinery (no rate/regret/eps-explicit constants). Numerics code/p17/spreij_cond.py: transient lengthens as eps falls (condition number of G_t ~3e4 at n=0.97).", append=True)

r = find_row(s4, 'Hawkes Identification with a Prescribed Causal Basis')
s4 = set_cell(s4, f'H{r}', ' CHECKED 25 Oct 2026 in the PDF: Lemma 4(b) Bartlett spectrum Lambda/|1-phi0^|^2; Lemma 13 and eq. (5.2)-(5.5) spectral integrals for R*, R*^(1,0); Thm 15 alpha*=R*^-1 R*^(1,0), c*=Lambda(1-Gamma*). So the spectral weighted-L2 form of the LS pseudo-truth and the surrogate rate matching are known; Gamma*<1 not proved; CLT without mixing relies on the linear resolvent martingale representation (no transfer to the KL score).', append=True)
r = find_row(s4, 'Goodness-of-Fit Test for Mismatched Self-Exciting Processes')
s4 = set_cell(s4, f'H{r}', ' Further 25 Oct 2026: consistency and CLT are taken from modifications of Ogata (1978) Thms 2, 4, 5; only equality constraints appear (no KKT/inequality constraints); the piecewise-constant model is linear in the bin heights (same linear-in-parameter structure as ours). Proof of Lemma 1 itself not retrieved.', append=True)
tmp = 'tracker_new_rn.xlsx'
with zipfile.ZipFile(tmp, 'w', zipfile.ZIP_DEFLATED) as zo:
    for i in zin.infolist():
        data = zin.read(i.filename)
        if i.filename == 'xl/worksheets/sheet1.xml':
            data = s1.encode()
        if i.filename == 'xl/worksheets/sheet4.xml':
            data = s4.encode()
        zo.writestr(i, data)
zin.close()
shutil.move(tmp, src)
print('ok', zipfile.ZipFile(src).testzip())
