"""One-off: record the reading of the four reference papers (Spreij 1986/1991, Piggott-Solo 2018, Pasha-Solo 2019)."""
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


# Problems sheet (ID n -> row n+1)
s1 = set_cell(s1, 'M18', " READ 25 Oct 2026 (refpapers_conclusions.md): the EF recursion is Spreij's AML recursion (Stochastics 1986; JMVA 1991) applied to the Hawkes regressor; his well-specified a.s. consistency (Lyapunov x'Q^-1 x, Lai-Wei growth condition, known coordinate box) and Cramer-Rao CLT are the A=I special case of Theorems 17A/17B, which are therefore NOT claimed as new. "
              "Still ours: misspecified target theta-dagger, drift matrix I_F^-1 H != I and the Lyapunov-equation limit law (+ Lemma 17.5), observed-information preconditioner, exact convex projection with a binding stability face / alpha_k=0, active-constraint limit law, regret. Overlap risk: High for the well-specified part, Low-Medium for the rest.", append=True)
s1 = set_cell(s1, 'M6', " ATTRIBUTION (25 Oct 2026): the Fisher-preconditioned recursion of this problem is Spreij's AML recursion (1986/1991) applied to Hawkes; novelty of C1 lies in the log-loss regret analysis, the exact F-norm projection, the observed-information variant and the misspecified analysis, not the recursion. See refpapers_conclusions.md.", append=True)
s1 = set_cell(s1, 'M23', " UPDATE 25 Oct 2026: Piggott-Solo (linear Hawkes, least squares, multiplicative non-negative update, no theorem, stability not enforced), Pasha-Solo (log-linear model, likelihood steepest ascent, heuristic averaging, unproved PD of E[phi phi' lambda]) and Spreij 1986/1991 now read in full; paragraphs corrected in literature_difference_statements.md. Main consequence: Spreij anticipates the well-specified recursion and its CLT.", append=True)
s1 = set_cell(s1, 'M24', " UPDATE 25 Oct 2026: Spreij 1986/1991, Piggott-Solo 2018 and Pasha-Solo 2019 read in full (refpapers_conclusions.md). Spreij confirmed as the key prior art for Problem 17 (well-specified part); the Solo papers have no theorems.", append=True)

# Literature sheet
r = find_row(s4, 'Recursive parameter estimation for counting processes with linear intensity')
s4 = set_cell(s4, f'G{r}', 'High')
s4 = set_cell(s4, f'H{r}', 'READ IN FULL 25 Oct 2026 (both papers; equations from page images). Intensity lambda=theta\'phi with observed predictable phi, well-specified. AML recursion = our EF Newton in continuous time; 1986 Thm 3.2 and 1991 Thm 3 give a.s. consistency (Lyapunov x\'Q^-1 x; Lai-Wei growth; known box [eps,1/eps]); 1986 Thm 5.4 / 1991 Thm 7 give the Cramer-Rao CLT. No Hawkes example, misspecification, regret, near-critical. Anticipates the well-specified parts of Thms 17A/17B and the algorithm of Problem 5. See refpapers_conclusions.md.')
r = find_row(s4, 'Non-Negative Online Estimation')
s4 = set_cell(s4, f'G{r}', 'Medium')
s4 = set_cell(s4, f'H{r}', 'READ IN FULL 25 Oct 2026. Linear Hawkes, Laguerre basis (one shared decay constant), exponentially weighted LEAST SQUARES via a one-step NMF multiplicative update, fixed gain, discrete bins; simulations only; stability not enforced; no theorem. Claims projections can stall. Pasha-Solo (CDC 2019) read in full too: log-linear model, likelihood steepest ascent, heuristic averaging analysis assuming well-specification. See refpapers_conclusions.md.')
r = find_row(s4, 'Adaptive Estimation of a Vector Point Process Model')
s4 = set_cell(s4, f'G{r}', 'Low')
s4 = set_cell(s4, f'H{r}', 'READ IN FULL 25 Oct 2026. Log-linear (nonlinear) Hawkes, no sign constraints, fixed-gain likelihood steepest ascent, heuristic averaging stability sketch (stationary, well-specified; PD of E[phi phi\' lambda] and the hovering theorem left unproved). Not the linear model of this project; no theorem on misspecification, regret or tracking.')
tmp = 'tracker_new_ref.xlsx'
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
