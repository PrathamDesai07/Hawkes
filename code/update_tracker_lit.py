"""One-off: Problems 22 and 23 status/notes (Problems sheet) and Literature sheet updates."""
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


# ---- Problems sheet: 22 -> row 23, 23 -> row 24
n22 = (" 25 Oct 2026 (second pass). literature_difference_statements.md: one 'difference from' paragraph for each core paper "
       "(Herrera-Marin I and II, Hall-Willett, Yang et al., Piggott-Solo + Pasha-Solo) and for the papers that bear on Problems 16-18 "
       "(Rong-Nair 2602.20795, Wei-Zhu-Zhang-Xie 2006.09439, Spreij 1986/1991), plus a related-work draft. Each paragraph states its access level. "
       "KEY FINDINGS: the generic KL pseudo-truth/sandwich content of Problem 16 is anticipated in batch form by Wei et al. 2021; the LS pseudo-truth and a CLT without mixing "
       "by Rong-Nair 2026 (possible overlap with Prop 18.2); well-specified recursive estimation by Spreij. "
       "STILL MISSING: full texts of Piggott-Solo and Pasha-Solo (abstract only), Spreij scans (image-only, not OCRed), Herrera-Marin II Thm 3.1, end of Rong-Nair. "
       "Statement text should be rechecked against those before submission.")
s1 = set_cell(s1, 'L23', 'Done')
s1 = set_cell(s1, 'M23', n22, append=True)
n23 = (" 25 Oct 2026. literature_unverified_items.md. Searched/fetched all listed items. Zhang 2023: no matching paper found (identity unresolved). "
       "Wang 2023 = Technometrics online CUSUM change-point detector (not an estimator). Online EM for Hawkes: none found; Zhou-Zha-Song is batch ADMM; "
       "Lemonnier-Vayatis attribution unconfirmed. Cartea-Cohen-Labyad: batch LS gradient sampling, no regret. Bacry et al.: no misspecification paper located. "
       "Citation trees screened by titles (Hall-Willett via Semantic Scholar: 64; Yang et al. via OpenAlex: 26; a second Semantic Scholar query hit HTTP 429). "
       "Found beyond the list: Wei-Zhu-Zhang-Xie AISTATS 2021 (KL pseudo-truth, QMLE sandwich), Rong-Nair 2026 (LS pseudo-truth, CLT), Spreij 1986/1991 "
       "(recursive estimators for counting processes), Pasha-Solo CDC 2019 (fixed-gain), Lima 2021 (eps-margin stabilisation). "
       "Overlap-risk column updated. UNRESOLVED: Zhang 2023, online EM, full texts of Spreij/Solo papers.")
s1 = set_cell(s1, 'L24', 'Done')
s1 = set_cell(s1, 'M24', n23, append=True)


# ---- Literature sheet: update existing rows by matching paper text
def find_row(x, needle):
    for m in re.finditer(r'<row r="(\d+)"[^>]*>.*?</row>', x, re.S):
        if needle in m.group(0):
            return int(m.group(1))
    raise KeyError(needle)


r = find_row(s4, 'Hawkes Identification with a Prescribed Causal Basis')
s4 = set_cell(s4, f'G{r}', 'High')
s4 = set_cell(s4, f'H{r}', 'Read 25 Oct 2026 (first ~200k of 213k characters of the HTML). Batch LS on a prescribed unit-mass basis; pseudo-true parameters via the LS contrast '
              '(Section 5.3, Thm 15) and CLTs under correct spec (Thm 17) and misspecification (Thm 21) by a martingale route with no mixing; numerics P=1..5 vs a '
              'sum-of-exponentials truth (no rate, no power law). No KL, no online, no proof of Gamma*<1. May anticipate Prop 18.2 and the surrogate rate matching; '
              'route for the CLT item. See literature_difference_statements.md.')
r = find_row(s4, 'Non-Negative Online Estimation')
s4 = set_cell(s4, f'H{r}', ' Also: Pasha and Solo, CDC 2019, Adaptive Estimation of a Vector Point Process Model (abstract: fixed-gain adaptive estimator on the '
              'instantaneous likelihood, non-negativity constraints, genomic application). Full texts of both still not obtained.', append=True)
for needle, new_h, new_g in [
    ('Zhang et al. 2023 online estimator',
     'Searched 25 Oct 2026 (4 queries): no matching 2023 paper found; identity of the cited work unresolved. Wang et al. 2023 (Technometrics 65(1):44-56) is an online '
     'recursive CUSUM change-point detector, not an estimator. See literature_unverified_items.md.', 'Unknown'),
    ('Online EM for Hawkes',
     'Searched 25 Oct 2026: no online EM for Hawkes found. Zhou-Zha-Song (AISTATS 2013, PMLR 31:641-649) is batch ADMM (ADM4). Lemonnier-Vayatis attribution NOT confirmed '
     '(search returned Lemonnier-Scaman-Kalogeratos, low-rank large-scale inference).', 'Low'),
    ('Bacry, Muzy, Delattre, Hoffmann',
     'Searched 25 Oct 2026: Cartea-Cohen-Labyad (2111.10637) is a batch gradient method for the least-squares contrast (adaptive stratified sampling of the gradient); '
     'no regret/online/misspecification in the abstract. No Bacry et al. misspecification paper located.', 'Low'),
    ('Citation trees of Yang',
     'Screened 25 Oct 2026 by titles: Hall-Willett tree 64 entries (Semantic Scholar), Yang et al. tree 26 entries (OpenAlex); second Semantic Scholar query rate-limited. '
     'Found: Pasha-Solo 2019, Lima 2021 (eps-margin stabilisation of parametric Hawkes MLE), Ghassemi et al. 2022 (online mixtures of MHP), Jaffard et al. 2024 '
     '(learning rule for Hawkes networks), Xu et al. 2018. None on misspecification or the near-critical ratio by title.', 'Low')]:
    r = find_row(s4, needle)
    s4 = set_cell(s4, f'H{r}', new_h)
    s4 = set_cell(s4, f'G{r}', new_g)

# ---- add new rows before the footer note
rows = list(re.finditer(r'<row r="(\d+)"[^>]*>.*?</row>', s4, re.S))
foot = rows[-1]
assert 'Details come from' in foot.group(0)
last_data = int(rows[-2].group(1))
nxt = last_data + 1


def mk(r_, order, paper, auth, link, probs, extract, risk, status):
    def c(col, st, t):
        return f'<c r="{col}{r_}" s="{st}" t="inlineStr"><is><t>{escape(t)}</t></is></c>'
    return (f'<row r="{r_}" ht="46.5" customHeight="1" s="27"><c r="A{r_}" s="29" t="inlineStr"><is><t>{escape(order)}</t></is></c>'
            + c('B', 31, paper) + c('C', 30, auth) + c('D', 30, link) + c('E', 29, probs) + c('F', 30, extract)
            + c('G', 29, risk) + c('H', 32, status) + '</row>')


new = [
    ('25', 'Goodness-of-Fit Test for Mismatched Self-Exciting Processes', 'Wei, Zhu, Zhang, Xie, AISTATS 2021 (PMLR 130:1243-1251)', 'arXiv 2006.09439', '16, 17',
     'QMLE theory under misspecification: KL pseudo-truth by concavity (Prop 2), a.s. convergence (Lemma 1), sandwich covariance.', 'Medium',
     'Read 25 Oct 2026 (about 100k of 136k characters of ar5iv). Piecewise-constant triggering, batch, near-critical not analysed. Anticipates the generic part of '
     'Theorem 16 / Prop 16.7; not the exponential dictionary with stability set, rate matching, active constraint, online.'),
    ('26', 'Recursive parameter estimation for counting processes with linear intensity; Recursive approximate maximum likelihood estimation for a class of counting process models',
     'Spreij, Stochastics 18 (1986) 277-312; J. Multivariate Anal. 39 (1991) 236-245', 'CWI R 8513 / ir.cwi.nl pub 2491, 6311; staff.fnwi.uva.nl/p.j.c.spreij', '17',
     'Recursive MLE-type estimators for intensities linear in theta: consistency and asymptotic normality by stochastic Lyapunov functions.', 'High',
     'Abstract-level only (25 Oct 2026); the scans are image-only and were not OCRed. Likely anticipates the well-specified recursive consistency/CLT of Problem 17; '
     'misspecification, active constraint, regret, near-critical not known to be covered. Read first.'),
    ('27', 'Adaptive Estimation of a Vector Point Process Model', 'Pasha and Solo, CDC 2019, pp. 2904-2909', 'IEEE Xplore (dblp conf/cdc/PashaS19)', '13, 15, 21',
     'Fixed-gain adaptive estimator on the instantaneous likelihood of a log-intensity model with non-negativity constraints.', 'Medium',
     'Abstract only (25 Oct 2026). Companion of the Piggott-Solo row.'),
    ('28', 'Fine-Grained epsilon-Margin Closed-Form Stabilization of Parametric Hawkes Processes', 'Lima, arXiv, May 2021', 'arXiv 2105.03800', '17',
     'Uses an epsilon margin on the stability condition to stabilise parametric MLE.', 'Low',
     'Abstract only (25 Oct 2026). Anticipates the device n<=1-eps, not the analysis.'),
    ('29', 'Gradient-based estimation of linear Hawkes processes with general kernels', 'Cartea, Cohen, Labyad, arXiv Nov 2021', 'arXiv 2111.10637', '5, 21',
     'Adaptive stratified sampling of the least-squares gradient for large datasets.', 'Low', 'Abstract only (25 Oct 2026). Batch LS, no regret.'),
    ('30', 'Online recursive CUSUM for Hawkes networks (Wang et al.)', 'Wang, Xie, Xie, Cuozzo, Mak, Technometrics 65(1):44-56 (2023)', 'DOI 10.1080/00401706.2022.2054862', '14, 20',
     'Online change-point detection for Hawkes networks.', 'Low', 'Abstract only (25 Oct 2026). Detection, not estimation.'),
    ('31', 'Provable local learning rule by expert aggregation for a Hawkes network', 'Jaffard, Vaiter, Muzy, Reynaud-Bouret, AISTATS 2024', 'arXiv 2304.08061', '-',
     'Hawkes network used as a classifier learned by expert aggregation; local regret bounds.', 'Low', 'Abstract only (25 Oct 2026). Not parametric estimation.'),
    ('32', 'Online Learning for Mixture of Multivariate Hawkes Processes', 'Ghassemi et al., ICAIF 2022', 'arXiv 2208.07961', '5',
     'Streaming mixtures of MHP.', 'Low', 'Title/abstract-level only (25 Oct 2026).'),
]
add = ''.join(mk(nxt + i, *t) for i, t in enumerate(new))
newfoot = nxt + len(new)
foot_new = foot.group(0).replace(f'<row r="{foot.group(1)}"', f'<row r="{newfoot}"').replace(f'r="B{foot.group(1)}"', f'r="B{newfoot}"')
foot_new = foot_new.replace('most new entries were read only in part (see Read status).',
                            'most new entries were read only in part (see Read status). Updated again 25 Oct 2026 (second pass): rows 25-32 added; Rong-Nair risk raised to High; '
                            'the four to-check rows filled in (literature_unverified_items.md).')
s4 = s4[:foot.start()] + add + foot_new + s4[foot.end():]
s4 = re.sub(r'<dimension ref="A1:H\d+"', f'<dimension ref="A1:H{newfoot}"', s4)
tmp = 'tracker_new_lit.xlsx'
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
