"""One-off: record the 25 Oct 2026 literature search in the tracker.
 - Literature sheet (sheet4): correct row 7 (Piggott and Solo) and row 17 (CIR paper title/venue); insert rows 24-29 (new sources) before the footer, which moves to row 31;
   extend merge/dimension/validation/conditional-formatting ranges accordingly.
 - Problems sheet (sheet1): append short literature notes to the Notes of problems 10 (row 11) and 11 (row 12)."""
import zipfile, re, shutil
from xml.sax.saxutils import escape

src = 'hawkes_research_problem_tracker.xlsx'
zin = zipfile.ZipFile(src)
lit = zin.read('xl/worksheets/sheet4.xml').decode()
pro = zin.read('xl/worksheets/sheet1.xml').decode()

def set_cell(x, ref, text):
    m = re.search(rf'<c r="{ref}"( s="\d+")[^>]*?>(?:<is><t[^>]*>.*?</t></is>)?</c>', x, re.S)
    assert m, ref
    return x[:m.start()] + f'<c r="{ref}"{m.group(1)} t="inlineStr"><is><t>{escape(text)}</t></is></c>' + x[m.end():]

assert 'Non-Negative Online Estimation' in lit and 'CIR drift estimation theory' in lit
lit = set_cell(lit, 'C7', 'Piggott and Solo, ICASSP 2018, pp. 6333-6337 (listed earlier as "Solo")')
lit = set_cell(lit, 'H7', 'Reading: abstract read 25 Oct 2026 (via OpenAlex); full text NOT obtained. Fixed-gain exponentially weighted least squares on a Laguerre basis, non-negativity preserving, two time scales, simulations; abstract states no regret bound or near-critical analysis. See literature_search_report.md.')
lit = set_cell(lit, 'B17', 'Local asymptotic properties for Cox-Ingersoll-Ross process with discrete observations')
lit = set_cell(lit, 'C17', 'Ben Alaya, Kebaier, Tran, Scand. J. Stat. 47(4), 2020 (arXiv 1708.07070v5)')
lit = set_cell(lit, 'D17', 'arXiv 1708.07070')
lit = set_cell(lit, 'H17', 'Reading: intro, Thm 2.1, Cor 2.3 read 25 Oct 2026; proofs NOT read. Fixed sub/critical/supercritical regimes only; no local-to-critical (b=c/T) case. Thm 2.1 information matrix and Cor 2.3(i) critical coefficient match my K=1 Hawkes Fisher matrix and cold-start coefficient as eps->0 (lit_check_cir_fisher.py).')

new_rows = [
 (19, 'Two methods of estimation of the drift parameters of the Cox-Ingersoll-Ross process: continuous observations', 'Dehtiar, Mishura, Ralchenko, arXiv, May 2020', 'arXiv 2005.05262', '10, 11',
  'MLE of (a,b) defined only for 2a > sigma^2; stationary Gamma law of shape 2a/sigma^2; MLE formula uses int dr/r and int dt/r (denominator int r * int 1/r - T^2 mirrors det I_T).', 'Medium',
  'Reading: abstract, Prop 1 and Section 3 read 25 Oct 2026; proofs and the alternative estimator not read.'),
 (20, 'Near-Unit-Root Theory for Affine Processes', 'Anne, Lu, Yu, Zhou, arXiv, 27 May 2026', 'arXiv 2605.27859', '10, 11, 18',
  'CIR limit of local-to-unity discrete-time affine processes with a Feller remark (2mu >= sigma^2); OLS intercept inconsistent in the local-to-unity regime. No Hawkes (0 hits), OLS only, no likelihood information or lower bound.', 'Medium',
  'Reading: abstract, introduction, Thm 1 discussion and parts of Thm 3-4 read 25 Oct 2026 (71 pages, mostly unread).'),
 (21, 'Asymptotics of nearly critical Galton-Watson processes with immigration', 'Kevei, Acta Sci. Math. (Szeged) 77 (2011) 681-702', 'arXiv 1107.2440', '11',
  'Limit laws: Poisson/compound Poisson when offspring second factorial moment / (1-rho_n) -> 0, negative binomial (discrete Gamma) when it is of the same order. Distribution only, no estimation.', 'Low',
  'Reading: abstract, introduction, Sections 2.1-2.2 discussion read 25 Oct 2026; proofs not read.'),
 (22, 'A fluctuation limit theorem of branching processes with immigration and statistical applications', 'Ma, arXiv, 2009', 'arXiv 0906.2586', '11',
  'Survey of the CLSE literature for nearly critical GWI (Sriram 1994: branching-diffusion limit for positive variance; Ispany et al. 2005: n^{3/2} normal for vanishing variance); Thm 3.2: CLSE of (m_n, omega_n).', 'Low',
  'Reading: abstract, introduction, Thm 3.2 read 25 Oct 2026; proofs not read.'),
 (23, 'Discrete CLSE papers (abstracts only): Wei-Winnicki 1990; Ispany-Pap-van Zuijlen 2005; Li 2011', 'Ann. Statist. 18 (1990); Adv. Appl. Probab. 37 (2005); Commun. Statist. Theory Methods (2011)', 'DOI 10.1214/aos/1176347876; 10.1017/s000186780000029x; 10.1080/03610926.2010.500109', '11',
  'Li 2011: with non-zero initial values the CLSE of the immigration mean is not consistent (nearly critical Jirina). Wei-Winnicki abstract states only the aim; Herrera-Marin I attributes to them inconsistency of an immigration-mean estimator in the critical case: UNVERIFIED and in conflict with another summary.', 'Low',
  'Reading: abstracts only (OpenAlex), 25 Oct 2026. Full texts not obtained.'),
 (24, 'Regret / redundancy lower bounds: Rissanen 1984 (IEEE IT 30), Rissanen 1986 (Ann. Statist. 14), Merhav-Feder 1998 (IEEE IT 44)', 'Rissanen; Merhav and Feder', 'DOI 10.1109/TIT.1984.1056936; 10.1214/aos/1176350051; 10.1109/18.720534', '9, 5, 6',
  'Classical (d/2) log n-type lower bounds are for fixed parameters; no Hawkes- or point-process-specific minimax regret lower bound found; none for a triangular array with eps -> 0.', 'Unknown',
  'Reading: abstracts only; Merhav-Feder PDF host failed TLS verification (not bypassed). Hypotheses of the Rissanen theorems not checked.'),
]
rows_xml = ''
for k, (order, title, auth, link, probs, extract, risk, status) in enumerate(new_rows):
    r = 24 + k
    cells = (f'<c r="A{r}" s="29" t="n"><v>{order}</v></c>'
             f'<c r="B{r}" s="31" t="inlineStr"><is><t>{escape(title)}</t></is></c>'
             f'<c r="C{r}" s="30" t="inlineStr"><is><t>{escape(auth)}</t></is></c>'
             f'<c r="D{r}" s="30" t="inlineStr"><is><t>{escape(link)}</t></is></c>'
             f'<c r="E{r}" s="29" t="inlineStr"><is><t>{escape(probs)}</t></is></c>'
             f'<c r="F{r}" s="30" t="inlineStr"><is><t>{escape(extract)}</t></is></c>'
             f'<c r="G{r}" s="29" t="inlineStr"><is><t>{escape(risk)}</t></is></c>'
             f'<c r="H{r}" s="32" t="inlineStr"><is><t>{escape(status)}</t></is></c>')
    rows_xml += f'<row r="{r}" ht="75" customHeight="1" s="27">{cells}</row>'

foot = re.search(r'<row r="25".*?</row>', lit, re.S).group(0)
foot_new = foot.replace('r="25"', 'r="31"').replace('B25', 'B31').replace('Overlap risk is my judgement.', 'Overlap risk is my judgement. Updated 25 Oct 2026: rows 7 and 17 corrected, rows 24-29 added from literature_search_report.md; most new entries were read only in part (see Read status).')
assert foot_new != foot
lit = lit.replace(foot, rows_xml + foot_new)
for a, b in (('<dimension ref="A1:H25" />', '<dimension ref="A1:H31" />'), ('<mergeCell ref="B25:F25" />', '<mergeCell ref="B31:F31" />'),
             ('sqref="G2:G23"', 'sqref="G2:G29"'), ('sqref="H2:H23"', 'sqref="H2:H29"')):
    assert a in lit, a
    lit = lit.replace(a, b)

def append(x, r, text):
    m = re.search(rf'<c r="M{r}"( s="\d+")[^>]*?(?:/>|><is><t[^>]*>(.*?)</t></is></c>)', x, re.S); assert m, r
    return x[:m.start()] + f'<c r="M{r}"{m.group(1)} t="inlineStr"><is><t>{m.group(2) or ""}{escape(text)}</t></is></c>' + x[m.end():]
assert re.search(r'<c r="A11"[^>]*><v>10</v>', pro) and re.search(r'<c r="A12"[^>]*><v>11</v>', pro)
pro = append(pro, 11, " LITERATURE 25 Oct (literature_search_report.md): my exact K=1 Hawkes Fisher matrix converges (relative diff ~ eps) to the CIR Fisher matrix of Ben Alaya-Kebaier-Tran Thm 2.1 under (a,b)=(beta mu, beta eps), 2 sigma=beta^2; the conjectured cold-start coefficient 2/(beta(2mu-beta)) equals the CIR critical-case information rate of their Cor 2.3(i). So the K=1 information structure is the Hawkes prelimit of classical CIR theory (fixed regimes); no local-to-critical CIR result found.")
pro = append(pro, 12, " LITERATURE 25 Oct (literature_search_report.md): CIR side: information about the drift constant finite iff a > sigma (Ben Alaya et al. Thm 2.1); MLE defined only for 2a > sigma^2 (Dehtiar et al.); critical case sqrt(log T) rate. Discrete side: all estimation results found are CLSE/OLS; Li 2011: nonzero initial values => CLSE of immigration mean inconsistent; Anne-Lu-Yu-Zhou 2026: OLS intercept inconsistent in local-to-unity affine processes (no Hawkes, no information bound; overlap risk Medium). Herrera-Marin's Wei-Winnicki attribution UNVERIFIED.")

tmp = 'tracker_new.xlsx'
with zipfile.ZipFile(tmp, 'w', zipfile.ZIP_DEFLATED) as zo:
    for i in zin.infolist():
        if i.filename == 'xl/worksheets/sheet4.xml': zo.writestr(i, lit.encode())
        elif i.filename == 'xl/worksheets/sheet1.xml': zo.writestr(i, pro.encode())
        else: zo.writestr(i, zin.read(i.filename))
zin.close()
shutil.move(tmp, src)
z = zipfile.ZipFile(src)
print('zip test:', z.testzip())
y = z.read('xl/worksheets/sheet4.xml').decode()
print(re.findall(r'<row r="(\d+)"', y))
print(re.search(r'<dimension[^>]*>', y).group(0), re.search(r'<mergeCell [^>]*>', y).group(0))
for rn in (7, 17, 24, 29, 31):
    r = re.search(rf'<row r="{rn}".*?</row>', y, re.S).group(0)
    print(rn, re.sub(r'<[^>]+>', ' | ', r)[:170])
