import zipfile, re, shutil
from xml.sax.saxutils import escape
src = 'hawkes_research_problem_tracker.xlsx'
zin = zipfile.ZipFile(src); x = zin.read('xl/worksheets/sheet1.xml').decode()
upd = {
 6: ("In progress", " 25 Oct 2026. PARTIAL (problem05_E5_expected_regret.md, theoremR_proof_complete.md, K2_regret_chain.md, box_cost_and_lower_bound.md): Theorem R gives E R_T <= (1/2 + O(delta+eta)) d ln T + C with C = O(mu d ln(1/(eps delta)) / (beta eps^2 delta^2)) (order only, K=1; general K with the tangent information level lambda_*/2); burn-in is ~1/eps^2 events (~1/(beta eps) in time units), not 1/(mu eps); no C log(1/eps) term visible in time units (leading term d/2 ln T + 1/2 ln det I with det I bounded for k0>1). MISSING for Done: every constant of the theorem written as an explicit function of (eps, mu, K, beta-grid); the eps-dependence of C is an order, not a proved sharp rate."),
 8: ("In progress", " 25 Oct 2026. PARTIAL: the notion used throughout is the log-loss (KL) regret of a predictive intensity, E R_n = KL(P_theta^(n) || Q^(n)) (box_cost_and_lower_bound.md B1), with pseudo-regret R_n = sum int lam* psi(lam_theta/lam*) <= X_n pathwise; parameter error enters via the Fisher-norm quantities D_n = ||theta_{n+1}-theta*||^2_{F_n} (Lemma C1) and the tracking theorem (problem14_tracking_theorem.md); the lower bound passes from KL to MSE through mutual information (Prop LB). MISSING for Done: one definition stated at the top of the paper and a single lemma connecting regret and Fisher-norm error (local quadratic equivalence KL ~ 1/2 ||.||^2_I with constants on a ratio set), not yet written."),
 9: ("In progress", " 25 Oct 2026. PARTIAL: (i) path-term two-point bound TV(Q_theta,Q_theta') <= TV(pi_theta,pi_theta') + sqrt(rho T u'I u/2) (thm B); (ii) minimax lower bound for ANY predictor sup_Theta0 E R_n >= (d/2) ln(n/ln n) - O(1) conditional on a uniformity hypothesis (U) (box_cost_and_lower_bound.md, Prop LB) matching the upper bound (d/2) ln T + C in the leading constant; (iii) non-consistency at eps=c/T for every k (problem11_observed_information.md, Thm NC, K=1). Upper-vs-lower comparison: leading constant d/2 matches (fixed eps); eps-dependence of the additive constant: upper O(eps^-2 ln 1/eps) (proof artefact), lower bound has no eps-explicit form; (U) unproved. MISSING for Done: the comparison table in eps, mu, T and an eps-explicit lower bound for eps = eps(n) -> 0."),
 19: ("In progress", " 25 Oct 2026. PARTIAL: the simulation scripts exist and have been run (code/p5_*.py regret slopes: E R_n - ln T flat in eps, slope ~0.9 vs d/2 = 1 for K=1, ~1.3 vs 1.5 for K=2; code/p14_*.py tracking; code/p11_obs_info.py observed information; code/p13*.py from the other session), with tables in the notes. MISSING for Done: plots of error and regret versus T for several eps with fitted slopes compared with theory (no figures have been produced; the numbers are in tables only)."),
 22: ("In progress", " 25 Oct 2026. PARTIAL: literature_crosscheck.md, literature_search_report.md, novelty_reassessment.md (Herrera-Marin I/II, Jaisson-Rosenbaum, CIR drift papers, Hall-Willett, Yang et al., Piggott-Solo abstract only, Anne-Lu-Yu-Zhou 2026) with overlap-risk updates on the Literature sheet. MISSING for Done: 'difference from' statements for each core paper written as final text, and full-text reading of Piggott-Solo (closed access)."),
}
for pid, (status, note) in upd.items():
    row = pid + 1
    m = re.search(rf'<c r="L{row}"( s="\d+")[^>]*?>.*?</c>', x, re.S); assert m, row
    x = x[:m.start()] + f'<c r="L{row}"{m.group(1)} t="inlineStr"><is><t>{status}</t></is></c>' + x[m.end():]
    m = re.search(rf'<c r="M{row}"( s="\d+")[^>]*?(?:/>|><is><t[^>]*>(.*?)</t></is></c>)', x, re.S); assert m, row
    x = x[:m.start()] + f'<c r="M{row}"{m.group(1)} t="inlineStr"><is><t>{m.group(2) or ""}{escape(note)}</t></is></c>' + x[m.end():]
tmp = 'tracker_new.xlsx'
with zipfile.ZipFile(tmp, 'w', zipfile.ZIP_DEFLATED) as zo:
    for i in zin.infolist(): zo.writestr(i, x.encode() if i.filename == 'xl/worksheets/sheet1.xml' else zin.read(i.filename))
zin.close(); shutil.move(tmp, src); print('ok', zipfile.ZipFile(src).testzip())
