"""List control sequences used in the given .tex files that are neither standard (LaTeX/amsmath/booktabs/graphicx) nor defined in hawkes_master.tex.
Usage: python scan_macros.py part6_epsilon_validation.tex [more.tex ...]   (no compiler needed; complements check_tex.py)"""
import re, sys, io
STD = set("""section subsection subsubsection paragraph label ref emph textbf textit texttt textsf tfrac dfrac frac sum prod int max min ln log det exp le ge leq geq ll gg approx cdot cdots dots ldots
operatorname mathrm mathbb mathcal mathbf boldsymbol bm begin end item hline Big big Bigg bigg left right in notin subset subseteq cup cap times to infty alpha beta gamma delta epsilon varepsilon
lambda mu nu rho sigma tau theta kappa omega Omega Gamma Delta Sigma phi psi chi pi eta xi zeta hat bar tilde widehat widetilde sqrt qquad quad centering small includegraphics linewidth
caption toprule midrule bottomrule textbackslash sup inf lim text forall exists partial nabla preceq succeq prec succ pm mp neq equiv sim propto emptyset setminus mid colon langle rangle
lceil rceil lfloor rfloor not circ ast star dagger top bot vert Vert lVert rVert nonumber hfill vspace noindent cite footnote cline multicolumn arraystretch tabular table figure center
proof theorem lemma proposition corollary definition remark assumption itemize enumerate abstract document emptyset overline underline stackrel substack limits nolimits Pr sin cos tanh
tanh sinh cosh arg argmin argmax newline linebreak par ensuremath mathop mathbin mathrel binom choose lesssim gtrsim ll dd tr""".split())
mac = io.open('hawkes_master.tex', encoding='utf8').read()
defined = set(re.findall(r'\\(?:new|renew|provide)command\*?\s*\{?\\([A-Za-z]+)', mac)) | set(re.findall(r'\\DeclareMathOperator\*?\s*\{\\([A-Za-z]+)', mac))
for f in sys.argv[1:]:
    t = io.open(f, encoding='utf8').read(); t = re.sub(r'(?<!\\)%.*', '', t)
    used = set(re.findall(r'\\([A-Za-z]+)', t))
    print(f, 'unknown control sequences:', sorted(u for u in used if u not in STD and u not in defined))
print('defined in master:', sorted(defined))
