import io, re

def rd(p): return io.open(p, encoding='utf8').read()
def wr(p, s): io.open(p, 'w', encoding='utf8', newline='\n').write(s)

# 1. insert the tangent subsection at the end of the spectrum section (before the projection section)
p1 = rd('part1_model_fisher.tex')
marker = "%==============================================================\n\\section{Projection onto the constraint set}"
assert marker in p1
if 'part1b_tangent' not in p1:
    p1 = p1.replace(marker, "\\input{part1b_tangent}\n\n" + marker, 1)
wr('part1_model_fisher.tex', p1)

# 2. Theorem FK: tangent-based information level (replace the truncation-based statement)
p2b = rd('part2b_extensions.tex')
a = p2b.index("\\begin{theorem}[concentration of the information without Assumption M]")
b = p2b.index("\\end{theorem}", a) + len("\\end{theorem}")
new = r"""\begin{theorem}[concentration of the information without Assumption M]\label{thm:FK}\Ps{}
Let $\lambda_*=\sup_{\lambda_0}\lambda_{\min}(N(\lambda_0))>0$ (Theorem~\ref{thm:tangent}). Pathwise $J_T\succeq\int_0^Taa\T(2/\lambda_0-\lambda^*/\lambda_0^2)\dd t$, whose entries are time integrals of the cubic polynomials $f_{ij}(S)=a_ia_j(2/\lambda_0-\lambda(S)/\lambda_0^2)$ with $\E f_{ij}=N(\lambda_0)_{ij}$. Clipping the state at $s_{\max}$ (componentwise) makes $f_{ij}$ Lipschitz in $d_\ell$ with constant $\Lambda_f(s_{\max})=O(s_{\max}^2)$; on $\{\max_{t\le T}S\le s_{\max}\}$ (probability $\ge1-T^{-\kappa}$ for $s_{\max}=O(\eps^{-1}\ln T)$, Lemma~\ref{lem:expmoment}) the clipping is inactive, and $|\E\tilde f-\E f|$ is exponentially small. The Poisson-equation/Freedman argument ($\hat u$ is $\Lambda_f/\gamma$-Lipschitz by Lemma~\ref{lem:ellcontr}, $|\Delta\mathcal M|\le\Lambda_f/\gamma$, $\langle\mathcal M\rangle\le(\Lambda_f/\gamma)^2\Lambda_T$) and a union over the $(K+1)(K+2)/2$ entries with $\tau=\lambda_*/2$ give, with probability $\ge1-3(K+1)(K+2)e^{-c_fT}-T^{-\kappa}$, $c_f\asymp\tau^2\gamma^2\eps/(\Lambda_f(s_{\max})^2\mu)$:
\[
J_T\succeq\frac{T\lambda_*}{2}\,I_{K+1}.
\]
\end{theorem}"""
p2b = p2b[:a] + new + p2b[b:]
p2b = p2b.replace(r"under a known polytope $\Theta_c\ni\theta^*$ and $\lambda_L>0$", r"under a known polytope $\Theta_c\ni\theta^*$ (the information level is $\lambda_*/2>0$ of Theorem~\ref{thm:FK})")
p2b = p2b.replace(r"Under a known polytope $\Theta_c\ni\theta^*$ and $\lambda_L>0$", r"Under a known polytope $\Theta_c\ni\theta^*$ (information level $\lambda_*/2$ of Theorem~\ref{thm:FK})")
wr('part2b_extensions.tex', p2b)
s = rd('part2b_extensions.tex')
s = s.replace(r"Constants are uninformative: for $\beta=(1,2)$, $\mu=1.5$, $\eps=0.4$ the explicit bound on $\lambda_{\min}(I)$ is $5.6\times10^{-4}$ against a true value $\approx6\times10^{-2}$ and $L\approx50\mu/\eps$, so the information level $\lambda_L/(4L)\approx10^{-6}$.",
              r"With the tangent level the constants are far less pessimistic than with truncation: for $\beta=(1,2)$, $\mu=1.5$, $\eps=0.4$ (0.3) the information level is $\lambda_*/2\approx0.0126$ ($0.0132$) against $\approx10^{-6}$ for the truncation device, a gain of about $10^4$ in $T_1$; they remain far from sharp (the rate $c_f$ contains $\Lambda_f(s_{\max})^{-2}$).")
wr('part2b_extensions.tex', s)

# 3. Theorem R: full proof (replace the outline)
p2 = rd('part2_regret.tex')
a = p2.index("\\begin{proof}[Proof outline (steps written in the notes)]")
b = p2.index("\\end{proof}", a) + len("\\end{proof}")
new = r"""\begin{proof}
Fix $\delta\in(0,\frac14]$, $\eta\in(0,1]$, a non-increasing $\bar\ell(\cdot)\to0$, a threshold $t_*$, and the \emph{predictable} good indicator $g_i=\1\{t_{i-1}\ge t_*\}\1\{\theta_i\in R(\delta)\}\1\{\bar\ell_i\le\bar\ell(t_{i-1})\}$. Inputs: (H3) $\P(\exists i:t_{i-1}\ge T,g_i=0)\le C_1e^{-c_1T}$ for $T\ge t_*$ (Theorem~\ref{thm:Cexp}, information concentration, Lemma~\ref{lem:E0}); (H4) $t_{i-1}\ge i\eps/(2\mu_{hi})$ except with probability $\le Ce^{-ci}$ (upper tail of $\Lambda_T$); (H5) uniformly bounded moments of the gap quantities. Then $\sum_i\P(g_i=0)^{1/2}\le C_b$, hence $\E\sum_i(1-g_i)\Phi_i\le C_B$ for every $\Phi_i\ge0$ with $\E\Phi_i^2\le C_\Phi$, independently of $n$.
\emph{(F1)} $R_i\le X_i:=v_i\T\Delta F_iv_i$ pathwise ($\psi(r)\le(r-1)^2/r$). \emph{(F2)} On $\{\theta_i\in R(\delta)\}$, $\psi(r)\le(\frac12+\delta)(r-1)^2/r$ so $\E[g_iR_i]\le(\frac12+\delta)\E[g_iX_i]$. \emph{(F3)} $\E[g_iH_i\mid\cG_{i-1}]\le g_i\frac{\bar\ell(t_{i-1})}{\mu_{lo}}\E[X_i\mid\cG_{i-1}]$ (Lemma~\ref{lem:E5d} and the compensator identity). \emph{(F4)} $\E[g_iQ_i\mid\cG_{i-1}]=g_iq_i$ and, with $m_i=2\ln(1/\bar\ell(t_{i-1}))$, $q_i\le(1+\delta)r(\bar x_i)\E[\log\det(F_i/F_{i-1})\mid\cG_{i-1}]+3\bar\ell^3\ln(e/\bar\ell)/\mu^*$ (Lemma~\ref{lem:logdet} with $\lambda_\theta\le(1+\delta)\lambda^*$), $\bar x_i=\bar\ell m_i/(c^-\mu^*)\to0$, so $r(\bar x_i)\le1+\delta$ for $t_{i-1}\ge t_*$.
Assembly: by Proposition~\ref{prop:reg}(c), $\E X_n\le P_0+\E G$ with $G=\sum g_i\T F_i^{-1}g_i\le\sum g_i\T F_{i-1}^{-1}g_i$. Bad gaps contribute $\le C_B/\lambda_{\min}(F_0)$. Good gaps: $g\T F^{-1}g\le(1+\eta)Q+(1+1/\eta)H$ and (F3)--(F4) give
\[
\E X_n\le P_0+(1+\eta)(1+\delta)^2\E\log\det\tfrac{F_n}{F_0}+C_\rho+\tfrac{C_B}{\lambda_{\min}(F_0)}+B\,\E X_n,\qquad B=(1+\tfrac1\eta)\tfrac{\bar\ell(t_*)}{\mu_{lo}} .
\]
Choose $t_*$ with $B\le\delta$. Then by (F1)--(F2), $\E R_T\le(\frac12+\delta)\E X_T+C_B$, and $\E\log\det(F_n/F_0)\le\log\det(F_0^{-1}(F_0+T\,I/c^-))+C=d\ln T+O(1)$ (concavity of $\log\det$, $\Delta F_i\preceq\Delta J_i/c^-$, deterministic horizon $T$). Hence $\E R_T\le(\frac12+\delta)(1-\delta)^{-1}(1+\eta)(1+\delta)^2\,d\ln T+C$. No exponential tail of $G$, restart or Ville inequality is used in the assembly; they enter only through (H3). Full details: \texttt{theoremR\_proof\_complete.md}.
\end{proof}"""
p2 = p2[:a] + new + p2[b:]
wr('part2_regret.tex', p2)

# 4. review section
p3 = rd('part3_tracking_review.tex')
p3 = p3.replace(r"$\lambda_{\min}(I)$ truncation bound & \Ps{} & Lemma~\ref{lem:moments} constants; $20$--$400\times$ loose\\",
                r"$\lambda_{\min}(I)$ truncation bound & \Ps{} & Lemma~\ref{lem:moments} constants; $20$--$400\times$ loose\\" "\n" r"Tangent minorant: lower bound for every eigenvalue (Thm.~\ref{thm:tangent}) & \Pt{} & exact moment systems solved in 60-digit arithmetic (not interval arithmetic); within $1.3$--$4\times$ of Monte Carlo\\")
p3 = p3.replace(r"Spectrum: three tiers & \Pt{}/\Nt{} & ridge and top bounds proved; middle tier heuristic\\",
                r"Spectrum: three tiers & \Pt{}/\Nt{} & all tiers bounded below (Thm.~\ref{thm:tangent}); upper bounds only for the ridge; eigenvector convergence numerical\\")
wr('part3_tracking_review.tex', p3)
print('patched')
