import io

def sub(path, a, b):
    s = io.open(path, encoding='utf8').read()
    assert a in s, (path, a[:70])
    io.open(path, 'w', encoding='utf8', newline='\n').write(s.replace(a, b, 1))

# --- K2 note: replace the loose Z_omega version of Lemma E1 by the exact Z_c identity
s = io.open('K2_regret_chain.md', encoding='utf8').read()
a = s.index("### 2.1 The identities of Lemma E1 for general K")
b = s.index("### 2.2 Theorem FK")
new = """### 2.1 The identities of Lemma E1 for general K  **[P]**
With c_k = α_k/β_k, Z_c = Σc_kS_k (so Σc_k = n): dZ_c = −(λ − μ)dt + n dN. Writing N = Λ + M′ (M′ = N − Λ) gives the **exact** identity
 εΛ_T = μT + n M′_T − (Z_c(T) − Z_c(0)),
the integrated form of the key identity A Z_c = μ − ελ (Problem 3), valid for every K and identical to the K = 1 formula (a) of Lemma E1. Hence the lower tail of Λ_T (Lemma E2) holds verbatim:
 P(Λ_T ≤ (1 − η)μT/ε, |Z_c(T) − Z_c(0)| ≤ y) ≤ exp(−(ημT/n − y)²ε/(2μT)),
and the upper tail is Bernstein, with Z_c having the exponential moment of Lemma `lem:expmoment`. (An earlier draft of this note used Z_ω = ωᵀS and the inequalities (b − γ)∫Z ≤ Σα_k∫S_k ≤ (B − γ)∫Z, which lose a constant factor B·ν̃ ≥ 1 in the upper tail for K ≥ 2; the exact Z_c identity makes that unnecessary.)

"""
s = s[:a] + new + s[b:]
io.open('K2_regret_chain.md', 'w', encoding='utf8', newline='\n').write(s)

# --- LaTeX
sub('latex/part2b_extensions.tex',
    r"and Lemma~\ref{lem:E1} by (with $Z=\omega\T S$, $\omega_k=\alpha_k/(\beta_k-\gamma)$): $\dd Z=-(\lambda-\mu+\gamma Z)\dd t+\dd N$, so $\gamma\int_0^TZ=\mu T+M'_T-(Z_T-Z_0)$ and $\Lambda_T\in\mu T+[b-\gamma,B-\gamma]\int Z$.",
    r"while Lemma~\ref{lem:E1}(a) holds verbatim for every $K$ with $Z_c=\sum_k(\alpha_k/\beta_k)S_k$: $\eps\Lambda_T=\mu T+nM'_T-(Z_c(T)-Z_c(0))$ (integrated key identity of Proposition~\ref{prop:keyid}).")
print('ok')
