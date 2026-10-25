import io

def sub(path, pairs):
    s = io.open(path, encoding='utf8').read()
    for a, b in pairs:
        assert a in s, (path, a[:60])
        s = s.replace(a, b)
    io.open(path, 'w', encoding='utf8', newline='\n').write(s)

sub('problem05_E5_expected_regret.md', [
 ("Let x̄ := c⁺ℓ̄m/μ* and", "Let x̄ := ℓ̄m/(c⁻μ*) and"),
 ("On {E_i ≤ m}: tr(F⁻¹ΔF) ≤ ℓ̄c⁺E_i/μ* ≤ x̄,",
  "On {E_i ≤ m}: tr(F⁻¹ΔF) = ∫φᵀF⁻¹φ/λ_θ du ≤ (ℓ̄/c⁻)τ ≤ ℓ̄E_i/(c⁻μ*) ≤ x̄ (λ_θ ≥ c⁻λ*, τ ≤ E_i/μ*),"),
])
sub('latex/part2_regret.tex', [
 (r"$\bar x=c^+\bar\ell m/\mu^*$", r"$\bar x=\bar\ell m/(c^-\mu^*)$"),
 (r"$\tr(F^{-1}\Delta F)\le\bar\ell c^+E_i/\mu^*\le\bar x$",
  r"$\tr(F^{-1}\Delta F)=\int\phi\T F^{-1}\phi/\lambda_\theta\dd u\le\bar\ell\tau/c^-\le\bar\ell E_i/(c^-\mu^*)\le\bar x$"),
])
print('fixed')
