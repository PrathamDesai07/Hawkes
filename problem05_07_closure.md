# Closing Problems 5 and 7 against their "Done when" criteria

Tracker criteria: **Problem 5** "A theorem with proof, in expectation or with high probability, under stated assumptions." **Problem 7** "A lemma bounding the joint deviation of the two gradient components."
Verdict after checking the notes: both criteria are met by results that already exist in pieces; what was missing was a single clean statement for each. Both are given here. Everything is a consolidation of lemmas proved in `problem05_online_newton.md`, `problem07_leverage.md`, `problem05_exponential_tails.md`, `problem05_E5_expected_regret.md`, `problem05_bootstrap.md` (no new estimates, apart from the elementary pathwise bound R_n ≤ X_n and the Gamma-sum tail used in the burn-in part).

Notation: gap i, F = F_{i−1}, v_i = θ_i − θ\*, E_i = Λ\*(τ_i) (conditionally Exp(1), independent of the past), ℓ̄_i the point leverage (sup_u φᵀF⁻¹φ/λ\*), s_i = b_i − a_i/λ\*(τ_i), h_i = g_i − s_i.

## Problem 7: Lemma JD (joint deviation of the two gradient components)  **[P]**

Let K be arbitrary, θ\* the truth, λ\* ≥ μ\*, and ℓ̄_i ≤ ℓ̄ on the gap. The two parts of the gradient at θ\* are b_i = ∫₀^τ φ du (compensator part) and a_i/λ\*(τ_i) (jump part); s_i is their difference.

1. **Joint zero mean and covariance.** E[b_i − a_i/λ\*(τ_i) | past] = 0, so the two parts have the same conditional mean, and E[s_is_iᵀ | past] = E[ΔJ_i | past], ΔJ_i = ∫₀^τ φφᵀ/λ\* du. Because b_i is not parallel to a_i (S decays strictly, so φ(u) = (1, S e^{−βu}) spans R² on any gap with S > 0), ΔJ_i ≻ 0 on such gaps: the joint deviation is non-degenerate and carries no hidden constraint.
2. **Pathwise bound on each part and on the joint deviation.** |a_i/λ\*(τ_i)|²_{F⁻¹} ≤ ℓ̄/μ\*, |b_i|²_{F⁻¹} ≤ ℓ̄E_i²/μ\*, and therefore
 **|s_i|²_{F⁻¹} ≤ (2ℓ̄/μ\*)(1 + E_i²).**
3. **Tail.** For every m > 0, P(|s_i|²_{F⁻¹} > (2ℓ̄/μ\*)(1 + m²) | past) ≤ e^{−m}.
4. **Cumulative (time-uniform) deviation.** Truncating at E_i ≤ m, for 0 < ε′ ≤ e − 2 and B̄ = (2ℓ̄/μ\*)(1 + m²), with probability ≥ 1 − ne^{−m} − e^{−x}, for all j: Σ_{i≤j} |s_i|²_{F⁻¹} ≤ (1 + ε′)Σ_{i≤j} q_i + (e − 2)B̄x/ε′, with q_i = tr(F⁻¹E[ΔJ_i | past]) (Lemma E5b).

*Proof.* (1) The stopped martingale M_t = ∫₀^t(φ/λ\*)(dN − λ\*du) satisfies M_τ = a/λ\*(τ) − ∫₀^τφ = −s_i and ⟨M⟩_τ = ΔJ_i; optional stopping gives mean zero and E[M_τM_τᵀ] = E⟨M⟩_τ. (2) λ\* ≥ μ\* gives the first bound; for b, |b|_{F⁻¹} ≤ ∫√(φᵀF⁻¹φ)du ≤ √ℓ̄∫√λ\* ≤ √(ℓ̄τE_i) (Cauchy–Schwarz) and τ ≤ E_i/μ\*; combine with |x − y|² ≤ 2|x|² + 2|y|². (3) P(E_i > m) = e^{−m} since E_i is Exp(1) given the past (time-change theorem). (4) Lemma E5b. ∎

*Numerical check.* The pathwise bound in (2) had 0 violations in 1.44·10⁶ simulated gaps (maximum |s|²_{F⁻¹}/bound = 0.26, ε = 0.3, 0.1, 0.03); Σ|s|²/Σq had median 0.95–1.02 (5–95 %: 0.81–1.23); the drift part is empirically negligible (ΣH ≈ 0.04–0.39 vs ΣQ ≈ 15–19).

**Problem 7: Done** (K = 1 numerics; the lemma and proof are dimension-free).

## Problem 5: Theorem HP (high-probability regret bound for the Fisher-preconditioned online Newton)  **[P\*]**

**Setting.** Data from the stationary (or any) Hawkes law with parameter θ\*; K = 1 (d = 2), or any K with d = K + 1 under the replacements of `K2_regret_chain.md`. Θ_c a known compact convex set (box, resp. polytope) containing θ\*, with constants c⁻ ≤ λ_θ/λ\* ≤ c⁺, λ_θ ≥ μ_lo; iteration θ_{i+1} = Π^{F_i}_{Θ_c}(θ_i − F_i⁻¹g_i), F_i = F_{i−1} + ΔF_i, F_0 ≻ 0; pseudo-regret R_n = Σ_i∫₀^{w_i}λ\*ψ(λ_{θ_i}/λ\*)du, ψ(r) = r − 1 − log r.

**Assumptions.** Fix η > 0, i_0 ≥ 1 and ℓ̄. Let 𝓔 be the event that ℓ̄_i ≤ ℓ̄ for all i ∈ [i_0, n] with B_X := (1 + 1/η)(e − 1)ℓ̄/μ_lo ≤ 1/4. (By Theorem F/E3 and Lemma E0 for K = 1, or Theorem FK for general K, P(𝓔^c) ≤ 3e^{−cγT_0} + e^{−x} once t_{i_0−1} ≥ T_0 = O((x + ln n)/γ).)

**Theorem HP.** Let x > 0, m = ln n + x, ψ_max = max(ψ(c⁻), ψ(c⁺)), κ₁ = max(1, 1/c⁻ − 1), θ_{1/2} = min(1/(2κ₁), c⁻/(8(e − 2))). Then with probability ≥ 1 − 6e^{−x} − P(𝓔^c),
 **R_n ≤ ψ_max(2i_0 + 3x) + 4[ D′_0 + x/θ_{1/2} + 𝒢_n(x) ],**
where D′_0 = ‖θ_{i_0} − θ\*‖²_{F_{i_0−1}} and 𝒢_n(x) is the X-free part of Proposition E5:
 𝒢_n(x) = (1+η)(1+ε′)c⁺r(x̄)(1−ε″/2)⁻¹(log det(F_n/F_{i_0−1}) + b′x/ε″) + (1+η)[(1+ε′)(ℓ̄/μ\*)(m+1)e^{−x} + (e−2)B̄x/ε′] + (1+1/η)(ℓ̄/μ_lo)b_Ax.
In particular, if ℓ̄(ln n + x)² → 0 and c⁺ → 1: R_n ≤ 4(1 + o(1))·d·ln T_n + O(x + i_0 + D′_0) with probability ≥ 1 − 6e^{−x} − P(𝓔^c).

*Proof.* (i) *Pathwise R ≤ X.* ψ(r) ≤ (r−1)²/r gives, gap by gap, ∫λ\*ψ ≤ ∫(φᵀv)²/λ_θ du = v_iᵀΔF_iv_i, so R_{[i_0,n]} ≤ X_n := Σ_{i≥i_0}v_iᵀΔF_iv_i pathwise. (ii) *Newton telescoping* (Lemma C1, valid for the projected iteration): 0 ≤ D_n ≤ D′_0 − X_n − M_n + G_n. (iii) *Martingale term* (Lemma E4 with ρ = ½): with probability ≥ 1 − e^{−x}, −M_n ≤ ½X_n + x/θ_{1/2}. (iv) *Gradient-norm term* (Proposition E5): with probability ≥ 1 − 4e^{−x}, G_n ≤ 𝒢_n(x) + B_XX_n. (v) Combining, X_n(1 − ½ − B_X) ≤ D′_0 + x/θ_{1/2} + 𝒢_n(x), and B_X ≤ ¼ gives X_n ≤ 4(D′_0 + x/θ_{1/2} + 𝒢_n(x)). (vi) *Burn-in gaps i < i_0:* ∫λ\*ψ(λ_θ/λ\*) ≤ ψ_max·E_i (ψ convex, ratio in [c⁻, c⁺]), and ΣE_i ≤ 2i_0 + 3x w.p. ≥ 1 − e^{−x} (Gamma/Bernstein tail: P(Σ(E_i − 1) ≥ 2√(i_0x) + 2x) ≤ e^{−x}). The failure probabilities add to 6e^{−x}. ∎

Remarks. (1) The factor 4 is not optimal: the sharp constant ½ (Theorem R / Reg-sharp) needs the ψ-expansion on R(δ) and the bootstrap; Theorem HP needs neither. (2) The assumptions are on the information level (leverage), certified by Theorem F-exp with Assumption M (K = 1; `assumptionM_exact.md`) or Theorem FK (any K). (3) The expectation version (E R_n ≤ P_0 + E G, Prop. 1) is the K-free statement with no leverage assumption on 𝓔 beyond those inside G.

**Problem 5: Done** — a theorem with proof, with high probability (Theorem HP) and in expectation with the sharp constant (Theorem R, `problem05_E5_expected_regret.md`), under stated assumptions (known compact convex Θ_c; leverage/information event; Assumption M certified for K = 1, replaced by λ_L for K ≥ 2). Residual items (constants 10⁵–10⁶ conservative; K ≥ 2 numerics at n ≤ 4·10⁴ only) are quality issues, not gaps in the criterion.
