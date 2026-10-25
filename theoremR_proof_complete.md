# Theorem R: a complete proof of the expected-regret bound (every step written)

Replaces the "proof outline" of Theorem R (`problem05_E5_expected_regret.md` §4) by a full argument. It is *simpler* than the outline: because the statement is in expectation and the good-gap indicators are predictable, no exponential tail (E5), no Ville inequality and no truncation of the horizon is needed in the final assembly. Those exponential tools are used only for one input, the per-time failure probability (H3). Tags: **[P]** proved here, **[P\*]** proved modulo the cited inputs. K arbitrary, d = K + 1.

## 0. Setting and inputs

Algorithm: Π^{F_i}-projected online Newton on a known compact convex Θ_c ∋ θ\* (F_0 ≻ 0), as in `problem05_online_newton.md`; ratio constants c⁻ ≤ λ_θ/λ\* ≤ c⁺ and λ_θ ≥ μ_lo on Θ_c (Prop. `prop:ratio`). Notation per gap i (history 𝒢_{i−1}): F = F_{i−1}, v_i = θ_i − θ\*, X_i = v_iᵀΔF_iv_i = ∫₀^{τ}(φᵀv)²/λ_θ du, R_i = ∫₀^{τ}λ\*ψ(λ_θ/λ\*)du, g_i the gradient, Q_i = s_iᵀF⁻¹s_i, H_i = h_iᵀF⁻¹h_i, E_i = Λ\*(τ_i), ℓ̄_i the point leverage (sup_u φᵀF⁻¹φ/λ\*).

Fix δ ∈ (0, ¼] and η ∈ (0, 1]. Let ℓ̄(·) be a non-increasing function with ℓ̄(t) → 0 (for K = 1, ℓ̄(t) = C_ℓ/(εt) up to logs; for general K, C_ℓ(ln t)²/(ε²t)), and define the **good indicator**
 g_i := 1{t_{i−1} ≥ t_\*}·1{θ_i ∈ R(δ)}·1{ℓ̄_i ≤ ℓ̄(t_{i−1})}, which is 𝒢_{i−1}-measurable (predictable).
(R(δ) = {|μ/μ\* − 1| ≤ δ, |α_k/α\*_k − 1| ≤ δ ∀k}; on R(δ), |λ_θ/λ\* − 1| ≤ δ.)
The threshold t_\* is chosen below (depends only on δ, η and the constants).

**(H3) per-time failure bound [P\*].** There are C₁, c₁ > 0 with P(∃ i: t_{i−1} ≥ T, g_i = 0) ≤ C₁e^{−c₁T} for all T ≥ t_\*. This is what Theorem Cexp (bootstrap with exponential tails) + the information-concentration theorem (Theorem F/E3 for K = 1, Theorem FK/§4 of `tangent_bound.md` for any K) + Lemma E0/`lem:expmoment` (bounded states) deliver: those statements are uniform in i ≥ i_0 and have failure probability ≤ 5e^{−x} + 3e^{−c γ T_0}, and taking x = c′γT gives the bound with c₁ ∝ γδ².
**(H4) event counting.** P(N_T > 2μ_hiT/ε′) ≤ C e^{−cT} (upper tail of Λ_T, Lemma E2/`lem:E1`), so t_{i−1} ≥ iε/(2μ_hi) except with probability ≤ Ce^{−ci}.
**(H5) moments.** E[X_i²], E|g_i|⁴ ≤ C_Φ for all i (the gap quantities are polynomial in (w_i, S, ∫S) with all moments bounded uniformly in time, Lemma `lem:moments`; θ_i ∈ Θ_c bounded; F_0 ≻ 0 bounds F_i⁻¹ ⪯ F_0⁻¹).

Consequence of (H3)–(H5) **[P]:** S_bad := Σ_i P(g_i = 0)^{1/2} ≤ C_b < ∞, independent of n. (P(g_i = 0) ≤ 1 for the O((μ/ε′)t_\*) gaps with t_{i−1} < t_\*, and for later gaps P(g_i = 0) ≤ P(t_{i−1} < iε/(2μ_hi)) + P(∃ later bad after time iε/(2μ_hi)) ≤ Ce^{−ci}.) Hence for any non-negative Φ_i with E Φ_i² ≤ C_Φ: E Σ_i(1 − g_i)Φ_i ≤ C_Φ^{1/2}·(N_{early} + C_b) =: C_B, independent of n.

## 1. Four pointwise-in-expectation facts

**F1 (R ≤ X).** R_i ≤ X_i pathwise, since ψ(r) ≤ (r − 1)²/r. [P]
**F2 (sharp factor on good gaps).** On {θ_i ∈ R(δ)}, λ_θ/λ\* ∈ [1−δ, 1+δ] and ψ(r) ≤ (½ + C₂δ)(r−1)²/r (ψ(1+x) = x²/2 − x³/3 + …, (r−1)²/r = x²(1 − x + …), ratio = ½ + x/6 + O(x²), C₂ = 1 for δ ≤ ¼). Hence E[R_i | 𝒢_{i−1}] ≤ (½ + C₂δ)E[X_i | 𝒢_{i−1}] on {θ_i ∈ R(δ)} and, g_i being predictable, E[g_iR_i] ≤ (½ + C₂δ)E[g_iX_i]. [P]
**F3 (drift part).** E[g_iH_i | 𝒢_{i−1}] ≤ g_i(ℓ̄(t_{i−1})/μ_lo)E[X_i | 𝒢_{i−1}]. Proof: H_i ≤ (ℓ̄_i/μ_lo)(a_iᵀv_i)²/(λ\*λ_θ)(τ_i) (Lemma `lem:E5d`) and E[(aᵀv)²/(λ\*λ_θ)(τ) | 𝒢_{i−1}] = ∫e^{−Λ\*}(φᵀv)²/λ_θ du = E[X_i | 𝒢_{i−1}] (compensator identity). [P]
**F4 (score part).** E[g_iQ_i | 𝒢_{i−1}] = g_iq_i exactly (Lemma `lem:score`: E[s sᵀ|past] = E ΔJ), and for any m_i > 0 that is 𝒢_{i−1}-measurable, on g_i = 1:
 q_i ≤ (1+δ)·r(x̄_i)·E[log det(F_i/F_{i−1}) | 𝒢_{i−1}] + (ℓ̄(t_{i−1})/μ\*)(m_i + 1)e^{−m_i}, x̄_i = ℓ̄(t_{i−1})m_i/(c⁻μ\*)
(Lemma `lem:logdet`, with c⁺ replaced by 1+δ because θ_i ∈ R(δ), and the corrected x̄). Choose m_i := 2 ln(1/ℓ̄(t_{i−1})) (so e^{−m_i} = ℓ̄²): then x̄_i = 2ℓ̄ ln(1/ℓ̄)/(c⁻μ\*) → 0 as t_{i−1} → ∞, so r(x̄_i) ≤ 1 + δ for t_{i−1} ≥ t_\*, and the remainder is ≤ (3ℓ̄³/μ\*)ln(e/ℓ̄) =: ρ_i. [P]

## 2. Assembly  **[P\*]**

Let P_0 = ‖v_1‖²_{F_0}, X_n = Σ_{i≤n}X_i, G = Σ_i g_iᵀF_i⁻¹g_i, and use F_i ⪰ F_{i−1}. By Proposition `prop:reg`(c): **E X_n ≤ P_0 + E G.** Split G over good and bad gaps:
- *Bad gaps:* g_iᵀF_i⁻¹g_i ≤ |g_i|²/λ_min(F_0), E|g_i|⁴ ≤ C_Φ ⇒ E Σ(1−g_i)gᵀF⁻¹g ≤ C_B/λ_min(F_0) by §0.
- *Good gaps:* gᵀF⁻¹g ≤ (1+η)Q_i + (1+1/η)H_i. By F4 and the tower property: E Σ g_i Q_i ≤ (1+δ)²E log det(F_n/F_0) + E Σ g_iρ_i, and E Σ g_iρ_i ≤ E∫(3ℓ̄(t)³/μ\*)ln(e/ℓ̄(t))λ_tdt ≤ C_ρ < ∞ (ℓ̄(t)³ ln(1/ℓ̄) integrable at infinity, λ_t ≤ μ_hi/ε′ in mean). By F3, E Σ g_iH_i ≤ B′·E X_n with B′ := ℓ̄(t_\*)/μ_lo.
Therefore, with B := (1+1/η)ℓ̄(t_\*)/μ_lo,
 E X_n ≤ P_0 + (1+η)(1+δ)²E log det(F_n/F_0) + (1+η)C_ρ + C_B/λ_min(F_0) + B·E X_n.
Choose t_\* so large that B ≤ δ (possible since ℓ̄(t) → 0, μ_lo is a box constant): then
 **E X_n ≤ (1−δ)⁻¹[P_0 + (1+η)(1+δ)²E log det(F_n/F_0) + C_3].**
By F1, F2, §0 and (1 − g_i)R_i ≤ (1 − g_i)X_i:
 E R_n ≤ (½ + C₂δ)E X_n + E Σ(1−g_i)X_i ≤ (½ + C₂δ)E X_n + C_B.
Finally log det is concave: E log det(F_n/F_0) ≤ log det(F_0⁻¹(F_0 + E[F_n − F_0])) and F_n − F_0 ⪯ J_n/c⁻ with E J_n = E∫₀^{T_n}aaᵀ/λ\* dt ⪯ T·I(θ\*) + C_cold (stationary or empty start; the transient excess is O(ln T) in the μ-entry only, Prop. `prop:cold`, and negligible against T·I_μμ for fixed ε), so
 E log det(F_n/F_0) ≤ d ln T + log det(I(θ\*)/c⁻) + C_4, where the horizon is a deterministic time T and F_n is the matrix after the last gap completed before T (so E J = E∫₀^T aaᵀ/λ\* dt, no random stopping time enters; n = N_T is then random, and the statement is E R_T ≤ …, which is the same statement as for n events up to the O(1) boundary gap).
Collecting,
 **E R_n ≤ (½ + C₂δ)(1−δ)⁻¹(1+η)(1+δ)²·d·ln T_n + C(δ, η, ε, μ, β, Θ_c, F_0),**
with C independent of n. The prefactor is ½·d·(1 + O(δ+η)), which gives limsup_n E R_n/ln T_n ≤ (½ + O(δ+η))d for all δ, η, hence ≤ d/2. ∎

## 3. What changed relative to the outline, and what is still an input

- The assembly needs no exponential bound for G, no restart, no Ville inequality; the E5 machinery is used only inside (H3) (to show the iterates are in R(δ) and the leverage is small at a *fixed late time* with exponentially small failure probability).
- The constant C collects: P_0, the early gaps (their number is O((μ_hi/ε′)t_\*)), C_B (bad gaps), C_ρ, C_3, C_4. Its dependence on ε is through t_\* ∝ the time for ℓ̄(t) ≤ δμ_lo/(1+1/η) (K = 1: ℓ̄ ≈ C/(εt) ⇒ t_\* ∝ 1/(εδμ_lo)), the event rate μ/ε, and the information level; order O(μ d ln(1/(εδ))/(βε²δ²)) as stated before, not improved.
- **Inputs not proved here:** (H3) rests on Theorem Cexp and the information-concentration theorem; (H4), (H5) are moment/tail facts of Section `dynamics`. The honest statement is therefore: *Theorem R is proved assuming (H3)*, and (H3) is proved modulo the already-listed conservative constants of Theorem F-exp (K = 1) or the polynomial-Lipschitz constants of Theorem FK (any K).
- A genuine subtlety: g_i uses ℓ̄_i, which depends on F_{i−1} and on the state S_{t_{i−1}}: both 𝒢_{i−1}-measurable, so g_i is predictable and F3, F4 are legitimate conditional-expectation steps.

**Status: Theorem R — all steps written (this file); remaining input (H3).**
