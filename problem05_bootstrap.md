# Problem 5: a true bootstrap (iterates enter a ratio set that shrinks to 1), K = 1

Status: a two-stage argument with explicit (very conservative) constants and polynomial tails, plus numerical checks of every ingredient. It proves that the iterates enter and stay in R(δ) = {|μ/μ* − 1| ≤ δ, |α/α* − 1| ≤ δ} for t ≥ T_1(δ, q), T_1 = O(1/(βε δ²))-type. The sharp constant (d/2) log T is proved only on the good event; as an expected-regret statement it needs exponential tails (open).
Builds on `problem05_circularity.md` (compact box Θ_c, exact F-norm projection, c₀^± deterministic), `problem05_fisher_concentration.md` (Theorem F), `problem07_leverage.md` (Theorem 7).
Code: `code/p5_bootstrap_check.py`.

## Stage 1: a high-probability bound on D_j = ‖θ_{j+1} − θ*‖²_{F_j}, using only the box constants c₀^±

**Lemma C1 (martingale structure).** For the projected iteration and every restart index i_0, for j ≥ i_0,
D_j ≤ D'_0 − X_j − M_j + G_j,
with D'_0 = ‖θ_{i_0} − θ*‖²_{F_{i_0−1}}, X_j = Σ_{i=i_0}^{j} v_iᵀΔF_iv_i ≥ 0, G_j = Σ g_iᵀF_i⁻¹g_i (nondecreasing in j), and M_j = Σ m_i with
m_i = −2(a_iᵀv_i)/λ_θ(τ_i) + 2∫_0^{τ_i}(φᵀv_i)λ*/λ_θ du = −2∫_0^{τ_i}(φᵀv_i/λ_θ)(dN − λ* du).
Moreover E[m_i | past] = 0 and E[m_i² | past] ≤ (4/c⁻) E[v_iᵀΔF_iv_i | past].
Proof. Non-expansiveness (Problem 5 circularity note) gives the Newton inequality 2⟨g_i,v_i⟩ ≤ ‖v_i‖²_{F_i} − ‖v_{i+1}‖²_{F_i} + g_iᵀF_i⁻¹g_i; telescope. Then ⟨g,v⟩ = ⟨s,v⟩ + (aᵀv)²/(λ*λ_θ) and −2⟨g,v⟩ + vᵀΔFv = −2⟨s,v⟩ − A − vᵀΔFv with A = 2(aᵀv)²/(λ*λ_θ) − 2vᵀΔFv; m = 2⟨s,v⟩ + A simplifies to the displayed form because (1 − φᵀv/λ_θ) = λ*/λ_θ. It is a stopped point-process martingale increment with integrand w = −2φᵀv/λ_θ; its predictable variation ∫_0^τ w²λ* du = 4∫(φᵀv)²λ*/λ_θ² du ≤ (4/c⁻)∫(φᵀv)²/λ_θ du, and the last integral has expectation E vᵀΔFv. ∎

**Lemma C2 (sup bound).** Pathwise D_j ≤ D'_0 + sup_{i≤j}(−M_i) + G_j. By Doob's L² maximal inequality and Lemma C1, E sup_{i≤n} M_i² ≤ 4E⟨M⟩_n ≤ (16/c₀⁻) E X_n; by Prop. 1 (restart form), E X_n ≤ E D'_0 + E G_n =: Γ. Since G is nondecreasing, Markov applies to its supremum. For x ≥ 1,
P(sup_{j≤n} D_j ≥ 3x) ≤ P(D'_0 ≥ x) + P(sup(−M) ≥ x) + P(G_n ≥ x) ≤ (1 + 16/c₀⁻)Γ/x.

**Quantities in Γ.** E D'_0 and E G are bounded by the Prop. 1 / Theorem 7 chain with box constants: E G_{[i_0,n]} ≤ A_0 E log det(F_n/F_{i_0}) + B_0 E X, hence Γ ≤ (E D'_0 + A_0 E log det(F_n/F_{i_0}))/(1 − B_0) ≈ O(d log T). (These use the leverage conditions of Theorem F for t ≥ T_0 and the stopping-time localisation, sketched but not written out.)

## Theorem C (bootstrap, K = 1)

Under the standing assumptions (Θ_c, exact projection, Assumption M, event E_F of Theorem F) and for t ≥ T_0 define λ_T := mTλ_Q/(4c₀⁺L) (so F_j ⪰ λ_T I when t_j ≥ T, by F_j ⪰ J_{t_j}/c₀⁺ ⪰ (m t_j/(4c₀⁺L))Q), ρ_* = min(μ*, α*) and Γ(T) := E D'_0 + E G up to time 2T. Then for δ ∈ (0,1) and q ∈ (0,1),

P( ∃ j with t_j ≥ T_1 : θ_{j+1} ∉ R(δ) ) ≤ q + P(E_F^c) whenever δ² ≥ δ_*²(T_1, q) := 24 (1 + 16/c₀⁻) c₀⁺ L Γ(T_1) / ( q m λ_Q ρ_*² T_1 ).

Proof. If θ_{j+1} ∉ R(δ) then |v_{j+1}|² ≥ δ²ρ_*², so D_j ≥ λ_min(F_j)|v_{j+1}|² ≥ λ_Tδ²ρ_*² (for t_j ≥ T). Put 3x_T = λ_Tδ²ρ_*². Lemma C2 on the block t_j ∈ [T, 2T] gives failure probability ≤ (1 + 16/c₀⁻)Γ(T)/x_T = 3(1 + 16/c₀⁻)Γ(T)/(λ_Tδ²ρ_*²). Sum over dyadic T_k = 2^kT_1: Γ(T) grows only logarithmically and 1/λ_T halves each block, so the sum is ≤ 2× the first term. ∎

Scaling: δ_* ∝ (L Γ/(m λ_Q T))^{1/2} with L ∝ μ/ε and Γ ∝ d log T, i.e. the iterates are in R(δ) after T_1 ≳ d log T_1/(εδ²) (time units), or T_1γ ≳ d log T_1/δ², independent of ε apart from the constants m, λ_Q, c₀. The rigorous constant is very large (c₀⁺L/(mλ_Q) is of order 10²–10³, and 1/q); the failure probability decays only like 1/T_1 (Doob/Markov), not exponentially.

## Numerical checks (K = 1; ε = 0.3, 0.1, 0.03; 16 reps; n = 4·10⁴–2·10⁵; box [0.5μ*, 2μ*] × [0.5α*, β(1−ε/2)]; start (1.4, 0.8α*))

- (A) Lemma C1's inequality D_j ≤ D'_0 − X − M + G: maximum violation over all runs and steps ≤ 9×10⁻¹⁵ (rounding).
- (B) E Σ m_i² / E X = 3.74, 3.97, 4.28 (the lemma gives ≤ 4/c⁻; with the box value c⁻ = 0.5 that is 8, and with c⁻ ≈ 1 it is 4, which is what is observed because the iterates are near θ*).
- (C) sup_j |M_j|: median 6–9, 95 % quantile 12–18, against a Doob scale of 80–95 at q = 0.05 (a factor of 5 conservative).
- (D) sup_j D_j: median 5.5–6.9, 95 % ≤ 14, never above 14 in any run; E X = 10.5–14.1, E G = 16–18, E log det(F_n/F_0) = 18–19 (so D = O(d log T) as the lemma predicts).
- (E) Last time outside R(δ), median of T·γ over runs (γ = βε):

| δ | ε = 0.3 | ε = 0.1 | ε = 0.03 |
|---|---|---|---|
| 0.5 | 21.8 | 7.1 | 8.0 |
| 0.3 | 28.9 | 19.7 | 14.0 |
| 0.2 | 84 | 21.5 | 22.9 |
| 0.1 | 317 | 177 | 98 |

Roughly ∝ δ⁻² (δ 0.2 → 0.1 multiplies by 3.8, 8.2, 4.3) and of order 1–3/δ² in units of 1/γ, compared with the rigorous requirement of order 10⁵/(qδ²). No run was outside R(0.3) at the end of the horizon.
- (F) Slope of E R_n against log T_n between consecutive checkpoints: ε = 0.3: 1.1, 0.87, 1.3, 1.01; ε = 0.1: 0.85, 0.99, 0.72, 0.79; ε = 0.03: 1.1, 0.56, 0.58, 0.89. Mean ≈ 0.9, close to d/2 = 1 and not to d = 2.

## Stage 2: constants tend to the ideal ones (conditional)

On the event of Theorem C (iterates in R(δ) for t ≥ T_1), re-run Theorem 7 with c^± = 1 ± δ instead of c₀^±: x̄ ≲ ℓ̄ ln(N/q)/(μ(1−δ)) → 0 and B → 0 as ℓ̄ = O(1/(εTm)) → 0, so A = (1+η)(1+δ)r(x̄) → 1 + η with η → 0; and ψ(r) ≤ (½ + O(δ))(r−1)²/r for r ∈ [1−δ, 1+δ], so the factor 2 lost in Prop. 1(b) becomes 1 + O(δ). On that event
regret on [T_1, T_n] ≤ (½ + O(δ + η + x̄ + B)) · (log det(F_n/F_{i_1}) + O(1)) = (d/2 + o(1)) log(T_n/T_1) + O(1),
since log det(F_n/F_{i_1}) ≈ log det(J_n/J_{i_1}) ≈ d log(T_n/T_1). The regret on [0, T_1] is O(d log T_1) from the loose-constant bound plus the crude burn-in on [0, T_0]. With T_1 = polylog(T_n) the loss is lower order, giving the sharp leading constant (d/2) log T on the good event. Empirically (F) agrees.

## What is not proved

1. **Expected regret with the sharp constant.** The bad event has probability ≈ 1/T_1 only (polynomial tails). To turn the good-event statement into E R_n ≤ (d/2 + o(1)) log T_n one needs failure probability ≲ 1/T_n^κ, i.e. exponential tails (Freedman/Bernstein for M and for the Fisher concentration in Theorem F) or a controlled bound for the regret on the bad event. Not done.
2. The localisation of the leverage bounds (stopping time ν) and of Γ(T) is sketched in `problem05_circularity.md`, not written out; Theorem C's Γ inherits that.
3. Assumption M (π-mass of the bump functions) is numerical; the box [0.5μ*, 2μ*] × [0.5α*, ·] requires a priori knowledge of the parameters within a factor of 2.
4. K ≥ 2 (QP projection; the contraction and two-interval arguments extend but were not tested).
5. The constants in δ_* are 10²–10⁴ times worse than observed; no attempt was made to tighten them.

Re-derive Lemma C1's formula for m_i and the Doob step yourself; the numerics confirm the inequality, the martingale variance bound (B) and the scalings, but not the bookkeeping of constants.

---

## Update (25 Oct 2026): item 1 (exponential tails) is addressed in `problem05_exponential_tails.md`

Lemmas E2–E4 there give exponential-tail versions of Theorem F and of the Newton-martingale step (with probability ≥ 1 − e^{−x}, −M_t ≤ X_t + x/θ for all t, so sup D_j ≤ D'_0 + x/θ + G_n), hence a bootstrap with T_1 = O((ln(1/q) + d log T)/(γδ²)). The control of G_n with exponential tails is sketched, not proved in detail, and the expected-regret statement with the sharp constant still needs an episode/restart argument.
