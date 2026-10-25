# Problem 5 (K = 1): E5 proved, and expected regret with the sharp constant

Closes the two items left open in `problem05_exponential_tails.md` (E5; "expected regret needs an episode/restart argument").
Status of each part is marked **[proved]**, **[proved modulo X]**, or **[checked numerically]**. Code: `code/p5_e5_check.py`, `code/p5_e0_expmoment.py`.
Notation as in `problem05_bootstrap.md` and `problem07_leverage.md`: gap i, F = F_{i−1} (predictable), φ(u) = (1, S(u)), λ* = φᵀθ* ≥ μ*, v = θ_i − θ*, s = b − a/λ*(τ) (score), h = g − s = a/λ*(τ) − a/λ_θ(τ), E_i := Λ*(τ_i) (compensator over the gap).

## 0. A correction to Theorem 7 first

Theorem 7 and Lemma 7c assumed "(L2) x_i ≤ x̄ almost surely", x_i = λ_max(F⁻¹ΔF). This is **not** implied by anything: tr(F⁻¹ΔF) ≤ c⁺ℓ̄E_i/μ* and E_i = Λ*(τ_i) is Exp(1) given the past, hence unbounded. The fix below (truncate at E_i ≤ m, pay (m+1)e^{−m}) replaces (L2) and costs nothing in the constant. Theorem 7 as stated is therefore only valid with the truncated form of Lemma 7c below.

## 1. Four lemmas **[proved]**

**Lemma E0 (exponential moments of the state).** Let θ_* > 0 solve (e^{θ}−1)/θ = β/α (θ_* ≈ 2ε/(1−ε)), and κ_θ := βθ − α(e^θ−1) > 0 for 0 < θ < θ_*. Then with V = e^{θS}: 𝓛V = V[μ(e^θ−1) − κ_θ S] (exact; 𝓛 is the generator, jump rate μ+αS). Hence 𝓛V ≤ −(κ_θ s_0/2)V + b_θ 1{S ≤ s_0}, s_0 = 2μ(e^θ−1)/κ_θ, and Foster–Lyapunov gives sup_t E e^{θS_t} ≤ C(θ, S_0) < ∞. By Markov and a union over the n gap-starts, max_{i≤n} S_{t_i} ≤ (ln C + 2 ln n + x)/θ with probability ≥ 1 − e^{−x}/n.
Numerics: E e^{θS_t} plateaus (ε = 0.3: ≈ 5.4 at θ = θ_*/2; ε = 0.1: ≈ 4.3). At θ = 0.8θ_* the Monte Carlo is noisy (10–60) but bounded.
Consequence: the scale of the state is 1/ε, the same as ES = μ/γ.

**Lemma E5a (pathwise bound on one score term).** If ℓ̄_i ≥ sup_{y∈[0,S_{t_{i−1}}]} (1,y)F_{i−1}⁻¹(1,y)ᵀ/(μ*+α*y), then
 Q_i := s_iᵀF_{i−1}⁻¹s_i ≤ (2ℓ̄_i/μ*)(1 + E_i²).
Proof. Q ≤ 2bᵀF⁻¹b + 2aᵀF⁻¹a/λ*². Since λ* ≥ μ*, aᵀF⁻¹a/λ*² ≤ ℓ̄/μ*. For b = ∫φ: |b|_{F⁻¹} ≤ ∫√(φᵀF⁻¹φ) ≤ √ℓ̄ ∫√λ* ≤ √ℓ̄ √(τ Λ*(τ)) (Cauchy–Schwarz), and τ ≤ E/μ* because λ* ≥ μ*. So bᵀF⁻¹b ≤ ℓ̄E²/μ*. ∎
(P1) checked: 0 violations in 1.44·10⁶ gaps, max Q/bound = 0.26.
Also E_i is Exp(1) and independent of the past (time-change theorem); this is what gives sub-Weibull control with no moment assumption on the state.

**Lemma E5b (score part, exponential tail).** Fix i_0, ℓ̄ with ℓ̄_i ≤ ℓ̄ for i_0 ≤ i ≤ n, truncation level m > 0, B̄ := (2ℓ̄/μ*)(1+m²), and 0 < ε' ≤ e−2. Put q_i := tr(F_{i−1}⁻¹E[ΔJ_i | past]) = E[Q_i | past] (Lemma 7a). Then with probability ≥ 1 − n e^{−m} − e^{−x}, simultaneously for all j ∈ [i_0, n],
 Σ_{i=i_0}^{j} Q_i ≤ (1+ε') Σ_{i=i_0}^{j} q_i + (e−2) B̄ x/ε'.
Proof. On {E_i ≤ m ∀i ≤ n} (probability ≥ 1 − ne^{−m}) we have Q_i = Q'_i := min(Q_i, B̄) by E5a. Q'_i ∈ [0, B̄], E[Q'_i | past] ≤ q_i, E[Q'_i² | past] ≤ B̄ q_i. For 0 < θ ≤ 1/B̄ use e^y ≤ 1+y+(e−2)y² (y ≤ 1): exp(θΣ(Q'_i − E[Q'_i|past]) − (e−2)θ²B̄Σq_i) is a supermartingale; Ville gives, with probability ≥ 1 − e^{−x} for all j, ΣQ' ≤ ΣE[Q'|past] + (e−2)θB̄Σq + x/θ. Take θ = ε'/((e−2)B̄) (≤ 1/B̄ iff ε' ≤ e−2). ∎
Interpretation: the fluctuation term is B̄ x/ε', which is small as soon as ℓ̄ (ln n + x)² x ≪ ε'. With ℓ̄_i ≈ C(ε)/i this is a *slowly growing* index requirement, not a burn-in on n.

**Lemma E5c (conversion q → log det, truncated 7c).** Suppose λ_θ ≤ c⁺λ* on the gap and ℓ̄_i ≤ ℓ̄. Let x̄ := c⁺ℓ̄m/μ* and ℓ'_i := log det(F_i/F_{i−1})·1{E_i ≤ m} ∈ [0, b], b := d log(1+x̄). Then
 q_i ≤ c⁺ r(x̄) E[ℓ'_i | past] + (ℓ̄/μ*)(m+1)e^{−m}.
Proof. On {E_i ≤ m}: tr(F⁻¹ΔF) ≤ ℓ̄c⁺E_i/μ* ≤ x̄, hence λ_max(F^{−1/2}ΔFF^{−1/2}) ≤ x̄ and Lemma 7c applies to ΔJ·1{E_i≤m} ⪯ c⁺ΔF·1{E≤m}. On {E_i > m}: tr(F⁻¹ΔJ_i) ≤ ℓ̄E_i/μ* so E[·; E>m] ≤ (ℓ̄/μ*)(m+1)e^{−m}. ∎
Conversion of conditional expectations to realised sums: for ℓ' ∈ [0,b] and 0 < ε'' ≤ 1, E[e^{−θℓ'}|past] ≤ exp(−θ(1−θb/2)E[ℓ'|past]); with θ = ε''/b, Ville gives with probability ≥ 1 − e^{−x}, for all j,
 Σ E[ℓ'_i | past] ≤ (1−ε''/2)⁻¹ ( Σ ℓ'_i + b x/ε'' ) ≤ (1−ε''/2)⁻¹ ( log det(F_j/F_{i_0−1}) + b x/ε'' ).

**Lemma E5d (drift part).** With λ_θ ≥ μ_lo, ℓ̄_i ≤ ℓ̄, and w_A := (φᵀv)²/(λ_θλ*) ≤ κ_1²c⁺ =: b_A (κ_1 = max(1, 1/c⁻−1), because φᵀv/λ_θ = 1 − λ*/λ_θ): h_iᵀF⁻¹h_i ≤ (ℓ̄/μ_lo)(a_iᵀv_i)²/(λ*λ_θ)(τ_i), and (a_iᵀv_i)²/(λ*λ_θ)(τ_i) = v_iᵀΔF_iv_i + ν_i, ν_i the increment of N^A_t = ∫w_A(dN − λ*du), with ⟨N^A⟩ ≤ b_A X (because w_A²λ* ≤ b_A(φᵀv)²/λ_θ). Same Ville argument with θ = 1/b_A: with probability ≥ 1 − e^{−x}, N^A_j ≤ (e−2)X_j + b_A x for all j. Hence
 Σ_{i≤j} h_iᵀF⁻¹h_i ≤ (ℓ̄/μ_lo) ((e−1)X_j + b_A x).
(P2) checked: 0 violations; the identity (P3) checked: N^A_n/X_n has median ≈ 0, max 0.05 (the drift term is negligible: E ΣH = 0.04–0.39 vs E ΣQ = 15–19).

## 2. Proposition E5 (exponential tail for G_n) **[proved]**

Assume for i ∈ [i_0, n]: c⁻λ* ≤ λ_θ ≤ c⁺λ* (guaranteed by the box), λ_θ ≥ μ_lo, ℓ̄_i ≤ ℓ̄ (predictable condition on F_{i−1}, S_{t_{i−1}}; holds with high probability by Theorem F and Lemma E0), and B_X := (1+1/η)(e−1)ℓ̄/μ_lo ≤ 1/4. Take m = ln n + x. Then with probability ≥ 1 − 4e^{−x}, for all j ∈ [i_0, n]:
 G_j ≤ (1+η)(1+ε')·c⁺r(x̄)/(1−ε''/2)·( log det(F_j/F_{i_0−1}) + b x/ε'' ) + (1+η)[(1+ε')(ℓ̄/μ*)(m+1)e^{−m}n + (e−2)B̄x/ε'] + B_X X_j + (1+1/η)(ℓ̄/μ_lo)b_A x.
Proof. G_j ≤ Σ g_iᵀF_{i−1}⁻¹g_i (F_i ⪰ F_{i−1}) and g = s + h gives gᵀF⁻¹g ≤ (1+η)Q + (1+1/η)hᵀF⁻¹h. Apply E5b, E5c (summing the (ℓ̄/μ*)(m+1)e^{−m} remainder over ≤ n gaps, which is (ℓ̄/μ*)(m+1)e^{−x}) and E5d; the four failure events are the truncation (≤ e^{−x} since ne^{−m} = e^{−x}), and one Ville event each in E5b, E5c, E5d. ∎

**Sharp form.** As ℓ̄·(ln n+x)² → 0 and c⁺ → 1: x̄ → 0, r(x̄) → 1, b → 0, B̄ → 0, and every x-term has a vanishing coefficient. Then G_j ≤ (1+o(1)) log det(F_j/F_{i_0−1}) + o(1)·x + B_X X_j. The ratio of the realised ΣQ to Σq in simulation: median 0.95–1.02 (5–95 %: 0.81–1.23), consistent with this.

## 3. Closing the bootstrap with E5 **[proved modulo Theorem F-exp and Assumption M]**

In Lemma E4 replace the constant by ρ ∈ (0,1): for 0 < θ ≤ min(1/(2κ_1), ρc⁻/(4(e−2))), with probability ≥ 1 − e^{−x}, −M_t ≤ ρX_t + x/θ_ρ for all t. With Lemma C1 and Proposition E5 (take ρ = 1/2 and B_X ≤ 1/4):
 sup_{i_0≤j≤n} D_j + ¼X_n ≤ D'_0 + x/θ_ρ + 𝒢_n(x), where 𝒢_n(x) := the X-free part of Proposition E5.
This is the statement left as "sketch" in `problem05_exponential_tails.md`; Theorem C-exp is now fully justified (given the stated inputs), with g_n = 𝒢_n(x) = (1+o(1))d·log T_n + O(x).

## 4. Expected regret with the sharp constant — episode argument **[proved modulo the inputs listed in §5]**

**Theorem R (K = 1).** Under the standing assumptions, for every δ, η > 0 there is C = C(δ, η, ε, μ, β, box, ρ_0), independent of n, with
 E R_n ≤ (½ + O(δ + η)) · d · ln T_n + C,   d = 2.
Hence limsup_n E R_n/ln T_n ≤ d/2 = 1.

Proof (the point: bad events are charged to *their own time t*, not to the horizon n, so there is no union bound over n and no ln n-length crude burn-in).

*Step 1: per-time failure event.* Let 𝓑_t = {θ_t ∉ R(δ)} ∪ {λ_min(F_t) < λ_t} ∪ {ℓ̄(t) > ℓ̄_t} ∪ {S_t > s_t}, where λ_t is the Theorem F threshold, ℓ̄_t = C_ℓ/(εt)·(1+ s_t-terms), s_t = (C + 2 ln t)/θ_*·… from Lemma E0. By Theorem F-exp (E3: failure ≤ 3e^{−cγt}), Lemma E0 (polynomially small, then use x = c γ t) and Theorem C-exp (with q = e^{−c'γt}, which is allowed once t ≥ t_* since the burn-in T_1(q) = O((ln(1/q) + d ln t)/(γδ²)) ≤ t for t ≥ t_* independent of n), P(𝓑_t) ≤ C_1e^{−c_1γδ²t} for t ≥ t_*.

*Step 2: the pathwise decomposition.* Per-gap regret r_i satisfies E[r_i | past] ≤ (½ + C_2δ) E[v_iᵀΔF_iv_i | past] if θ_i ∈ R(δ) (ψ-expansion, `problem05_bootstrap.md`, Stage 2) and ≤ C_3 E[v_iᵀΔF_iv_i | past] always (Prop. 1(b), box constants). Hence
 E R_n ≤ E R_{[0,t_*]} + (½ + C_2δ) E X_{(t_*,n]} + C_3 Σ_i E[v_iᵀΔF_iv_i 1_{𝓑_{t_i}}].
The last sum: v_iᵀΔF_iv_i is bounded by box-diameter² times (μ/ε)-order constants, and Σ_i P(𝓑_{t_i}) ≤ (μ_max/ε)∫_{t_*}^∞ C_1e^{−c_1γδ²t}dt < ∞, independent of n. The first term is a constant.

*Step 3: restart.* Let ν_1 < ν_2 < … be the successive stopping times at which a predictable condition (λ_min(F) ≥ λ_t, ℓ̄_i ≤ ℓ̄_t) fails after having held, and restart the telescoping of Lemma C1 at each ν_k and at each recovery. Lemma C1 holds for any stopping-time restart index i_0 with D'_0 = ‖θ_{i_0} − θ*‖²_{F_{i_0−1}}. Summing the episodes, the log-determinants add to at most log det(F_n/F_0) = d ln T_n + O(1) (one *global* telescoping, because the episodes are disjoint intervals), while each restart costs E D'_0(ν_k) ≤ diam(box)²·E λ_max(F_{ν_k}) ≤ C_4 (μ_max/ε)E[ν_k ; ν_k<n]. The failure times satisfy, as in Step 1, Σ_k E[t_{ν_k}1{ν_k<n}] ≤ ∫ t P(𝓑_t) dt < ∞ independent of n.

*Step 4: expected value of the good part.* Inside an episode where the conditions hold for all i, Proposition E5 gives G ≤ A·Ld + o(1)x + B_XX with 4e^{−x} tails; integrating over x, E[G 1_good] ≤ A E Ld + O(1) with A = (1+η)(1+ε')c⁺r(x̄)/(1−ε''/2) → 1 + O(η) once t ≥ T_1(η) (where c⁺ ≤ 1+δ, x̄ → 0). With E X ≤ (E D'_0 + E G)/(1 − B_X) (Prop. 1 with B_X ≤ 1/4 → 0), and E Ld ≤ d ln T_n + O(1) by the E2 concentration of J_n/T. The bad part of G (events of probability e^{−cγt}) is bounded by Cauchy–Schwarz with E G_i² ≤ poly (g is bounded by box constants over λ_min(F_0)), giving the same convergent integral. Collecting: E R_n ≤ (½ + C_2δ)(1+O(η)) d ln T_n + C. ∎

## 5. What this does and does not establish

Proved here: Lemmas E0, E5a–d; Proposition E5; the closed form of the exponential bootstrap; Theorem R *given* the following inputs, none of which is new:
1. **Theorem F-exp (E3)** (lower tail of the accumulated information with rate cγt): proved in `problem05_exponential_tails.md` modulo integrability, with a rate c ≈ 10⁵× smaller than observed. This only affects C, not the leading constant.
2. **Assumption M** (π-mass of the bump functions): numerical only.
3. **Compact box with known ratio** [0.5μ*, 2μ*] × [0.5α*, β(1−ε/2)]: an a-priori knowledge assumption (circularity note).
4. The exact F-norm projection and the non-expansiveness inequality (proved, `problem04_projection.md`, `problem05_circularity.md`).
5. In Step 3 I treated the restart as "the telescoping restarts at a stopping time". I have checked that Lemma C1 and Prop. 1 are stated for arbitrary predictable restart indices, but I did **not** write out the optional-stopping details for the martingale M across *random* restart indices; they are standard (M_j − M_{ν} has zero conditional mean) but unwritten.
6. E5 and Theorem R are K = 1 only. For K ≥ 2 the exact projection is a QP, and Theorem F, Assumption M and the Fisher spectrum (three tiers, Problem 2) are open.

Not claimed: any improvement of C (it is as conservative as the constants in Theorem F-exp and Theorem C); the sharp constant when the box is not known a priori; ε-uniformity of C (C blows up as ε → 0 like (μ/ε)/(γδ²)² at least, from Step 2, as expected: it is the cost of the critical slowing-down).
Numerical status: E5's ingredients were checked (`p5_e5_check.py`); Theorem R's statement itself is consistent with the earlier slope measurement of E R_n against ln T_n (mean ≈ 0.9 vs d/2 = 1) but that was not rerun here.

---
## Update (25 Oct 2026): dependencies completed in `K1_completion.md`
Item 5 (optional stopping) is closed there and the restarts of §4 Step 3 are replaced by predictable good-gap indicators (one global telescoping). Assumption M (item 2) is closed for ε ≤ ε_M via the Gamma comparator of Theorem T (E_Γ h = 0.164 and 0.044 vs simulated 0.15 and 0.05). The integrability in Theorem F-exp is closed. Still assumed: the known box, K ≥ 2 open, constants conservative.
