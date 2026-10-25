# Problem 5 / 7 follow-up: high-probability bounds on the accumulated Fisher information (K = 1)

Goal: turn the a.s. leverage assumptions (L2), (L3) of Theorem 7 (`problem07_leverage.md`) into high-probability statements. Scope of what is proved here: the **oracle** accumulated information J_T = ∫_0^T a aᵀ/λ* dt (data from the true law, evaluated at θ*), with Chebyshev-level tails, under one numerically certified assumption. Not proved: (L1) for the iterates (see "Circularity").
Code: `code/p5_fisher_conc.py`, `p5_variance_check.py`, `p5_leverage_from_occupation.py`.

## Notation (K = 1, stationary start)

γ = βε (mixing rate, Problem 3), ES = μ/(βε), σ_S = (E S² + 3(E S)²)^{1/2} (both moments are closed form: E λ = μ/ε, E λ² = μ(2βμ + α²)/(2βε²), S = (λ − μ)/α).
Two disjoint intervals I_j = [c_j − d, c_j + d], j = 1,2, with centres c_1 < c_2 and h_j: R_+ → [0,1] a trapezoid supported in I_j with ramps of width ρ, so Lip(h_j) = 1/ρ. m_j := E_π h_j.

## Lemma F1 (variance of an occupation functional) — rigorous given Problem 3's contraction

For h with values in [0,1] and Lipschitz constant Λ_h, and S stationary,
Var(∫_0^T h(S_t) dt) ≤ 2 T Λ_h √m σ_S / γ, m = E_π h.
Proof. Cov(h(S_0), h(S_t)) = E[(h(S_0) − m)(P_t h(S_0) − m)] and |P_t h(s) − m| = |E h(S^s_t) − E h(S̃_t)| ≤ Λ_h E|S^s_t − S̃_t| ≤ Λ_h e^{−γt}(s + E S), using the monotone-coupling contraction E|S¹_t − S²_t| ≤ e^{−γt}|s¹ − s²| of Problem 3 (Lemma C; for K = 1, ℓ = 1) and conditioning on S̃_0. Cauchy–Schwarz with Var h ≤ m(1−m) ≤ m gives |Cov(t)| ≤ Λ_h √m σ_S e^{−γt}; and Var(∫_0^T h) ≤ 2T ∫_0^∞ |Cov(t)| dt. ∎
Numerical check (120 reps; trapezoid in I_1 and I_2, ε = 0.3, 0.1; Tγ = 5, 20): Var/bound = 0.0006–0.002, so the lemma holds with a 500–1500× margin; Var grows linearly in T (×3.5 from Tγ = 5 to 20) as the lemma predicts.

## Lemma F2 (two-interval lower bound on J_T) — rigorous, deterministic

Let O_j = ∫_0^T 1{S_t ∈ I_j} dt and L = max of λ* over I_1 ∪ I_2 = μ + α(c_2 + d). If Δ = c_2 − c_1 > 2√2 d then for every v ∈ R²,
vᵀJ_T v ≥ (min_j O_j / L) · vᵀ Q v, Q = [[1, (c_1+c_2)/2], [(c_1+c_2)/2, (c_1² + c_2²)/2 − 2d²]] ≻ 0.
Proof. On I_j, vᵀa = (v_0 + v_1 c_j) + v_1(s − c_j), so (vᵀa)² ≥ ½(v_0 + v_1 c_j)² − v_1² d²; divide by λ* ≤ L, sum over j with weights O_j ≥ min O. The resulting 2 × 2 form has determinant Δ²/4 − 2d² > 0. ∎
(An earlier, cruder version, with the interval half-width d treated as an additive error, fails at scale 1/ε because d² ≫ 1; the exact form above is what makes the scale-1/ε intervals usable.)

## Theorem F (K = 1): J_T concentrates at the mixing time scale — rigorous modulo Assumption M

Assumption M (certified numerically): m := min_j E_π h_j > 0. (Computable to any accuracy from the stationary law; here from long simulations: m ≈ 0.15 for I_1 and 0.05 for I_2 at ε = 0.3, 0.1, 0.03 with centres 0.5·ES and 1.8·ES, half-width d = 0.15·ES, ramp ρ = d/2.) A rigorous lower bound on E_π h_j (anti-concentration of the stationary law) is not proved.

Statement. For q ∈ (0,1) and T_0 ≥ 16 σ_S / (ρ γ m^{3/2} q), with probability ≥ 1 − 2q (simultaneously for all T ≥ T_0),
J_T ⪰ (m T / (4L)) · Q.
Proof. Chebyshev with Lemma F1: P(∫_0^T h_j ≤ Tm/2) ≤ 8Λ_hσ_S/(γ T m^{3/2}) =: p(T). Take T_k = 2^k T_0; Σ_k p(T_k) ≤ 2p(T_0) (union over j gives the factor 2, so total ≤ 4p(T_0) ≤ 2q with the stated T_0 up to constants). O_j ≥ ∫h_j and J_T is nondecreasing in T, so for T ∈ [T_k, 2T_k], J_T ⪰ J_{T_k} ⪰ (mT_k/(2L))Q ⪰ (mT/(4L))Q by Lemma F2. ∎

Numerics. T_0 γ for q = 0.1, m = 0.05 is 4×10⁵, the same for ε = 0.3, 0.1, 0.03 (everything scales in units of 1/γ). Empirically λ_min(J_T)/T is within a factor ≈ 2 of its limit (0.13, 0.067, 0.024 at ε = 0.3, 0.1, 0.03) by Tγ ≈ 20–50. So the rigorous threshold is about 10⁴× the empirical one; the dependence on ε (T_0 ∝ 1/γ = 1/(βε)) is the right one. The resulting λ_min bound mT λ_Q/(4L) is ≈ 100× below the truth and scales as ε^{1.1} (truth ε^{0.75} over this range).

## Corollary (leverage under (L1)): the a.s. assumptions of Theorem 7 hold w.h.p.

On the event of Theorem F, if in addition the iterates satisfy λ_{θ_i} ≤ c⁺λ* and ≥ c⁻λ* on [0, t_{i−1}] (so F_{i−1} ⪰ J_{t_{i−1}}/c⁺ ⪰ (m t_{i−1}/(4c⁺L)) Q), then for t_{i−1} ≥ T_0:
- (L3, point leverage) ℓ̄_i ≤ (4c⁺L/(m t_{i−1})) · sup_{s ≥ 0} φ(s)ᵀQ⁻¹φ(s)/(μ + αs), φ(s) = (1, s), with sup over s ∈ [0, s_max] finite. The dominant term is the quiet limit s → 0: ℓ ≈ (4c⁺L/(m t)) · (c_1² + c_2² − 4d²)/(2 det Q · μ) = O(1/(ε t m)) since L ∼ μ/ε and the ratio is O(1/μ) in units where c ∼ 1/ε. The case s → ∞ is controlled by s_max ≤ ES · log(n/q) w.h.p. (exponential tail of π).
- (L2, gap leverage) x_i ≤ tr(F_{i−1}⁻¹ΔF_i) ≤ ℓ̄_i τ_i/c⁻ and P(τ_i > t) ≤ e^{−μt}, so with prob ≥ 1 − q, τ_i ≤ ln(N/q)/μ for all i ≤ N events: x̄ ≤ ℓ̄ ln(N/q)/(μ c⁻).
Numerical check of the pathwise linear algebra (oracle J_T from simulation; 60 reps; ε = 0.3, 0.1, 0.03; Tγ ∈ {2, 5, 10, 20, 50}): ℓ_orc(T) ≤ ℓ_bd(T) in all ≈ 1300 checks (0 violations). ℓ_bd is 26–115× larger than ℓ_orc. ℓ_orc decays like 1/T (medians at Tγ = 50: 0.058, 0.031, 0.025), consistent with the (k−1)/(εT) heuristic of `problem07_leverage.md`.

## What this does and does not settle

- Settled (K = 1, oracle information): J_T grows linearly with the correct ε-scaling at the mixing time scale 1/(βε), with explicit constants; the point and gap leverage then satisfy the assumptions of Theorem 7 w.h.p. after T ≳ T_0.
- Not settled — **circularity**. (L1) (iterates stay in a ratio set around θ*) is needed to pass from J_T to the actual F_{i−1}, and (L1) itself is a consistency statement that uses the information bound. A proof needs a stopping time ν = first index where (L1) fails, an induction showing the information bound ⇒ error bound ⇒ (L1) persists, and control after ν. Not attempted.
- Not settled: a rigorous lower bound on m (Assumption M); exponential rather than Chebyshev tails (would reduce T_0's dependence on q to log(1/q) and probably the constants by orders of magnitude, via a Bernstein-type inequality for the occupation functional); the projected algorithm; K ≥ 2 (Lemma F2 needs d + 1 well-separated intervals in R^K, and the contraction used in F1 is the K-dimensional Lemma C, which does hold for K ≥ 2).
- The rigorous constants are far from tight (10⁴× in T_0, ~100× in the information level), so this should be read as a proof of the right scaling, not a usable bound.

Re-derive Lemmas F1, F2 and the proof of Theorem F yourself. I wrote each in one pass; the numerics support F1 and the pathwise step of the corollary but do not independently verify the Chebyshev/union-bound bookkeeping.
