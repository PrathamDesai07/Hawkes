# Completing the dependencies of Theorem R (K = 1)

Request: "complete that thing" = the list "What Theorem R still rests on" in `problem05_E5_expected_regret.md` §5. Result per item: **closed**, **reduced**, or **not closable here**. Code: `code/p5_assumption_M_gamma.py`.

## Item 4 (optional stopping across random restarts): closed, and the restarts are unnecessary

Work with the gap filtration 𝒢_i = σ(first i gaps). Every "good gap" condition used in §4 is 𝒢_{i−1}-measurable:
 good_i := {θ_i ∈ R(δ)} ∩ {λ_min(F_{i−1}) ≥ λ_{t_{i−1}}} ∩ {ℓ̄_i ≤ ℓ̄_{t_{i−1}}} ∩ {S_{t_{i−1}} ≤ s_{t_{i−1}}}.
(θ_i and F_{i−1} are 𝒢_{i−1}-measurable; ℓ̄_i is a function of F_{i−1} and S_{t_{i−1}}.) So 1_{good_i} is a **predictable** indicator and no stopping time or restart is needed:

1. *Telescoping, once, globally.* Lemma C1 from i_0 = 1: D_n ≤ D_0 − X_n − M_n + G_n with all terms defined over all gaps. Because E[m_i | 𝒢_{i−1}] = 0 and E[m_i² | 𝒢_{i−1}] ≤ (4/c⁻)E[v_iᵀΔF_iv_i | 𝒢_{i−1}], with E X_n < ∞ (the box and E ΔF_i < ∞ bound it by n·const), M_n is a square-integrable martingale and E M_n = 0 without any optional-stopping theorem.
2. *Split G.* G_n = Σ_i g_iᵀF_i⁻¹g_i 1_{good_i} + Σ_i g_iᵀF_i⁻¹g_i 1_{bad_i}. The good sum is a sum of the increments Q_i, h-terms and conversions of E5b–d multiplied by the predictable indicator 1_{good_i}. Each supermartingale used there (Freedman for Q', Ville for the log-det conversion, Ville for N^A) stays a supermartingale when its increments are multiplied by a predictable indicator in [0,1] (the conditional MGF bound is applied to the product, and the bound on the truncated variable is unchanged). The hypotheses of E5 are used only through the bounds ℓ̄_i ≤ ℓ̄_{t_{i−1}}, λ_θ ∈ [c⁻λ*, c⁺λ*] and c⁺ ≤ 1+δ on good gaps, which is exactly what good_i supplies.
3. *Time-varying ℓ̄.* E5 was stated with a constant ℓ̄ on [i_0, n]. Replace it by ℓ̄_i ≤ C_ℓ ln t_{i−1}/(εt_{i−1}) (Theorem F: λ_min(F) ≳ mλ_Q t/(4c⁺L); Lemma E0: S ≤ O(ln t/ε); the point leverage grows linearly in S). Use truncation levels m_i = 2 ln i + x (Σ_i e^{−m_i} ≤ 2e^{−x}). Then B̄_i = (2ℓ̄_i/μ*)(1+m_i²) = O((ln i)³/i) → 0, and in E5b one takes B̄ := sup_{good i ≥ i_1}B̄_i, which is a number that tends to 0 as the index i_1 from which the sharp constants are used tends to infinity. Likewise x̄_i, b → 0.
4. *Bad part.* g_iᵀF_i⁻¹g_i ≤ |g_i|²/λ_min(F_0) with E|g_i|⁴ < ∞ (|g_i| ≤ (w_i, b_{1,i}, 1/μ_lo) and moments of w_i, b_{1,i} are bounded by E0), so E[g_iᵀF_i⁻¹g_i 1_{bad_i}] ≤ C(P(bad_i))^{1/2} ≤ C'e^{−c_1γδ²t_{i−1}/2}, and Σ_i of that is ≲ (μ/ε)·C'/(c_1γδ²) =: C_bad, independent of n.
5. *Regret versus X.* E[r_i | 𝒢_{i−1}] ≤ (½ + C_2δ)E[vᵀΔFv | 𝒢_{i−1}] on {θ_i ∈ R(δ)} and ≤ C_3E[vᵀΔFv | 𝒢_{i−1}] otherwise (box). The bad gaps contribute at most C_3Σ_i E[vᵀΔFv 1_{bad_i}] ≤ C_3·diam²·Σ_i (E‖ΔF_i‖²)^{1/2}P(bad_i)^{1/2} = O(C_bad).

Consequently the argument of §4 runs as one global telescoping with predictable indicators; Step 3 (restarts) of Theorem R is replaced by items 1–5 above. The log-determinants need no splitting: Σ_i ℓ'_i ≤ log det(F_n/F_0).

## Item 2 (Assumption M): closed for small ε (given Theorem T)

Theorem T, Step 7 gives TV(π, Γ) ≤ Cε^{min(1,k−½)} for the Gamma comparator Γ = Gamma(k, scale a/ε) in λ-space (k₀ > 1). The trapezoid h_j of `problem05_fisher_concentration.md` (support [c_j−d, c_j+d] in units of ES, plateau [c_j−d/2, c_j+d/2], ramps d/2; c = 0.5, 1.8; d = 0.15) satisfies 0 ≤ h ≤ 1, so
 E_π h_j ≥ E_Γ h_j − TV(π,Γ) ≥ E_Γ h_j − Cε^η.
E_Γ h_j is an explicit incomplete-gamma integral. Computed (`p5_assumption_M_gamma.py`, μ = β = 1, k₀ = 2): I_1: 0.1365 (ε = 0.3), 0.1535 (0.1), 0.1607 (0.03), 0.164 (limit); I_2: 0.050, 0.0462, 0.045, 0.0444 (limit). Simulation gave 0.15 and 0.05, consistent. Hence **Assumption M holds with m = ½·min(E_Γ h_j) ≈ 0.022 for all ε ≤ ε_M(β, μ)**, where ε_M is determined by Cε^η ≤ ½·0.044, with C from Theorem T (not explicit).
Limits of this closure: (i) it inherits Theorem T, which is a proof outline in which every analytic estimate was verified numerically but the chain of constants was not formalised; (ii) it is asymptotic in ε (C unspecified), so for a given finite ε it is not a certificate; the finite-ε values above are numerical; (iii) k₀ ≤ 1 is not covered (E_π[1/λ] no longer O(ε), different regime). For k₀ ≤ 1 the same Gamma mass is positive but the TV comparison in Theorem T needs q > ½.

## Item 1 (integrability in Theorem F-exp): closed

The only integrability needed in E3 is that 𝓜_T = ∫(û(S_{t−}+1) − û(S_{t−}))(dN − λdt) is a martingale and Dynkin's formula holds for û.
- û(s) = ∫_0^∞(P_th(s) − m)dt is well defined: by the monotone-coupling contraction (Lemma C of Problem 3, K = 1) |P_th(s) − P_th(s')| ≤ Λ_he^{−γt}|s−s'|, so |P_th(s) − m| ≤ Λ_he^{−γt}(s + ES) (integrate s' against π; E_π S = ES), hence |û(s)| ≤ Λ_h(s+ES)/γ and û is Λ_h/γ-Lipschitz.
- 𝓜 has predictable integrand bounded by b_0 = Λ_h/γ, and E Λ_T = ∫E λ_t dt < ∞ (E0, or the linear moment equation), so ⟨𝓜⟩_T ≤ b_0²Λ_T is integrable and 𝓜 is a square-integrable martingale.
- Dynkin: û is Lipschitz and 𝓛û = −(h−m) holds in the sense of the resolvent (Fubini on ∫_0^∞P_t𝓛... is justified by the same bound), and E|û(S_t)| ≤ Λ_h(E S_t + ES)/γ < ∞; so ∫_0^T(h−m)dt = û(S_0) − û(S_T) + 𝓜_T holds almost surely (û ∈ domain of the extended generator because it is Lipschitz and the jump part is controlled by E Λ_T < ∞).
The numerical rate c of E3 (≈ 10⁵ below observation) is a constants issue, not an integrability issue, and is untouched.

## Item 6 (ε-dependence of C): computed, order only

From Step 2 and item 4 above: C_bad ≈ C_3 diam²·(μ/ε)·(C_1/(c_1γδ²)) = O(μ/(βε²δ²)); the pre-t_* period costs at most (μ/ε)t_* with t_* = O(d ln(1/(εδ))/(γδ²)), i.e. O(μ d ln(1/ε)/(βε²δ²)); the sharp phase contributes ½d·ln(T_n/T_1) plus o(1) terms. Thus
 **C(ε) = O( (μ/(βε²δ²)) · d·ln(1/(εδ)) )**
(c, c_1, C_1 are ε-independent in the K = 1 scaling: the bump masses m and the normalised Lipschitz constant κ = Λ_h·ES do not depend on ε). The ε⁻² is the same order as the leverage burn-in (≈ 1/ε² events) and as the mixing time in events, so it is the cost of the critical slowing-down and cannot be avoided by this method. The log-determinant coefficient d/2 is ε-free, so the regret is (d/2)ln T_n + O(ε⁻²ln(1/ε)); this matches the E R_n measurements (slope ≈ 0.9) at fixed ε but the ε⁻² law for C was not measured.

## Item 3 (known box): not closable here

The box [0.5μ*, 2μ*] × [0.5α*, β(1−ε/2)] enters (a) through c⁻ ≥ 0.5, c⁺ ≤ 2, which make the Newton regret and the martingale bounds finite, and (b) through the a-priori bound used in the circularity and first-gap arguments. Removing it needs a self-bounding argument (a first stage that localises θ* inside a data-dependent region), which is not available without extra assumptions. What can be said: only the *ratios* μ_max/μ*, α_min/α* matter, i.e. the assumption is "a constant-factor prior guess of (μ*, α*)", and (L1)-(L3) never use the box beyond these constants. Stated as an assumption of the theorem, not removed.

## Item 5 (K ≥ 2): not closable here

For K ≥ 2 the projection is a QP (not O(1) exact), Theorem F needs d+1 well-separated intervals in R^K and a K-dimensional Assumption M, and the Fisher spectrum has three tiers (Problem 2), so ℓ̄ involves the ridge direction with eigenvalue ≈ ε/(μ(k−1)). The ingredients of the present note that do carry over verbatim: Lemmas E5a–E5d (they use only λ* ≥ μ* and the predictable structure), the predictable-indicator argument (item 4), and the monotone-coupling contraction. What does not: Theorem T (no Gamma comparator), exact projection, Assumption M. Open.

## Net status of Theorem R (K = 1, k₀ > 1)

| Dependency | Before | Now |
|---|---|---|
| Optional stopping / restarts | unwritten | closed; restarts replaced by predictable indicators |
| Integrability in Theorem F-exp | "modulo integrability" | closed |
| Assumption M | numerical | closed for ε ≤ ε_M, via Theorem T (itself an outline with numerically verified estimates) |
| C's ε-dependence | not computed | O(μ d ln(1/(εδ))/(βε²δ²)) (order of magnitude, not a proved sharp rate) |
| Known box | assumed | still assumed |
| Constants (rate c, ε_M, C) | conservative | unchanged |
| K ≥ 2 | open | open |

## Update (25 Oct 2026, later): Assumption M with explicit ε-range
Theorem T's constant chain is now explicit (`thmT_explicit_constants.md`). Assumption M is certified for k₀ = 2, μ = β = 1 and ε ≤ 10⁻⁴ (m ≥ 0.0248). The ε-range is far smaller than the simulation range (ε = 0.3…0.03); for those ε the assumption remains numerical.
