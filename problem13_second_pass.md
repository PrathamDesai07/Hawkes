# Problem 13, second pass: resolving the four caveats of the first pass

Extends [problem13_tracking_bound.md](problem13_tracking_bound.md). The first pass proved Theorem TR with four weaknesses: (1) the assumptions A1–A4 were measured, not derived; (2) constants were 12–100× loose; (3) a jump cost of D²/η² instead of D²/η; (4) step-hold drift, K = 1 numerics. This note resolves 1 (for K = 1), 3 and 4 completely, and 2 partly, with the code checks in §7. Re-derive the key steps yourself; all numerical checks are K = 1 or K = 2.

## 0. Summary

| caveat | outcome |
|---|---|
| (3) jump rate | **resolved:** the bound is now ≍ J·D²/(η n) (sharp rate), via an L²-decomposition of the error into a noise part and a drift part and a Minkowski argument (§3). Measured scaling η^{−0.94}. |
| (4) step-hold / dimension | **resolved:** continuous drift inside gaps is handled exactly (§1), no step-hold assumption; the proof is dimension-free; checked numerically at K = 1 and K = 2 (§7). |
| (1) assumptions | **resolved for K = 1:** A1, A2 follow from A4 with explicit dependence on η/ε² (§5); A4 itself is proved for K = 1 (stationary parameters) by a window argument, and the window condition η ≲ ε²/(κ μ C) comes out of the proof (§6). **Not done for K ≥ 2.** |
| (2) constants | **partly resolved:** noise term now within 2× of its sharp value (that factor is the dropped −E w'ΔF w term); jump term sharp up to a factor ≈ 2 (measured ratio 0.30–0.42); c_F → 1/(c⁺κ) as η/ε² → 0 instead of being a run minimum. **Smooth-drift term remains about 25× loose** (§8: factor 4 from contracting the squared norm by 1−η instead of 1−2η, plus the triangle inequality, plus F-versus-Ī and 1/c_F). |

## 1. Exact recursion with continuous drift

Let θ*(t) be the true parameter (any measurable path), λ*(u) = φ(u)ᵀθ*(t_{i−1}+u) the true intensity in gap i, v_i = θ_i − θ*(t_{i−1}), d_i = θ*(t_i) − θ*(t_{i−1}) (the parameter change over gap i), s_i = b_i − a_i/λ*(w_i) (score increment under the true time-varying intensity), A_i = a_ia_iᵀ/(λ*(w_i)λ_{θ_i}(w_i)), G_i = F_i⁻¹A_i. Since λ_{θ_i}(w_i) − λ*(w_i) = a_iᵀ(θ_i − θ*(t_i)) = a_iᵀ(v_i − d_i):

g_i = s_i + A_i(v_i − d_i), and **v_{i+1} = (I − G_i)v_i − (I − G_i)d_i − F_i⁻¹s_i** (exact).

Two facts hold for *time-varying* truth and any F_{t_{i−1}}-measurable w, because P(τ_i > u) = e^{−Λ*(u)} with Λ*(u) = ∫₀^u λ*:
- (F1) E[w'A_iw | past] = ∫e^{−Λ*}(φ'w)²/λ_θ du = E[w'ΔF_iw | past] (Fubini on the survival function; this is Problem 5's Proposition 1(a) for arbitrary w, not only w = v_i);
- (F2) E[s_i | past] = 0 and E[s_is_i' | past] = E[ΔJ_i | past] with ΔJ_i = ∫φφ'/λ* (optional stopping for the true compensated martingale).

Pathwise: ‖G_id_i‖_{F_i} = |â'd||â'F_i⁻¹â|^{1/2} ≤ (â'F_i⁻¹â)‖d_i‖_{F_i} ≤ (ℓ̄/c⁻)‖d_i‖_{F_i}, with â = a/√(λ*λ_θ). So the forcing (I − G_i)d_i has F_i-norm at most (1 + ℓ̄/c⁻)‖d_i‖_{F_i} **even though d_i is not predictable** (it depends on w_i). This is why no step-hold is needed.

**Assumption (A1′, point leverage).** â_i'F_i⁻¹â_i ≤ ℓ̄ < 1 a.s.

## 2. Contraction and noise lemmas

**Lemma K (contraction).** For every F_{t_{i−1}}-measurable w: E[‖(I − G_i)w‖²_{F_i} | past] ≤ (1−η)‖w‖²_{F_{i−1}} − (1−ℓ̄)E[w'ΔF_iw | past] ≤ (1−η)‖w‖²_{F_{i−1}}.
*Proof.* ‖w − F⁻¹Aw‖²_{F_i} = ‖w‖²_{F_i} − 2w'Aw + w'AF_i⁻¹Aw; ‖w‖²_{F_i} = (1−η)‖w‖²_{F_{i−1}} + w'ΔFw; by (F1) E w'Aw = E w'ΔFw; and w'AF_i⁻¹Aw = (â'w)²(â'F_i⁻¹â) ≤ ℓ̄(â'w)², whose conditional expectation is ℓ̄E w'ΔFw by (F1). ∎

**Lemma N (noise).** Let w_i be predictable and define u_{i+1} = (I − G_i)u_i − F_i⁻¹s_i. Then E‖u_{i+1}‖²_{F_i} ≤ (1−η)E‖u_i‖²_{F_{i−1}} + (1−ℓ̄)⁻¹E[s_i'F_i⁻¹s_i], and E s_i'F_i⁻¹s_i ≤ q_i/(1−η) with q_i = E tr(F_{i−1}⁻¹ΔJ_i).
*Proof.* Expand ‖(I−G)u − F⁻¹s‖²_F. The cross term is −2E[u's] + 2E[u'AF_i⁻¹s]; the first vanishes by (F2), and |u'AF⁻¹s| = |â'u||â'F⁻¹s| ≤ |â'u|(ℓ̄)^{1/2}(s'F⁻¹s)^{1/2}, so 2|·| ≤ ρ″(â'u)² + (ℓ̄/ρ″)s'F⁻¹s with ρ″ = 1−ℓ̄; E(â'u)² = E u'ΔFu by (F1) is absorbed by the −(1−ℓ̄)E u'ΔFu term of Lemma K. Finally F_i ⪰ (1−η)F_{i−1} gives s'F_i⁻¹s ≤ s'F_{i−1}⁻¹s/(1−η), and E s's-expectation is q_i by (F2). ∎

## 3. Theorem TR2 (any K; sharp noise constant; sharp jump rate; continuous drift)

Split the actual error with the **same realised matrices**: v^ξ_{i+1} = (I−G_i)v^ξ_i − F_i⁻¹s_i (v^ξ_1 = v_1) and v^d_{i+1} = (I−G_i)v^d_i − (I−G_i)d_i (v^d_1 = 0). By linearity given the matrices, **v_i = v^ξ_i + v^d_i exactly**, and both are F_{t_{i−1}}-measurable. Write U^ξ_i = ‖v^ξ_i‖²_{F_{i−1}} etc.

**Theorem TR2.** Under (A0), (A1′) and the discounted-trace hypothesis (A2) of the first pass:
1. (noise) (1/n)ΣE U^ξ_i ≤ (1/(ηn))E U_1 + (1/((1−ℓ̄)η(1−η)))·(1/n)ΣE q_i, and by the discounted-trace lemma (1/n)ΣE q_i ≤ c⁺κ_x[d(−ln(1−η)) + Λ_n/n], so the noise part is ≤ (1+o(1))·c⁺κ_x d/((1−ℓ̄)(1−η)²) in F-units.
2. (drift) (E U^d_n)^{1/2} ≤ Σ_{i<n}(1−η)^{(n−1−i)/2}δ_i, δ_i = (1+ℓ̄/c⁻)(E‖d_i‖²_{F_i})^{1/2}. Consequently (1/n)ΣE U^d_i ≤ (1−√(1−η))⁻²(1/n)Σδ_i² for smooth drift, ≤ (1/(ηn))(Σδ_i)² in path-length form, and ≤ J δ²_max/(ηn)(1+o(1)) for J isolated jumps.
3. (total) For every ϑ > 0: (1/n)ΣE U_i ≤ (1+ϑ)·(noise bound) + (1+1/ϑ)·(drift bound).
With E‖d_i‖²_{F_i} = d_i'E[F_i]d_i ≤ (c_I/(ηc⁻))‖d_i‖²_Ī + (1−η)^i d_i'F_0d_i (A3) and (A4) ηF ⪰ c_FĪ, the Fisher-metric average is ≤ (η/c_F)×[total]:
**MSE_Ī ≲ (1+ϑ)(c⁺κ_x d/((1−ℓ̄)c_F))η + (1+1/ϑ)(c_I(1+ℓ̄/c⁻)²/(c_Fc⁻))·{4Δ_Ī²/η² (smooth) | JD_Ī²/(ηn)·(1+o(1)) (J jumps) | V_Ī²/(ηn) (path length V)}.**

*Proof.* (1) Lemma N with the recursion summed as in the first pass (E U^ξ_{i+1} ≤ (1−η)E U^ξ_i + c q_i), plus the discounted-trace lemma. (2) From v^d_{i+1} = (I−G_i)v^d_i − (I−G_i)d_i: pathwise ‖v^d_{i+1}‖_{F_i} ≤ ‖(I−G_i)v^d_i‖_{F_i} + (1+ℓ̄/c⁻)‖d_i‖_{F_i}; take L²(Ω) norms (Minkowski) and use Lemma K: u^d_{i+1} ≤ (1−η)^{1/2}u^d_i + δ_i with u^d_i = (E U^d_i)^{1/2}. Unrolling gives the convolution bound; Young's inequality with ℓ¹–ℓ² gives the three forms (ℓ¹ kernel for smooth, ℓ² kernel for path length: Σ_k(1−η)^k = 1/η). (3) The F_{i−1}-norm is a norm pathwise, so Minkowski in L²(Ω) gives (E U_i)^{1/2} ≤ (E U^ξ_i)^{1/2} + (E U^d_i)^{1/2}; square and use (a+b)² ≤ (1+ϑ)a² + (1+1/ϑ)b². ∎

**What changed relative to the first pass.** There is no (1+ζ) Young step on the Newton recursion, so the noise constant is not doubled; there is no 1/η Young factor on the drift, so an isolated jump of size D costs ≍ D²/η (the sharp rate; the first pass gave D²/η²); there is no step-hold assumption.

## 4. Why the jump rate is sharp and what remains for smooth drift

For one jump of size D let δ = (1+ℓ̄/c⁻)(E‖d‖²_F)^{1/2}, so δ² ≈ (1+ℓ̄/c⁻)²(c_I/(ηc⁻))D_Ī² (D_Ī = jump size in the Ī-norm). The drift part decays like (1−η)^{k/2} in norm, so Σ_k E U^d_k ≤ δ²/η, and the Fisher-metric cost per jump is at most (η/c_F)·δ²/η = δ²/c_F ≈ c_I(1+ℓ̄/c⁻)²D_Ī²/(c⁻c_Fη). The truth is D²/(2η) (contraction 1−2η in squared norm), so the bound is within a factor 2c_I/(c⁻c_F) of sharp. For smooth drift the same contraction loss is a factor 4 in the squared lag (u^d ≈ 2δ/η instead of δ/η).

## 5. A1′ and A2 follow from A4 (K = 1, explicit)

Assume (A4) ηF_{i−1} ⪰ c_FĪ with Ī = I/r, Θ := mE[λ] = μm/ε > 1 (true for k > 1: Θ → k/(k−1)). In (p, q)-coordinates u'Iu = p²m + 2pq + q²Eλ with p = u_μ − μu_α/α, q = u_α/α and φ = (1, s) ↔ vector (1, λ_φ), λ_φ = μ + αs; hence φ'I⁻¹φ = (Eλ − 2λ_φ + mλ_φ²)/(mEλ − 1). Since F_i⁻¹ ⪯ (η/((1−η)c_F))Ī⁻¹ = (ηr/((1−η)c_F))I⁻¹:

- **(A1′)** ℓ̄ ≤ [ηr/((1−η)c_Fc⁻(Θ−1))]·(Eλ/λ*² − 2/λ* + m) ≤ **[η/(ε²c_Fc⁻(Θ−1))]·(1 + O(ε))** (the bracket is decreasing in λ* on λ* < Eλ, so its supremum is at λ* = μ; and ηr/(εμ) = η/ε²).
- **(A2)** x̄ ≤ [η/((1−η)c_F)]·max_i tr(Ī⁻¹ΔF_i), tr(Ī⁻¹ΔF_i) ≤ (1/c⁻)∫φ'Ī⁻¹φ/λ* du ≤ r E_i(Eλ/μ² + m)/(c⁻(Θ−1)) with E_i = Λ*(τ_i) ∼ Exp(1) independent of the past, so **x̄ ≲ [η/(ε²c_Fc⁻(Θ−1))]·(ln n + x)** with probability ≥ 1 − e^{−x}.

Both are small exactly when η ≪ ε²: the window condition of Problem 14/15 is the condition under which the proof's leverage assumptions hold. (Numerically: ℓ̄ formula 0.062, 0.049, 0.233, 0.306 against measured maxima 0.065, 0.058, 0.232, 0.289 at (ε, η) = (0.4, 0.004), (0.2, 0.001), (0.2, 0.004), (0.4, 0.016); the gap-leverage bound with the 5 % quantile of c_F covers the measured maxima 0.72, 0.21, 1.10, 2.95 in all four cases.)

## 6. A4 for K = 1 (stationary parameters): the window argument

*Level of detail:* Steps 1, 2 and 4 are complete arguments. In Step 3 the Poisson-equation reduction, the Freedman step, the truncation and the small-ball bound are stated with their key inequalities, but the explicit constant C_* is not derived line by line (the numerical checks in §7 test the identity of Step 2, the Poisson-equation bound, and the scaling of the window information in T_0ε, not the final constant).

**Step 1 (discounted sum ⪰ time window).** Since ΔF_j ⪰ ΔJ_j/c⁺ (λ_θ ≤ c⁺λ*) and the last m_e events carry weights ≥ (1−η)^{m_e}, if the last m_e events cover the time window [t−T_0, t] then F_{i−1} ⪰ ((1−η)^{m_e}/c⁺)𝒥[t−T_0, t], 𝒥 = ∫aa'/λ* ds (PSD integrand: extra gaps only help).
**Step 2 (reduction to two window means).** For any u: u'𝒥u = T_0(p²x + 2pq + q²ℓ), x = ⟨1/λ⟩, ℓ = ⟨λ⟩ over the window (exact identity, checked to 10⁻⁸). So with Q_∞ = [[m,1],[1,Eλ]]: Q_w ⪰ cQ_∞ iff x − cm ≥ 0, ℓ − cEλ ≥ 0 and (x−cm)(ℓ−cEλ) ≥ (1−c)². If x ≥ (1−δ)m and ℓ ≥ (1−δ)Eλ this holds with **c = 1 − δ/(1 − Θ^{−1/2})** (positive iff δ < 1 − Θ^{−1/2}).
**Step 3 (concentration of window means).** *Poisson equation.* For f with |f(S)−f(S′)| ≤ L|S−S′| (K = 1: d_ℓ = |ΔS|, γ = βε), u(s) = ∫₀^∞(P_tf(s) − π f)dt is (L/γ)-Lipschitz by the exact contraction E|S_t − S′_t| ≤ e^{−γt}|s − s′| (the other session's Lemma C, also my Theorem 2 with rate γ), and ∫₀^T(f − πf)ds = u(S_0) − u(S_T) + M_T with jumps of M bounded by L/γ (|S+1−S| = 1) and ⟨M⟩_T ≤ (L/γ)²Λ_T, Λ_T = ∫λ. Freedman then gives |∫(f−πf)| ≤ (L/γ)[√(2vx) + x/3 + S_0 + S_T + E_πS] with probability ≥ 1 − 2e^{−x} on {Λ_T ≤ v}. *For ℓ:* f = λ, L = α. *For x:* 1/λ is not Lipschitz uniformly, so use g_z = 1/max(λ, z), z = c₁/ε, L_g = α/z² ∝ ε², and x ≥ ⟨g_z⟩; E_πg_z ≥ m − E[1/λ; λ < z], and the Gamma majorant of Theorem A (φ_ε ≤ P_ε) gives P(λ < z) ≤ e(c₁/(a + c₁))^{k_ε} and, integrating min(φ, P(λ<z)) over u, E[1/λ; λ < z] ≤ (ε/a)(k/(k−1))p_z^{1−1/k} (a = α²/(2β)), which is a vanishing fraction of m ≈ ε·2/(2μ−β) as c₁ → 0. The requirement is T_0 ≥ C_*(δ, x)/ε with C_* ≍ max(8x/(δ²μ), 32xμ/(δ²c₁⁴κ_m²)) (m = κ_mε), plus boundary terms of the same order. *Covering:* N(window) ≤ m_e := ⌈(1+δ_N)rT_0⌉ with probability ≥ 1 − e^{−x} (the integral ∫λ over the window is a sum of Exp(1) variables).
**Step 4 (the window condition).** With θ_w = ηT_0r: c_F ≥ θ_we^{−κ_Nθ_w}c_win/c⁺ ≥ c_win/(eκ_Nc⁺) at θ_w = 1/κ_N, i.e. T_0 = 1/(κ_Nηr). This needs T_0 ≥ C_*/ε, that is
**η ≤ ε²/(κ_N μ C_*).**
In the limit η/ε² → 0 with the window length growing, c_win → 1 and c_F → 1/(c⁺κ) (the discount time 1/(ηr) then far exceeds the mixing time 1/γ and the exponentially weighted version of Step 3 applies directly).

## 7. Numerical verification (K = 1 and K = 2, α unknown, continuous drift; code in `code/p13b_*.py`)

| check | result |
|---|---|
| exact split v = v^ξ + v^d (after resynchronising 40–437 clipped steps of ≈ 6·10⁴) | max error 3·10⁻¹⁵ (K = 1), 5·10⁻¹⁵ and 8·10⁻¹⁵ (K = 2) |
| Lemma K, w = v^d / v^ξ / a random predictable w, K = 1 (ε = 0.4, η = 0.004) | LHS − RHS = −0.00004 ± 0.00002, −0.00013 ± 0.00009, −0.00008 ± 0.00033 |
| Lemma K, K = 2 (η = 0.002) | −0.00005 ± 0.00004, −0.00027 ± 0.00014, −0.00084 ± 0.00216; (η = 0.001) −0.00012 ± 0.00007, −0.00007 ± 0.00003 |
| point leverage (A1′) | mean 0.008, max 0.079 (K = 1); mean 0.006, max 0.059 (K = 2, η = 0.002) |
| Minkowski, K = 1 | √(mean U) = 1.264 ≤ √U^ξ + √U^d = 1.749 |
| noise bound (Lemma N), K = 1 | mean U^ξ = 1.108 vs bound 2.140 (ratio 0.52); K = 2, η = 0.001: 1.698 vs 3.083 (0.55). The remaining factor ≈ 2 is the dropped −(1−ℓ̄)E w'ΔF w ≈ −ηU term. (At η = 0.002, K = 2, U^ξ = 3.77 exceeded the bound 3.11; this run had 1 % clipped steps, and the diagonal clip is not the F-norm projection assumed by the theory, so the excess is a diagnostic artefact; at η = 0.001, 0.06 % clipped, the bound holds.) |
| sharp jump rate (D = 0.3, 800 time units per regime, ε = 0.4), η = 0.002, 0.004, 0.008 | per-jump Σ U^d = 1379, 369, 94 vs bound 3299, 974, 313 (ratios 0.42, 0.38, 0.30); per-jump Fisher-metric cost 3.55, 1.87, 0.96, log-log slope **−0.94** (sharp rate −1; the first-pass bound has −2) |
| window information c_win(T_0 = C/ε) (K = 1) | ε = 0.4: median / 1 % / min = 0.66 / 0.26 / 0.15 at C = 16 and 0.84 / 0.55 / 0.45 at C = 64; ε = 0.2: 0.66 / 0.25 / 0.15 and 0.84 / 0.53 / 0.43: **the same distribution at equal C = T_0ε**, confirming T_0 ∝ 1/ε (m ∝ μ/ε² events) |
| quadratic-form reduction | u'𝒥u = 16347.237177 vs p²X + 2pqT + q²L = 16347.237293 |
| Poisson-equation bound (K = 1, ε = 0.4) | sd of window integrals of λ: 7.3, 15.0, 30.4 vs bound 18.8, 26.3, 41.3 at T_0ε = 4, 16, 64 |
| leverage formulas (§5) | see §5 (within ≈ 20 % of the measured maxima) |

## 8. What is still not established

- **A4 for K ≥ 2.** The reduction of Step 2 is specific to K = 1 (a 2×2 matrix depending on two window means). For K ≥ 2 the same Poisson-equation argument applies entrywise to the window information, but whitening by Ī requires relative accuracy in the tiny-eigenvalue directions (ridge and shape modes of Problem 12), and I did not do that perturbation analysis. A1′ and A2 then also remain assumptions for K ≥ 2 (measured: K = 2, leverage max 0.06, c_F median 0.8–0.88, min 0.19–0.39).
- **A4 with drifting parameters.** Step 3 uses the time-homogeneous dynamics (constant θ*); for a drift box I only have the pathwise comparison λ*/λ_ref ≤ c_box (which transfers the lower bound with a factor 1/c_box) and not the concentration of the window means under non-stationary dynamics.
- **Constants.** C_* contains 1/c₁⁴ and is large; the proved c_F ≥ c_win/(eκ_Nc⁺) is conservative at finite η/ε² (measured median c_F is 0.8). The smooth-drift term remains ≈ 25× loose at the tested settings (4 from squared-norm contraction, ≈ 2 from (1+1/ϑ), ≈ 1.3 from F vs Ī and 1/c_F, ≈ 2.6 from the triangle inequality ignoring sign cancellation of an oscillating drift); sharpening it needs a lower bound on E w'ΔF_iw in terms of ηU_i, i.e. a decoupling between the current gap's information and the error direction, which the measured correlation structure (both depend on the slow state) does not give for free.
- **Projection.** The theory assumes the unprojected iteration or the exact F-norm projection of Problem 4; the diagonal clip used in these simulations is neither, and K = 2 runs at η = 0.002 clipped 1 % of steps.
- **Almost-sure assumptions.** (A1′), (A2) are used as a.s. bounds; turning them into high-probability statements uses the stopping argument of Problem 5's Theorem HP, not repeated here.
- The Gamma majorant (Theorem A) used in Step 3 exists only for K = 1.
