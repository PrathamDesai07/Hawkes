# Problem 13: tracking bound for a constant-gain Fisher-preconditioned estimator

> Later passes: `problem13_second_pass.md`, `problem13_third_pass.md`.

**Update (second pass):** the four weaknesses listed under "What is established and what is not" below are addressed in [problem13_second_pass.md](problem13_second_pass.md): sharp jump rate D²/η (Minkowski decomposition), continuous drift without step-hold, dimension-free proof checked at K = 1 and 2, A1/A2 derived from A4 for K = 1, A4 proved for K = 1 (stationary parameters) by a window argument that yields the condition η ≲ ε²/(κμC), and the noise constant improved about 2×. Still open: A4 for K ≥ 2, the smooth-drift constant (≈ 25× loose), and a.s. versus high-probability versions.

**Status: first pass.** A theorem (Theorem TR) under stated assumptions, proved in the realised-information norm, with every step of the chain checked numerically (K = 1, α unknown). The bound has the right form and the right ε-dependence but its constants are conservative by one to two orders of magnitude; the assumptions are on realised leverage and information and are *measured*, not derived, here. Built on Problem 5 (Proposition 1, Lemma JD, leverage lemmas) and Problem 14 (drift class), and it is what Problem 15 was waiting for.

## Setting

Event-indexed. Gap i has length w_i, a_i = φ(w_i), b_i = ∫₀^{w_i}φ, with φ(u) = (1, S(t_{i−1}+u)). Recursion with forgetting factor 1 − η:

F_i = (1−η)F_{i−1} + ΔF_i(θ_i), ΔF_i = ∫₀^{w_i} φφᵀ/λ_{θ_i} du, g_i = b_i − a_i/(a_iᵀθ_i), θ_{i+1} = θ_i − F_i⁻¹g_i

(unprojected; for the projected version use the exact F_i-norm projection of Problem 4, which is non-expansive in the F_i-norm, so the first step below is unchanged, as in Problem 5's Lemma C1).

**Drift.** The true parameter in force on gap i is θ*_i, F_{t_{i−1}}-measurable ("step-hold": the parameter is constant on each gap and changes at event times; for a continuous drift this approximates the true path by an error of order Δ·w_i within a gap). d_i = θ*_{i+1} − θ*_i. This contains Problem 14's class P(Δ, J, D) in discretised form.

Notation: v_i = θ_i − θ*_i, U_i = ‖v_i‖²_{F_{i−1}}, s_i = b_i − a_i/λ*_i (score increment), h_i = g_i − s_i = a_i a_iᵀv_i/(λ*_iλ_{θ_i}), ΔJ_i = ∫φφᵀ/λ* du, and Ī a fixed reference information matrix per event (for K = 1, Ī = I(θ*_ref)/r).

## Assumptions

- **(A0)** θ_i is F_{t_{i−1}}-measurable; c⁻ ≤ λ_{θ_i}/λ*_i ≤ c⁺ along the path.
- **(A1) point leverage.** E[|h_i|²_{F_i⁻¹} | past] ≤ B_h·E[v_iᵀΔF_iv_i | past] with (1 + 1/ρ)B_h ≤ 1 for a chosen ρ > 0. (Problem 5's Lemma E5c gives B_h in terms of the point leverage ℓ̄; here it is measured.)
- **(A2) gap leverage.** λ_max(((1−η)F_{i−1})^{−1/2}ΔF_i((1−η)F_{i−1})^{−1/2}) ≤ x̄ a.s. (In practice a high-probability event, handled by stopping as in Problem 5's Theorem HP.)
- **(A3) information envelope.** E[ΔJ_j] ⪯ c_I Ī for all j (unconditional; stationary-type bound, with c_I absorbing the box constant κ_I of Problem 14).
- **(A4) information lower bound (only for the Ī-norm corollary).** ηF_{i−1} ⪰ c_F Ī.

## Theorem TR (expectation, F-norm form)

Let 0 < η ≤ ½, ζ = η/2, C_q = (1+ζ)(1+ρ)/(1−η), C_d = 1 + 1/ζ = 1 + 2/η, κ_x = (1−η)x̄/log(1+x̄). Under (A0)–(A3), for every n,

**(1/n) Σ_{i=1}^n E‖v_i‖²_{F_{i−1}} ≤ (2/(ηn)) E U_1 + (2/η)·{ C_q c⁺ κ_x [ d(−ln(1−η)) + Λ_n/n ] + C_d [ (c_I/(η c⁻))·(1/n)Σ‖d_i‖²_Ī + (1/n)Σ(1−η)^i d_iᵀF_0 d_i ] }**,

with Λ_n = E[log det F_n − log det F_0]⁺ ≤ d ln( (c_I trĪ/(c⁻η) + tr F_0)/d ) − log det F_0 = O(d) when F_0 ≍ Ī/η.

**Corollary (Ī-norm, with A4).** (1/n)ΣE‖v_i‖²_Ī ≤ (η/c_F)·(right-hand side). For η small the leading terms are
**(2(1+ρ)c⁺κ_x d/c_F)·η + (4c_I/(c_F c⁻))·(1/n)Σ‖d_i‖²_Ī/η²** + transients.
This is the C5 target η σ² + Δ²/η² (σ² = 2(1+ρ)c⁺κ_x d/c_F in the Ī-metric); for smooth drift (1/n)Σ‖d_i‖²_Ī ≈ Δ_Ī² (per-event drift in the Ī-norm); for J jumps of size D_Ī it is JD_Ī²/n, giving the conservative term 4c_IJD_Ī²/(c_Fc⁻η²n) (the sharp rate is conjecturally JD²/(ηn), see below).

## Proof

1. **Newton algebra (exact).** v_{i+1} = v_i − F_i⁻¹g_i − d_i, and ‖v_i − F_i⁻¹g_i‖²_{F_i} = ‖v_i‖²_{F_i} − 2⟨g_i,v_i⟩ + g_iᵀF_i⁻¹g_i, with ‖v_i‖²_{F_i} = (1−η)U_i + v_iᵀΔF_iv_i (because F_i = (1−η)F_{i−1} + ΔF_i).
2. **Moving target (Young).** U_{i+1} = ‖v_{i+1}‖²_{F_i} ≤ (1+ζ)‖v_i − F_i⁻¹g_i‖²_{F_i} + (1+1/ζ)‖d_i‖²_{F_i}.
3. **Conditional expectation.** By Proposition 1(a) of Problem 5 (gap i has constant true parameter θ*_i and v_i is measurable), E[⟨g_i,v_i⟩ | past] = E[v_iᵀΔF_iv_i | past]. Hence E[U_{i+1}|past] ≤ (1+ζ)[(1−η)U_i − E[vᵀΔFv] + E[gᵀF_i⁻¹g]] + C_dE[‖d_i‖²_{F_i}].
4. **Gradient norm.** |g|² ≤ (1+ρ)|s|² + (1+1/ρ)|h|² in the F_i⁻¹ norm. Since F_i ⪰ (1−η)F_{i−1}: E|s|²_{F_i⁻¹} ≤ E|s|²_{F_{i−1}⁻¹}/(1−η) = q_i/(1−η) with q_i = tr(F_{i−1}⁻¹ E[ΔJ_i|past]) (information identity E[ss^T|past] = E[ΔJ_i|past], Lemma JD(1)); by (A1) and (1+ζ)(1+1/ρ)B_h ≤ 1 the −E[vᵀΔFv] term absorbs the h part, giving E[U_{i+1}|past] ≤ (1−η/2)U_i + C_q q_i + C_dE[‖d_i‖²_{F_i}|past], using (1+η/2)(1−η) ≤ 1 − η/2.
5. **Summing.** Taking expectations and summing: (η/2)Σ_{i≤n}E U_i ≤ E U_1 + C_qΣE q_i + C_dΣE‖d_i‖²_{F_i}.
6. **Discounted trace lemma.** Since ΔJ ⪯ c⁺ΔF, E q_i ≤ c⁺E tr(F_{i−1}⁻¹ΔF_i). Put M_i = ((1−η)F_{i−1})^{−1/2}ΔF_i((1−η)F_{i−1})^{−1/2} ⪰ 0 with eigenvalues ≤ x̄ (A2). Then log det F_i − log det((1−η)F_{i−1}) ≥ log det(I + M_i) = Σ log(1+μ_j) ≥ (log(1+x̄)/x̄) tr M_i (concavity, μ_j ≤ x̄), so tr(F_{i−1}⁻¹ΔF_i) ≤ (1−η)(x̄/log(1+x̄))[log det F_i − log det F_{i−1} + d(−ln(1−η))]. Summing telescopes to log det(F_n/F_0) + nd(−ln(1−η)); take expectations and use E log det F_n ≤ d ln(E tr F_n/d) (Jensen) with E tr F_n ≤ (c_I trĪ)/(c⁻η) + tr F_0 (since ΔF ⪯ ΔJ/c⁻ and (A3)).
7. **Drift term.** d_i is deterministic given the past, so E‖d_i‖²_{F_i} = d_iᵀE[F_i]d_i, and E F_i = Σ_j(1−η)^{i−j}E ΔF_j + (1−η)^iF_0 ⪯ (c_I/(ηc⁻))Ī + (1−η)^iF_0. Insert steps 6 and 7 into step 5 and divide by n. The corollary multiplies by η/c_F using (A4). ∎

Steps 1, 2, 5, 6, 7 are elementary and exact; step 3 is Proposition 1(a); step 4 uses the information identity and the assumption (A1).

## Numerical check of every step (K = 1, α unknown, step-hold sinusoidal drift in μ, A = 0.3, P = 1200; `code/p13_check.py`)

| setting | η | clip fraction | E⟨g,v⟩ − E vᵀΔFv | mean U | steps 1–5 bound | E q_i vs d(−ln(1−η)) | x̄ (99.9 %, max) | B_h (measured) | c_F median / 1 % / min |
|---|---|---|---|---|---|---|---|---|---|
| ε = 0.4 | 0.004 | 0.0000 | +0.00017 ± 0.00030 | 1.56 | 17.0 | 0.00851 vs 0.00802 | 0.20, 0.61 | 0.014 | 0.80 / 0.23 / 0.145 |
| ε = 0.2 | 0.001 | 0.0001 | (not printed) | 2.37 | 49.6 | 0.00210 vs 0.00200 | 0.094, 0.32 | | 0.85 / 0.29 / 0.22 |
| ε = 0.2 | 0.004 (η/ε² = 0.1, outside the window) | 0.0065 | | 1.03 | 9.47 | 0.00860 vs 0.00802 | 0.36, 1.21 | | 0.65 / 0.10 / 0.067 |
| ε = 0.4 | 0.016 (η/ε² = 0.1, outside) | 0.0129 | | 0.94 | 9.46 | 0.0364 vs 0.0323 | 0.82, 2.75 | | 0.66 / 0.10 / 0.048 |

- **Identity (step 3):** holds within one standard error (first row; 67 494 gaps).
- **Chain (steps 1–5):** the measured mean of U is below the bound in all four settings, by factors 9–21 (the bound is conservative mainly through the Young constants and the dropped negative term).
- **Score norm (step 4):** E|s|²_{F_i⁻¹} = 0.00804 ≤ E q_i/(1−η) = 0.00854 (first row).
- **Discounted trace lemma (step 6):** measured mean tr(F_{i−1}⁻¹ΔF_i) = 0.00823, 0.00205, 0.00834, 0.03396 against the lemma's right-hand side 0.01024, 0.00230, 0.01217, 0.06609 (it holds, with room 1.1×–1.9×; with x̄ at its maximum, which is what the lemma needs, the room is smallest where x̄ is smallest).
- **Drift envelope (step 7):** the ratio η·E‖d‖²_F/‖d‖²_Ī, which stands in for c_I/c⁻, is 1.13–1.42.
- **Fisher-metric MSE** E vᵀĪv = 0.0071 (ε = 0.4, η = 0.004), 0.0025, 0.0086, 0.0311 in the four rows; the corollary's bound with c_F = minimum is 0.47, 0.22, 0.57, 3.2 (loose by 65–100×) and with the median c_F 0.085, 0.059, 0.059, 0.23 (loose by 12–24×).

**What the measured assumptions say about the window.** Inside Problem 14's window (η/ε² ≲ 0.03) the gap leverage satisfies x̄ ≲ 0.2 (max 0.61) and the information lower bound is c_F ≥ 0.15 over the whole run; outside it (η/ε² = 0.1 in both of the last two rows), the 99.9 % quantile of x̄ rises to 0.36 and 0.82 with a maximum above 1 (2.75), the minimum of c_F falls to 0.07 and 0.05, and the clip fraction rises to 0.7–1.3 %. So (A2) and (A4) are exactly the conditions that fail where the empirical window rule says tracking becomes noise dominated, which is the link to Problem 15.

## What is established and what is not

- **Established:** the F-norm bound of Theorem TR under (A0)–(A3), and the Ī-norm corollary under (A4); the chain is verified numerically at K = 1.
- **The assumptions are not derived:** (A1)–(A4) are measured quantities (a.s. versions are stronger than what one expects; a high-probability version needs the stopping argument of Problem 5's Theorem HP together with a concentration bound for the discounted information F_i). A lower bound on ηF_i in terms of Ī is the genuinely near-critical step: it needs the discount window 1/η (events) to exceed the mixing time μ/(βε²) events of Problem 3, which is the same condition as η ≲ ε². I did not prove it; it is measured in the last column above.
- **Conservatism:** constants are one to two orders of magnitude loose (Young step, c_F taken as the run minimum, B_h absorption). The *form* and the dependence on ε (σ² ∝ (r/eff) ∝ ε^{−2} in μ-coordinates) are right; the constants should not be used for tuning (use Problem 15's formulas, which are fitted to the measured MSE).
- **Jumps:** the Young step makes a jump of size D cost D²/η² per jump in the average, whereas the filter analysis and the simulations of Problem 14 give D²/η. Obtaining the sharp rate needs a non-quadratic argument (for example a pathwise bound for ‖v_{i}‖_F along the jump) which I have not done.
- **Drift model:** step-hold per gap, not continuous drift within a gap. The α-unknown case in μ-coordinates carries an extra factor m/eff = k (for k > 1) in the Ī-to-μ conversion for pure-μ drift; the simulations show the actual excess is far below that factor.
- Only K = 1 was simulated; the proof is dimension-free (d = K + 1).
- Re-derive steps 2–6 yourself; the discounted trace lemma (step 6) is the new ingredient relative to Problem 5.
