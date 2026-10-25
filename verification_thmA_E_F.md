# Verification of Theorems A, E, F (and the impossibility corollary), 25 Oct 2026

Scope: the proofs in `problem10_rigorous.md` (Theorem A) and `problem11_dichotomy.md` §4b (Lemma G, Lemma H, Theorems E and F) and the use made of them in Corollary E′. Method: I re-derived every algebraic step by hand from the stated definitions, checked every numerical constant or inequality claimed with my own code (independent of the other session's scripts), and tested the stationary-measure steps by simulation. Code: `code/verify_thmA.py`, `code/verify_lemmaH_EF.py`, `code/verify_EF_mc.py`, `code/verify_stationarity.py`, `code/lit_check_kl_lemma.py`. "Verified" below means: re-derived and numerically consistent. I did **not** formalise the measure-theoretic justifications (Dynkin's formula, true-martingale property, integrability); the notes sketch them and I read them and found nothing wrong, but that part is a read, not a check.

## Verdict

| Item | Verdict | Details |
|---|---|---|
| Derivation of the reduction (Problem 10 §1–2: Fisher entries, Laplace ODE) | **correct** | generator computation re-derived; matches Herrera-Marín Thm 4.1 |
| **Theorem A, Steps 1 and 3** | **correct** | identity and dominated-convergence arguments re-derived; all three cases |
| **Theorem A, Step 2** | **correct statements, two inadequate justifications** | see below; both repairable |
| Theorem A limits | **confirmed numerically** | k₀=2: m/ε = 1.899, 1.984, 1.998 → 2; k₀=0.6: m/ε^{k₀} = 5.25, 5.73, 5.93 → C = 6.052; k₀=1: m/(ε ln 1/ε) = 1.897, 1.926, 1.945 → 2 |
| Lemma G | **correct** | |
| **Lemma H** | **correct** | every step re-derived; pathwise inequalities hold in all 30 (parameter set, z) cases |
| **Theorem E** | **correct given Lemma H** | |
| **Theorem F** | **correct given Lemma H** | uniformity of C₄ in z needs one sentence (below) |
| KL lemma, Theorem B / Corollary E′ | **lemma correct; the stationary-record application has a gap** | initial-state law, below |

## Theorem A

**Re-derivation.** With c(u) = εβu + h(αu), h(x) = e^{−x} − 1 + x: (1) φ′ = −uβμφ/c from the generator A e^{−uλ} = e^{−uλ}[uβ(λ−μ) + λ(e^{−uα} − 1)] and stationarity; (2) s/c(s) = 1/(εβ + α²s/2) + D_ε(s) with D_ε = (α²s²/2 − h)/((εβs + h)(εβ + α²s/2)) (algebraic identity, residual 5×10⁻¹⁶ numerically); (3) βμ∫₀^u ds/(εβ + α²s/2) = k_ε ln(1 + au/ε), a = α²/(2β), giving φ_ε = P_εF_ε, P_ε = (1+au/ε)^{−k_ε}, F_ε ∈ (0,1] decreasing; (4) dominated convergence: k₀>1 with domination (1+βv/4)^{−(1+k₀)/2} (valid since a ≥ β/4 for ε ≤ 0.29); k₀<1 with domination (au)^{−k_ε}F_ε and the exponential decay of F_ε; k₀=1 with the explicit integral (ε/(a(k_ε−1)))[1 − (1+aδ/ε)^{−(k_ε−1)}] ≈ (ε/a)ln(1/ε), the tail on [δ,∞) being O(ε^{k_ε}) = O(ε) and the factor F_ε on [0,δ] lying in [e^{−4βμδ/(3α)}, 1] with δ → 0 afterwards. All as in the note.

**Step 2, claim by claim** (β = 1; grids up to x = 10⁷):
- h bounds, x²/2 − h ≤ x³/6, x²/2 − h ≥ 0.132x² for x ≥ 1: true.
- "x ≤ 1.5: D_ε ≤ 4/(3α)": re-derived (h ≥ x²/4 at x ≤ 1.5; second factor ≥ αx/2) and true numerically (max of D_εα·3/4 = 0.556).
- **"x ≥ 1.5: D_ε bounded by a constant (numerator ≤ x²/2, denominator ≥ h(1.5)αx/2)": the stated denominator bound only gives D_ε ≤ x/(αh(1.5)), which grows linearly.** The statement is true; a valid argument: h is convex with h(0)=0, so h(x) ≥ x·h(1.5)/1.5 for x ≥ 1.5, hence D_ε ≤ 3/(2h(1.5)α) = 2.07/α. Numerically sup D_εα = 1.000 over x ≥ 1.5 (attained as x → ∞).
- **"for ε ≤ 1/4 and x ≥ 1 the denominator is ≤ αx²": the written reasoning does not prove it** (h ≤ x, εβs ≤ x/3, εβ ≤ α/3 give 20/18·αx², which exceeds αx² for x ≳ 4.75). The statement is true: sup of denominator/(αx²) over x ≥ 1, ε ≤ 1/4 is 0.667. A correct proof needs the exact form (x/3 + h(x))(1/3 + x/2)/x² at ε = 1/4 (and monotonicity in ε).
- "D_ε ≥ 0.13/β for s ≥ 1/α": true, with room (min D_εβ = 0.30). F_ε decay bound exp(−0.13μ(u − 1/α)): true on the grid tested (3 values of μ, 3 of ε, 7 of u).

**Impact:** the theorem is correct; the proof text needs the two repairs above. Nothing downstream changes.

## Lemma H, Theorems E and F

**Lemma H, steps re-derived:** (i) A log W = (μ−γZ)/W + λ[log(1+1/W) − 1/W] (W jumps by Σω_k = 1; the drift of Z equals μ − γZ − λ); (ii) ψ does not jump since Σκ_k = c₀ν̃ − n/2 = 0; A(ψ/W²) = −Σκ_kβ_kS_k/W² + 2ψ(Σβ_kω_kS_k − λ)/W³ + λψ(3W+2)/(W³(W+1)²), using 1/W² − 1/(W+1)² = 2/W³ − (3W+2)/(W³(W+1)²) (checked); (iii) Σβ_kω_kS_k − λ = γZ − μ, |ψ| ≤ κ̄Z ≤ κ̄W, λ ≤ B_zW, (3W+2)/(W+1)² ≤ 3/W (checked); (iv) λ/2 + Σκ_kβ_kS_k = μ/2 + c₀Z; (v) −x²/2 ≤ log(1+x) − x ≤ −x²/2 + x³/3. The stated upper bound follows with C₄ = B_z/3 + C₂. For the lower bound the terms γz/W and c₀z/W² that appear are non-negative and can be dropped, which gives the stated form.
**Pathwise test (own generator, exact A f_z, 2×10⁵ random states per case):** 6 parameter sets (K = 1, 2, 3; k_K from 0.61 to 3.34) × 5 values of z (0.1 … 1000): both inequalities hold in all 30 cases, minimum slack +3×10⁻⁷ (always positive).

**Theorem E, re-derived:** with 2κ̄γ ≤ δ₀/4 and 1/W² ≤ 1/(zW), z ≥ 4(μ/2 + C₂)/δ₀: A f_z ≥ (δ₀/2)/W − γ, so E[1/W] ≤ 2γ/δ₀; and λ/W ≥ min(μ/z, b−γ) gives the bound on E[1/λ]. C₂ is non-increasing in z (via B_z), so the z-condition is satisfiable. **Theorem F, re-derived:** taking expectations in the upper bound gives (c₀z + C₄)E[1/W²] ≥ γ + (c₀ − μ − γz − 2κ̄γ)E[1/W]; integrate −a′/a ≥ θ/(c₀z+C₄) over [1, η/γ]; Jensen a(η/γ) ≥ γ/(μ+η); exponent 1 − θ/c₀ = k_K + (η + 2κ̄γ)/c₀; λ ≤ max(μ,B)(1+Z). One sentence missing in the text: for z ≥ 1, B_z ≤ max(μ, B−γ), so C₄ can be taken independent of z on [1, η/γ].

**Stationary-measure steps by simulation** (K = 2, β = (1,2), α_k = (1−ε)β_k/2):
- E_π[A f_z] = 0: −4.6×10⁻⁴ ± 9.2×10⁻⁴ (z = 1) and −9.4×10⁻⁵ ± 4.7×10⁻⁴ (z = 10) over 16 independent replicates (z-scores −0.5, −0.2). An earlier run with a single sample set looked 4–5σ off because the error bar ignored autocorrelation; with independent replicates it is clean.
- Theorem E (μ = 1.5, k_K = 2.38, ε = 0.05 and 0.02, z₀ ≈ 10.3–10.6): E[1/W] = 0.0377, 0.0204 against the bound 2γ/δ₀ = 0.152, 0.063; E[1/λ] = 0.0564, 0.0229 against explicit bounds 1.04, 0.44 (loose by 18–19×); E[1/λ]/ε = 1.13, 1.145, approaching the conjectured limit 2/(2μ − 1/ν₁) = 1.2.
- Theorem F inequality (μ = 0.4, k_K = 0.635): left side ≥ right side at z = 1, 5, 10 (0.52 vs 0.12; 0.087 vs 0.051; 0.046 vs 0.033).

## Corollary E′ and Theorem B (the use made of E)

The Le Cam step is standard and the KL lemma KL ≤ ρ·uᵀI_Tu is proved by the two inequalities x log(x/y) − x + y ≤ (x−y)²/y and 1/λ′ ≤ ρ/λ; I re-derived the ridge value vᵀIv = (1+ε)²m − (ε/μ)(1+2ε) (Proposition 3) and tested the lemma numerically (8/8 configurations, `lit_check_kl_lemma.py`). **The gap:** the lemma compares the two path laws *with the same initial state*, but the experiment the note calls "stationary observation" has an initial state drawn from π_θ versus π_θ′, which differ. Two statements are available:
1. **Rigorous (given Theorem E):** if the statistician is given the pre-sample state S₀ (equivalently the full pre-sample history), no estimator of μ is consistent when k_K > 1. This is the conditional experiment.
2. **Stationary record with unobserved past:** the record is a marginal of (S₀, events), so by data processing KL(record) ≤ KL(π_θ‖π_θ′) + (conditional KL). It suffices to show KL(π_θ‖π_θ′) = O(1) for the ridge perturbation, which the note does not do. Under a Gamma approximation (not a proof) it is small and T-independent: 0.0040 and 0.0336 (ε = 0.1; δ = 0.1, 0.3), 0.0028 and 0.0226 (ε = 0.01). So the gap is probably closable, but it is the same difficulty that Herrera-Marín I flags for stationary records ("would require controlling the information carried by the unobserved past and remains open"). Do not state the impossibility for the unconditional stationary record as proved.

**Observation model matters for the headline.** For the streaming estimator the natural start is empty (S = 0), where the μ-information diverges at least like ln T for every k (Problem 10 §4, a short argument I read and find plausible but did not re-derive). So "μ is not consistently estimable for k > 1" is a statement about the stationary/conditional experiment, whereas in the empty-start experiment the information diverges logarithmically. Herrera-Marín's heavy-tailed contrast (μ is learned) is for empty-start records and for windows long after the start; the comparison with the light-tailed case should be made within one observation model.

## Not verified
- The proof of Theorem C (superseded by E) and of Problems 11's §3–4c claims on the MLE and the borderline.
- The Problem 12 claims (middle tier, three-tier spectrum) and Problem 10's cold-start coefficients.
- Measure-theoretic justifications (Dynkin, true-martingale, existence and moments of π) — read, not formalised.
- Anything about k_K = 1 for K ≥ 2.

---

## Update (25 Oct 2026): both issues addressed in `theoremA_step2_and_stationary_tv.md`

- **Theorem A, Step 2:** the two inadequate justifications are replaced by Lemma 2 there (D_ε ≤ 2.08/α for all s via convexity of h; D_ε ≥ 0.166/α for s ≥ 1/α when ε ≤ 1/4, via the exact factor r ≤ 1/3 and the monotonicity of h(x)/x²). All claimed constants hold numerically. Theorem A is now fully justified.
- **Corollary E′ gap (initial-state law), K = 1, k₀ > 1:** reduced to bounding TV(π_θ, π_θ′) (not KL) and bounded by C|δ| + Cε^{min(1,k₀−1/2)} via the exact characteristic function, a Gamma comparator and a Sobolev-norm bound on total variation. Verified numerically (TV(π_θ,π_θ′) ≈ 0.026 / 0.070 for δ = 0.1 / 0.3, flat in ε from 0.1 to 0.003). **My proof; not independently reviewed.** K ≥ 2 and explicit constants remain open.
