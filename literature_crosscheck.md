# Literature cross-check (25 Oct 2026)

What was done: for the papers that bear most directly on results I rely on, I obtained the primary text, extracted the exact statements, re-derived the relevant result from the paper's own setup (not from my earlier notes), and tested it numerically with estimators independent of the ones used before. Papers were fetched from arXiv as PDFs, text-extracted with pypdf in a throw-away virtual environment (outside the project), and treated as data. Code: `code/lit_check_herrera.py`, `code/lit_check_jr.py`.

## Read (full text unless noted) and not read

| Sheet # | Paper | Status |
|---|---|---|
| 1 | Herrera-Marín, arXiv 2609.38268 (29 Sep 2026, unrefereed) | read in full text; Thm 2.2, Cor 2.3, Thm 4.1, Remark 4.8 checked |
| 2 | Herrera-Marín, arXiv 2609.38267 | only the abstract/introduction, Section 5 and Section 6.2 read; Theorem 3.1 and its proof not checked |
| 3 | Hall & Willett, arXiv 1409.0031 | Theorem 1 statement, loss parametrisation (Appendix C) read; proof not checked |
| 4 | Yang et al., NeurIPS 2017 | not read (its journal-style version, #5, was read instead) |
| 5 | Yang et al., arXiv 1801.08273 | read: Assumptions 3.1–3.2, eq. (8)–(9), Theorem 2, Corollary 1; proof of Theorem 2 only skimmed |
| 9 | Jaisson & Rosenbaum, arXiv 1310.2033 | read: assumptions, Prop. 2.1–2.2, Theorem 2.2 and its heuristic derivation; proof not checked |
| 10 | Hazan, Agarwal, Kale (2007) | the original returned HTTP 403; the same material was read in Hazan's monograph (arXiv 1909.05207, §4.2 and §4.4, Thm 4.5) |
| 6, 7, 8, 11–18, and the four "to check" rows | Solo; Rong–Nair; Deschatre et al.; Horst–Xu; Kresin et al.; Khabou–Talbi; Hardiman et al.; Lallouache–Challet; Ben Alaya et al.; Filimonov–Sornette; change-point paper; Zhang/Wang; online EM; Bacry/Cartea | **not read**. In particular Solo (the closest prior art for constant step size) was not obtained, and the "unverified" overlap rows remain unverified. |

## Claims checked

### Herrera-Marín I (2609.38268)
- **Theorem 2.2** (kernel-agnostic KL bound): KL(P₀ᵀ‖P₁ᵀ) ≤ δ²(1+c)I_T/(ρ₀ρ₁), I_T = Tμ(1−ρ₀), 1−ρ₁ = (1+δ)(1−ρ₀), empty start, same μ and kernel shape.
 *Re-derivation:* on a common filtration KL = E₀∫φ(λ₀,λ₁)dt with φ(x,y) = x log(x/y) − x + y ≤ (x−y)²/y; pathwise λ₁ = μ + (ρ₁/ρ₀)(λ₀−μ), so |λ₁−λ₀| ≤ (m/ρ₀)λ₀ with m = δ(1−ρ₀) and λ₁ ≥ (ρ₁/ρ₀)λ₀, giving φ ≤ m²λ₀/(ρ₀ρ₁); E₀λ ≤ λ̄ = μ/(1−ρ₀) for an empty start. I obtain the paper's bound step by step.
 *Test:* KL estimated by Monte Carlo for my K=1 model (μ=β=1; ε₀ = 0.3, 0.1, 0.03; δ = 0.5, 0.1; I_T = 3 and 30; 200 reps each). The bound held in all 12 configurations; bound/KL between 2.2 and 6.9. The local quadratic approximation ½J_ηη(δη)²T was within 0.6–1.4 of the MC KL.
- **Theorem 4.1** (Fisher information of (log μ, log(1−ρ))): J_T = I_T[[A_T, −(1−aA_T)/ρ], [·, (D_T − 2a + a²A_T)/ρ²]], A_T = E[λ̄/λ], D_T = E[λ/λ̄], a = 1−ρ.
 *Re-derivation:* with α = β(1−e^η) the three entries equal my I_μμ = E[1/λ], I_μα = (1−μE[1/λ])/α, I_αα = (Eλ−2μ+μ²E[1/λ])/α² after the Jacobian (μ, −βa); the efficient information for log μ equals my Schur complement det/I_αα because A_T = E[1/λ]·E[λ]. Agreement is exact algebraically.
 *Test:* the score covariance (an estimator independent of E[aaᵀ/λ]; 1200 stationary paths, T = 200/ε, ε = 0.3, bootstrap errors): I_μμ 0.383 ± 0.016 vs exact 0.3766 (z = +0.39), I_μα 0.921 ± 0.043 vs 0.8906 (z = +0.71), I_αα 3.66 ± 0.15 vs 3.49 (z = +1.18). Mean scores ≈ 0.
- **Remark 4.8:** in the light-tailed regime the limit is CIR with a Gamma stationary law of shape κ; the inverse moment is infinite exactly when κ ≤ 1 (Feller-type), so information about the baseline diverges on one side of the threshold. The author states explicitly that the converse direction (boundedness for κ > 1) "would require uniform integrability of 1/V^T, which we do not prove".
 *Consistency with my Problem 10:* κ = 2μ/β·(1−ε)⁻² = my k = 2βμ/α², the same threshold; my numerics gave A_T → k/(k−1) for k > 1 and A_T ∼ ε^{−(1−k)} for k < 1. These are **numerical** findings, not proofs; they point at the direction the paper leaves open and give a rate it does not state.

### Jaisson & Rosenbaum (1310.2033)
- Setup: kernel a_Tφ, ∫φ = 1, m = ∫sφ, T(1−a_T) → λ_J. Theorem 2.2: C_t^T = (1−a_T)λ^T_{tT} → CIR dX = (λ_J/m)(μ−X)dt + (√λ_J/m)√X dB, and V^T = (1−a_T)N_{tT}/T → ∫X ds.
- *Caveat on reading the PDF:* the extracted text writes the volatility as "√ λ m". I resolved it as √λ_J/m (not √(λ_J/m)) from the paper's own derivation (eq. (4): ψ^T(Tx) ≈ e^{−xλ_J/m}/m times √λ_J) and from dimensional consistency.
- *Mapping to my model* (φ = βe^{−βs}, m = 1/β, a_T = 1−ε): in real time the mean-reversion rate is εβ = γ (my spectral gap) and the stationary law is Gamma with shape 2mμ = 2μ/β (my k as ε → 0).
- *Test* (β = 2, μ = 0.5, chosen so that the two readings differ: 0.5 vs 1.0): measured mean²/var of λ = 0.567 (ε = 0.05) and 0.780 (ε = 0.2), against exact 2βμ/α² = 0.554 and 0.781; the √(λ_J/m) reading (shape 1.0) is rejected. Autocorrelation of λ(t) matches e^{−γt}: 0.612/0.607, 0.372/0.368, 0.139/0.135 at t = 0.5, 1, 2 over γ (ε = 0.2); 0.600/0.607, 0.361/0.368, 0.128/0.135 (ε = 0.05). E[V₁^T] for an empty start (exact formula) → the CIR integral as ε → 0 (0.4691 vs 0.4688 at ε = 0.01, λ_J = 8).

### Hazan (monograph) / Hazan–Agarwal–Kale
- §4.2: for the portfolio loss f_t(x) = −log(r_tᵀx) the Hessian is r_tr_tᵀ/(rᵀx)², "a rank one matrix", far from strongly convex, but ∇²f = ∇f∇fᵀ so f is 1-exp-concave (Def. 4.1, Lemma 4.2: ∇²f ⪰ α∇f∇fᵀ). Theorem 4.5 (ONS): Regret ≤ 2(1/α + GD)·n·log T with γ = ½min{1/(GD), α}, which needs a bound G on the gradients and a diameter D.
- Consequence for the Hawkes loss l = −log(aᵀθ) + bᵀθ: its Hessian is the same rank-one aaᵀ/λ², but ∇l = b − a/λ is not in the span of a (b is the compensator integral), so ∇²l ⪰ α∇l∇lᵀ fails for every α > 0 whenever b ∦ a (my Problem 1, Lemma 1, now matched to the textbook criterion), and G is unbounded since b grows with the gap length. I did not re-derive the textbook proof of Theorem 4.5, only checked the statement and the criterion.

### Yang et al. (arXiv 1801.08273)
- The O(log T) regret comes from **Tikhonov regularisation** (eq. (8): ½ω_iμ_i² + Σζ_ij‖f_ij‖²) which makes the per-step loss strongly convex, plus online gradient descent (Theorem 2: R ≤ C₁(1+log M(t)), C₁ = 2(1+pκ_z²)ζ⁻¹|δ−μ_min⁻¹|², step 1/(ζk+b)). The comparator is the *regularised, discretised, truncated* loss; Corollary 1 adds a term C₂N(t) for discretisation/truncation. The authors state the trade-off: large ζ → stable but "deviates from the negative log-likelihood".
- Assumption 3.1 asks for N(t)−N(t−z) ≤ κ_z for all t (bounded increments), plus μ ∈ [μ_min, μ_max].
- **My inference (not a statement in the paper):** for a Hawkes process near criticality the increment bound κ_z should scale like λ̄z = μz/ε, which would make C₁ scale like ε⁻²/ζ. The paper does not discuss the branching ratio at all.

### Hall & Willett (1409.0031)
- Tracks the *intensity vector* λ_t ∈ [λ_min, λ_max]^p with a known (or learned W) dynamical model Φ_t(λ, W) by mirror descent in θ = log(δλ), where the per-bin Poisson loss −⟨x,θ⟩ + ⟨exp θ, 1⟩ is strongly convex with parameter δλ_min. Theorem 1: tracking regret ≤ C(1 + Σ‖λ_{t+1} − Φ_t(λ_t)‖²)√T with C depending on δ, p, x_max, λ_max, λ_min. This is a different object (rates, √T) from estimating fixed parameters with log T regret.

### Herrera-Marín II (2609.38267), Section 5 only
- For fits with K exponentials to a long-memory kernel the quadratic surrogate predicts the misspecification bias to ~10 %; "for larger K the criterion develops an almost flat direction" in which slow components trade mass with the baseline (from K ≈ 5 in his experiment, which he calls a property of the experiment, not a universal constant, and not proved). This is qualitatively the mechanism behind the fast drop of λ_min(I) with K in my Problem 2 (a slow exponential is nearly constant over the window, hence nearly collinear with the baseline). I did not test this paper's claims.

## Consequences for the project

1. **C3 (two-point lower bound KL ≲ δ²Tμ(1−n))** is Theorem 2.2 + Corollary 2.3 of Herrera-Marín I, for any normalised kernel shape including finite sums of exponentials (empty start, shared μ and shape). It is prior art as stated. What is not covered there: lower bounds for *online* estimators' regret, bounds in directions other than the branching ratio (e.g. individual α_k), and non-empty-start/stationary records (the paper says this "remains open").
2. **C4:** the divergence half (information about μ diverges for κ ≤ 1) is Remark 4.8 (via Fatou). The bounded half (κ > 1) is explicitly unproved there; I have numerical evidence and a rate, but no proof of my own that I would stand behind.
3. **C1:** the textbook criterion confirms that exp-concavity cannot be used for the unregularised Hawkes loss; Yang et al. circumvent it by regularising (changing the objective), Hall–Willett by tracking intensities. My Fisher-preconditioned analysis (Proposition 1 in `problem05_online_newton.md`) avoids both, but its novelty relative to Solo's fixed-gain recursion is unchecked because I did not obtain Solo.
4. **Nothing I re-derived contradicted a paper.** The only trap was a PDF-extraction ambiguity (J–R's volatility coefficient), resolved above. Several statements are only checked in the K=1 exponential case.

## Not done
- Papers not read (see table), notably Solo (#6), Horst–Xu (#11), Ben Alaya et al. (#16, exact title still unconfirmed) and the unverified "overlap" rows.
- Proof-level checking of Herrera I Theorem 3.x, Herrera II Theorem 3.1, J–R's proof, Yang et al.'s Theorem 2 proof, Hall–Willett's Theorem 1 proof.
- The NeurIPS 2017 version of Yang et al. (#4).
