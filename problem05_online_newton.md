# Problem 5: stochastic analysis of the Fisher-preconditioned online Newton (recursive MLE)

Status: first pass. One short rigorous reduction (Proposition 1), numerical evidence for the rates for K=1, and an explicit list of what is not proven. Code: `code/p5_online_newton_k1.py`, `p5_sweep.py`, `p5_prior_test.py`, `p5_logT.py`, `p5_identity_check.py`.

## Algorithm (event-indexed)

Gap i has length w_i, state s just after event i−1, a_i = (1, S(t_i−)), b_i = ∫_gap (1, S) dt, φ(u) = (1, S(t_{i−1}+u)) on the gap, so a_i = φ(w_i) and b_i = ∫_0^{w_i} φ. Loss l_i(θ) = −log(a_iᵀθ) + b_iᵀθ, gradient g_i = b_i − a_i/(a_iᵀθ_i).
Fisher increment ΔF_i = ∫_0^{w_i} φ φᵀ / λ_{θ_i} du (closed form for K=1). F_i = F_{i−1} + ΔF_i, and θ_{i+1} = Π(θ_i − F_i⁻¹ g_i), with F_0 a prior weight.
Pseudo-regret R_n = Σ_i ∫_gap λ*ψ(λ_{θ_i}/λ*) du, ψ(r) = r − 1 − log r. By the identity below, E R_n equals the expected log-loss regret against θ*.

## Proposition 1 (rigorous; unprojected iteration)

Let θ_i be F_{t_{i−1}}-measurable, v_i = θ_i − θ*, and P_0 = ‖θ_1 − θ*‖²_{F_0}. Then

(a) E[⟨g_i, v_i⟩ | F_{t_{i−1}}] = E[v_iᵀ ΔF_i v_i | F_{t_{i−1}}];
(b) E[R-increment_i | past] ≤ E[⟨g_i, v_i⟩ | past] (since ψ(r) ≤ (r−1)²/r);
(c) with the unprojected update, E R_n ≤ P_0 + E Σ_i g_iᵀ F_i⁻¹ g_i.

Proof. (a) Given the past, the gap has no event before τ = w_i, so φ, λ*(u) = φᵀθ*, λ_θ(u) = φᵀθ_i are deterministic functions of u, with P(τ > u) = e^{−Λ*(u)} and density λ*(u)e^{−Λ*(u)}. Since φᵀv = λ_θ − λ*, ⟨g,v⟩ = ∫_0^τ (λ_θ − λ*) du − (λ_θ − λ*)(τ)/λ_θ(τ). By Fubini its expectation is ∫ e^{−Λ*}(λ_θ − λ*)(1 − λ*/λ_θ) du = ∫ e^{−Λ*}(λ_θ − λ*)²/λ_θ du, which is exactly E ∫_0^τ (vᵀφ)²/λ_θ du = E vᵀΔF v.
(b) λ*ψ(λ_θ/λ*) ≤ (λ_θ − λ*)²/λ_θ, by log r ≥ 1 − 1/r.
(c) Without projection, 2⟨g_i, v_i⟩ = ‖v_i‖²_{F_i} − ‖v_{i+1}‖²_{F_i} + g_iᵀF_i⁻¹g_i. Summing, with ΔF_i = F_i − F_{i−1}: 2Σ⟨g,v⟩ = P_0 + Σ v_iᵀΔF_i v_i + Σ g_iᵀF_i⁻¹g_i − ‖v_{n+1}‖²_{F_n}. Take expectations and use (a): X := E Σ⟨g,v⟩ = E Σ vᵀΔF v satisfies 2X = P_0 + X + G − E‖v_{n+1}‖²_{F_n}, where G = E Σ gᵀF⁻¹g. So X ≤ P_0 + G, and E R_n ≤ X by (b). ∎

Why it matters: exp-concavity is never used. The Fisher increment on the left cancels one copy of itself in the telescoping sum, so only the gradient-norm term G remains to be controlled. This replaces the "adversarial ONS with exp-concavity" step that Problem 1 closed.

## Numerical check of Proposition 1 (K=1, n=2000, unprojected runs that never left Θ_ε)

- Pathwise Newton identity: residual ≤ 7×10⁻¹³ in all settings.
- Expectation identity (a): E Σ⟨g,v⟩ vs E Σ vᵀΔF v = 11.46 ± 0.25 vs 11.40 ± 0.47 (ε = 0.3, start (2, 0.5)); 10.45 ± 0.29 vs 10.78 ± 0.46 (ε = 0.1); 12.49 ± 1.02 vs 12.99 ± 1.04 and 15.86 ± 2.37 vs 16.34 ± 2.26 for the near-truth starts.
- (c): E R_breg = 5.98, 5.58, 7.46, 11.09 against bound P_0 + G = 13.72, 12.54, 14.78, 18.69. R ≤ P_0 + G in 97–98 % of paths. The bound is about 2× loose, since R_breg/⟨g,v⟩ ≈ 0.52–0.7.
- **Caveat: selection.** The unprojected iteration exits Θ_ε in 55–85 % of runs within 2000 events (only 121, 165, 44, 67 of 300 runs stayed inside). These statistics are conditional on never clipping, which can bias martingale identities. They agree within error, but this does not test the projected algorithm.

## Empirical rates (K=1, d=2, projected iteration, μ* = β = 1)

1. **Slope.** E R_n grows like (d/2) log n, i.e. slope 1 (slope between 1.0 and 1.1 for ε ≥ 0.1; noisy at smaller ε).
2. **ε in time units.** E[R_n − log T_n] is nearly flat in ε: −1.2, −2.3, −2.4, −2.8 for ε = 0.3, 0.1, 0.03, 0.01 (starting at θ*, ρ₀ = 1, 16 reps, standard errors 0.3–1.1). In events, R_n − log n drifts as −2.4, −4.6, −6.2, −7.3, i.e. roughly −(d/2) log(1/ε) because T ≈ nε/μ. So no positive C log(1/ε) term is visible for K=1; if one exists in time units its coefficient is ≲ 0.5 and of the favourable sign. The problem-list target "(d/2) log T + C log(1/ε) + burn-in" is therefore not contradicted but not confirmed either; the evidence is for (d/2) log T + O(1).
3. **Normalised error.** E[(θ_n − θ*)ᵀ J_n (θ_n − θ*)] ≈ 1.3–2.7 vs d = 2, as for an efficient estimator.
4. **Burn-in cost is the prior.** Starting far away with a strong prior gives an apparent excess growing like 0.5/ε (6, 14, 24, 55 for ε = 0.05, 0.03, 0.02, 0.01). This is the prior term P_0 = ‖θ_1 − θ*‖²_{F_0}: it scales linearly with the prior weight ρ₀ (≈ +4.9ρ₀ at ε = 0.02, α₀ = 0.5) and disappears for a weak prior (ρ₀ = 0.2: c = −3.2, same as starting at the truth). μ₀ errors of 100 % cost almost nothing. So burn-in cost is a design parameter of F_0, not an intrinsic 1/ε.
5. **Transient length.** With a strong prior, e_n is large for n of order 1/ε² events (ε = 0.03: ≈ 1100), consistent with the Problem 3 mixing time.

## What is not proven (in decreasing order of importance)

1. **Control of G = E Σ gᵀF_i⁻¹g.** Writing g_i = s_i + a_i a_iᵀ v_i/(λ*_i λ_θ,i) with s_i the score, E[s_i s_iᵀ | past] = E ΔJ_i (information identity), so the score part gives Σ tr(F_i⁻¹ΔF_i) ≤ log det(F_n/F_0) ≈ d log T. The drift part needs a leverage bound a_iᵀF_i⁻¹a_i/λ_i = O(1) uniformly (this is where ε and the Problem 2 conditioning enter) and local control λ_θ/λ* ∈ [1 − δ, 1 + δ]. This is the main missing lemma, close to Problem 7.
2. **Sharp constant.** Proposition 1 gives d log T, not (d/2) log T, because ψ ≤ (r−1)²/r loses a factor 2. Recovering ½ needs ψ(r) ≤ (½ + O(δ))(r−1)²/r on a localisation event, i.e. consistency of the iterates. Numerically the ratio is 0.5–0.7.
3. **Projection.** The diagonal projection from Problem 4 is not non-expansive in the F-norm, so (c) fails as written for the projected algorithm; an exact F-norm projection (a QP) would restore it. Empirically the unprojected iteration leaves Θ_ε often, so this is not a technicality.
4. **K ≥ 2.** All numerics are K = 1 (d = 2). Problem 4 shows the K ≥ 2 Fisher matrix is much worse conditioned, and the leverage bound in item 1 will degrade accordingly.
5. **High probability / adaptivity.** Everything above is in expectation with θ_i adapted to the past; no Freedman-type tail bound was attempted.
6. **Dependencies.** Items 1 and 5 rely on the Problem 3 moment and mixing lemmas (Lemma B is itself a draft) and on Problem 2's smallest-eigenvalue bound.

Re-derive Proposition 1 yourself before relying on it. The numerical agreement supports it, but it is a five-line argument I wrote in one pass.

---

## Update (25 Oct 2026): item 1 above is now addressed conditionally

`problem07_leverage.md` splits the gradient into the score and a drift term and proves, for K = 1 and the unprojected iteration, E R_n ≤ [P_0 + (1+η)c⁺r(x̄)·E log det(F_n/F_0)]/(1−B) under a.s. assumptions on three quantities: iterate/truth intensity ratio c⁺, gap leverage x̄, and point leverage ℓ̄ (with B = (1+1/η)ℓ̄/μ_min < 1). The remaining gap is turning those assumptions into stopping-time / high-probability statements. Empirically the point leverage needs about 4(k−1)/(ε² μ_min) events to fall below the threshold, so the effective burn-in is of order 1/ε² events, not the 1/(με) guessed in the problem list. Items 2–5 of the list above are unchanged.

---

## Update (25 Oct 2026): projection (item 3) is resolved for K = 1

Use a compact box Θ_c known a priori and the exact F-norm projection. Then Proposition 1 holds for the projected iterates (non-expansiveness in the F-norm), the ratio-set assumption (L1) holds by construction, and the circularity noted in `problem05_fisher_concentration.md` disappears. Details and numerics in `problem05_circularity.md`. The diagonal clip used in some earlier experiments here is not adequate for tight boxes (it can inflate regret ~10×), so figures above for small ε with a tight region should be read with that in mind.
