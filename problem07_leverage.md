# Problem 7 / the missing piece of Problem 5: score, drift and leverage

Status: first pass for K = 1. Three short rigorous lemmas, a conditional theorem closing Problem 5's gradient-norm term under explicit leverage assumptions, and numerics. Whether those assumptions hold a.s./w.h.p. is not proven; it is argued heuristically and measured.
Code: `code/p5_leverage.py`, `code/p5_theorem_check.py`, `code/p5_component_check.py` (uses `code/p5_online_newton_k1.py`).

## Setup and decomposition of the gradient

As in `problem05_online_newton.md`: gap i has φ(u) = (1, S) on the gap, τ = w_i the first-event time, a = φ(τ), b = ∫_0^τ φ, λ*(u) = φᵀθ*, λ_θ(u) = φᵀθ_i, v = θ_i − θ*, g_i = b − a/λ_θ(τ).
Split g_i = s_i + h_i with
- s_i = b − a/λ*(τ), the score at θ* (this is the "compensator term" of Problem 7: the two parts b = ∫φ and a/λ*(τ) have the compensator-minus-jump structure of a martingale),
- h_i = a (aᵀv)/(λ*(τ) λ_θ(τ)), the drift caused by the iterate not being at θ*.

## Lemma 7a (the two parts of the score are jointly centred, with covariance the Fisher increment) — rigorous

Given the past (θ_i, F_{i−1} predictable): E[s_i | past] = 0 and E[s_i s_iᵀ | past] = E[ΔJ_i | past], ΔJ_i = ∫_0^τ φφᵀ/λ* du.
Proof. M_t = ∫_0^t (φ/λ*)(dN − λ* du) is a martingale, stopped at the first jump τ: M_τ = φ(τ)/λ*(τ) − ∫_0^τ φ du = −s_i. Optional stopping gives mean zero; the predictable quadratic variation of the stopped martingale is ⟨M⟩_τ = ∫_0^τ φφᵀ/λ* du, so E[M_τ M_τᵀ] = E⟨M⟩_τ (square-integrable since ΔJ has finite mean). ∎
This is the "both parts concentrate jointly with zero mean at θ*" statement in expectation, with b and a/λ* not parallel entering through the covariance ΔJ being non-degenerate. Concentration beyond second moments is not done.

## Lemma 7b (drift bound via point leverage) — rigorous

Let ℓ̄_i ≥ sup_{u ≥ 0} φ(u)ᵀF_{i−1}⁻¹φ(u)/λ*(u) (the supremum over every point the gap could end at, including the quiet limit φ → (1,0)), and λ_θ ≥ μ_min. Then
E[h_iᵀF_{i−1}⁻¹h_i | past] ≤ (ℓ̄_i/μ_min) · E[vᵀΔF_i v | past].
Proof. hᵀF⁻¹h = [φᵀF⁻¹φ/λ*](τ) · (φᵀv)²/(λ*λ_θ²)(τ). With τ having density λ*e^{−Λ*}, the expectation is ∫ e^{−Λ*}[φᵀF⁻¹φ/λ*](φᵀv)²/λ_θ² du ≤ (ℓ̄/μ_min)∫ e^{−Λ*}(φᵀv)²/λ_θ du = (ℓ̄/μ_min) E vᵀΔFv (the last equality is the one in Prop. 1(a) of Problem 5). ∎

## Lemma 7c (score part → log-determinant) — rigorous

Suppose λ_{θ_i} ≤ c⁺λ* on the gap (so ΔJ_i ⪯ c⁺ΔF_i) and x_i := λ_max(F_{i−1}^{−1/2}ΔF_iF_{i−1}^{−1/2}) ≤ x̄ ("gap leverage"). Then
tr(F_{i−1}⁻¹ E[ΔJ_i | past]) ≤ c⁺ r(x̄) · E[log det(F_i/F_{i−1}) | past], r(x) = x/log(1+x).
Proof. If y_j are the eigenvalues of F_{i−1}^{−1/2}ΔF_iF_{i−1}^{−1/2} ∈ [0, x̄], then Σy_j = tr(F_{i−1}⁻¹ΔF_i) and log det(F_i/F_{i−1}) = Σ log(1+y_j) ≥ Σ y_j · log(1+x̄)/x̄ by concavity. ∎

## Theorem 7 (conditional): the gradient-norm term of Problem 5

Assume the unprojected iteration, P_0 = ‖θ_1 − θ*‖²_{F_0}, and that a.s. for all i ≤ n: (L1) μ_min ≤ λ_{θ_i} ≤ c⁺λ*; (L2) x_i ≤ x̄; (L3) the point leverage satisfies sup_u φᵀF_{i−1}⁻¹φ/λ* ≤ ℓ̄. For any η > 0 with B := (1+1/η)ℓ̄/μ_min < 1,

**E R_n ≤ [ P_0 + (1+η) c⁺ r(x̄) · E log det(F_n/F_0) ] / (1 − B).**

Proof. Prop. 1(c) of Problem 5: E R_n ≤ X ≤ P_0 + G with X = E Σ vᵀΔFv. Since F_i ⪰ F_{i−1}, G ≤ G' := E Σ g_iᵀF_{i−1}⁻¹g_i, and g = s + h gives gᵀF⁻¹g ≤ (1+η)sᵀF⁻¹s + (1+1/η)hᵀF⁻¹h. Lemma 7a gives E[sᵀF_{i−1}⁻¹s | past] = tr(F_{i−1}⁻¹E ΔJ_i), then Lemma 7c bounds the sum by (1+η)c⁺r(x̄)·E log det(F_n/F_0); Lemma 7b bounds the h-sum by B·X. So X ≤ P_0 + A D + B X with A = (1+η)c⁺r(x̄), D = E log det(F_n/F_0); solve for X. ∎

For K = 1, E log det(F_n/F_0) ≈ 2 log T + log det I − log det F_0, so the theorem gives a bound of order A · d log T/(1−B) with A ≥ 1.44, i.e. the right growth in T and a constant roughly 3–6× the optimal one, plus the Problem 5 caveat that the factor ½ from ψ is lost (so another 2×).

## What the leverage assumptions mean, and how big they are (K = 1)

- Point leverage. The supremum over u is attained either at the quiet limit φ = (1,0), where ℓ = [F⁻¹]_{μμ}/μ*, or at large s, where ℓ ≈ I_μμ s/(α det I · T) grows linearly in the state s. With F ≈ T·I and the K=1 values of Problem 10 (I_αα/det I = (k−1)μ/ε), the quiet-limit value is ℓ_quiet ≈ (k−1)/(εT). This is the dominant term, and it is controlled by the weakly informed baseline direction (Fisher information about μ is only ∼ ε per unit time).
 Check: at the start, ℓ_max = 0.500 = 1/(F_{0,μμ}μ*) exactly (ε = 0.1, ρ₀ = 20, F_{0,μμ} = ρ₀ε/μ = 2).
- Hence ℓ̄ < μ_min/(2(1+1/η)) requires T ≳ 2(1+1/η)(k−1)/(ε μ_min), i.e. n_0 ≈ (μ/ε)T ≈ 4(k−1)/(ε² μ_min) events for η = 1 (≈ 20/ε² events at μ_min = 0.1, k = 2). Observed (mean runs, sampled ℓ): median first index with 2ℓ/μ_min < 1 for good = 1520 (ε = 0.1; max 4191) and 10086 (ε = 0.03; max 19486), against 2000 and 22000 predicted. So the leverage burn-in is of order 1/ε² events, the same scale as the Problem 3 mixing time. The Problem 6 target "burn-in(με)" should therefore be read as ∼ 1/ε² events (∼ 1/ε in time units).
- Mean leverage decays like 1/i (1.46, 0.32, 0.057, 0.0063, 0.0017 over index windows 1–10, …, 10⁴–2·10⁴ at ε = 0.1, ρ₀ = 1), as the identity E[aᵀI⁻¹a/λ] = d predicts. The max over reps decays much more slowly (0.44 in the window 10³–10⁴), driven by quiet gaps and rare large states.
- Gap leverage x_i: mean ≲ 10⁻⁴ after 10⁴ events, max over reps 0.06–0.8 in the transient; c⁺ ≤ 2 once the iterates are near θ* (up to 5–10 before).

## Numerical checks (K = 1; 20 reps, n = 20000; μ_min = 0.5)

Component steps, iterates off the truth (start (1.2, 0.8) at ε = 0.1, (1.2, 0.93) at ε = 0.03):
1. Lemma 7a: E Σ sᵀF⁻¹s = 6.91 / 9.18 / 2.23 vs E Σ tr(F⁻¹ΔJ) = 6.77 / 8.99 / 2.20 (agree within about 2 %).
2. Lemma 7b: E Σ hᵀF⁻¹h = 0.024 / 0.038 / 0.003 ≤ bound 3.4 / 8.9 / 1.8. True but extremely loose: the drift term is negligible in practice and G is essentially the score part.
3. Lemma 7c: max over runs of Σ tr(F⁻¹ΔF)/log det = 1.008 / 1.019 / 1.008 ≤ r(x_max) = 1.10 / 1.25 / 1.09.
Full bound (start at θ*, so P_0 = 0; assumptions taken as the worst values seen): bound 63.7 (ε = 0.1, ρ₀ = 60), 25.0 (ε = 0.03, ρ₀ = 200), 4.3 (ε = 0.03, ρ₀ = 600) vs measured E R_n = 1.87, 0.67, 0.27. Valid but 15–40× loose. ε = 0.3 (ρ₀ = 5) and ε = 0.1 (ρ₀ = 20) violate B < 1 once the u = ∞ supremum is included. These runs start at θ* with a strong prior, which pins the iterates, so E R_n is small there by construction; this checks the inequalities, not the sharpness.

## What is not done

1. (L1)–(L3) are assumed a.s. The honest statement needs a stopping time ν at which one fails, a bound on P(ν < n), and a worst-case bound on regret after ν. Showing ℓ̄ small requires λ_min(F_{i−1}) ≳ T·λ_min(I) with high probability: a concentration result for the accumulated Fisher information, which needs the Problem 3 mixing lemma and the Problem 2 eigenvalue bound. Not attempted.
2. The unprojected iteration leaves Θ_ε in 55–85 % of runs (Problem 5). Projection onto Θ_ε in the F-norm would be needed; not done.
3. K ≥ 2. The quiet-limit leverage generalises to [F⁻¹]_{μμ}/μ*, which by Problem 2 is much worse conditioned; nothing measured.
4. Sharp constant. A ≥ 1.44 and the ψ-factor 2 are both artefacts of the argument.
5. Lemma 7b's bound is 10²–10³× loose; a sharper treatment of the drift term would remove assumption (L3) from the theorem, since the drift term is empirically negligible.

Re-derive Lemmas 7a–7c and the theorem yourself; I wrote each in one pass and checked them only numerically.
