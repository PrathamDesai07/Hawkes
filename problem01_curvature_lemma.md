# Problem 1: curvature of the per-event loss

Setup: θ = (μ, α_1..α_K) ∈ Θ_ε. On event interval i, a_i = (1, S_1(t_i−), …, S_K(t_i−)) and b_i = (Δ_i, ∫S_1, …, ∫S_K) (the compensator increments).
Per-event loss: l_i(θ) = −log(a_iᵀθ) + b_iᵀθ. Write λ_i = a_iᵀθ.

## Lemma 1 (rank-one Hessian, not exp-concave)

1. ∇l_i = b_i − a_i/λ_i and ∇²l_i = a_i a_iᵀ / λ_i². This is PSD and has rank one, so l_i is convex, but strongly convex in no direction orthogonal to a_i.
2. l_i is η-exp-concave iff ∇²l_i ⪰ η ∇l_i ∇l_iᵀ (equivalently exp(−η l_i) is concave).
   Here the left side has range span(a_i), while the right side has range span(b_i − a_i/λ_i).
   If b_i is not parallel to a_i then b_i − a_i/λ_i ∉ span(a_i), so the inequality fails for every η > 0.
   If b_i = c a_i it holds with η ≤ 1/(cλ_i − 1)².
3. b_i ∥ a_i would need ∫S_k dt = S_k(t_i−)·Δ_i for every k, i.e. S constant on the interval. S decays strictly, so this fails whenever some S_k > 0. Exp-concavity fails on essentially every event.

Counterexample (K=1, checked numerically in `code/p1_check.py`): μ=1, α=0.5, β=1, S=0.7, Δ=0.8 gives det[a b] = −0.1745 ≠ 0, and min-eig(H − η g gᵀ) = −0.0204·η < 0 for every η tried (1e-1, 1e-3, 1e-6).

Consequence: the "adversarial ONS with exp-concavity" route is closed. Curvature has to come from the expected loss.

## Lemma 2 (corrected curvature: the expected loss)

Take data from the stationary Hawkes law with true θ*, λ*_t = a_tᵀθ*. The expected per-unit-time loss is L(θ) = E[−λ* log λ_θ + λ_θ]. Then

- ∇L(θ*) = 0 and ∇²L(θ) = E[ a aᵀ λ* / λ_θ² ].
- At θ = θ*, ∇²L = I(θ*) = E[ a aᵀ / λ* ], the Fisher information per unit time.
- For any θ, λ_θ/λ*_t ≤ c(θ) := max(μ/μ*, max_k α_k/α*_k), which needs α*_k > 0. Hence
  ∇²L(θ) ⪰ c(θ)^{-2} I(θ*) on all of Θ_ε.

So the loss is strongly convex in expectation, with a constant controlled by λ_min(I(θ*)).

## Lemma 3 (smallest eigenvalue of I(θ*))

For any unit v, by Cauchy–Schwarz,
(E|vᵀa|)² = (E[ |vᵀa|/√λ* · √λ* ])² ≤ E[(vᵀa)²/λ*] · E[λ*].
Since E[λ*] = μ*/ε (ε = 1 − n*),

**λ_min(I(θ*)) ≥ (ε/μ*) · inf_{|v|=1} (E|vᵀa|)².**

The infimum is positive exactly when a has no support on a hyperplane (identifiability). Quantifying it in K and the β-grid is Problem 2. The linear dependence on ε is the point: curvature degrades like ε, not faster, from this argument.

## K = 1 check (exact, no simulation)

Using the stationary Laplace transform of λ (see `code/p10_fisher_k1.py`), I_μμ = E[1/λ], I_μα = (1 − μ I_μμ)/α, I_αα = (E λ − 2μ + μ² I_μμ)/α², det = (E[1/λ]·E λ − 1)/α².
With μ = β = 1:

| ε | λ_min(I) | λ_min / ε | I_μμ | I_αα |
|---|---|---|---|---|
| 0.5 | 0.178 | 0.357 | 0.558 | 2.23 |
| 0.1 | 0.0669 | 0.669 | 0.155 | 10.07 |
| 0.02 | 0.0174 | 0.868 | 0.0367 | 50.0 |
| 0.01 | 0.00917 | 0.917 | 0.0190 | 100.0 |

λ_min ≈ ε/μ* as ε → 0, so the bound of Lemma 3 is tight in order in ε for K=1. The two diagonal entries separate: I_μμ ~ ε and I_αα ~ 1/ε, with det bounded. The ill-conditioning (condition number ~ 1/ε²) is therefore genuine, and it matters for the preconditioner in Problem 4/5.

## Status and open points

- Done: Lemma 1 (with counterexample), Lemma 2, Lemma 3, and the K=1 check. This is a derivation plus a numerical check, not a formal proof of the stationary-law claims.
- Not done: the explicit lower bound on inf_v (E|vᵀa|)² for general K (Problem 2); the case α*_k = 0 for some k (c(θ) is infinite, so restrict to α*_k ≥ α_min or handle separately); the assumption of stationarity, which needs the burn-in discussion of Problem 3.
- Re-derive Lemma 1 and Lemma 3 yourself, as the problem list asks. I have not independently checked them beyond the numerics above.
