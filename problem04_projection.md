# Problem 4: projection onto Θ_ε

Θ_ε = {μ ≥ μ_min, α_k ≥ 0, Σ_k α_k/β_k ≤ 1 − ε}. Code: `code/p4_projection.py`.

## Algorithm (diagonal metric D = diag(d_μ, d_1..d_K)), O(K log K)

- μ is separable: μ ← max(y_μ, μ_min).
- α: minimise Σ d_k(a_k − y_k)² s.t. a ≥ 0, Σ w_k a_k ≤ r with w_k = 1/β_k and r = 1 − ε.
  - If the clipped point max(y,0) already satisfies the budget, return it.
  - Otherwise KKT gives a_k = max(0, y_k − τ w_k/d_k) with τ > 0 chosen so the budget is tight.
  - g(τ) = Σ w_k max(0, y_k − τ w_k/d_k) is piecewise linear and decreasing with breakpoints t_k = y_k d_k/w_k. Sort the breakpoints, form suffix sums, and take the unique suffix whose root τ_j lies in [t_{j−1}, t_j].
- Cost: one sort, O(K log K). A median-of-breakpoints search would give O(K) expected (not implemented).

## Checks

- 300 random instances (K ≤ 7, random β, ε, μ_min, random diagonal metric) against SLSQP: worst coordinate difference 4e-6, which is SLSQP's tolerance. Feasibility asserted on every instance.
- Time per projection for violating inputs: about 0.3–0.9 ms for K between 100 and 10⁴, 16 ms at K = 10⁵ (numpy; no tuning).

## Effect on the regret constants

- Projection onto a convex set in the same metric as the update is non-expansive in that metric, so it adds no extra constant to the standard ONS-type analysis, provided the update and the projection use the same matrix.
- This is an exact result only for a *diagonal* preconditioner. With a dense Fisher/ONS matrix the projection is a QP with a coupled objective and no closed form; that case is not solved here (a warm-started active-set or dual method would cost about O(K²) per event).
- So the price of using the cheap projection is the condition number of the diagonally scaled Fisher matrix, D^{-1/2} I D^{-1/2}:

| case (cond at ε = 0.3 → 0.1 → 0.05) | raw cond | diag-scaled cond |
|---|---|---|
| K=1 (exact, ε = 0.5 → 0.001) | 15 → … → 1.0×10⁶ | 8.6 → … → 5.88 |
| K=1 MC | 27 → 156 → 521 | 8.0 → 7.2 → 6.7 |
| K=2, r=4 (MC, 80k events) | 79 → 163 → 343 | 27 → 44 → 72 |
| K=2, r=2 | 138 → 280 → 574 | 118 → 215 → 368 |
| K=3, r=2 | 1503 → 2967 → 5261 | 765 → 1423 → 2477 |

- K=1: scaling removes the ε-dependence completely. The scaled cond tends to (1+ρ)/(1−ρ) with ρ² → 1/2 (from det/(I_μμ I_αα) → 1/2 for every k), i.e. 5.83. Observed 5.88 at ε = 0.001.
- K ≥ 2: diagonal scaling helps only by a factor of 1.2 to 3 and the scaled cond still grows as ε falls. It is strongly K- and r-dependent (K=3, r=2: thousands). This looks like near-collinearity of the exponentials (Problem 2), not the μ-versus-α scale split. MC is 80k events, so treat these as indicative.
- **Correction (25 Oct 2026):** an earlier version of this table was computed with `p2_identifiability.py`, whose thinning loop redrew the waiting time after each rejection instead of accumulating it. That biased the simulated state and made the K ≥ 2 numbers wrong. The bug is fixed and the table above uses the corrected simulator. The qualitative conclusion (a diagonal metric is not enough for K ≥ 2) is unchanged; the numbers are not.

## Status

- Done: projection algorithm, implementation, and a correctness/timing check.
- Not done: projection in a dense (full-matrix) ONS metric; an O(K) implementation; a proof-level statement of the regret-constant effect (above is the standard non-expansiveness argument plus a numerical condition-number study, not a theorem); a better preconditioner for K ≥ 2 (this depends on Problem 2).
- Decision for Problem 5: with a diagonal preconditioner the projection is exact and cheap, but the curvature constant carries the scaled condition number. That is a constant for K=1 and unresolved for K ≥ 2. If K ≥ 2 constants matter, the dense metric (or a low-rank-plus-diagonal one) has to be handled.

---

## Update (25 Oct 2026): for K = 1 the diagonal clip is not adequate; use the exact F-norm projection

The Problem 5 experiments (`problem05_circularity.md`) show the diagonal clip violates non-expansiveness in the F-norm in 0.7–4.6 % of steps and inflates the regret by about 10× when the box is tight (E R_n = 16.0 vs 3.5 at ε = 0.1; 78.9 vs 5.4 at ε = 0.03). The exact F-norm projection onto a box in R² is O(1) (the optimum is among the unconstrained point, four edge minimisers and four corners) and had zero violations in 24 configurations. The conclusion above that "a diagonal metric is enough for K = 1" concerned the conditioning of the scaled Fisher matrix, which is unchanged; it does not make the diagonal clip a safe projection.
For K ≥ 2, projecting in the F-norm onto {μ ≥ μ_min, α ≥ 0, Σα_k/β_k ≤ 1−ε} is a general convex QP (active-set or dual method, roughly O(K³) per event naively); not implemented.


---
## Update (25 Oct 2026): exact projection in the dense F-norm for any K (closes the Fisher/ONS-metric part of the criterion)

Code: `code/p4_dense_projection.py`. For the regret analysis the projection must be the metric projection in the Newton matrix F_i (Prop. `prop:reg` needs ‖Π^F(y) − θ\*‖_F ≤ ‖y − θ\*‖_F), and F_i is dense for K ≥ 2. Algorithm for Θ = {lo ≤ x ≤ hi, cᵀx ≤ r} (c = (0, 1/β), r = 1 − ε): if y is feasible return y; else solve the box-constrained QP (inner problem) for τ = 0; if the budget holds, done; otherwise find the multiplier τ ≥ 0 of the budget by bracketing and bisection of the monotone function g(τ) = cᵀx(τ), x(τ) = argmin_box ½(x−y)ᵀF(x−y) + τcᵀx. The inner box-QP is a bounded least-squares problem ‖Lᵀ(x − F⁻¹q)‖², F = LLᵀ, solved exactly by BVLS (Stark–Parker active set).
- **Tests** (400 random instances, K = 1..4, F with condition numbers up to 10⁸, 393 with an active projection): objective excess vs brute-force enumeration of all active sets ≤ 4.9·10⁻¹⁴ (relative), vs SLSQP ≤ 4·10⁻⁷ (SLSQP's tolerance); **non-expansiveness: max over (instance, feasible z) of ‖Π(y) − z‖_F − ‖y − z‖_F = 0.0** over 2 000 pairs.
- **Cost:** O(K³) per inner solve (dense F), about 50 inner solves for the bisection: 19 ms (K = 2), 97 ms (K = 5), 187 ms (K = 10), 588 ms (K = 50) in unoptimised numpy+SciPy. The projection is active in 0.2 % (K = 2, ε = 0.3–0.4) to ~5 % (K = 1, ε = 0.03) of steps, so the amortised cost is small for small K. The O(K log K) bound of the diagonal-metric algorithm above does **not** hold in the dense metric; I do not claim it.
- **Effect on the regret constants:** none beyond the usual (exact metric projection is non-expansive in F_i, Prop. `prop:reg`); the diagonal clip must not be used inside the Newton recursion.

**Verdict on the tracker criterion** ("an algorithm with O(K) or O(K log K) cost, a note on its effect on the regret constants, and a working implementation"): met by the diagonal-metric algorithm (O(K log K), 300-instance check, note) and, for the Fisher/ONS metric that the regret analysis really needs, by the exact dense algorithm above (O(K³), tested, note). **Problem 4: Done**, with the caveat that the dense case is not O(K log K).
