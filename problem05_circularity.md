# Problem 5: removing the circularity (K = 1)

The circularity: Theorem 7 / the Fisher-concentration corollary need the iterates in a ratio set R = {θ : c⁻λ* ≤ λ_θ ≤ c⁺λ*} (assumption L1), but showing the iterates stay near θ* needs the information bound that (L1) feeds. Status: **I did not prove a consistency bootstrap.** What I prove instead is that, for K = 1, the circularity can be avoided by design: make Θ a compact box known a priori and project exactly in the F-norm. Then (L1) holds deterministically. The price is that c⁺, c⁻ are the constants of the box, not of a shrinking neighbourhood of θ*.
Code: `code/p5_projected_exact.py`, `p5_projected_interior.py`, `p5_F_vs_J.py` (uses `p5_online_newton_k1.py`).

## Setup

Θ_c = [μ_lo, μ_hi] × [α_lo, α_hi] ∋ θ*, with α_hi ≤ β(1 − ε/2) (stability margin), all bounds known a priori. Constants
c⁺ = max(μ_hi/μ*, α_hi/α*), c⁻ = min(μ_lo/μ*, α_lo/α*).
Since (μ + αs)/(μ* + α*s) lies between μ/μ* and α/α* for every s ≥ 0, every θ ∈ Θ_c satisfies c⁻λ* ≤ λ_θ ≤ c⁺λ* on every path. Algorithm: θ_{i+1} = Π^F_{Θ_c}(θ_i − F_i⁻¹g_i), where Π^F is the projection in the norm ‖x‖²_{F_i} = xᵀF_ix. For a 2-dimensional box this convex QP has its optimum among {the unconstrained point, the four edge minimisers (fix one coordinate at a bound, minimise over the other, clip), the four corners}, so the exact projection costs O(1).

## Theorem B (K = 1)

1. **(L1) holds deterministically.** For all i, c⁻λ* ≤ λ_{θ_i} ≤ c⁺λ*, with c^± as above. In particular F_{i−1} ⪰ F_0 + J_{t_{i−1}}/c⁺ pathwise, where J is the oracle information.
2. **Proposition 1 survives projection.** Since Θ_c is convex, θ* ∈ Θ_c and F_i ≻ 0, ‖θ_{i+1} − θ*‖_{F_i} ≤ ‖θ_i − F_i⁻¹g_i − θ*‖_{F_i}. Hence 2⟨g_i, v_i⟩ ≤ ‖v_i‖²_{F_i} − ‖v_{i+1}‖²_{F_i} + g_iᵀF_i⁻¹g_i (an inequality in the right direction), and the telescoping/expectation argument of `problem05_online_newton.md` gives 2X ≤ P_0 + X + G, i.e. E R_n ≤ X ≤ P_0 + E Σ g_iᵀF_i⁻¹g_i with the same G, now for the projected algorithm.
3. **Leverage w.h.p.** On the event of Theorem F (`problem05_fisher_concentration.md`), for t_{i−1} ≥ T_0 = O(1/(βε)): F_{i−1} ⪰ (m t_{i−1}/(4c⁺L)) Q, so the point leverage ℓ̄ ≲ 4c⁺L/(m t_{i−1}) · sup_s φᵀQ⁻¹φ/λ* = O(c⁺/(ε m t)) and the gap leverage satisfies x̄ ≤ ℓ̄ ln(N/q)/(μ_lo)·(1/c⁻)-type bound w.h.p. Combining with Theorem 7 gives its conclusion for the steps after T_0.
4. **Burn-in.** Steps before T_0 contribute at most P_0 + E Σ g_iᵀF_0⁻¹g_i (using F_i ⪰ F_0 and Prop. 1), a finite but crude term of order (number of events in [0,T_0]) = O(μT_0/ε).

Proof: (1) is the ratio computation above. (2) is non-expansiveness of the metric projection onto a convex set, applied in the norm of F_i; this is the exact property the diagonal clip lacks. (3) and (4) are Theorem F and Theorem 7 with c^± fixed by Θ_c.

**What is only sketched.** Theorem 7 was proved in expectation under a.s. assumptions. The assumptions of step 3 hold on an event, so the final bound needs a localisation: let ν be the first i > i_0 at which the F_{t_{i−1}}-measurable condition J_{t_{i−1}} ⪰ (m t/(4L))Q fails (a stopping time); run the argument for the stopped process and add a term √(P(ν ≤ n) · E R_n²) ≤ √(2q E R_n²). The truncation of the gap length τ_i ≤ ln(N/q)/μ is handled similarly. I have not written out this bookkeeping.

## Numerical evidence (K = 1; μ* = β = 1; 16–24 reps; n = 20000; stationary-warmed start; ρ_0 = 1)

Box [0.5μ*, 2μ*] × [0.5α*, β(1−ε/2)] (so c⁺ = 2, c⁻ = 0.5):

| ε | projection | E R_n | E R_n − log T_n | steps with projection active | steps violating non-expansiveness | P_0 + E G ≥ E R_n |
|---|---|---|---|---|---|---|
| 0.3 | exact | 5.64 | −3.06 | 0.13 % | 0 | yes |
| 0.3 | diag clip | 5.56 | −3.14 | 0.14 % | 0.010 % | yes |
| 0.1 | exact | 3.46 | −4.13 | 1.0 % | 0 | yes |
| 0.1 | diag clip | 15.99 | +8.40 | 1.1 % | 0.73 % | **no** |
| 0.03 | exact | 5.35 | −1.02 | 4.9 % | 0 | yes |
| 0.03 | diag clip | 78.88 | +72.5 | 6.3 % | 4.6 % | **no** |

(These first runs started at a corner of the box. Interior starts, the truth itself, and a wider box [0.2, 5] × [0.1, hi] give the same picture: the exact projection has 0 violations in every one of 24 configurations and R_n − log T_n between −4.4 and −0.5; the diagonal clip has violations in all of them and inflates the regret to 12–55 when the box is tight, with R_n − log T_n up to +48; with the wide box the damage is mostly gone but violations remain, 0.004–1 %.)
Also: F − F_0 − J/c⁺ was PSD at every step of all 60 runs (ε = 0.3, 0.1, 0.03; c⁺ = 2), as step 1 requires.

## Corrections to earlier notes

- **Problem 4 / Problem 5.** I recommended the diagonal-metric clip as the cheap projection and noted it is non-expansive only in the diagonal metric. Here it is clearly worse than just "not covered by the theory": at ε = 0.03 with a tight box it multiplies the regret by 10 or more, because the Fisher matrix is strongly correlated (ρ² → ½) and the diagonal metric misjudges the geometry. For K = 1 the exact F-norm projection is O(1) and should be used. Some of the "burn-in cost" and exit-from-Θ_ε phenomena reported in the Problem 5 note used the diagonal clip with a wide box (μ_min = 0.1, α_min = 0.02); the conclusion that burn-in cost is mostly the prior term P_0 is unaffected, but the figures there for tight regions are not representative of the exact projection.
- The projection is active in 1–5 % of steps even when the iterate is at the truth, because α* = β(1−ε) is only ε/2 (relative) below the stability bound α_hi. Projection is therefore part of the typical behaviour for small ε, not an edge case.

## Not done

1. A true bootstrap (iterates in a ratio set that shrinks around θ*, giving c^± → 1). This would improve constants but not the scaling in T and ε.
2. The localisation bookkeeping above; exponential tails in Theorem F.
3. A priori bounds: Θ_c must contain θ* with known scale (here a factor 2 in μ and α).
4. K ≥ 2: the exact F-norm projection onto the polytope {μ ≥ μ_lo, α_k ≥ α_lo,k, Σα_k/β_k ≤ 1−ε} is a general QP, not O(1); the contraction (Lemma C) and the two-interval argument extend, but nothing was measured.

Re-derive steps 1 and 2 yourself. They are short, but I wrote them in one pass and the numerics, while consistent, only test them through the final inequalities.

---

## Update (25 Oct 2026): item 1 under "Not done" (a true bootstrap) is now addressed for K = 1

`problem05_bootstrap.md` proves, with polynomial tails and conservative constants, that the iterates enter and stay in R(δ) for t ≥ T_1(δ, q) ≈ d log T_1/(βε δ²) (rigorous constants 10²–10⁴ × the observed ones), and derives the sharp constant (d/2) log T on the good event. The expected-regret version of the sharp constant still needs exponential tails.
