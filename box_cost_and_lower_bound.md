# (1) What the a-priori box costs, and (3) a minimax lower bound for any algorithm

Code: `code/p5_box_cost.py`, `code/p5_box_prior.py` (K = 1, exact F-projection). Tags: **[P]** proved on paper, **[P\*]** modulo listed inputs, **[N]** numerical, **[O]** open.

## Part A. The box

### A1. Where the box enters the proof (traced line by line)

The box Θ_c = [μ_lo, μ_hi] × [α_lo, α_hi] ∋ θ\* enters only through five quantities, and the leading constant (d/2) ln T enters through none of them:

| use | quantity | where it appears |
|---|---|---|
| ratio set (L1) | c⁺ = max(μ_hi/μ\*, α_hi/α\*), c⁻ = min(μ_lo/μ\*, α_lo/α\*) | information level λ_T = mTλ_Q/(4c⁺L); constant A = (1+η)(1+ε′)c⁺r(x̄) in Prop. E5 until iterates enter R(δ) |
| martingale bounds | κ₁ = max(1, 1/c⁻ − 1), b_M = 2κ₁, b_A = κ₁²c⁺ | Lemma E4: x/θ with θ ∝ min(1/(2κ₁), ρc⁻/(4(e−2))), so the deviation term scales like x/c⁻ |
| drift lemma | μ_lo | leverage condition B_X = (1+1/η)(e−1)ℓ̄/μ_lo ≤ ¼, i.e. ℓ̄ ≲ μ_lo: leverage burn-in T₀ ∝ c⁺/μ_lo |
| start | P₀ = ‖θ₁ − θ\*‖²_{F₀} | additive; ∝ ρ₀·(μ_hi − μ_lo)² |
| stability | α_hi ≤ β(1 − ε₁/2) for a known ε₁ ≤ ε | only needed to keep Θ_c inside the stable set |

After the bootstrap (iterates in R(δ)) the box constants are replaced by 1 ± δ (Stage 2), so the leading term is box-free. Hence:

**Proposition A1 [P\*].** Theorem R holds for any compact box containing θ\*, with the same leading constant. The additive constant C(box) depends on the box polynomially: from T₁(δ, q) ∝ (c⁺/(mλ_Qδ²ρ\*²))·(P₀ + x/θ_M + d ln T) it grows at most like c⁺·max(1, 1/c⁻) (times the dependence of P₀ on the prior weight), and from the leverage burn-in like c⁺/μ_lo. No other dependence. (The exponents are read off from the proof's formulas, not optimised.)

So "a ratio-2 box" was a convenience; what is really assumed is "θ\* ∈ Θ_c for a known compact Θ_c".

### A2. Empirical cost (K = 1, n = 30 000 events, 6 replicates, exact F-projection, box [μ\*/s, sμ\*] × [α\*/s, β(1−ε/2)])

E R_n (and E R_n − ln T_n) for a start at the truth and at the worst corner (sμ\*, α\*/s), prior weight ρ₀ = 1:

| ε | s | truth start | corner start | extra (corner − truth) |
|---|---|---|---|---|
| 0.3 | 2 | 6.60 | 6.82 | 0.2 |
| 0.3 | 5 | 7.06 | 11.66 | 4.6 |
| 0.3 | 20 | 7.06 | 87.1 | 80 |
| 0.3 | 100 | 7.06 | 2455.9 | 2449 |
| 0.1 | 2 | 3.64 | 4.00 | 0.4 |
| 0.1 | 5 | 6.46 | 7.54 | 1.1 |
| 0.1 | 20 | 6.98 | 29.8 | 22.8 |
| 0.1 | 100 | 6.97 | 596.7 | 589.7 |

The extra regret grows like s^{2.0–2.1} (ε = 0.3: ×17 from s = 5 to 20, ×30 from 20 to 100), the signature of the start term P₀ ∝ (sμ\* − μ\*)²ρ₀. It is a prior-weight effect, not an intrinsic cost of the box. With a weak prior (ε = 0.1, corner start; `p5_box_prior.py`):

| s | ρ₀ = 1 | ρ₀ = 0.1 | ρ₀ = 0.01 |
|---|---|---|---|
| 5 | 0.9 | 0.01 | 0.02 |
| 20 | 21.1 | 1.7 | 0.35 |
| 100 | 643.9 | 69.5 | 11.9 |

(extra regret; the truth-start baseline itself rises from 7.9 to 10.6–13.5 when ρ₀ is small, because a weak prior lets the first iterates wander.) So the practical cost of not knowing the scale of θ\* within a factor s is roughly s² times the prior weight ρ₀, and ρ₀ can be traded against the early-wandering cost; I measured E R_n at one horizon only (n = 30 000), so the table does not by itself show that the log T slope is unaffected; the theory says the slope (the leading constant) is box-free.

### A3. What cannot be removed

1. **In the regime ε = c/T the scale of μ\* cannot be learnt.** By Corollary `cor:Eprime` / Theorem `thm:V` the stationary record carries bounded information about μ when k_K > 1, so no data-driven rule can certify μ_lo, μ_hi with ratio O(1) with probability → 1. A guarantee that is uniform in θ\* in that regime needs prior knowledge of the scale, and the box is a genuine requirement, not an artifact of the proof.
2. **For fixed ε** the baseline is identifiable and the box can in principle be learnt: run the same algorithm in a very wide box (cost A1/A2: polynomial in the width), then restart in a box 1 + o(1) around the first-stage estimate; the second stage has the sharp constant. I did not write this two-stage argument; it needs a consistency statement for the first stage that holds with probability 1 − O(n^{-κ}) (available from Theorem C-exp with the wide-box constants).
3. A *lower* bound μ_lo > 0 and α_lo > 0 is always needed (λ_θ ≥ μ_lo is used in the drift lemma and in the boundedness of ψ-arguments); the *upper* bounds can be taken generously (μ_hi arbitrary, α_hi just below β).

**Status: [P\*] for Proposition A1; [N] for the empirical cost; the two-stage removal for fixed ε is [O].**

## Part B. A minimax lower bound for any algorithm (K = 1, fixed ε)

### B1. Regret is a KL divergence
An algorithm is a rule producing, before each gap, a predictive intensity λ̂_i(u) ≥ μ̂_i > 0 on the gap. Its hazard-form law of the next event time is the density λ̂_i(u)e^{−Λ̂_i(u)}, a proper probability law, so the sequence defines a probability law Q^{(n)} on records of n events. For the data law P_θ^{(n)}, the cumulative log-loss regret against θ is
 E_θ Σ_i [ℓ_i(θ̂_i) − ℓ_i(θ)] = E_θ log(dP_θ^{(n)}/dQ^{(n)}) = KL(P_θ^{(n)} ‖ Q^{(n)}),
and this equals E R_n. (The log-ratio of the two densities is log(λ\*/λ̂)(τ_i) − Λ\* + Λ̂.) This applies to every predictor, not only to online Newton.

### B2. Proposition B (lower bound) **[P\*]**
Let Θ₀ ⊂ int Θ_ε be a compact box of positive Lebesgue measure |Θ₀| (K = 1, d = 2; or any K with d = K + 1), fixed independent of n. Suppose there is an estimator θ̂_n (a function of the record of n events) and a constant C₁ with
 (U)  ∫_{Θ₀} E_θ|θ̂_n − θ|² dθ / |Θ₀| ≤ C₁ ln n / n  for all n.
Then for every predictive law Q^{(n)},
 **sup_{θ∈Θ₀} KL(P_θ^{(n)}‖Q^{(n)}) ≥ ln|Θ₀| + (d/2) ln( n/(ln n) ) − (d/2) ln(2πe C₁/d).**
If (U) holds with C₁/n in place of C₁ ln n/n (an efficient estimator), the bound is (d/2) ln n − O(1).
Proof. (i) Sup ≥ average over Θ ~ Unif(Θ₀). (ii) For any Q, E_Θ KL(P_Θ‖Q) = I(Θ; X) + KL(P̄‖Q) ≥ I(Θ; X), P̄ the mixture (the compensation identity). (iii) Data processing: I(Θ; X) ≥ I(Θ; θ̂_n). (iv) I(Θ; θ̂_n) = h(Θ) − h(Θ|θ̂_n), h(Θ) = ln|Θ₀|, and h(Θ|θ̂_n) ≤ h(Θ − θ̂_n) ≤ (d/2) ln(2πe E|Θ − θ̂_n|²/d) (the Gaussian maximises entropy for a given covariance; AM–GM on the determinant). Insert (U). ∎

### B3. When (U) holds
For K = 1 and a compact Θ₀ in the interior of {ε ≥ ε₀ > 0}, Theorem C-exp and Theorem FK/F (information level λ_n ≥ c n w.h.p.) give, for our own estimator θ̂_n = θ_{n+1}, D_n = ‖θ_{n+1} − θ\*‖²_{F_n} ≤ C(d ln n + x) w.p. ≥ 1 − 5e^{−x} and λ_min(F_n) ≥ cn w.p. ≥ 1 − 3e^{−c′n}, hence E|θ̂_n − θ\*|² ≤ C₁ ln n/n + diam²·e^{−c′n}. The constants in those theorems are explicit continuous functions of (μ, α, β) and of the certified quantities m or λ_L, so they are bounded on a compact Θ₀ — I state this uniformity as the hypothesis (U) rather than proving it. It is a statement about the already-proved chain, not a new estimate.

### B4. Consequences
1. Theorem R and Proposition B together give, for K = 1 and fixed ε, **inf_Q sup_{θ∈Θ₀} E R_n = (d/2) ln n (1 + o(1))**: the leading constant d/2 is minimax optimal; the algorithm is not merely having a good upper bound. The matching is at the level of the leading constant (the lower bound loses (d/2) ln ln n from the ln n in (U)).
2. ε-dependence: in time units the first-order term is (d/2) ln T + ½ ln det I(θ) + (prior term), and det I is bounded as ε → 0 for k₀ > 1 (β²det I → β/(2μ−β), Theorem `thm:A`). Hence no log(1/ε) term appears in time units, consistent with the flat E R_n − ln T in the experiments. An ε-explicit lower bound for ε = ε(n) → 0 would need (U) uniformly down to ε ≈ 1/n and is **[O]**: there the information in the ridge direction is bounded and the geometry is not LAN.
3. K ≥ 2: the same statement with d = K + 1, given (U) from `K2_regret_chain.md` (Theorem FK).

### B5. Limits
- Proposition B is conditional on (U) (uniformity of constants, not written out).
- It lower-bounds the *average* over a fixed box, hence the sup over the box; it says nothing about a single θ (for a single θ a trivial estimator that outputs θ has zero regret — this is why the sup/average over a set is the right formulation).
- The leading-constant-only match leaves open whether C(ε) ~ ε^{-2} ln(1/ε) in Theorem R is necessary; the experiments (R_n − ln T flat or slightly negative in ε) suggest it is a proof artefact.

**Status: Proposition B [P\*] (conditional on (U)); no Hawkes-specific regret lower bound was found in the literature, so this appears to be new at this level of generality; the argument itself is the standard Bayes-mixture/mutual-information route (Rissanen, Clarke–Barron style) with an explicit estimator in place of a LAN theorem.**
