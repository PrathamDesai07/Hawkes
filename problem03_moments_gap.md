# Problem 3: uniform moments and mixing rate as n → 1 (first pass, K=1 plus the K-dimensional gap)

Notation: ε = 1 − n, k = 2βμ/α² (CIR/Gamma shape), λ = μ + Σ α_k S_k. Code: `code/p3_moments_gap.py`.

## Moments (K=1, stationary)

Apply the generator L f = −β(λ−μ) f' + λ(f(λ+α) − f(λ)) to f = λ^p and set E[Lf] = 0:

p β ε M_p = p β μ M_{p−1} + Σ_{j=2..p} C(p,j) α^j M_{p−j+1},  M_p = E λ^p.

- M_1 = μ/ε, M_2 = μ(2βμ + α²)/(2βε²). Var λ = μα²/(2βε²), which is exactly the Gamma(k) value.
- Scaling: if M_q ~ c_q (μ/ε)^q for q < p, the recursion gives M_p ~ c_p (μ/ε)^p, with each c_p an explicit function of (k, p) only. So E λ^p ≲ C_p (μ/ε)^p uniformly in ε.
- Numerics (μ=β=1): M_p (ε/μ)^p → 1, 3/2, 3, 15/2, 45/2, 315/4 for p=1..6 (k=2), which are the Gamma(2) moments (k)_p/k^p. The recursion agrees with simulation to about 1–3 %.

So the moment half of the lemma is explicit: the p-th moment is (μ/ε)^p times a constant that does not blow up with ε. It depends on k, and it grows quickly with p when k is small.

## Mixing rate

- K=1: L λ^j = −jβε λ^j + (lower order). L is triangular on polynomials, so −jβε, j = 0,1,2,…, are eigenvalues, with eigenfunctions polynomial in λ. The slowest nonzero one is λ − μ/ε at rate βε. Spectral gap ≤ βε, and equal to it on the polynomial eigenspace.
- General K: E[S] solves dE[S] = A E[S] dt with A = −diag(β) + 1 αᵀ. The slowest rate is γ = −max Re eig(A), the root of 1 = Σ α_k/(β_k − γ) in (0, β_1), and γ ≈ ε / Σ_k(α_k/β_k²) as ε → 0. Numerically γ matches this to within 0.3–4 % for ε ≤ 0.05 (K=1,2,3).
- So the mixing time is of order 1/γ ~ 1/(βε) in time units.

## Burn-in in events

Event rate is μ/ε, so a burn-in of 1/γ in time is about μ/(βε²) events. This is 1/ε², not the 1/(με) guessed in Problem 6. The Problem 6 guess is probably too optimistic if burn-in is measured in events; that should be revisited there.

## Status

- Done (derivation plus numerical check): the K=1 moment recursion and its ε-scaling; the mean-relaxation rate γ for general K.
- Not done: a proof that the L² spectral gap equals γ (needs polynomial density in L²(π)); an ε-explicit bound on c_p for all p; general-K moments (the recursion couples S_k, S_l and needs the multivariate analogue); a Wasserstein/Lyapunov contraction statement, which is probably what Problem 5 actually needs.
- The mean-relaxation rate is an upper bound on the gap. The "gap = γ" claim is an expectation, not a theorem.

---

## Second pass: general K (code: `code/p3_generalK.py`)

Setup: A = −diag(β) + 1αᵀ with α_k ≥ 0, −γ its rightmost eigenvalue (γ ∈ (0, β_1), root of 1 = Σα_k/(β_k − γ)), ℓ the left eigenvector, ℓ_k ∝ α_k/(β_k − γ) > 0, normalised so c = ℓᵀ1 = 1. Z = ℓᵀS.

**Lemma A (scalar Lyapunov variable).** L Z = −γ Z + μ exactly (using ℓᵀ(−B) + αᵀ = −γℓᵀ). So E Z_t relaxes at exactly rate γ and E_π Z = μ/γ.

**Lemma B (ε-uniform moments, K general).** With λ ≤ μ + κZ, κ = max_k α_k/ℓ_k, applying L to Z^p gives
pγ M_p ≤ pμ M_{p−1} + Σ_{j=2..p} C(p,j)(μ M_{p−j} + κ M_{p−j+1}), M_p = E_π Z^p.
Solving upward gives M_p ≲ (μ/γ)^p · (constant depending on p and κ/μ). For K=1 it is an equality and reduces to the Section 1 recursion. Because γ ≍ ε/Σ(α_k/β_k²), E λ^p ≲ C_p (μ/ε)^p holds uniformly in ε for general K.
Numerics: for K=1 the bound is tight (MC within 1–5 % for p ≤ 2; p=3,4 at ε=0.05 are noisy, up to 18 % above the bound, but the recursion is exact there so this is MC error). For K=2 the bound holds with margin, 1.3× at p=2 up to about 5× at p=4. It is valid but conservative.

**Lemma C (exact contraction in a weighted L¹ metric).** Couple two chains with the same Poisson random measure on R_+ × R_+ (an event fires when the mark is ≤ λ). Decay preserves componentwise order, and so do jumps, since α ≥ 0 gives λ¹ ≥ λ². For an ordered pair (s¹ ≥ s²), D = ℓᵀ(S¹ − S²) ≥ 0 satisfies L D = −γ D, so E D_t = e^{−γt} D_0 exactly. For a general pair, let s∨s' be the componentwise maximum; the monotone coupling of three chains gives |S¹_t − S²_t| ≤ (S^∨_t − S¹_t) + (S^∨_t − S²_t), hence
**E ℓᵀ|S¹_t − S²_t| ≤ e^{−γt} ℓᵀ|s¹ − s²|.**
This is a clean W¹-type contraction at rate γ, with no ε-dependent constant, and the rate is sharp (it equals the mean-relaxation rate). It implies existence and uniqueness of the stationary law π and
W_ℓ(δ_s P_t, π) ≤ e^{−γt}(ℓᵀs + μ/γ).
So burn-in is about (1/γ) log(initial distance / tolerance), with 1/γ ≍ Σ(α_k/β_k²)/ε in time units.
Numerics (β = (1,4), ε = 0.2): for an ordered start, E d_ℓ(t) matches d_0 e^{−γt} within MC error at all t ≤ 6; for an unordered start it lies below the bound.

**What this settles.** The Lyapunov/contraction statement Problem 5 needs is now explicit and has no ε-dependent constants beyond γ and μ/γ. The "L² gap = γ" question is no longer needed: use the W_ℓ contraction instead.

**Still open.**
- Lemma B's constants as an explicit closed form in (p, κ/μ), not just a recursion.
- Contraction in a metric suitable for martingale concentration (Freedman needs bounds on predictable quantities along the path; Lemma C gives path-wise decay of the *expected* distance, which should be enough for a Poisson-equation / resolvent argument, but I haven't written that).
- The assumption α_k ≥ 0. Without it the coupling's monotonicity fails (this affects Θ_ε only if the constraint set is relaxed).
- Lemma C is a derivation plus numerical check; re-derive the coupling step yourself.
