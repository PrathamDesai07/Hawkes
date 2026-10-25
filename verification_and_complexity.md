# Verification record and complexity analysis (25 Oct 2026)

Request: verify each claim, lemma and theorem against "proper maths", and say where time or space complexity can be improved.
Scripts: `code/verify_symbolic.py`, `code/verify_inequalities.py`, `code/complexity_bench.py`, `code/tangent_fast.py` (sympy/mpmath/scipy in the scratchpad venv).
Honest scope: results of type (A) and (B) below are verified *as mathematics* (exact symbolic identities; randomised tests of inequalities with zero violations); type (C) are formulas checked against independent simulation; type (D) are probabilistic theorems whose proofs are chains of (A)–(C) with standard tools — they can be re-derived step by step (done in `final_audit.md`) but not machine-checked; type (E) are heuristic or numerical-only and labelled as such everywhere.

## 1. Verification matrix

### (A) Algebraic identities — 34/34 verified symbolically (sympy, exact derivatives of the generator, 12 random rational points each, |diff| < 10⁻²⁵, most exactly 0)
- generator identities: A Z_c = μ − ελ; A Z^p (p = 2, 3, 4); A e^{θZ_c}; A e^{θS} (K = 1); A(S_iS_j) and A(S_iS_jS_l) (the moment systems of the tangent bound); A e^{−uλ} (Laplace ODE);
- Perron objects: A Z_ω = μ − γZ_ω; 1 = n + γν̃; Σκ_k = 0; Lemma H (i)–(iv) and the elementary identity 1/W² − 1/(W+1)² = 2/W³ − (3W+2)/(W³(W+1)²); in particular Lemma H (i) and (ii), the two generator computations on which Theorems E and F rest, are exact;
- Theorem A step 1 (Gamma factorisation of s/c(s)); Thm T: N(y) and E₁(y) as the real/imaginary parts of h(−iy) and c(−iτ);
- Fisher reduction (det I), efficient information (m − ε/μ)/(1 − 2ε + εμm), ridge Rayleigh quotient (1+ε)²m − (ε/μ)(1+2ε), Hessian of the per-event loss;
- the Lemma C1 identity 2⟨s,v⟩ + A = −2aᵀv/λ_θ + 2∫(φᵀv)λ*/λ_θ; the Newton step identity; the tanh jump law and its periodicity condition; the sinusoid response |1 − κ/(κ+iω)|² = ω²/(κ²+ω²).

### (B) Inequalities — 18 families, 0 violations over 10⁴–10⁵ random cases each
Lemma 1′ (pathwise, p ≤ 6, K ≤ 4); Lemma 4 (exponential moment); Lemma G sandwich and γ = ε/ν̃; Theorem T elementary bounds (sup|d| = 1.1429 ≤ 19, |E₁| lower bounds, |E₁|² upper bound); Lemma 4.3 (α D_ε ≤ 0.9974 ≤ 2.08; ≥ 0.226 ≥ 0.166 for s ≥ 1/α); the positive-part trick and positive definiteness of Q in Lemma F2; e^y ≤ 1 + y + (e−2)y² (y ≤ 1) and e^{−t} − 1 + t ≤ (e−2)t²; the MGF inequality of Lemma E5c; ψ(r) ≤ (½ + δ)(r−1)²/r (max ratio 0.508, 0.516, 0.537 for δ = 0.05, 0.1, 0.25 — note the true constant is ½ + δ/6-ish, so ½ + δ is safe); Σy_j ≤ r(x̄)Σlog(1+y_j); the tangent minorant; the Gamma tail used in Theorem HP; the path-length bounds; the entropy step of the lower bound.

### (C) Formulas validated against independent simulation or exact solvers
Second and third moments of S (0.4–0.9 %); Fisher matrix entries and cold-start transient formula (earlier notes, 0.3–2 %); cluster-survival flow (MC 0.5–1 %); characteristic function and Φ-bound (exact function vs majorants, 10 values of ε); Assumption-M Fourier formula (agreement with simulation and with the Gamma limit); dense F-projection (brute-force enumeration 5·10⁻¹⁴, SLSQP 4·10⁻⁷, non-expansive on 2 000 pairs); K = 2 Newton/martingale chain (per-step inequality to 4·10⁻¹³).

### (D) Probabilistic theorems — proofs re-derived step by step; consequences simulated
Proposition `prop:reg`, Lemmas C1, E4, E5a–d, Proposition E5, Theorems Cexp, HP, R; Theorems E, F; Theorems T′, V; Theorem FK; Proposition LB. All steps are reductions to (A)/(B) plus Freedman/Ville, optional stopping for stopped martingales, Poisson-equation, Le Cam, Pinsker, Plancherel and the chain rule for TV. No step is non-standard, but these are proofs written by the author; independent human review remains advisable for Theorem V (superposition and cluster survival) and Theorem R's input (H3).

### (E) Heuristic or numerical only
Middle-tier upper bounds and eigenvector convergence (Problem 12); cold-start coefficients (Problem 10); matching estimator rate for k_K < 1 (Problem 11); amplitude effect of the data-driven preconditioner (Problem 14); K = 2 regret slope 1.3 vs 1.5.

**Result of this round: no new defect found by the symbolic and randomised checks.** (The three defects found earlier — E5c constant, overstated impossibility, loose general-K identity — were all in type (D) reasoning that no identity or inequality test could have caught; this is why (D) is the category to scrutinise.)

## 2. Complexity: where it was wasteful, what was improved, what remains

| component | before | after / remark | measured |
|---|---|---|---|
| K = 1 online Newton, per event | O(1) time and memory (closed-form ΔF, 2×2 solve, 8-candidate projection) | optimal | – |
| general K, per event: ΔF | quadrature, q nodes: O(qK²) | same (no closed form for K ≥ 2); q = 8–16 suffices | – |
| general K, per event: Newton direction | dense solve O(K³) | O(qK²) by q Sherman–Morrison updates of F⁻¹ (stable: ‖F⁻¹F − I‖ ≤ 1.4·10⁻¹³ after 2 000 events); pays only for K ≳ 100 in NumPy (measured: dense faster for d = 3–51, Sherman–Morrison 1.14× faster at d = 101); in compiled code the crossover should be near d ≈ q | see below |
| general K, per event: state, memory | O(K), O(K²) | optimal for a full-matrix preconditioner; no data are stored (streaming) | – |
| projection, diagonal metric | O(K log K), O(K) | optimal, but invalid inside the Newton recursion (not F-non-expansive) | – |
| projection, dense F-norm (needed by the regret proof) | ~55–90 inner box-QP solves (bisection on the budget multiplier), each O(K³) | **safeguarded Newton on the multiplier: 1.9–6.6 inner solves**, same solution (F-norm difference 1.7·10⁻¹³) | **×15–30 faster** (K = 2: 0.8 ms; K = 20: 21 ms); further gain possible by warm-started active sets with rank-one Cholesky updates (O(K²) per change) — not implemented |
| moments of S up to order 3 (tangent bound, Theorem T) | dense K³×K³ solve: O(K⁹) time, O(K⁶) memory | **Kronecker-sum diagonalisation of A = diag β − 1αᵀ: O(K⁴) time, O(K³) memory** (A⊕A⊕A is the operator of the third-moment system, A⊕A of the second; verified equal to 10⁻¹⁴) | K = 6: ×43; K = 40: 21 ms, where the dense method would need a 32.8 GB matrix |
| tangent bound λ₀-search | 75-point grid, 60-digit eigendecompositions | **concave in t = 1/λ₀** (N(t) = 2tM − t²M₃ is operator-concave since M₃ ⪰ 0, λ_min is concave and monotone) ⇒ Brent search, 13–26 evaluations; double precision suffices for K ≤ 3 | ~3 ms per bound; reproduces the reference values (even marginally better, 0.07913 vs 0.07910) |
| Assumption M by Fourier quadrature | O(T_max/panel · n), T_max up to 3·10⁵ for k₀ < 1 (minutes) | **no longer needed**: the tangent information level replaces it (and works for every K) | – |
| double precision at tight grids | – | K ≥ 10 with grid ratio 1.2 returns ~10⁻¹³ or negative: the true λ_min is below double precision (Cauchy-type conditioning); use mpmath there | – |

Space summary: the estimator needs O(K²) memory (F or F⁻¹) and O(K) for the state, independent of the number of events; the analysis tools needed O(K⁶) memory before and need O(K³) now.

Not improved (and why): (i) the O(K³) dense solve for K ≤ 50 is already cheaper than the Python overhead of the rank-update alternative; (ii) a low-rank-plus-diagonal representation of F could reduce memory to O(K) for large K because the spectrum has one dominant direction and K−1 middle ones, but this would need a proof that the truncation does not destroy the regret bound (open); (iii) the 1/ε² burn-in in events is a statistical, not a computational, cost.

## 3. Open points flagged by this round
1. Type (D) proofs remain author-written; the most valuable independent review targets are Theorem V and the per-time failure bound (H3) of Theorem R.
2. The multiprecision/double-precision evaluations (tangent bound, Assumption M, Φ) are not interval arithmetic.
3. The fast routines (`moments_fast`, `project_newton`, `tangent_fast`) are research code; `project_newton` was tested on 120 random instances for K = 2–20 against the bisection version, and the K⁴ moment solver against the dense one for K ≤ 6.
