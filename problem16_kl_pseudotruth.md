# Problem 16: the exact KL pseudo-truth θ† (contribution C6)

Date: 25 Oct 2026. Status tags as in `LABELS.md` / `K1_main_results.md`: **[P]** proved on paper, **[P*]** proved modulo stated inputs, **[N]** numerical, **[O]** open. Nothing here is refereed. Numerical evidence is in section 6.

## 0. What has to be settled

Problem 16 asks for a theorem on existence and uniqueness of θ†, the minimiser of the exact Hawkes KL (not the quadratic surrogate), with the first-order condition written out. Problem 17 then needs two extra facts from this note: how the stability constraint n ≤ 1 − ε enters θ† (active or not), and a computable characterisation of θ† to test convergence against. So the target theorem has four parts: (a) the objective is well defined, convex, smooth; (b) a unique minimiser exists on Θ_ε and on the unconstrained cone; (c) first-order (KKT) conditions; (d) an explicit criterion for when the constraint n ≤ 1 − ε binds.

## 1. Setting (what "misspecified" means here)

Data: a stationary ergodic point process N on ℝ with finite rate Λ* = E[dN]/dt and predictable intensity λ*_t with respect to its own history (no model assumed: it may be a Hawkes process with a power-law kernel, a different dictionary, a Cox process, ...).

Model: the dictionary {β_k}, β_1 < … < β_K distinct and known, and
λ_θ(t) = μ + Σ_k α_k S_k(t−) = θᵀa_t,  S_k(t) = ∫_{−∞}^{t} e^{−β_k(t−s)} dN_s,  a_t = (1, S_1(t−), …, S_K(t−))ᵀ.

**Key structural point.** The regressor a_t is computed from the *observed* path, not from a simulation of the model. So L(θ) below is defined for every θ with μ > 0 and α ≥ 0, including θ with n(θ) = Σα_k/β_k ≥ 1; the constraint n ≤ 1 − ε is about the stability of the model as a generative object and about the online algorithm (box/projection, Problem 5), not about the existence of the objective. This is what lets us prove n(θ†) < 1 *without* imposing it (Corollary 2).

Pseudo-KL rate:
L(θ) = E[ D(λ*_t, λ_θ(t)) ],  D(x, y) = x log(x/y) − x + y ≥ 0 (Bregman divergence of x log x − x).

Equivalent forms: (i) L(θ) = L_0 + L̃(θ), L̃(θ) = E[λ_θ − λ* log λ_θ], L_0 = E[λ* log λ* − λ*] constant; (ii) by the compensator identity E∫ log λ_θ dN = E∫ λ* log λ_θ dt, **L̃(θ) is the a.s. limit of the per-unit-time negative log-likelihood** −(1/T)ℓ_T(θ) (ergodic theorem). Hence θ† is the almost-sure limit of the batch quasi-MLE. This is the reason the exact KL (and not the surrogate E(λ*−λ_θ)²) is the right object for the online estimator of Problem 5, whose per-event loss is −log(a_iᵀθ) + b_iᵀθ.

Parameter sets: Θ_ε = {μ ≥ μ_min, α_k ≥ 0, Σα_k/β_k ≤ 1 − ε} (μ_min > 0, as in Part 1), and the cone Θ_∞ = {μ ≥ μ_min, α_k ≥ 0}.

Assumptions:
- **(A1)** Λ* ∈ (0, ∞); λ*_t > 0 a.s. (true for any Hawkes truth with μ* > 0).
- **(A2)** E[λ* (1 + |a|²)] < ∞ and E[λ* log⁺ λ*] < ∞. (For a subcritical Hawkes truth these hold: cluster sizes are Galton–Watson with Poisson offspring, which has all moments.)

## 2. Results

**Lemma 16.1 (derivatives) [P].** Under (A2), L is finite and convex on Θ_∞, and C² in the relative interior of each face, with
∇L(θ) = E[ a (1 − λ*/λ_θ) ],  ∇²L(θ) = E[ λ* a aᵀ / λ_θ² ].
*Proof.* λ_θ ≥ μ_min on Θ_∞, so |a (1 − λ*/λ_θ)| ≤ |a|(1 + λ*/μ_min) and |λ* a aᵀ/λ_θ²| ≤ λ*|a|²/μ_min²; dominated convergence under (A2). Convexity: x ↦ −log x is convex and λ_θ is linear in θ. ∎

**Lemma 16.2 (strict convexity) [P].** Under (A1), (A2) and distinct β_k, ∇²L(θ) ≻ 0 for every θ ∈ Θ_∞. No assumption on the truth beyond a nonzero rate.
*Proof.* vᵀ∇²L v = E[λ*(vᵀa)²/λ_θ²] = 0 forces vᵀa_t = 0 P-a.s. since λ* > 0. Take a gap [τ_i, τ_{i+1}) of length > δ (exists with positive probability: the Palm mean gap is 1/Λ* > 0 and Λ* < ∞). On it S_k(τ_i+) = s_k ≥ 1 (the jump at τ_i adds 1 to every coordinate, and S_k ≥ 0 before) and S_k(τ_i + u) = s_k e^{−β_k u}, so v_0 + Σ v_k s_k e^{−β_k u} = 0 on [0, δ]. The left side is real-analytic, so it vanishes for all u ≥ 0. The functions 1, e^{−β_1 u}, …, e^{−β_K u} are linearly independent for distinct β_k > 0, so v_0 = 0 and v_k s_k = 0, hence v = 0. ∎
(Quantitative form: λ_min(∇²L) is continuous and positive on the compact Θ_ε ∩ {μ ≤ M}, so a uniform strong-convexity constant exists; it is *not* explicit and it degrades with K and with grid compression as in Problem 2.)

**Lemma 16.3 (coercivity on the cone) [P].** Under (A1)–(A2), L(θ) → ∞ as |θ| → ∞ in Θ_∞.
*Proof.* Write θ = rω, |ω| = 1, ω ≥ 0. Since E a_k = Λ*/β_k, E λ_θ = r (ω_0 + Σ ω_k Λ*/β_k) ≥ c_0 r with c_0 = min(1, Λ*/β_K) (using |ω|_1 ≥ |ω|_2 = 1). Since log x ≤ x, E[λ* log λ_θ] ≤ Λ* log r + E[λ*|a|] . Hence L̃(θ) ≥ c_0 r − Λ* log r − E[λ*|a|] → ∞. ∎

**Theorem 16 (existence, uniqueness, first-order conditions) [P].** Assume (A1)–(A2), distinct β_k.
(a) L has a unique minimiser θ†_∞ on Θ_∞ and a unique minimiser θ†_ε on Θ_ε (for any ε ∈ (0, 1) with Θ_ε ≠ ∅, i.e. μ_min feasible and α = 0 allowed).
(b) (KKT) θ† = θ†_ε is characterised by: there exist ρ_0 ≥ 0, ρ_k ≥ 0, κ ≥ 0 with
 ∂_μ L(θ†) = ρ_0,  ∂_{α_k} L(θ†) = ρ_k − κ/β_k,
 ρ_0 (μ† − μ_min) = 0, ρ_k α†_k = 0, κ (1 − ε − n†) = 0,
 where ∂_μ L = E[1 − λ*/λ_θ†] and ∂_{α_k} L = E[S_k (1 − λ*/λ_θ†)]. For θ†_∞ put κ = 0.
 In particular, if all α†_k > 0, μ† > μ_min and n† < 1 − ε: **E[a_t (1 − λ*_t/λ_θ†(t))] = 0**, i.e. the *prediction error ratio* is orthogonal to every dictionary feature, in the same sense as the likelihood equations E[a(dN − λ_θ dt)]/dt = 0.
(c) θ†_∞ is the global minimiser of the convex L; if n(θ†_∞) ≤ 1 − ε then θ†_ε = θ†_∞. Otherwise n(θ†_ε) = 1 − ε exactly and κ > 0.
*Proof.* (a) L is continuous (Lemma 16.1), strictly convex (16.2); Θ_ε is compact in α and closed in μ, and L is coercive in μ (Lemma 16.3), so a minimiser exists on both sets; strict convexity gives uniqueness. (b) The constraints are affine, so KKT holds at the minimiser (Slater/linear constraint qualification, Θ_ε ≠ ∅). (c) If θ†_ε ≠ θ†_∞ and n(θ†_ε) < 1 − ε, move along the segment towards θ†_∞: it stays in Θ_ε for a short time and strictly decreases L (convexity), contradicting minimality. The multiplier is positive because the constraint binds and θ†_ε is not a global minimiser (otherwise ∇L = 0 on the free coordinates and the segment argument forces κ > 0 or a contradiction). ∎

**Proposition 16.4 (rate matching identity) [P].** For every θ ∈ Θ_∞,
E[λ_θ(t)] = μ + n(θ)Λ*  (because E S_k = Λ*/β_k by stationarity of dS_k = −β_k S_k dt + dN).
At the constrained minimiser, multiplying the KKT conditions by θ† and summing (complementary slackness kills the ρ_k α_k terms),
 **μ† + n†Λ* − Λ* = μ_min ρ_0 − κ n†.**
*Proof.* θᵀ∇L(θ) = E[λ_θ − λ*] = E λ_θ − Λ*, and by KKT θ†ᵀ∇L(θ†) = μ†ρ_0 + Σα†_k(ρ_k − κ/β_k) = μ_min ρ_0 − κ n†. ∎

**Corollary 2 (the stability constraint is a statement about the immigrant share) [P].**
(i) With ρ_0 = κ = 0 (interior in μ, constraint slack): **μ† = (1 − n†) Λ*.** In words, the pseudo-true Hawkes model reproduces the true mean rate: μ†/(1 − n†) = Λ*. In particular n(θ†_∞) = 1 − μ†_∞/Λ* < 1 **automatically**, even though L does not know about stability.
(ii) For μ† > μ_min: θ†_ε = θ†_∞ iff μ†_∞ ≥ εΛ* (equivalently n†_∞ ≤ 1 − ε). So the constraint is active iff the unconstrained pseudo-true immigrant share 1 − n†_∞ would fall below ε.
(iii) When active (n† = 1 − ε, μ† > μ_min): **μ† = εΛ* − κ(1 − ε) < εΛ*** with κ > 0 the multiplier from (b). The constrained model *under-predicts the mean rate*: E λ_θ† = Λ* − κ(1 − ε) < Λ*.
(iv) **[known for the surrogate: Rong–Nair 2026, eq. (5.7) c* = Λ(1 − Γ*), checked in their PDF 25 Oct 2026; the KL version (i)–(iii) and its consequence n† < 1, which they cannot prove for the LS contrast, are ours]** The same row (E[λ*] = E[λ_θ] at the optimum) is shared by the quadratic surrogate (its normal equations are E[aaᵀ]θ = E[λ* a], whose first row is E λ_θ = Λ*), so rate matching does not distinguish the two criteria. They differ in the *shape* conditions (weights 1/λ_θ² vs 1; section 3).

**Proposition 16.5 (local structure and Newton decrement) [P*].** For any θ_0 in the relative interior with residual r_t = λ*_t − λ_θ0(t),
∇L(θ_0) = −E[a r/λ_0],  H_0 = ∇²L(θ_0) = E[λ* aaᵀ/λ_0²] = I(θ_0) + O(‖r‖) (Fisher information when r ≡ 0, Lemma `expcurv`),
and θ† − θ_0 ≈ H_0⁻¹ E[a r/λ_0] (Fisher-weighted projection of the residual onto span(a)). Quantitatively: f(θ) = L̃(θ) is a sum of terms λ*·(−log λ_θ) with λ* ≥ λ_min* > 0 when the truth has μ* > 0, so f is self-concordant with constant 2/√μ*; if the Newton decrement δ_0 = (∇Lᵀ H_0⁻¹ ∇L)^{1/2} < √μ* and the Dikin ellipsoid lies in Θ, then ‖θ† − θ_0‖_{H_0} ≤ δ_0/(1 − δ_0/√μ*). Tag [P*] because the closedness hypothesis of the self-concordance theorem is verified only on the face where L → ∞ at the boundary (μ → μ_min barrier supplied by the constraint set); stated here as a diagnostic that is also checked numerically.
Contrast with the quadratic surrogate: D(λ*, λ) = (λ* − λ)²/(2λ) + O(|λ* − λ|³), so the KL criterion is a **1/λ_θ-weighted** least squares, the quadratic surrogate is unweighted. The weighting matters exactly when λ_t fluctuates by a large factor, i.e. near criticality (Gamma shape k small, Problem 10).

**Corollary 3 (link to the online estimator, Problem 17) [P*].** (i) The batch QMLE θ̂_T = argmin Σℓ_i converges a.s. to θ† on Θ_ε and on Θ_∞: pointwise a.s. convergence of the convex functions L̂_T(θ) → L(θ) (ergodic theorem for the stationary marked process, E|ℓ| < ∞ by (A2)), hence uniform on compacts (convexity lemma) and argmin convergence by uniqueness. (ii) If the pseudo-truth is interior, the online Newton estimator has the right limit point by Prop. 1 of `problem05_online_newton.md` (E⟨g_i, v_i⟩ = E vᵀΔF_i v_i is a statement about the *mean* of the gradient, which vanishes exactly at θ†). (iii) When the constraint binds the algorithm's projection onto {n ≤ 1 − ε} is a non-expansive F-norm projection (problem04_projection.md), and the fixed point is the KKT point of Theorem 16(b) (stationary point of the projected dynamics = θ† by convexity). The convergence *rate* and the active-constraint dynamics are Problem 17.

## 3. Why this differs from the quadratic surrogate (the "difference from" statement)

The prior-work criterion is Q(θ) = E(λ* − λ_θ)² (least squares on the intensity); its minimiser solves the linear system E[aaᵀ]θ_Q = E[λ*a] (no positivity, no constraint handling, no 1/λ weight, and not the limit of any likelihood-based online method). The exact KL pseudo-truth θ†:
1. is the limit of the (quasi-)likelihood estimators, which Q is not;
2. weights squared residuals by 1/λ_θ (second-order expansion above), so in principle it differs from θ_Q by an amount governed by how much λ_t varies. **Numerically (V4, section 6) the difference is small:** 0.1–0.3 % relative in θ for a power-law truth with n* = 0.8 and 0.95, K = 1..4. The honest claim is therefore a difference in principle and in what is proved (items 1, 3, 4), not a large numerical gap; I did not search for examples where the gap is large;
3. is automatically in the stability region (Corollary 2(i)) and has the exact rate-matching identity;
4. has a nonlinear first-order condition (not a linear system), so θ† is computed by Newton on the empirical objective, with the same Fisher-type curvature E[λ* aaᵀ/λ_θ²] that the online algorithm uses as preconditioner.

## 4. What was left open, and its resolution (sections 8 and 9)

The first pass left three items open: (1) an explicit strong-convexity constant, (2) the asymptotic variance of the QMLE under misspecification, (3) an a priori proof that μ† > μ_min. Section 8 gave (1) and (2), and an a-posteriori form of (3). Section 9 closes (3) a priori (Prop. 16.9) and adds the limit law with active constraints (Theorem 16.10). Still not covered: any bound on n(θ†) − n* (Problem 18), and the CLT itself without a mixing assumption.

## 5. What to code (pre-registered predictions)

Code in `code/p16/` (C simulator `sim16.c`, numpy estimator `p16lib.py`, drivers `t1.py`..`t7.py`; run with msys2 Python `C:\msys64\ucrt64\bin\python.exe`). Truth simulators: K'-exponential Hawkes by Ogata thinning with state vector; power-law kernel approximated by M exponentials (Laplace quadrature of (1 + t/τ)^{−(1+η)}, rescaled so ∫φ = n*). Estimator of θ†: streaming Newton on the empirical objective L̂_T(θ) = (1/T)[∫λ_θ dt − Σ log λ_θ(t_i)] (path regenerated from a seed each iteration, no storage), with a log-barrier for the constraints; exact gap formulas a_i = (1, s_k e^{−β_k w}), b_i = (w, s_k(1 − e^{−β_k w})/β_k).

| # | Experiment | Prediction | Pass criterion |
|---|---|---|---|
| V1 | Well-specified truth, K = 1, 2: θ̂ → θ* | θ† = θ* (Lemma `expcurv`) | error shrinks like T^{−1/2} |
| V2 | Rate matching in the sample: ∫₀ᵀ λ_θ̂ dt = N_T exactly at an interior optimum | μ̂ = (1 − n̂)N_T/T + O(1/T) | identity to machine precision (code check) |
| V3 | Misspecified: truth = 2-exponential Hawkes, model K = 1 (β on grid); truth = power-law, model K = 1..4 geometric grid | θ†∞ well-defined, μ†/(1 − n†) = Λ*, n† < 1; H ≻ 0 | multiple starts give one minimiser; rate identity holds to MC error |
| V4 | KL pseudo-truth vs quadratic surrogate θ_Q (same data) | θ† ≠ θ_Q; same μ(1−n)-row; θ_Q may violate α ≥ 0 | report ‖θ† − θ_Q‖ and n(θ†), n(θ_Q) |
| V5 | Active constraint: choose 1 − ε below n(θ†_∞) | n† = 1 − ε, κ > 0, μ† = εΛ* − κ(1 − ε), ∇L(θ†) = (0, −κ/β_k) up to the ρ's | KKT identity holds to MC error |
| V6 | Convergence of θ̂_T to θ† with T (misspecified) | error ~ T^{−1/2} (sandwich) | slope ≈ −1/2 on log-log |
| V7 | Newton decrement bound of Prop 16.5 | actual ‖θ† − θ_0‖_H ≤ δ_0/(1 − δ_0/√μ*) | inequality in all trials where δ_0 < √μ* |

## 6. Results of the coding study (all [N]; truncated power law, τ = 1, η = 0.5, M = 24 exponentials, μ* = 1; one seed per cell unless stated)

- **V1/V2 (code check).** Well specified K = 1 (μ, α, β) = (1, 0.7, 1), 4·10⁵ events: (μ̂, α̂) = (0.994, 0.704). At the QMLE, (∫λ_θ̂ dt)/N = 1.0000000059, so μ̂/(1 − n̂) = N/T to 10⁻⁵ relative (Prop. 16.4, exact up to the boundary term).
- **V3 (existence, uniqueness, rate matching), n* = 0.8, Λ* = 5.** 1.5·10⁶ events per K. Different starting points give the same minimiser (difference ≤ 7·10⁻¹⁴). μ†/(1 − n†) = 4.9967 for K = 1..4, equal to the sample rate (true 5.000). The Hessian is positive definite but its spread grows fast with K (K = 4: eigenvalues 0.006 … 2170), as Problem 2 predicts. n(θ†) = 0.463, 0.621, 0.697, 0.733 for K = 1..4 against n* = 0.8; at n* = 0.95, K = 3: n(θ†) = 0.868. So **n(θ†) < n*, increasing in K** (input to Problem 18). K = 3, n* = 0.8 was run three times with different seeds/lengths: 0.697, 0.706, 0.702.
- **V4 (KL vs quadratic surrogate).** ‖θ_KL − θ_Q‖/‖θ_KL‖ = 0.003 (K = 1) and 0.001–0.003 (K = 3, n* = 0.8 and 0.95). The surrogate also satisfies μ/(1 − n) = Λ (4 digits), as in Corollary 2(iv), and stayed in α ≥ 0. The advantage of the exact KL is structural, not numerical, in these examples.
- **V5 (active constraint), K = 3, n(θ†_∞) = 0.706.** For 1 − ε = 0.6, 0.5, 0.3: n† = 1 − ε; μ† vs εΛ − κ(1 − ε): 1.98582 vs 1.98564, 2.47682 vs 2.47652, 3.45641 vs 3.45585 (≤ 2·10⁻⁴ relative); ∂_μ L ≈ 10⁻¹²; −β_k ∂_{α_k} L equals one κ for all k (0.02469, 0.04797, 0.14904). At 1 − ε = 0.3 the slowest α hits 0 and its entry (0.107) falls below κ, consistent with ρ_1 ≥ 0. E λ_θ† = Λ − κ n†. All as in Corollary 2(iii).
- **V6 (convergence to θ†), K = 2, n* = 0.8, 16 seeds per size.** RMS error of μ̂: 0.186, 0.094, 0.048, 0.037 at 10⁴, 4·10⁴, 1.6·10⁵, 6.4·10⁵ events: halves per 4× events until the floor set by the reference error (0.01). α entries show the same behaviour. Mean errors are within noise.
- **V7 (Newton-decrement bound).** Four perturbed starts: actual H-distance 0.0091, 0.0498, 0.1887, 0.7106 vs bound δ/(1 − δ/√μ*) = 0.0092, 0.0538, 0.2105, 1.1649. Never violated; tight for small δ.

## 7. Verdict

Theorem 16 and Props. 16.4/16.5 are proved on paper (Lemma 16.2 via linear independence of exponentials on a gap, Lemma 16.3 coercivity, KKT with rate matching) and each quantitative prediction was reproduced numerically. One expectation failed: the size of the KL-versus-quadratic difference (section 3, item 2), now corrected. Open items are in section 4. Not independently reviewed.

## 8. Second pass: the three open items

**Proposition 16.6 (explicit strong convexity) [P].** Assume additionally λ*_t ≥ λ_lo > 0 a.s. (any Hawkes truth: λ_lo = μ*). Let Y = Σ_k β_k S_k and, for M ≥ μ_min,
G_M = E[ a aᵀ/(M + Y)² ].
Then G_M ≻ 0 (same argument as Lemma 16.2) and for every θ ∈ Θ_ε with μ ≤ M:
**∇²L(θ) ≽ λ_lo · G_M**, hence λ_min(∇²L(θ)) ≥ m := λ_lo λ_min(G_M) > 0, and L(θ) − L(θ†) ≥ (m/2)‖θ − θ†‖².
*Proof.* On Θ_ε, α_k ≤ β_k·(Σα_j/β_j) ≤ β_k, so λ_θ = μ + Σα_k S_k ≤ M + Y; and λ* ≥ λ_lo, so E[λ*(vᵀa)²/λ_θ²] ≥ λ_lo E[(vᵀa)²/(M+Y)²]. ∎
**The relevant M.** Take M = Λ*: by Prop. 16.4, μ† ≤ Λ*(1 − n†) ≤ Λ* whenever ρ_0 = 0 (μ† > μ_min), and μ† = μ_min ≤ Λ* otherwise. So every θ in a Λ*-ball in μ, in particular θ†, satisfies the bound. G_{Λ*} is a functional of the stationary law of S only (a time average over the data, no knowledge of θ† or of the truth's kernel). Consequence for the estimator: ‖θ̂_T − θ†‖ ≤ ‖∇L̂_T(θ†)‖/m_T on the event that L̂_T inherits the constant.

**Proposition 16.7 (consistency and sandwich asymptotics) [P for consistency, P* for the CLT].**
(i) *Consistency.* Let the truth be ergodic and (A1)–(A2) hold. L̂_T(θ) = (1/T)(∫₀ᵀ λ_θ dt − ∫₀ᵀ log λ_θ dN) is convex in θ. For each fixed θ, L̂_T(θ) → L̃(θ) a.s.: (1/T)∫λ_θ dt → μ + Σα_k Λ*/β_k by the ergodic theorem (the initial state s(0) contributes O(1/T)), and (1/T)Σ_i log λ_θ(t_i) → E[λ* log λ_θ] by the ergodic theorem for the stationary marked point process, using E[λ*|log λ_θ|] ≤ E[λ*(|log μ_min| + |a|)] < ∞. Pointwise a.s. convergence of convex functions on a countable dense set implies uniform convergence on compacts (Rockafellar Thm 10.8), and strict convexity gives unique minimisers, so argmin L̂_T → θ† a.s. on Θ_ε and on Θ_∞ (Hjort–Pollard). [P]
(ii) *Sandwich.* Let u_t = (a_t/λ†_t)(dN_t − λ*_t dt) (martingale part) and R_t = a_t(λ*_t/λ†_t − 1) dt. The score at θ† is S_T = M_T + R_T with M_T = ∫u, R_T = ∫R, and E R_t = 0 by the first-order condition (Theorem 16(b), interior case). Then
√T(θ̂_T − θ†) ⇒ N(0, H⁻¹ J H⁻¹),  H = E[λ* aaᵀ/λ†²] = ⟨M⟩-rate,  **J = lim Var(S_T)/T = H + Γ_R + Γ_{MR} + Γ_{MR}ᵀ,**
where Γ_R is the long-run covariance of R and Γ_{MR} the long-run cross covariance. Well specified ⇒ R ≡ 0, J = H and the covariance is H⁻¹ = I⁻¹ (the Cramér–Rao value). Under misspecification J − H can have either sign and the *naive* H⁻¹/T is wrong. [P*]: requires a CLT for the stationary sequence of gap scores; it holds if the truth is a subcritical Hawkes process with exponentially decaying kernel (geometric ergodicity, Problem 3, Lemma `ellcontr`, gives strong mixing) with moments from (A2) and one more moment; stated here under that mixing assumption. The active-constraint case replaces the Gaussian limit by the projection of N(0, ·) onto the cone of feasible directions (Geyer 1994); not worked out.
Estimating J: split the path into B blocks, take block sums of the centred gap scores (a_i/λ_i − b_i) and divide their variance by the block length (batch means); the plug-in covariance is Ĥ⁻¹ Ĵ Ĥ⁻¹/T.

**Proposition 16.8 (μ† > μ_min is certified a posteriori, not a priori) [P].** Let μ̄ = argmin over {μ ≥ μ_min, ...}. If the minimiser has μ̄ > μ_min (equivalently ∂_μ L = 0 there), it is also the minimiser for every smaller μ_min (a local minimum of a convex function on a convex set is global on any larger convex set in which it is still a local minimum). So the choice of μ_min is immaterial whenever the optimiser stays off the edge; this is directly checkable from the estimator (and was true in every run). An *a priori* bound would follow from μ† ≥ sup_x(E[λ* 1{λ_θ† − μ† ≤ x}] − x), obtained from 1 = E[λ*/λ_θ†] ≥ E[λ* 1{Z ≤ x}]/(μ† + x) with Z = Σα†_k S_k; evaluated on the data this bound is *negative* in every example I tried (below), i.e. vacuous: the quiet-path mass of the dictionary state is too small to force μ† away from 0 through this inequality. Closed a priori in Prop. 16.9 below by a different argument.

### Numerical results (second pass)

- **Prop. 16.6 (certified constant).** Power-law truth, grid as in section 6, λ_lo = μ* = 1, G_{Λ*} computed by 10-point Gauss–Legendre on each gap, 1.5·10⁶ events. Certified m vs actual λ_min(H(θ†)): K = 1: 2.0·10⁻³ vs 3.3·10⁻² (ratio 0.059); K = 2: 4.4·10⁻⁴ vs 1.6·10⁻² (0.027); K = 3: 1.4·10⁻⁴ vs 9.2·10⁻³ (0.015); K = 4: 5.7·10⁻⁵ vs 6.0·10⁻³ (0.010); K = 3, n* = 0.95: 4.8·10⁻⁶ vs 1.4·10⁻³ (0.003). The certificate is valid in every case (m ≤ λ_min) and non-vacuous, but conservative by a factor 17 to 300, worsening with K and with n*, which is the same degradation as in Problem 2.
- **Prop. 16.8 (a priori μ bound).** The crude form (using only α_k ≤ β_k) gave −0.61, −2.2, −4.3, −6.5 for K = 1..4 and −28 at n* = 0.95; the sharper a-posteriori form with Z = α†ᵀS gave −0.17 (K=2), −0.56 (K=3), −4.3 (n*=0.95, K=3), against μ† = 1.90, 1.52, 2.81. Both vacuous. So item (3) stays open a priori.
- **Prop. 16.7 (sandwich), K = 2, n* = 0.8 power-law truth, 300 replicates of 4·10⁴ events, B = 20 blocks.** Variance of (μ̂, α̂_1, α̂_2): empirical (1.26·10⁻², 1·10⁻⁴, 7·10⁻⁴), sandwich (1.05·10⁻², 1·10⁻⁴, 7·10⁻⁴), naive H⁻¹/T (0.75·10⁻², 1·10⁻⁴, 6·10⁻⁴). For μ the naive formula underestimates the variance by 40 %; the sandwich by 17 %; the correlation of (μ̂, α̂_1) is −0.812 empirical vs −0.809 sandwich. Control, well-specified K = 2 (α = (0.12, 0.45), n = 0.7): empirical (1.3·10⁻³, ·, 4·10⁻⁴), sandwich (1.4·10⁻³, ·, 4·10⁻⁴), naive (1.5·10⁻³, ·, 4·10⁻⁴), all within MC error, as J = H predicts. The remaining 17 % underestimate for μ is within the batch-means bias at B = 20 blocks of 2·10³ events (long memory in the residual); not tuned.

### Verdict (second pass)

Resolved: explicit and computable strong-convexity constant (certified on the data, tested); consistency written out; sandwich variance derived and reproduced numerically (naive formula shown to be wrong for μ). Closed in section 9: the a priori μ† > μ_min and the active-constraint limit law. Not resolved: a CLT without the stated mixing assumption.

## 9. Third pass: a priori μ† bound and the limit law with active constraints

**Proposition 16.9 (a priori lower bound on μ†) [P].** Let A_ε = {α ≥ 0, Σα_k/β_k ≤ 1 − ε} (use ε = 0 for the cone problem, where n† < 1 by Corollary 2(i)), and
Φ(μ) = min_{α ∈ A_ε} E[ λ*_t / (μ + Σ_k α_k S_k(t−)) ] = min_{α ∈ A_ε} Λ* E⁰[ 1/(μ + αᵀS(0−)) ],
where E⁰ is the Palm expectation (average over the event times). Φ is continuous and strictly decreasing, Φ(μ) → 0 as μ → ∞. Put μ_0 = sup{μ ≥ 0 : Φ(μ) ≥ 1} (μ_0 = 0 if Φ(0) < 1). Then **μ† ≥ μ_0**, whatever μ_min is.
*Proof.* f(μ, α) = E[λ*/(μ + αᵀS)] is convex in α (1/z is convex) and strictly decreasing in μ, A_ε is compact, so Φ is well defined, continuous and strictly decreasing. Suppose μ† < μ_0. Then ∂_μ L(θ†) = 1 − f(μ†, α†) ≤ 1 − Φ(μ†) < 0 (strict by the strict decrease and Φ(μ_0) ≥ 1), and increasing μ keeps every constraint of Θ_ε satisfied, so L would decrease: contradiction. ∎
Remarks. (a) μ_0 depends only on the law of the data (a functional of the stationary law of (λ*, S)), not on θ†. (b) Empirically E[λ*/λ_θ] = (1/T)Σ_i 1/λ_θ(t_i) (compensator identity), so Φ is a convex program over α ≤ K variables on the event-time sample. (c) The crude Cauchy–Schwarz route, E[λ*/Z] ≥ Λ*²/E[λ*Z], is not enough because E[λ*S_k] > Λ*E S_k for a clustered process (positive correlation between λ* and S_k); this is why the first attempt in section 8 failed, and why Φ must be minimised over α rather than bounded in closed form.

**Theorem 16.10 (limit law with active constraints) [P*].** Assume θ† ∈ Θ_ε is a KKT point with the active set 𝒜 (constraints g_jᵀθ = b_j, j ∈ 𝒜: α_k ≥ 0 for k in some set, n ≤ 1 − ε if it binds, μ ≥ μ_min if it binds), **strict complementarity** (all active multipliers > 0), linearly independent normals {g_j}_{j∈𝒜}, and the CLT of Prop. 16.7: Z_T := √T(∇L̂_T(θ†) − ∇L(θ†)) ⇒ Z ∼ N(0, J). Let V be a matrix whose columns span 𝒩 = {d : g_jᵀd = 0, j ∈ 𝒜}. Then
√T(θ̂_T − θ†) ⇒ V (VᵀHV)⁻¹ Vᵀ Z̃, Z̃ ∼ N(0, J)  (sign irrelevant),
i.e. a Gaussian supported on 𝒩 with covariance **V(VᵀHV)⁻¹ VᵀJV (VᵀHV)⁻¹ Vᵀ**. In particular the active constraints hold exactly in the limit (n̂ = 1 − ε in every replicate), the covariance has rank d − |𝒜|, and the unconstrained sandwich H⁻¹JH⁻¹ overstates the variance.
*Proof (epi-convergence).* Put h_T(d) = T[L̂_T(θ† + d/√T) − L̂_T(θ†)] on the feasible set √T(Θ_ε − θ†). A second-order expansion (L̂_T convex, C² with Hessian converging uniformly on compacts) gives h_T(d) = d·√T∇L(θ†) + d·Z_T + ½dᵀHd + o_p(1). By KKT, ∇L(θ†) = Σ_{j∈𝒜} (multiplier_j) g_j with strictly positive multipliers. For feasible d, g_jᵀd ≥ 0 (with the sign convention of the constraints), so d·√T∇L(θ†) = √T Σ multiplier_j g_jᵀd ≥ 0, and it equals 0 exactly on 𝒩 and tends to +∞ for every fixed d ∉ 𝒩 (strict complementarity). Hence h_T epi-converges to d·Z + ½dᵀHd on 𝒩 and +∞ elsewhere; the minimiser of a strictly convex quadratic on a subspace is V(VᵀHV)⁻¹VᵀZ̃·(−1). By the argmin theorem for epi-convergent convex functions with a unique, tight minimiser (Geyer 1994; Davis, Knight and Liu 1992), √T(θ̂_T − θ†) converges in law to it. ∎ Without strict complementarity the limit is the projection onto a cone (not Gaussian); not treated.

### Numerical results (third pass)

- **Prop. 16.9 on the power-law truth** (600 000 events, Φ by SLSQP). Φ(0) and the crossing μ_0 with Φ(μ_0) = 1, against the pseudo-truth μ†:
 | n* | K | μ† | ε | Φ(0) | μ_0 (where Φ crosses 1) |
 |---|---|---|---|---|---|
 | 0.8 | 1 | 2.71 | 0.30 / 0.05 | 1.63 / 1.20 | > 1 / ≈ 0.67 |
 | 0.8 | 2 | 1.91 | 0.30 / 0.05 | 1.48 / 1.09 | > 1 / ≈ 0.38 |
 | 0.8 | 3 | 1.55 | 0.30 / 0.05 | 1.44 / 1.06 | > 1 / ≈ 0.28 |
 | 0.95 | 3 | 2.99 | 0.30 / 0.05 | 1.43 / 1.06 | > 1 / ≈ 1.0 |
 Φ(0) > 1 in every case, so the certificate is non-vacuous: μ† is strictly positive and bounded below by 0.3 to 1.0 for ε = 0.05. It is conservative by a factor 2 to 7 (μ_0 ≈ 0.28 vs μ† = 1.55 at K = 3), tighter at larger ε because the feasible α set shrinks. (Crossings read off the tabulated Φ values by interpolation; Φ was evaluated at μ = 0, 10⁻³, 0.05, 0.2, 0.5, 1.)
- **Theorem 16.10, K = 2, n* = 0.8 power-law truth, 1 − ε = 0.5 (active; multiplier κ ≈ 0.05), 300 replicates of 4·10⁴ events.**
 - n̂ = 0.5 to 10 digits in all replicates; empirical covariance has smallest eigenvalue 0 (rank 2 of 3), as predicted.
 - Variance of (μ̂, α̂_1, α̂_2): empirical (3.41·10⁻³, 2·10⁻⁵, 5.4·10⁻⁴); formula with the population J (estimated from the replicate scores at a reference θ†) (3.45·10⁻³, 2·10⁻⁵, 5.4·10⁻⁴): agreement to 1 %; corr(μ̂, α̂_1) −0.209 empirical vs −0.210 formula.
 - The unconstrained sandwich H⁻¹JH⁻¹ gives (18.5·10⁻³, 12·10⁻⁵, 6.8·10⁻⁴) for the same data, i.e. overstates var(μ̂) 5×.
 - With J estimated by batch means (B = 20) on a single path the formula gave 2.5·10⁻³ for var(μ̂) and corr −0.085: the discrepancy is the batch-means estimator of J (the same 17 to 27 % shortfall as in section 8), not the theorem.

### Verdict (third pass)

Closed: the a priori lower bound on μ† (Prop. 16.9, a computable certificate that is valid and non-vacuous in all runs) and the limit law under active constraints (Theorem 16.10, tested to 1 %). Remaining caveats: the CLT underlying Theorem 16.10 is itself [P*] (mixing), strict complementarity is assumed (the degenerate case gives a non-Gaussian cone-projection limit), and the batch-means estimator of J is biased low by 17 to 27 % at B = 20 on 4·10⁴ events; a better long-run variance estimator (HAC, larger B with longer paths) was not tried.
