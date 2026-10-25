# Conclusions from the four papers in `refPapers/` (25 Oct 2026)

Papers read: Spreij 1986 (Stochastics 18:277–312) in full, with the equations read from page images because the text layer is garbled; Spreij 1991 (J. Multivariate Anal. 39:236–245) in full; Piggott and Solo 2018 (ICASSP) and Pasha and Solo 2019 (CDC) in full (text layer clean). Statements below are tagged **[read]** when I saw the statement on the page, **[inferred]** when it is my reading of what a result implies for our work. Extracted text and page images are in the scratchpad, not in the project. All papers are well specified; none treats a Hawkes process with a misspecified kernel.

## 1. What each paper actually proves

### Spreij 1986, "Recursive parameter estimation for counting processes with linear intensity"
- **Model [read].** A counting process with intensity λ_t = θᵀφ_t, where φ is an observed predictable ℝ^d_+-valued process and the "true" θ₀ is a point of the model. No Hawkes example is given (examples: φ = (1, 1 + sin t), a two-state Markov switch), but a Hawkes dictionary φ_t = (1, S_1(t−), …, S_K(t−)) is exactly of this form.
- **Algorithm (2.18)–(2.21) [read].** dx_t = Q_tφ_t/(φ_tᵀθ̂_{t−}) (dn_t − φ_tᵀx_t dt), dQ_t = −Q_tφ_tφ_tᵀQ_t/(φ_tᵀθ̂_t) dt, so Q_t⁻¹ = Q_0⁻¹ + ∫φφᵀ/(φᵀθ̂) ds; and θ̂_t = x_t if φ_tᵀx_t > φ_tᵀε, else ε, with a known vector ε whose every component is below the corresponding component of θ₀. He calls it the approximate-maximum-likelihood (AML) algorithm, and derives it three ways (implicit-function argument, a second-order filter, the LAN expansion).
- **Theorem 3.2 [read].** If Φ_tᵀθ₀ → ∞ a.s. (Φ_t = ∫₀ᵗφ) and liminf (1/Φ_tᵀθ₀)∫φφᵀ/(φᵀθ₀) ds = C > 0, then θ̂_t → θ₀ a.s. and the weighted prediction error (1/θ₀ᵀΦ_t)∫(φᵀ(θ̂_s − θ₀))²/(φᵀθ₀) ds → 0 a.s.
- **Proof device [read].** Stochastic Lyapunov function x̃ᵀQ_t⁻¹x̃ with x̃ = x − θ₀ (equations (3.1)–(3.3)), and a nonnegative almost-supermartingale lemma (Lemma 3.1).
- **Theorem 5.4 and the remarks [read].** Under a deterministic-normalisation condition (P-lim μ(t)⁻¹θ₀ᵀΦ_t = 1) one has Q̄_t^{1/2}M_t → N(0, I); the covariance is the Cramér–Rao bound (6.2). Remark 6.1 handles a random limit of the covariation by a random time change τ = Φ_tᵀθ₀.

### Spreij 1991, "Recursive approximate maximum likelihood estimation for a class of counting process models"
- **Setting [read].** The same model; Assumption 1: θ₀ lies in a known compact set, ε < θ₀i < 1/ε for all i. The algorithm (2a)–(2e) is the 1986 recursion with a two-sided safeguard (indicators on φᵀx > εφᵀ1 and φᵀx < ε⁻¹φᵀ1).
- **Theorem 3 [read].** If Φ_t → ∞ a.s. and there is f with f(x)/x → ∞ and sup_t f(log λ_max(Q_t⁻¹))/λ_min(Q_t⁻¹) < ∞ a.s. (the Lai–Wei condition of stochastic-regression least squares), then θ̂_t → θ₀ a.s. The proof is the Lyapunov computation of Lemma 4 (P_t = x̃ᵀQ_t⁻¹x̃ = O(log λ_max)) plus martingale convergence lemmas.
- **Theorem 7 [read].** Under a random-time-change condition, Q_t^{−1/2}(θ̂_t − θ₀) → N(0, I). The closing remark says that for deterministic φ the asymptotic variance equals the Cramér–Rao bound. The discretised implementation is given in the appendix.

### Piggott and Solo 2018, "Non-Negative Online Estimation for Hawkes Process Networks"
- **Model [read].** Multivariate **linear** Hawkes with the kernel expanded in a Laguerre (Erlang) basis, all components sharing one decay constant; discrete time bins.
- **Algorithm [read].** Exponentially weighted **least squares** (the quadratic criterion) solved by a one-step non-negative-matrix-factorisation multiplicative update, with a max(·, ε) to avoid absorbing zeros; fixed gain; a two-time-scale recursion for the Laguerre states.
- **Claims [read].** Simulations on a 4-node network; the authors say projections "can lead to slow convergence and even stalling". Non-negativity is enforced; **the stability constraint is not**: "Extensions enforcing the stability constraint … will be pursued elsewhere as will be the challenging task of algorithm stability analysis". No theorem.

### Pasha and Solo 2019, "Adaptive Estimation of a Vector Point Process Model"
- **Model [read].** A **log-linear** (nonlinear) Hawkes model, λ = exp(cᵀ + Σ gk,j ⋆ dN) with a Laguerre basis; positivity is automatic, so there are no sign constraints and inhibition is allowed. Not the linear model of this project.
- **Algorithm [read].** Fixed-gain instantaneous steepest ascent on the point-process likelihood.
- **Theory [read].** A heuristic averaging analysis under stationarity (A1) and well-specification (A2: true parameters exist and the compensated process is a martingale): the averaged system dθ/dτ = μ[g(θ_e) − g(θ)] with g(θ) = E[φλ(θ)] is locally stable if R_{θ_e} = E[φφᵀλ] is positive definite. They write that positive definiteness of R is "plausible but … non-trivial and will be pursued elsewhere" and that a rigorous analysis "would require a full paper". The applied part is a genomic data set.

## 2. What this changes for our results

**A. The online algorithm itself is not new (Problems 5–7, C1).**
The EF recursion of Problems 5 and 17 (F_i = F_{i−1} + ∫φφᵀ/λ_{θ_i}, θ_{i+1} = Π(θ_i − F_i⁻¹g_i)) is the Spreij AML recursion applied to the Hawkes regressor, in event-indexed form with the gap integral done exactly (Spreij lets x_t drift continuously between events; his discretised version in the appendix of the 1991 paper is the time-step form [inferred]). The Fisher-preconditioned online Newton must therefore be presented as **Spreij's recursion applied to Hawkes processes**, and C1 reworded: what is ours is (i) the log-loss regret analysis with (d/2) ln T and explicit ε-dependence, (ii) the exact convex (F-norm) projection onto Θ_ε in place of a componentwise lower clip, (iii) the observed-information variant, (iv) the misspecified analysis.
He never evaluates his conditions for a Hawkes process, and gives no regret bound, no rate in ε and no near-critical analysis [read: none appear].

**B. The well-specified part of Theorem 17A and 17B is anticipated (Problem 17).**
- *Convergence.* Spreij 1986 Thm 3.2 and Spreij 1991 Thm 3 give a.s. convergence of the EF-type recursion in the well-specified case, under conditions that are **weaker than ours** (no ergodicity, no exponential moments, no Markov structure; only growth of the extreme eigenvalues of Q_t⁻¹) [read]. Our Lemma 17.2 is the misspecified generalisation of his Lyapunov identity (3.1)–(3.3): the F-norm Lyapunov function is the same, and the drift term that he bounds by a summable quantity becomes, under misspecification, the extra cross term 2∫(λ_θ − λ_†)(λ_† − λ*)/λ_θ du that cannot be controlled pathwise [inferred].
- *Limit law.* Spreij 1991 Thm 7 and the closing remark give asymptotic normality with the Cramér–Rao covariance. That is our Theorem 17B at A = I_F⁻¹H = I, J = H, Σ = Λ*I⁻¹.
- *Consequence for the well-specified near-critical case [inferred].* Spreij's Theorem 3 needs only λ_min(Q_t⁻¹) → ∞ with the Lai–Wei growth condition. Our Problems 2 and 12 give lower bounds on λ_min of the observed/expected Fisher information, so a.s. consistency of the EF recursion for each fixed ε would follow from his theorem plus those bounds, without the concentration machinery of Problem 5 (Theorems F, E5). This does not give regret or ε-explicit rates.

**C. What remains ours in Problem 17 (not in any of the four papers).**
1. Misspecification: the target is the KL pseudo-truth θ†, not a point of the model; the noise includes a non-martingale term a(λ*/λ† − 1), which is why Spreij's martingale route fails and the Poisson-equation route is needed.
2. The drift matrix A = Ī⁻¹H(θ†) with A ≠ I in general, the Lyapunov equation (A − ½I)Σ + Σ(A − ½I)ᵀ = Λ*Ī⁻¹JĪ⁻¹, the exact eigenvalue-1 fact (Lemma 17.5), and the threshold λ_min(A) > ½.
3. The observed-information preconditioner and its sandwich limit Λ*H⁻¹JH⁻¹ with no side condition.
4. A convex constraint set that can bind: Spreij's safeguard needs θ₀ strictly inside a coordinate box (all components, including every α_k, bounded away from 0 and from above); ours handles α_k = 0 at θ† and the stability face n = 1 − ε, with the active-constraint limit law and the O(1/n) distance to the face (Theorem 17C).
5. Log-loss regret, which is not in his papers at all.
**Overlap-risk for Problem 17: High for the well-specified Theorems 17A/17B, Low-to-Medium for the rest.** In the paper, state 17A/17B in the misspecified form and cite Spreij for the well-specified special case.

**D. Problem 16.** None of the four papers has a pseudo-truth. Spreij's identifiability condition, liminf (1/Φᵀθ₀)∫φφᵀ/(φᵀθ₀) > 0, is the well-specified analogue of our Lemma 16.2/Prop. 16.6, but our statement is for a non-model truth. Pasha–Solo's open item (positive definiteness of R_{θ_e}) is analogous to Lemma 16.2; in the linear model it is proved there (linear independence of exponentials on a gap).

**E. Problems 13–15 (fixed gain, tracking).** Piggott–Solo and Pasha–Solo are the only fixed-gain Hawkes recursions found and **neither has a theorem**: their stability analysis is stated as future work, and the Pasha–Solo averaging sketch assumes well-specification and unverified positive definiteness. Our constant-gain tracking bound (Problems 13–14) is not anticipated by them [read: they contain no tracking bound]. Their remark that projection "can lead to slow convergence and even stalling" matches what I found for observed-information with μ_min = 10⁻³ (1/300 runs stuck at the corner), but my exact F-norm projection never stalled for the expected-Fisher recursion in K ≤ 2 tests; their claim concerns an unspecified projection and I cannot say it is wrong.

**F. Problem 18.** Not touched by any of the four papers.

## 3. Corrections to what I wrote earlier
1. `literature_difference_statements.md` treated Piggott–Solo and Pasha–Solo together as "fixed-gain recursions on a Laguerre/log-intensity parametrisation, non-negativity preserved": wrong in detail. Piggott–Solo is **linear** Hawkes with a **least-squares** criterion and a multiplicative non-negative update; Pasha–Solo is a **log-linear** model with a **likelihood** gradient and no constraints. Neither enforces stability.
2. `problem17_online_convergence.md` and the difference statements did not credit the algorithm to Spreij. The EF recursion should be attributed to Spreij (1986, 1991).
3. `problem05_online_newton.md` and the topic statement say "Fisher-preconditioned online Newton or recursive MLE" as a new construction; this should be reworded as in §2A.

## 4. Still unverified
Whether Spreij's conditions hold for a near-critical Hawkes dictionary (never checked by him, and I did not check them either); his remaining sections (identifiability, the proof of Theorem 5.4's lemmas) read only for the statements; whether Spreij's later papers (e.g. *Stochastic Processes and their Applications* 26 (1987) 199, "Recursive estimation for a class of counting processes", DOI 10.1016/0304-4149(87)90099-8) extend the framework to misspecification or Hawkes; Herrera-Marín II Theorem 3.1; Rong–Nair's last 7 %.

## 5. Follow-up checks (25 Oct 2026, later)

### 5.1 Rong–Nair (arXiv 2602.20795): does it contain Prop. 18.2? **Yes (checked in the PDF).**
- Lemma 4(b): the Bartlett spectrum of a stationary Hawkes process is C̃(ω) = Λ/|1 − φ̃₀(ω)|².
- Lemma 13 and eq. (5.2)–(5.5): R* = Var[χ(0)] and R*^(1,0) = Cov[χ(0), φ₀(0)] are the double integrals ∫∫q(u)C(u−v)q(v) and ∫∫q(u)C(u−v)φ₀(v), with the Parseval forms (1/2π)∫q̃ C̃ q̃(−ω) dω and (1/2π)∫q̃ C̃ φ̃₀(−ω) dω.
- Theorem 15 (eq. (5.6)–(5.7)): α* = R*⁻¹R*^(1,0), c* = Λ − μᵀα* = Λ(1 − Γ*), Γ* = 𝟙ᵀα*.
- Consequences: my **Prop. 18.2** (weighted L² with the Hawkes spectrum) and the **surrogate rate matching** (Cor. 2(iv) of Problem 16) are **known** and are no longer claimed (notes edited). Their Γ* < 1 is not proved (c* may be negative); for the KL pseudo-truth n† < 1 follows from rate matching because μ > 0 is part of the domain (Cor. 2(i)): that distinction stays ours.
- Their CLT without mixing uses the linear martingale representation of λ₀ via the Hawkes resolvent; the LS score is bilinear in path functionals of N, so this works. The KL score has the factor 1/λ†, so the route does not transfer: the mixing/Poisson-equation caveat of Prop. 16.7, Theorem 16.10 and Theorem 17A stays.
- Wei et al. (the remaining ~25 % of the ar5iv text): consistency and the CLT are taken from "modifications of Ogata (1978) Theorems 2, 4, 5" (the classical MLE theory for stationary point processes), not from a mixing or Poisson-equation argument; the identifiability proof (Appendix D.1) and the KL interpretation were in the part I had already read; only equality constraints (for the test) appear, no KKT/inequality constraints; the piecewise-constant model is linear in the bin heights, the same linear-in-parameter structure as ours. Their proof of Lemma 1 itself was not in the part I could retrieve.

### 5.2 Spreij's Theorem 3 for a Hawkes dictionary: **the conditions hold** (checked, numerically and by a short argument).
**Corollary S (Spreij 1991 Thm 3 applied to Hawkes) [P modulo his theorem].** Let N be a Hawkes process with distinct β_k, true θ₀ = (μ, α) with every component positive, n(θ₀) < 1, any initial state, and a known ε with ε < θ₀i < 1/ε. Then Spreij's AML recursion with the box safeguard converges to θ₀ a.s.
*Proof.* φ_t = a_t = (1, S(t−)) is predictable and dN = θ₀ᵀφ dt + dm. Φ_t ≥ t → ∞. With θ̂ ∈ [ε, 1/ε]^d, φᵀθ̂ ∈ [εφᵀ𝟙, ε⁻¹φᵀ𝟙], so (his Lemma 2) ε G_t ≼ Q_t⁻¹ − Q₀⁻¹ ≼ ε⁻¹ G_t with G_t = ∫₀ᵗ φφᵀ/(φᵀ𝟙) ds. The map a ↦ aaᵀ/(𝟙ᵀa) is Lipschitz (gradient bounded by 3 on 𝟙ᵀa ≥ 1), the stationary law exists and the chain contracts in W_d at rate e^{−εbt} (Part 1, Theorem `contraction`), so G_t/t → Ḡ = E_π[aaᵀ/(𝟙ᵀa)] a.s. from any initial state (the coupling error is an a.s. finite integral). Ḡ ≻ 0 by the argument of Lemma 16.2 (exponentials on a gap are linearly independent; the weight 1/(𝟙ᵀa) is positive). Hence λ_min(Q_t⁻¹) ≥ ct and λ_max(Q_t⁻¹) ≤ Ct eventually, and with f(x) = x² the quantity f(log λ_max)/λ_min ≤ C′(log t)²/t → 0. Theorem 3 applies. ∎
*What it gives and does not give.* It gives a.s. consistency of the well-specified recursion at every fixed n < 1 without the concentration machinery of Problem 5, but no rate, no regret, no ε-explicit constants, and it is for Spreij's coordinate safeguard, not our F-norm projection. It does not cover ε = ε_T → 0 (Problem 11).
*Numerics* (`code/p17/spreij_cond.py`; K = 2, β = (0.3, 1.5), μ = 1; G_t by 10-point Gauss–Legendre per gap, 1.5–2·10⁶ events): λ_min(G_t)/t and λ_max(G_t)/t converge: 0.0146 and 9.0 (n = 0.70, t = 6·10⁵), 0.0082 and 28 (n = 0.90, t = 2·10⁵), 0.0029 and 94 (n = 0.97, t = 4.6·10⁴). The Spreij ratio f(log λ_max)/λ_min with f(x) = x^{1.5} falls below 1 at t ≈ 3·10³ (n = 0.70), 10⁴ (n = 0.90) and 3·10⁴ (n = 0.97), reaching 7·10⁻³, 4·10⁻² and 0.46. So the hypothesis is met, but the transient before it "kicks in" lengthens as ε falls (the condition number of G_t is ≈ 3·10⁴ at n = 0.97). (A first run with α = (0.12, 0.9) was supercritical, n = 1.27, and is not reported.)
