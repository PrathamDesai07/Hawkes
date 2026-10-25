# Problem 17: online convergence to θ† under misspecification, with the stability constraint active (contribution C6)

Date: 25 Oct 2026. Tags: **[P]** proved on paper, **[P*]** proved modulo stated inputs, **[N]** numerical, **[O]** open. Builds on `problem16_kl_pseudotruth.md` (Theorem 16, Props 16.4–16.10), `problem05_online_newton.md`, `problem05_circularity.md` (exact F-norm projection, non-expansive), `problem04_projection.md`. Not refereed.

## Attribution (added after reading Spreij 1986/1991, see `refpapers_conclusions.md`)

The EF recursion below is **Spreij's AML recursion** (Stochastics 18 (1986) 277–312; J. Multivariate Anal. 39 (1991) 236–245) applied to the Hawkes regressor φ = a_t, in event-indexed form with the gap integral done exactly. For a well-specified model Spreij proves a.s. consistency (1986 Thm 3.2; 1991 Thm 3, under the Lai–Wei growth condition on the eigenvalues of Q_t⁻¹, with a known coordinate box) and the Cramér–Rao CLT (1986 Thm 5.4; 1991 Thm 7). Hence **Theorems 17A and 17B at A = I (θ† = θ*, J = H = I_F) are his results and are not claimed here**. What is new in this note: the misspecified target θ†, the drift matrix A = Ī⁻¹H and the Lyapunov-equation limit law with the threshold λ_min(A) > ½ (and Lemma 17.5), the observed-information variant, a constraint set that can bind (the stability face n = 1 − ε, and α_k = 0 at θ†) handled by the exact F-norm projection, and the active-constraint limit law. Lemma 17.2 is the misspecified generalisation of his Lyapunov identity (3.1)–(3.3).

## 0. What has to be settled

The tracker asks: show the online estimator converges to θ† under misspecification, with n ≤ 1 − ε active, as an almost-sure or in-probability theorem. Three things are new relative to the well-specified analysis of Problem 5:
1. **The telescoping trick breaks.** Prop. 1 of Problem 5 uses E⟨g,v⟩ = E vᵀΔFv, which holds only when the data law is the model with θ*. Under misspecification the mean gradient does not vanish at the target, and in the constrained case it does not vanish even at θ†.
2. **The preconditioner is no longer the right curvature.** With the expected Fisher increment ΔF = ∫φφᵀ/λ_θ, F_n/n converges to I_F(θ†) = E[aaᵀ/λ†] (per time), whereas the curvature of the limit loss is H(θ†) = E[λ* aaᵀ/λ†²]. They coincide only if λ* = λ†. Their mismatch decides the rate.
3. **The constraint.** θ† can sit on the face n = 1 − ε (Theorem 16, Cor. 2). The iterates must reach and stay near it.

## 1. Setting and algorithm

Data, model and Θ_ε as in Problem 16 (stationary ergodic truth with predictable intensity λ*, dictionary {β_k}, a_t = (1, S_t), Θ_ε = {μ ≥ μ_min, α ≥ 0, cᵀα ≤ 1 − ε}, c_k = 1/β_k). Add a compact convex box Θ_c ⊇ Θ_ε ∩ {μ ≤ μ_hi} with μ_hi ≥ Λ* (Prop. 16.4 shows μ† ≤ Λ*), so θ† ∈ Θ := Θ_ε ∩ Θ_c.

Gap i: length w_i, a_i = (1, S(t_i−)), b_i = ∫_gap (1, S) dt, φ(u) = (1, S(t_{i−1}+) e^{−βu}), g_i = b_i − a_i/(a_iᵀθ_i). Online projected Newton:
 F_i = F_{i−1} + ΔF_i, θ_{i+1} = Π^{F_i}_Θ(θ_i − F_i⁻¹g_i), F_0 ≻ 0,
with the exact metric projection Π^F in the F-norm (dense QP, `problem04_projection.md`; for d ≤ 4 it is an enumeration of active sets). Two choices of the increment:
- **Expected Fisher (EF):** ΔF_i = ∫₀^{w_i} φφᵀ/(φᵀθ_i) du (the Problem 5 algorithm);
- **Observed information (OI):** ΔF_i = a_ia_iᵀ/λ_i² with λ_i = a_iᵀθ_i (the Hessian of the per-event loss, Hazan–ONS style).

Notation (per unit time, in the stationary law): Λ* the true rate; H(θ) = E[λ* aaᵀ/λ_θ²] = ∇²L(θ); I_F(θ) = E[aaᵀ/λ_θ]; J = long-run covariance of the score at θ† (Prop. 16.7). Per-event means are these divided by Λ* (n events take time n/Λ*).

## 2. Results

**Lemma 17.1 (one-gap identities) [P].** Suppose that between events the true intensity λ*(u) is a deterministic function of the past (true for any multi-exponential Hawkes truth; the general case needs the compensator and only changes the notation). Then for F_{t_{i−1}}-measurable θ_i and v = θ_i − θ†, with Λ*(u) = ∫₀^u λ*:
 E[⟨g_i, v⟩ | past] = ∫₀^∞ e^{−Λ*(u)} (λ_θ − λ_†)(λ_θ − λ*)/λ_θ du,  E[vᵀΔF_i^{EF} v | past] = ∫₀^∞ e^{−Λ*(u)} (λ_θ − λ_†)²/λ_θ du,
where λ_θ(u) = φᵀθ_i, λ_†(u) = φᵀθ†. In the well-specified case λ_† = λ* and the two coincide (Problem 5, Prop. 1(a)). In general
 E[vᵀΔF^{EF}v − 2⟨g, v⟩ | past] = −∫e^{−Λ*}(λ_θ − λ_†)²/λ_θ du + 2∫e^{−Λ*}(λ_θ − λ_†)(λ_† − λ*)/λ_θ du.
*Proof.* Fubini with P(τ > u) = e^{−Λ*(u)}, density λ*(u)e^{−Λ*(u)}, as in Prop. 1. ∎

**Lemma 17.2 (stationary mean drift of the Lyapunov function) [P].** Let W_n = ‖θ_n − θ†‖²_{F_{n−1}}. For a projected step, non-expansiveness of Π^F with θ† ∈ Θ gives the pathwise inequality
 W_{n+1} ≤ W_n + vᵀΔF_n v − 2⟨g_n, v_n⟩ + g_nᵀF_n⁻¹g_n  (v_n = θ_n − θ†),
valid also when θ† is on the boundary. In the stationary regime the time-average of the middle terms for a frozen θ is
 E[(vᵀa)²/λ_θ] − 2vᵀ∇L(θ) ≤ −vᵀ(2H̄ − I_F(θ))v,  H̄ = ∫₀¹ H(θ† + sv) ds,
because vᵀ∇L(θ) = vᵀ∇L(θ†) + vᵀH̄v and vᵀ∇L(θ†) ≥ 0 for every feasible θ (first-order condition on Θ, Theorem 16(b)), i.e. the constraint multiplier term has the *favourable* sign. So W_n is a stochastic Lyapunov function exactly where
 **(C½)  2H̄(θ, θ†) ≻ I_F(θ)  along segments to θ†** (EF);  for OI the middle terms give E[(vᵀa)²λ*/λ²] − 2vᵀ∇L ≤ −vᵀ(2H̄ − H(θ))v, and (C½) holds automatically near θ† (H̄ → H) and in general by convexity of s ↦ H along the segment is not claimed.
*Proof.* Expand ‖Π(y) − θ†‖²_{F_n} ≤ ‖y − θ†‖²_{F_n}, y = θ_n − F_n⁻¹g_n, and ‖v_n‖²_{F_n} = ‖v_n‖²_{F_{n−1}} + vᵀΔF_n v. The mean-field step is Lemma 17.1 averaged over the stationary law (Palm = time average). ∎

**Theorem 17A (almost-sure convergence, interior or active constraint) [P*; a self-contained proof under the Lyapunov condition (D_κ) is Theorem 17A* in section 7].** Assume (A1)–(A2) of Problem 16, geometric ergodicity of the (truth ⊗ dictionary) Markov state with exponential moments (true for subcritical multi-exponential Hawkes truth: Problem 3, Theorem `contraction`, Lemma `expmoment`), F_0 ≻ 0, Θ compact convex. Then for EF and OI, θ_n → θ†_Θ a.s., where θ†_Θ is the unique minimiser of L on Θ (Theorem 16), whether or not the constraint cᵀα ≤ 1 − ε binds.
*Proof route (ODE method with Markov noise and matrix gain).* (i) *Bounded gain.* On the compact set Θ, λ_θ ≥ μ_min and ΔF_n(θ) is Lipschitz in θ, so by the ergodic theorem uniformly over Θ, F_n/n converges to a matrix F̄(θ_n)/Λ*-type limit with F̄ ∈ [m_1 I, m_2 I] eventually a.s. (lower bound: the argument of Prop. 16.6 with E[aaᵀ/(M+Y)] ≻ 0). So the gain Γ_n = n F_n⁻¹ is bounded and bounded away from 0, and uniformly positive definite. (ii) *Poisson equation.* The mean field h(θ) = ∇L(θ)/Λ*, and the gradient noise g_n(θ) − h(θ) is a functional of the geometrically ergodic Markov chain; the Poisson equation has a bounded solution with Lipschitz dependence on θ (Benveniste–Métivier–Priouret, Part II, Thm 17), giving a martingale-difference noise plus a summable remainder. (iii) *Projected ODE.* The limit ODE is θ̇ = −Π_{T_Θ(θ)}(F̄(θ)⁻¹∇L(θ)) (projected dynamics); L is a strict Lyapunov function (d/ds L = −∇Lᵀ Γ ∇L ≤ 0 with equality only at KKT points), and L strictly convex on convex Θ has a unique KKT point θ†_Θ. (iv) By Kushner–Yin Thm 5.2.1 / BMP the iterates converge a.s. to the globally asymptotically stable set {θ†_Θ}. ∎ [P*] because (i) and the Poisson-equation regularity are cited, not re-proved; the algebra (Lemma 17.2, rates below) is derived.

**Theorem 17B (rate and limit law, free case) [P*].** Suppose θ† is in the relative interior of Θ, and let Ī be the limit of the preconditioner (Ī = I_F(θ†) for EF; Ī = H(θ†) for OI), A = Ī⁻¹H(θ†). If **λ_min(A) > ½** then
 √n(θ_n − θ†) ⇒ N(0, Σ),  (A − ½I)Σ + Σ(A − ½I)ᵀ = Λ* Ī⁻¹ J Ī⁻¹  (continuous Lyapunov equation),
and E‖θ_n − θ†‖² = O(1/n). If λ_min(A) < ½ the rate degrades to n^{−λ_min(A)} (Chung–Fabian). Consequences:
- *OI:* A = I and **Σ = Λ* H⁻¹JH⁻¹**, i.e. Cov(θ_n) ≈ H⁻¹JH⁻¹/T with T = n/Λ*: the online estimator attains the sandwich covariance of the batch QMLE (Prop. 16.7), asymptotically efficient *within* the misspecified class, with no condition.
- *EF:* Σ coincides with the OI value only if A = I; (C½) is λ_min(I_F⁻¹H) > ½, the threshold for the 1/√n rate. Well specified: A = I, J = H = I_F, Σ = Λ* I⁻¹ (Cramér–Rao).
*Proof.* Linearise the recursion at θ†: v_{n+1} = (I − (1/n)A)v_n + (Λ*/n)Ī⁻¹ξ_n where ξ_n are the (event-indexed) noises with covariance J/Λ*; standard SA CLT (Fabian 1968; Kushner–Yin Thm 10.1.1) with the Markov-noise long-run covariance. The noise scaling: Σ_{i≤n} g_i ≈ T_n ∇L̂ has covariance T_n J = (n/Λ*)J, and F_n ≈ n Ī/Λ*. ∎

**Theorem 17C (active constraint) [P*].** Assume θ† on a face with active set 𝒜, strict complementarity and independent active normals (as in Theorem 16.10). Let V span 𝒩 = {d : g_jᵀd = 0, j ∈ 𝒜}, and let A_V = (VᵀĪV)⁻¹VᵀHV (restricted drift). Then
(a) *the constraint is hit and respected:* the distance d_n = (1 − ε) − cᵀα_n ≥ 0 satisfies n·d_n = O_p(1) and is recurrent at 0 (the projection is active a positive fraction of the steps, e.g. Prop. in the circularity note: 0.2–5 %); in particular √n·d_n → 0 in probability;
(b) if λ_min(A_V) > ½ then √n(θ_n − θ†) ⇒ V ζ, ζ ∼ N(0, Σ_V) with (A_V − ½I)Σ_V + Σ_V(A_V − ½I)ᵀ = Λ*(VᵀĪV)⁻¹VᵀJV(VᵀĪV)⁻¹. For OI this is exactly Λ*·(VᵀHV)⁻¹VᵀJV(VᵀHV)⁻¹ = the constrained sandwich of Theorem 16.10 (with T = n/Λ*), which is again the batch constrained-QMLE covariance.
*Sketch.* Along the face the unconstrained step has component −Γ∇L = Γκ g_𝒜 pushing outward at rate κ/n against noise of the same order 1/n per step, so d_n is a reflected walk on the scale 1/n; (a) follows. On the face the iteration is the F-metric projected dynamics, whose linearisation in V coordinates is (VᵀĪV)⁻¹VᵀH V; (b) is then 17B in V coordinates. Strict complementarity is what keeps the face attracting; without it the limit is a cone-projected non-Gaussian law. [P*]: the reflection argument is a sketch.

**Proposition 17.4 (regret remains (d/2) ln T for the pseudo-regret against θ†, not against a truth) [O].** Under misspecification the regret against the best fixed θ† of the log-loss, Σℓ_i(θ_i) − Σℓ_i(θ†), has expectation (1/2)·tr(H⁻¹J)·ln n + O(1) for OI when the free sandwich applies (it equals d/2 ln n only when J = H, e.g. well specified). Not proved here; stated as the natural conjecture because the Lyapunov/Newton telescoping now produces (1/2)·tr(Ī⁻¹ J_ev) ln n in place of (d/2) ln n. Checked numerically below only through the covariance, not the regret itself.

**Lemma 17.5 (the drift matrix has eigenvalue 1 at θ†, interior case) [P].** If θ† is interior, then H(θ†)θ† = I_F(θ†)θ† (both equal E[a]); hence A = I_F⁻¹H has the eigenvector θ† with eigenvalue exactly 1, and H − I_F = E[(λ* − λ†)aaᵀ/λ†²] annihilates θ†.
*Proof.* Hθ† = E[λ* a (aᵀθ†)/λ†²] = E[λ*a/λ†] and I_Fθ† = E[a aᵀθ†/λ†] = E[a]; they are equal by the first-order condition E[a(1 − λ*/λ†)] = 0 (Theorem 16(b), interior). ∎ So at most d − 1 eigenvalues of A can differ from 1, and the rate condition (C½) can fail only through the other directions; numerically (below) every eigenvalue of A lies in [0.97, 1.01] even for severe misspecification.

## 3. What is not proved

1. **Partly closed (section 7).** Theorem 17A* is a self-contained proof (ergodic theorem, martingale strong law, non-expansive projection; no Poisson equation) under the drift condition (D_κ), which holds for κ small enough (EF) whenever sup λ†/λ* < ∞, and was checked numerically at κ ≤ 1 (EF). The statement for every gain constant without (D_κ) still relies on the cited ODE/Poisson-equation theory (Theorem 17A, [P*]); the numerics converge for κ = 2, 4 where (D_κ) fails, so the condition is a limitation of the proof. For OI, (D_κ) needs a ratio-restricted Θ (λ_θ ≥ κλ†/2). (D_κ) was verified on a random sample of Θ, not proved for all Θ in any example.
2. Theorem 17C(a): the excursion bound for d_n is a heuristic reflected-walk argument.
3. Strict complementarity (non-Gaussian cone limit otherwise), and the case λ_min(A) ≤ ½ beyond the quoted exponent.
4. The regret statement (Prop. 17.4).
5. K ≥ 2 beyond d = 3 numerics.

## 4. What to code and pre-registered predictions

C simulator with the algorithm in the loop (`code/p17/online17.c`): Ogata truth (power law as in Problem 16), model states from the observed path, EF/OI preconditioners, exact F-norm projection onto {μ ≥ μ_min, α ≥ 0, cᵀα ≤ 1 − ε} by active-set enumeration, checkpoints at n = 10², …, with R replicates. Python analysis (`code/p17/`): reference θ†, H, I_F, J from long batch runs (as in Problem 16), Lyapunov predictions.

| # | Test | Prediction | Pass |
|---|---|---|---|
| W1 | Well specified (code check): both algorithms | Cov ≈ Λ*I⁻¹/n | within MC error |
| M1 | Misspecified, free (ε large), OI | θ_n → θ†; n·Cov(θ_n) → Λ*H⁻¹JH⁻¹ | agreement to ≲ 15 % |
| M2 | Misspecified, free, EF | n·Cov → Σ from Lyapunov with A = I_F⁻¹H; report λ_min(A) vs ½ | agreement; λ_min(A) value |
| M3 | Misspecified, active (n ≤ 1 − ε below n(θ†_∞)), OI and EF | θ_n → constrained θ†; n·d_n = O(1), projection active in a positive fraction of steps; n·Cov → Σ_V (rank d − 1) | agreement ≲ 15 %; fraction active |
| M4 | Rate | E‖θ_n − θ†‖² ∝ 1/n | slope −1 |
| M5 | Wrong-preconditioner check | EF vs OI covariances differ by predicted amount | |

## 5. Results of the coding study (all [N])

Code: `code/p17/online17.c` (algorithm in the loop, exact F-norm projection by active-set enumeration, EF and OI), `lib17.py`, `w1.py`, `w1b.py`, `trap.py`, `ev.py`, `m1.py`, `m2.py`, `m3.py`. Truth: truncated power law (n* = 0.8, η = 0.5, τ = 1, μ* = 1), model K = 2 with β = (0.3, 1.5), Λ* = 5.00. Start θ_1 = (2, 0.05, 0.2), F_0 = I, μ_min = 0.3, warm-up 5 000 events before the algorithm starts. θ†, H, I_F from batch solutions on 2×1.5·10⁶ events; J from 600 independent paths of 10⁵ events (a first version of J from 200 paths of 4·10⁴ events was biased low by long memory and gave 15–25 % errors in the predictions).

- **W1 (well specified, code check).** K = 2, α = (0.12, 0.45), n = 0.7, 150 replicates, n up to 3·10⁴. Both EF and OI: n·Var(θ_n) plateaus at (52–55, 1.2–1.4, 14.7–17.0) against Λ*I⁻¹ = (58.5, 1.33, 15.2) (±12 % MC); I_F and H agree to 0.3 % (as they must). Mean errors → 0 like 1/n.
- **Lemma 17.5 and (C½).** λ(I_F⁻¹H) at the batch θ† for single- and two-exponential models of the power-law truth (n* = 0.8 and 0.95; β ∈ {0.05, 1, 10, (0.05, 10)}): all eigenvalues in [0.970, 1.002], with one equal to 1 as predicted. So (C½) holds with huge margin in every example, including grossly misspecified ones (n(θ†) = 0.109 vs n* = 0.8 for β = 10).
- **M1/M2 (misspecified, free, R = 600).** n·Var of (μ, α_1, α_2) for OI at n = 10⁴, 3·10⁴, 10⁵: (455, 3.49, 27.6), (451, 3.56, 27.5), (456, 3.68, 28.5), flat from 10⁴. Predicted from Theorem 17B with the population J: (481, 3.47, 24.6); corr(μ, α_1) −0.799 vs −0.795. Agreement 5 %, 6 %, 14 % within MC (≈ 6 % per entry). EF is indistinguishable from OI (R = 200: 420, 2.64, 23.2 vs 423, 2.66, 23.3 at n = 10⁵), as predicted because A ≈ I.
- **M4 (rate).** n·E‖θ_n − θ†‖² = 497, 482, 501 at n = 10⁴, 3·10⁴, 10⁵ (free) and 176, 176, 175 (active): E‖e‖² ∝ 1/n. Mean errors times √n shrink to the MC noise.
- **M3 (active constraint, n ≤ 0.5, θ† = (2.47, 0.092, 0.288) with n(θ†_∞) = 0.62, R = 600).**
 - n·Var: μ 149, α_1 0.99, α_2 24.5 at n = 10⁵ vs the predicted constrained Lyapunov law (Theorem 17C(b)) 148, 0.88, 21.9 (1 %, 12 %, 12 %; MC ≈ 6 %); corr −0.108 vs −0.053. The unconstrained prediction would be 5× too large for μ (Problem 16: 18.5 vs 3.4 per 4·10⁴ events).
 - Theorem 17C(a): n·d_n → 76 and stays there (n = 10³ … 10⁵: 58, 80, 73, 74, 77), where d_n = 0.5 − n(α_n) ≥ 0. So d_n = O(1/n), √n·d_n → 0, as stated.
 - The projection is active in a positive fraction of steps that stabilises: 8 %, 4.6 %, 3.6 %, 3.2 %, 3.1 % at n = 10³, 3·10³, 10⁴, 3·10⁴, 10⁵ (recurrence at the face, not "eventually always on it").
- **M5 (EF vs OI).** No visible difference in any of the settings above (the drift matrices coincide to 3 %). I did not find an example where they differ; the theory predicts a difference only when λ(I_F⁻¹H) departs from 1, which was never observed.
- **Failure mode found (OI, tiny μ_min).** With μ_min = 10⁻³ one run in 300 was trapped at the vertex (μ_min, 0, 0) for the whole run (error −0.999 in μ), because after an early overshoot λ_i ≈ μ_min makes a_ia_iᵀ/λ_i² huge and F_n locks the iterate in the corner. With μ_min = 0.1 no trapping; with 0.3 none. EF was never trapped (0/300 at all three μ_min). Also supports the Problem 5 design choice of a box known a priori that is not close to the singular boundary.

## 6. Verdict

Convergence to θ†_Θ is established as Theorem 17A [P*] in general and **as Theorem 17A* [P] without any cited stochastic-approximation theorem under the drift condition (D_κ)** (section 7; the active constraint is handled; the proof needs only o(t) for every non-drift term, so the ergodic theorem suffices); the limit laws are derived (17B, 17C [P*]) and reproduced numerically to 5–14 % in the interior and 1–12 % with an active constraint (MC ≈ 6 % per entry). Two exact facts came out of the algebra: the Lyapunov drift of ‖θ − θ†‖²_F equals −vᵀ(2H̄ − I_F)v with the multiplier term of favourable sign (Lemma 17.2), and A = I_F⁻¹H has eigenvalue 1 at θ† (Lemma 17.5). The observed-information preconditioner attains the batch sandwich covariance with no side condition and the Fisher one does too in every example tried. Not done (section 3): self-contained proof of the Poisson-equation step, rigorous excursion bound for d_n, degenerate (no strict complementarity) limit, regret statement 17.4, K > 2.

## 7. Theorem 17A*: a self-contained proof of almost-sure convergence (no stochastic-approximation theorem cited)

**Idea.** Spreij's Lyapunov function ‖θ − θ†‖²_F needs only that the *stationary mean* of its drift is ≤ 0, and the conclusion needs only that every other term is o(t), not O(log t) and not a rate. The o(t) statements follow from the ergodic theorem alone, so no Poisson equation is needed. The price is a drift condition (D_κ) on the pair (truth, θ†), which is checked below and which can always be met by scaling the gain.

**Recursion with gain constant κ > 0.** θ_{n+1} = Π^{κF_n}_Θ(θ_n − (κF_n)⁻¹g_n), F_n = F_0 + Σ_{j≤n}ΔF_j, with ΔF_j = ∫_{gap j}φφᵀ/λ_{θ_j} (EF) or a_ja_jᵀ/λ_j² (OI). κ = 1 is the recursion of the earlier sections. Notation: v = θ − θ†, λ†_t = θ†ᵀa_t.

**Assumptions.**
- (S1) N is a stationary ergodic point process with predictable intensity λ*_t ≥ λ_lo > 0, distinct β_k, and (1 + |S|)^6 (1 + λ*) integrable (moments as in Part 1; these hold for a subcritical Hawkes truth). The online path is stationary, or started from any state that couples to the stationary path (the effect of S(0) decays like e^{−βt}, which is summable against Lipschitz functionals).
- (S2) Θ ⊂ {μ ≥ μ_min > 0, α ≥ 0} is compact and convex (e.g. Θ_ε ∩ {μ ≤ μ_hi}); θ† ∈ Θ minimises L on Θ (exists by Theorem 16), hence vᵀ∇L(θ†) ≥ 0 for all θ ∈ Θ.
- (S3) **Drift condition (D_κ).** For all θ ∈ Θ,
 Ψ̄_κ(θ) := κ E[(vᵀa)²/λ_θ] − 2vᵀ∇L(θ) ≤ 0 (EF);  κ E[λ*(vᵀa)²/λ_θ²] − 2vᵀ∇L(θ) ≤ 0 (OI).

**Lemma 17.7 (closed form of the drift; sufficient conditions) [P].** With λ_θ = θᵀa,
 Ψ̄^{EF}_κ(θ) = −2vᵀ∇L(θ†) − E[(vᵀa)²(2λ* − κλ†)/(λ†λ_θ)],
 Ψ̄^{OI}_κ(θ) = −2vᵀ∇L(θ†) − E[λ*(vᵀa)²(2λ_θ − κλ†)/(λ†λ_θ²)].
Hence (D_κ) holds if **2λ*_t ≥ κλ†_t a.s. (EF)**, respectively **2λ_θ(t) ≥ κλ†_t for all θ ∈ Θ (OI)**. In particular for the well-specified model (λ* = λ†) it holds for every κ ≤ 2, and for a misspecified truth with sup λ†/λ* < ∞ it holds for every κ ≤ 2 inf λ*/λ†: **a convergent gain always exists**.
*Proof.* vᵀa = λ_θ − λ†, so vᵀ∇L(θ) = E[(λ_θ − λ†)(λ_θ − λ*)/λ_θ] and ∇L(θ) − ∇L(θ†) = H̄v with H̄ = ∫₀¹H(θ† + sv)ds = E[λ*aaᵀ/(λ†λ_θ)] (since ∫₀¹ ds/(λ† + s vᵀa)² = 1/(λ†λ_θ)). Insert, and use 1/λ_θ = 1/λ† − (λ_θ − λ†)/(λ†λ_θ) together with E[a(λ* − λ†)/λ†] = −∇L(θ†). ∎

**Theorem 17A\* [P].** Under (S1)–(S3), for any F_0 ≻ 0 and θ_1 ∈ Θ, the recursion (EF or OI) satisfies θ_n → θ† almost surely, and ‖θ_{n+1} − θ†‖²_{κF_n} = o(n). The constraint may bind (the first-order condition is used only through vᵀ∇L(θ†) ≥ 0).

*Proof.* Put v_n = θ_n − θ†, W_n = ‖v_n‖²_{κF_{n−1}}, g_n = b_n − a_n/λ_n.
**Step 0 (pathwise).** Π^{κF_n} is the metric projection onto the convex set Θ ∋ θ†, hence non-expansive in ‖·‖_{κF_n}. With ‖v_n‖²_{κF_n} = W_n + κv_nᵀΔF_nv_n,
 W_{n+1} ≤ W_n + κv_nᵀΔF_nv_n − 2⟨g_n, v_n⟩ + g_nᵀ(κF_n)⁻¹g_n.
**Step 1 (time-integral form).** Let θ(s) = θ_n during gap n (predictable) and v(s) = θ(s) − θ†. Since b_n = ∫_gap φ ds, a_n = φ(t_n), and for EF κv_nᵀΔF_nv_n = ∫_gap κ(vᵀφ)²/λ_θ ds,
 Σ_{n≤N}[κv_nᵀΔF_nv_n − 2⟨g_n, v_n⟩] = ∫₀^{t_N}[κ(vᵀφ)²/λ_θ − 2vᵀφ]ds + 2∫₀^{t_N}(vᵀφ/λ_θ)dN_s.
With dN = λ*ds + dM this equals ∫₀^{t_N}Ψ(Z_s, θ(s))ds + 2M^{(1)}_{t_N}, where Ψ = κ(vᵀφ)²/λ_θ − 2vᵀφ(1 − λ*/λ_θ) and M^{(1)} = ∫(vᵀφ/λ_θ)dM. (OI: the jump κ(vᵀa_n)²/λ_n² is ∫κ(vᵀφ)²/λ_θ² dN = ∫κλ*(vᵀφ)²/λ_θ²ds plus a second martingale.) The stationary mean of Ψ(·, θ) is Ψ̄_κ(θ).
**Step 2 (growth of the metric).** On Θ, λ_θ ≤ c₁|a|₁. Then F_n ≽ (1/c₁)∫φφᵀ/|φ|₁ ds (EF), or Σ_ja_ja_jᵀ/(c₁|a_j|₁)² (OI). By the ergodic theorem (time, resp. Palm) these dominate c·t_n (resp. c·n) eventually, with c the smallest eigenvalue of E[aaᵀ/|a|₁] ≻ 0, resp. E[λ*aaᵀ/|a|₁²]/Λ* ≻ 0 (positive definite by the gap argument of Lemma 16.2: exponentials with distinct rates are linearly independent on a gap, and the weights are positive).
**Step 3 (martingales).** ⟨M^{(1)}⟩_t ≤ (diam Θ)²μ_min⁻²∫|φ|²λ* ds = O(t) a.s. by the ergodic theorem, so M^{(1)}_t = o(t) by the martingale strong law (e.g. Liptser–Shiryaev, Thm 2.6.10); likewise for the OI martingale (quadratic variation O(t) by (S1)).
**Step 4 (the quadratic term of the Lyapunov inequality).** g_nᵀ(κF_n)⁻¹g_n ≤ κ⁻¹|g_n|²/λ_min(F_n) ≤ C|g_n|²/n for n ≥ n₀, and |g_n|² ≤ 2w_n²|φ(t_{n−1}+)|² + 2|a_n|²/μ_min², whose Palm average is finite by (S1). Abel summation gives Σ_{n≤N}g_nᵀ(κF_n)⁻¹g_n = O(log N).
**Step 5 (fluctuation of Ψ along the adapted path).** Let f(z, θ) = Ψ(z, θ) − Ψ̄_κ(θ). It is Lipschitz in θ on Θ with a constant K(z), E K < ∞ (|∂_θ(1/λ_θ)| ≤ |a|/μ_min²). For L > 0 put e_L(s) = sup_θ|L⁻¹∫_s^{s+L}f(Z_u, θ)du|. (i) E e_L(0) → 0 as L → ∞: by the L¹ ergodic theorem at each point of a finite ε-net of Θ and the Lipschitz bound, E e_L(0) ≤ o_L(1) + εE K. (ii) Fix L. Work on dyadic blocks of the event index, n ∈ [2^m, 2^{m+1}). On a block cut the time axis into windows of length L with an offset u_m ∈ [0, L) (chosen below). The iterate moves slowly: within a window whose events have index ≥ 2^m, Σ|θ_{j+1} − θ_j| ≤ Σ|(κF_j)⁻¹g_j| ≤ (C/2^m)Σ_{window}|g_j|. Replacing θ(s) by its value at the window start inside f therefore costs at most (C/2^m)(∫_{window}K)(Σ_{window}|g_j|) ≤ (C/2^{m+1})(A_k² + B_k²), with A_k = ∫_{window k}K and B_k = Σ_{window k}|g_j|, and with θ frozen a window contributes at most L·e_L(window start). Both bounds are sums of non-negative window functionals of the stationary process, so averaging over the offset u (pigeonhole) gives an offset u_m for which the total cost on the block is at most L⁻¹∫_{block}[L e_L(s) + C 2^{−m}(A_s² + B_s²)]ds, with A_s, B_s the same functionals for the window starting at s. The second part is O(1) per block (A_s, B_s have finite second moments, and the block has length O(2^m)), hence O(log N) in total, and the first is at most ∫₀ᵗe_L(s)ds. Therefore
 limsup_t t⁻¹|∫₀ᵗ f(Z_s, θ(s))ds| ≤ limsup_t t⁻¹∫₀ᵗ e_L(s)ds = E e_L(0)
by the continuous-time ergodic theorem. Letting L → ∞ gives ∫₀ᵗ f(Z_s, θ(s))ds = o(t) a.s.
**Step 6 (conclusion).** Summing Step 0 and inserting Steps 1, 3, 4, 5: W_{N+1} ≤ W_1 + ∫₀^{t_N}Ψ̄_κ(θ(s))ds + o(N) ≤ W_1 + o(N), because Ψ̄_κ ≤ 0 on Θ by (D_κ) and t_N = O(N). By Step 2, W_{N+1} ≥ κcN|v_{N+1}|², hence |v_{N+1}|² ≤ o(N)/(κcN) → 0. ∎

**What the proof does and does not use.** It uses: the ergodic theorem (time and Palm), compactness of Θ, the martingale strong law, non-expansiveness of the F-norm projection, and (S3). It does *not* use a Poisson equation, mixing rates, the Markov property of the truth, or any ODE/stochastic-approximation theorem. It gives no rate and requires (S3). For κ = 1 and the well-specified model it recovers Spreij's a.s. consistency for the Hawkes dictionary with the F-norm convex projection in place of his coordinate safeguard, and for any stationary ergodic truth. The general result without (S3) is Theorem 17A [P*] (ODE method).

### Numerical checks of (D_κ) and of the κ-scaled recursion (power-law truth n* = 0.8; `code/p17/t_drift.py`, `t_kappa.py`; `online17.c` now takes a kappa argument)
- *Drift condition on Θ.* 400 000 events, 500 random θ ∈ Θ plus the vertices (μ ∈ {0.3, 6}, α supported on one coordinate at the stability boundary); maximum over the sample of Ψ̄_κ(θ)/|θ − θ†|²:

 | case | EF, κ = 0.25 / 0.5 / 1 / 2 | OI, κ = 0.25 / 1 |
 |---|---|---|
 | K = 1 free (μ†/μ* = 2.70) | −0.048 / −0.041 / −0.027 / **+0.015** (99 % of θ violate) | **+34 / +215** |
 | K = 2 free (μ†/μ* = 1.89) | −0.028 / −0.024 / −0.016 / **+0.020** (58 %) | **+80 / +503** |
 | K = 2, n ≤ 0.5 active | −0.039 / −0.035 / −0.026 / **+0.011** (14 %) | **+40 / +253** |

 So (D_κ) holds for EF at κ ≤ 1 in all three cases, including the strongly misspecified K = 1 case and the active-constraint case, and fails at κ = 2. For OI it fails at extreme θ (μ near μ_min and α ≈ 0, where λ_θ ≪ λ†): for OI the proof needs a ratio-restricted Θ (λ_θ ≥ κλ†/2), the box design of Problem 5. This restricts the theorem, not the algorithm in practice: all earlier OI runs with μ_min = 0.3 converged. The sample is finite (not a proof that (D_κ) holds on all of Θ).
- *The pointwise sufficient condition is not necessary.* λ†/λ* at the events: maximum 1.35 (K = 1), 1.12 (K = 2), 1.31 (active); in the quiet limit μ†/μ* = 2.70, 1.89, 2.50, so 2λ* ≥ λ† fails (K = 1, active) while the integrated condition (D_1) holds.
- *Convergence of the κ-scaled recursion* (40 replicates, 10⁵ events, start far from θ†): median |θ_n − θ†| at n = 3·10², 10³, 3·10³, 10⁴, 3·10⁴, 10⁵. EF, K = 1: κ = 1: 0.58, 0.22, 0.25, 0.13, 0.065, 0.035; κ = 0.5: 0.41, 0.36, 0.25, 0.12, 0.084, 0.046; κ = 2: 1.18, 0.70, 0.42, 0.24, 0.15, 0.10; κ = 4: 1.65, 1.26, 0.97, 0.73, 0.57, 0.41. K = 2 and OI behave the same (κ = 1 at n = 10⁵: 0.038 and 0.043). **Every κ converges, including κ = 2 and 4 where (D_κ) fails**: (D_κ) is sufficient, not necessary, and the restriction on κ is an artefact of the Lyapunov function. The slowdown for large κ matches the deterministic prediction: the error decays like n^{−λ_min(A)/κ} with λ_min(A) ≈ 1 (κ = 4: slope −0.25 between n = 3·10³ and 10⁵).
