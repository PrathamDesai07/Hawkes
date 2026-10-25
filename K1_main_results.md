# The K = 1 story: consolidated statement, proof status, and gaps

Model: λ_t = μ + αS_t, dS = −βS dt + dN, ε = 1 − α/β ∈ (0,1), γ = βε, k = 2βμ/α², stationary mean ES = μ/γ. d = 2 parameters (μ, α), β known.
Each result points to the note holding the proof. Tags: **[P]** proved on paper, **[P*]** proved modulo the inputs in the last column, **[N]** numerical only, **[O]** open.
Date: 25 Oct 2026. Nothing here has been refereed; the proofs were written in one pass each and re-derived only numerically.

## A. Information

| # | Statement | Where | Tag |
|---|---|---|---|
| A1 | Exact stationary Fisher matrix: I_μμ = m = E[1/λ], I_μα = (1−μm)/α, I_αα = (Eλ − 2μ + μ²m)/α², det I = (mEλ − 1)/α². m = ∫φ with φ' = −uβμφ/(εβu + h(αu)), h(x) = e^{−x}−1+x. | `problem10_fisher_k1.md`, `problem10_rigorous.md` | [P] |
| A2 | Dichotomy at k₀ = 2μ/β = 1: E[1/λ] ≈ 2ε/(2μ−β) if k₀ > 1; (2/β) ε ln(1/ε) if k₀ = 1; C ε^{k₀} if k₀ < 1. So the baseline direction is learnable at rate ≈ ε per unit time when k₀ > 1 and slower otherwise. | `problem11_dichotomy.md` | [P] (Theorem A repaired in `theoremA_step2_and_stationary_tv.md`) |
| A3 | Matches the CIR information matrix (Ben Alaya–Kebaier–Tran Thm 2.1) as ε → 0 (relative difference ≈ ε), and the critical coefficient 2/(β(2μ−β)). | `literature_search_report.md`, `code/lit_check_cir_fisher.py` | [N] (limit identification), classical CIR side from literature |

## B. Dynamics and learning

| # | Statement | Where | Tag |
|---|---|---|---|
| B1 | Moments and contraction, mixing rate γ = βε; exponential moments of S up to θ_* ≈ 2ε/(1−ε) | `problem03_moments_gap.md`, E0 in `problem05_E5_expected_regret.md` | [P] |
| B2 | Concentration of the accumulated information (two-interval bound, Poisson equation, Freedman): failure ≤ 3e^{−cγT} | `problem05_fisher_concentration.md`, E3 in `problem05_exponential_tails.md` | [P*] Assumption M (numerical), integrability |
| B3 | Exact F-norm projection onto a box (2-D, O(1) work); non-expansive, so Newton regret telescoping survives projection; the diagonal clip does not | `problem04_projection.md`, `problem05_circularity.md` | [P] |
| B4 | Prop. 1: E⟨g_i, v_i⟩ = E vᵀΔF_iv_i, so E R_n ≤ P_0 + E Σ gᵀF⁻¹g without exp-concavity | `problem05_online_newton.md` | [P] |
| B5 | Leverage lemmas 7a–7c and conditional Theorem 7 (with the corrected truncated 7c, see below) | `problem07_leverage.md`, §0 of `problem05_E5_expected_regret.md` | [P] |
| B6 | Two-stage bootstrap: iterates enter R(δ) for t ≥ T_1(δ, q) = O((ln(1/q) + d ln T)/(γδ²)) | `problem05_bootstrap.md`, `problem05_exponential_tails.md` | [P*] B2, box |
| B7 | **E5 (new):** exponential tail for G_n: G_n ≤ (1+o(1)) log det + o(1)·x + B_X X with probability ≥ 1 − 4e^{−x} | `problem05_E5_expected_regret.md` §2 | [P] (given leverage conditions, which hold w.h.p. by B1, B2) |
| B8 | **Theorem R (new):** E R_n ≤ (½ + O(δ+η)) d ln T_n + C, C independent of n; limsup E R_n/ln T_n ≤ d/2 = 1 | `problem05_E5_expected_regret.md` §4 | [P*] B2, Assumption M, known box, unwritten optional-stopping details |

Observed: slope of E R_n against ln T_n ≈ 0.9 (ε = 0.3, 0.1, 0.03), against d/2 = 1 (`code/p5_bootstrap_check.py`, item F).

## C. Limits

| # | Statement | Where | Tag |
|---|---|---|---|
| C1 | Stationary-record impossibility: with only stationary data of length T, the two-point/Le Cam argument gives a lower bound of order 1/(uᵀI_Tu) in the weak (ridge) direction; the reduction needs TV(π_θ, π_θ′) = O(δ) | `theoremA_step2_and_stationary_tv.md` | [P] for K = 1 (characteristic-function TV bound) |
| C2 | The same for any K via ghost immigrants and the cluster-survival flow (Theorem V) | `general_K_impossibility.md` | [P] |
| C3 | Expected-regret lower bound (d/2) ln T for every algorithm, with explicit ε-dependence | none | [O] (no Hawkes-specific bound found in the literature either; classical Rissanen-type bounds are for fixed θ) |

## D. What changed in this session

1. E5 proved (Lemmas E5a–E5d, Proposition E5), with the pathwise inequalities checked in 1.4·10⁶ gaps and no violation.
2. A flaw in Theorem 7 / Lemma 7c was found and fixed: the assumption "gap leverage x_i ≤ x̄ almost surely" cannot hold (it depends on the Exp(1) variable E_i); truncation at E_i ≤ m costs (m+1)e^{−m} per gap.
3. The "expected regret needs an episode argument" gap was closed by charging bad events to their own time t rather than to the horizon: the cost of a failure at time t is ≲ (μ/ε)t and its probability e^{−cγt}, so the total is a constant independent of n.

## E. Honest list of what the K = 1 theorem still rests on

1. Assumption M: numerical.
2. A box containing θ*, known a priori (ratio ≤ 2 in μ, α ≥ α*/2): this is an assumption on the problem, not a limitation of the proof method alone.
3. Constants (rate c in B2) are ≈ 10⁵ worse than observed; the leading constant d/2 is unaffected.
4. Optional-stopping across random restart indices (B8, Step 3): standard, not written.
5. C (the additive constant) grows as ε → 0; its ε-dependence is not computed.
6. Lower bounds for algorithms (C3) and the k₀ ≤ 1 regime for the regret theorem: B2's rate and Assumption M were only tested for k₀ > 1 values in the main runs; k₀ ≤ 1 is not covered by the numerics.
7. K ≥ 2: projection (QP), the information spectrum (three tiers), TV for the stationary law (k_K = 1 borderline), and Theorem F are open for the regret side.

---
## Update (25 Oct 2026, later): dependencies of Theorem R
See `K1_completion.md`. E.1 (Assumption M) closed for ε ≤ ε_M via Theorem T; E.4 (optional stopping) closed (predictable indicators, no restarts); integrability in B2 closed; C(ε) = O(μ d ln(1/(εδ))/(βε²δ²)). Remaining: known box (E.2), constants (E.3), k₀ ≤ 1 (E.6), K ≥ 2 (E.7), and the fact that Theorem T is an outline with numerically verified estimates.
