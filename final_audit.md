# Final audit of everything marked Done (Problems 1, 2, 3, 4, 5, 7, 10, 11, 12, 14), 25 Oct 2026

Method: for every Done item I re-derived the proof steps that carry the claim (not just re-ran the authors' scripts), searched for places where a statement is stronger than its proof, and ran fresh numerical tests of the inequalities that the proofs use. Complements `audit_done_problems.md` (another session; Problems 1, 3, 10, 11, 12, 14).
**This is not a proof of perfection and cannot be one:** nothing has been refereed or machine-checked, and several Done items are Done in the sense of the tracker's "Done when" line (a draft with stated scope), not in the sense of a finished paper. What the audit does give is a list of defects that were found and repaired, and a list of what remains weaker than it sounds.

## Defects found in this audit and repaired

1. **Overstated impossibility (Problem 11, Theorem B / Corollaries D, E′; my LaTeX; verification note).** The notes say: "given the pre-sample state S₀, no estimator of μ is consistent when k_K > 1." The KL lemma compares path laws that start from the *same* initial state, and the Fisher bound I_T = T·I(θ) is for a stationary start. But if S₀ is drawn from π_θ (stationary experiment), its law depends on θ, so the correct two-point bound is TV(Q_θ, Q_θ′) ≤ TV(π_θ, π_θ′) + √(ρT uᵀIu/2); if S₀ is instead treated as ancillary (θ-free law, e.g. fixed), the path information is not T·I (from S₀ = 0 it diverges like ln T). So the "S₀ given" impossibility is **not** implied by the lemma. What *is* proved: K = 1, stationary record: Theorem T (TV(π_θ, π_θ′) = O(δ)); any K, stationary record with unobserved past: Theorem V. Statuses corrected in the LaTeX (Theorem B and Corollary E′ restated; review table updated). The headline "no consistent estimator of μ when k_K > 1 for every K" remains true **for the stationary record with unobserved past** (Theorem V); for K ≥ 2 with S₀ *observed* it is open (needs TV(π_θ, π_θ′) for K ≥ 2).
2. **Wrong constant in Lemma E5c (mine).** I bounded tr(F⁻¹ΔF) by ℓ̄c⁺E_i/μ*; since λ_θ ≥ c⁻λ* the correct bound is ℓ̄E_i/(c⁻μ*) (and for ΔJ the bound ℓ̄E_i/μ* was already right). Consequence: x̄ = ℓ̄m/(c⁻μ*). No effect on the leading constant (x̄ → 0), but the statement was false as written whenever 1/c⁻ > c⁺. Fixed in `problem05_E5_expected_regret.md` and the LaTeX. Not detectable by simulation (the K = 2 pathwise audit, `code/p5_k2_pathwise_check.py`, has 0 violations for both versions because the iterates never reach the box corners).
3. **Loose general-K Lemma E1 (mine).** `K2_regret_chain.md` used Z_ω = ωᵀS with inequalities losing a factor B·ν̃ ≥ 1 in the upper tail of Λ_T. The exact identity εΛ_T = μT + nM′_T − (Z_c(T) − Z_c(0)), Z_c = Σ(α_k/β_k)S_k, holds for every K (integrated key identity of Problem 3) and gives Lemma E2 verbatim. Replaced.

(Corrections already recorded earlier: thinning-loop bug, Problem 2 range artefact, Theorem A step-2 justifications, Hall–Willett norm vs squared norm, peek-ahead sampling and retracted c_lag, V²/(2ηT) factor, diagonal clip, (L2) "a.s." in Theorem 7.)

## Independent re-checks that passed in this audit

- **Newton/martingale chain for K = 2 with a polytope F-projection** (`code/p5_k2_pathwise_check.py`, 2 ε values, 3 replicates × 4 000 events, projection active in 1–4 % of steps): per-step Lemma C1 inequality D_i ≤ D_{i−1} − X_i − m_i + G_i holds to 4·10⁻¹³; R_i ≤ X_i holds pathwise; Lemma E5a bound holds; the corrected trace bound holds; the mean of ν_i is 0 within 1.7σ. (So the structural K ≥ 2 extension of the regret chain is numerically consistent.)
- **Dense F-projection**: matches brute-force enumeration to 5·10⁻¹⁴, non-expansive on 2 000 pairs.
- **Closed-form Problem 2 pieces**: C ⪰ (Eλ/4)G, λ_min(M) ≥ λ_min(C)/(1+2|m|²), Cauchy-inverse formula — all hold in 18/18 cases.
- Re-derivations by hand (no new numerics needed): Proposition 1 (a)–(c); Lemma C1 algebra (m_i simplification); Lemma E5d identity and the bound w_A ≤ κ₁²c⁺; Lemma E4; Lemma E5b/E5c Ville steps; the tanh jump law; the noise/efficient-information formula; Lemma 1′ recursion; Lemma 4 (exponential moment); Theorem 2 (contraction); Theorem A's three cases (including that the k₀ = 1 tail on [δ, ∞) is O(ε) because of the decay of F_ε, not of the Gamma tail); Lemma F2 (positive part trick); the lower-bound proposition (data processing, Gaussian maximum entropy).

## What remains weaker than "Done" suggests (honest list)

| Problem | Done in the sense of | Not established / weaker |
|---|---|---|
| 1 | Lemmas 1–3 and the K = 1 check | stationarity assumed; α* > 0 needed |
| 2 | explicit bound in K, r exists | closed-form chain is 10³–10⁸ below the truth; sharp scaling only numerical (three tiers) |
| 3 | explicit ε-dependence of moments and rate | constants carry β_max/β_min; "L² gap = γ" not claimed; sup over horizon of λ_t not bounded |
| 4 | O(K log K) diagonal algorithm + dense exact O(K³) one, both tested | dense case is not O(K log K); diagonal clip invalid for the regret recursion |
| 5 | Theorem HP (high probability) and Theorem R (expectation, sharp constant) | conditional on: known compact convex Θ_c; Assumption M (certified for K = 1 by quadrature, not interval arithmetic) or λ_L (K ≥ 2); Theorem R's proof is at the level of a detailed outline (Step 3 integrals not written line by line); constants 10⁵–10⁶ conservative; K ≥ 2 numerics only to 4·10⁴ events |
| 7 | Lemma JD (joint zero mean, covariance, tail, cumulative bound) | dimension-free proof; numerics K = 1, 2 only |
| 10 | Theorem A proved; transient formula validated | cold-start coefficients conjectural; K ≥ 2 stationary asymptotics via Problem 11 |
| 11 | K = 1 full dichotomy; any K: Theorems E, F | no matching upper bound (estimator rate) for k_K < 1; k_K = 1 open for K ≥ 2; efficient information for K ≥ 2, k_K < 1 not proved; **impossibility with S₀ observed for K ≥ 2 open** (item 1 above) |
| 12 | identification of the three tiers; ridge and top bounds proved | middle tier is a heuristic (spectral approximation, numerically within 2–15 %); eigenvector convergence only numerical |
| 14 | a decision (model class) + numerical study of a linear filter and the estimators | no theorem about the estimators' tracking error; amplitude effect not derived; K ≥ 4, ε < 0.4 (K ≥ 2) not run |

## Bottom line
After this audit I found three defects in the written results (listed above) and fixed them; I found no further error in the steps I re-derived. I cannot certify that no error remains: the proofs have been checked by their author and by re-derivation/numerics, not by an independent mathematician or a proof assistant. The items most worth independent human scrutiny are Theorem V (cluster/superposition argument), Theorem R's Step 3, and the claim list for Problem 12.


---
## Addendum (25 Oct 2026): strengthening round

Items of the "weaker than it sounds" table that were improved (details in `tangent_bound.md` and `theoremR_proof_complete.md`):
- **Problem 2:** rigorous lower bound within 1.3-3.5x of the Monte Carlo lambda_min (K = 1, 2, 3), correct eps-scaling, from the tangent minorant with exact moment systems; the 10^3-10^8 looseness is gone.
- **Problem 12:** every eigenvalue of I is bounded below rigorously (0.23-0.85 of Monte Carlo): the three tiers are now a theorem as regards lower bounds. Still open: matching upper bounds and eigenvector convergence.
- **Problem 5:** (a) every step of Theorem R's expected-regret proof is written (only the per-time failure bound H3 remains an input); (b) Assumption M is no longer needed: the information level comes from the tangent minorant (about 10^4 better than the truncation device for the K = 2 example).
- **Not strengthened (honestly):** Problem 11 (a matching estimator rate for k_K < 1 needs LAMN/weak-convergence arguments; the observed information is random at O(1) relative size when eps = c/T); Problem 14 (a tracking theorem for the estimators is the content of the not-started Problem 13; the amplitude effect of the data-driven preconditioner is not derived); the multiprecision evaluations are not interval arithmetic.
