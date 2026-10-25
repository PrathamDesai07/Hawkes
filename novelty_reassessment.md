# Novelty reassessment of C3 and C4 given Herrera-Marín (25 Oct 2026)

**Scope and limits.** Compared against: Herrera-Marín I (arXiv 2609.38268; Sections 1, 2, 4 read closely) and II (arXiv 2609.38267; abstract, introduction, Sections 5–6 read), and the project's own notes `problem10_rigorous.md`, `problem11_dichotomy.md`, `problem12_degenerate_directions.md` (written by another session; I read their statements, tables and structure, **not their proofs** except where stated). Not a systematic literature search: Solo, Ben Alaya et al., Horst–Xu, the four "to check" rows, and the discrete-time papers Herrera-Marín cites (Wei–Winnicki; Ispány et al.) were not read. "Not in the preprint" below never means "not known".

## 1. What Herrera-Marín covers (heavy-tailed regime; kernel shape known)

- **Thm 2.2 + Cor 2.3** (any normalised shape, empty start, shared μ): KL(P₀ᵀ‖P₁ᵀ) ≤ δ²(1+c)I_T/(ρ₀ρ₁), I_T = Tμ(1−ρ₀); hence the relative margin is not uniformly locally consistent when I_T is bounded. *I re-derived the proof and tested it (12/12 configurations, bound 2.2–6.9× above the Monte Carlo KL).*
- **Thm 4.1**: Fisher information of (log μ, log(1−ρ)) in terms of A_T = E[λ̄/λ], D_T = E[λ/λ̄]. *Re-derived; matches my K=1 entries exactly; score-covariance test z = +0.39, +0.71, +1.18.*
- **Introduction**: the margin information is also bounded in the light-tailed Jaisson–Rosenbaum regime with 1−ρ_T ≍ T⁻¹.
- **Remark 4.8**: in the light-tailed regime the baseline information diverges when the CIR shape κ ≤ 1 (Fatou); the converse for κ > 1 "would require uniform integrability of 1/V^T, which we do not prove". His ref. [22] (CIR drift estimation) assumes sup_t E X⁻¹⁻ᵉ < ∞.
- Everything else in I is about the rough heavy-tailed limit (MLE at a random rate; μ is *learned* there). II treats fixed ρ and long memory (minimax T^{−γ/(2γ+1)}); its Section 5 reports numerically, for fits of K exponentials to a long-memory kernel, an "almost flat direction" from K ≈ 5. Neither paper treats streaming estimators, regret, or individual α_k of a well-specified dictionary.

## 2. What the project has beyond that, and how much of it I have checked

| Claim (where) | Status of the claim | My verification |
|---|---|---|
| K=1, stationary: m_ε = E[1/λ] ∼ 2ε/(2μ−β) (k₀>1), ∼ (2/β)ε ln(1/ε) (k₀=1), ∼ Cε^{k₀} (k₀<1) (P10 Thm A) | proof written, "needs re-derivation" | **re-derived step by step and tested (`verification_thmA_E_F.md`): correct; Step 2 has two inadequate justifications (both statements true, repairable); limits confirmed numerically for k₀ = 2, 0.6, 1** |
| Any K, k_K = 2μΣα_k/β_k² > 1: E[1/λ] = O(ε) (P11 Thm E); k_K < 1: E[1/λ] ≥ cε^{k_K+η} (Thm F); borderline k_K = 1 open for K ≥ 2 | proofs written | my E3 (K=2, β=(1,4)): A_T equals the K=1 value at equal k (1.29/1.53/1.79 vs 1.26/1.55/1.78); k=0.5 grows ×2 per decade. **Lemma H, Theorems E and F re-derived: correct; Lemma H holds pathwise in 30/30 cases; E and F integral inequalities and E[A f]=0 confirmed by stationary MC (`verification_thmA_E_F.md`)** |
| Impossibility: no consistent estimator of μ when k_K > 1 (stationary), by Le Cam along the ridge u = (δ, −εαδ/μ) using KL ≤ ρ·uᵀI_T u (P11 Thm B, Cor E′) | proof written | **I tested the KL lemma independently: KL ≤ ρ·uᵀI_Tu in 8/8 configurations** (k₀ = 2, 0.6; ridge and μ-only; δ = 0.1, 0.3), bound 1.8–3.5× above the MC KL. **The lemma assumes the same initial state; for the stationary record with unobserved past one needs TV(π_θ,π_θ′) = O(δ), now shown for K=1, k₀>1 (`theoremA_step2_and_stationary_tv.md`; numerically 0.026/0.070 for δ=0.1/0.3, flat in ε); K ≥ 2 open.** |
| Three-tier Fisher spectrum: top ∼1/ε along ∇ₐn ∝ (0, 1/β_k); K−1 "middle" shape modes of size O(1) set by the grid; ridge ∼ ε/(μ(k−1)) along (1, −εα/μ) (P12) | ridge Rayleigh bound rigorous; middle tier heuristic | consistent with my Problem 2 data once the crossover is accounted for (see §5) |
| Observation model: stationary ⇒ bounded/log/polynomial; empty start ⇒ μ-information ≥ ln(1+β_maxT)/(μβ_max) for every k₀ (P10 §4) | short rigorous argument (Jensen + E λ_t ≤ μ(1+β_maxt)) | argument read and plausible; not re-derived |

## 3. Verdict on the original contributions

| Contribution as written in the Topic sheet | Verdict | Reason |
|---|---|---|
| **C3**: "two-point argument, KL ≲ δ²Tμ(1−n), which directions are unlearnable" | **Not novel as stated.** | Herrera-Marín Thm 2.2 + Cor 2.3, which already covers sums of exponentials. |
| **C4 divergence half** (κ ≤ 1) | **Not novel.** | Remark 4.8. |
| **C4 bounded half** (κ > 1) | **Open in the preprint; claimed in P10–P11** (K=1 Thm A; any K Thm E) | unverified by me (proofs). |
| **C4 rates / critical case / sums of exponentials / observation-model dependence** | **Not in the preprint; claimed in P10–P11** | as above. |
| **Impossibility of consistent μ-estimation for k_K > 1 (light-tailed) along the μ–n ridge** | **Not in the preprint; claimed in P11; proved for the experiment where S₀ is given; **for the stationary record with unobserved past: K=1, k₀ > 1 via a TV bound on the stationary laws (`theoremA_step2_and_stationary_tv.md`); any K, k_K > 1 via a different route (run θ′-dynamics from π_θ for L = K₀/γ, couple with ghost immigrants, cluster-survival bound q(s) ≤ C′γe^{−γs}), `general_K_impossibility.md`, Theorem V, TV ≤ Cδ√(ln 1/δ). Both are my proofs, unreviewed; TV(π_θ,π_θ′) itself for K ≥ 2 is still open** | His heavy-tailed result is the opposite (μ is *learned*); the light/heavy contrast for the exogenous rate is a clean, citable statement. |
| **C2** ((d/2) log T + C log(1/ε) + burn-in) | **No overlap found.** | Neither preprint concerns streaming estimators or regret. |

## 4. Proposed repositioning

**What to demote:** the Theorem-2.2-type margin bound and the κ ≤ 1 divergence become cited lemmas. Do not present either as a contribution, and avoid "first lower bound near criticality" or "first information dichotomy" wording.

**What to promote (headline):** *the light-tailed counterpart of Herrera-Marín's asymmetric learnability, and its consequence for streaming estimation.*
1. **C3′ — the unlearnable direction is the μ–n ridge, not the margin.** In the stationary light-tailed regime ε ≍ 1/T, for k_K = 2μΣα_k/β_k² > 1 the Fisher matrix has a ridge ∝ (1, −εα/μ) of size ε/(μ(k−1)); information along it stays bounded, so μ is not consistently estimable (opposite of the heavy-tailed case, where μ is learned and the margin is not). Needs: verification of P11 Theorems E/B and the P12 ridge bound (the KL lemma is verified numerically).
2. **C4′ — the Feller dichotomy with rates and for sums of exponentials**, in both observation models (stationary: bounded / log / polynomial; empty start: log vs polynomial), closing Remark 4.8's open converse. Needs: verification of P10 Thm A, P11 Thm E/F; open: k_K = 1 for K ≥ 2, matching estimator-side upper bound (P11 §4c).
3. **C2′ — streaming estimation with the ε-dependence made explicit.** I proposed the conjecture regret = (d/2) log T + ½ log det I(θ*;ε) + O(1), with ½ log det I ≈ ((1−k)⁺/2) log(1/ε) (E1 confirms that information side exactly). **The high-power experiment (E2, §5/§8) does not support it for this algorithm**: with a start at the truth, E[R − log T] is flat in ε (−3.44 ± 0.10 at ε = 0.3, −3.19 ± 0.32 at ε = 0.01, k₀ = 0.2) while ½ log det I rises by 1.30. So **drop the conjecture from the headline.** What remains for C2′ is the proven part (K=1: the reduction, oracle-information concentration, bootstrap; conservative constants) and the observation that, in time units, the measured regret has no visible log(1/ε) term for either k₀ = 2 or k₀ = 0.2 *in the start-at-truth design*. A lower bound with explicit ε-dependence is the real open question, and E2 cannot address it (see §5).
4. **C3″ (optional) — identifiability map for the K-dictionary.** P12 gives the three tiers, the ridge Rayleigh bound (an *upper* bound λ_min ≲ ε/(μ(k−1))), and heuristic middle-tier values; my Theorem T gives a rigorous *lower* bound. A rigorous upper bound on the middle-tier eigenvalues as K or the grid compression grows would turn Herrera-Marín II's numerical "flat direction" (K ≈ 5, misspecified setting) into a theorem for the well-specified dictionary. Not proved by anyone here.

**Safer framing sentence:** "We complete the light-tailed exponential case that Herrera-Marín (2026) left open, show that, in contrast with the rough case, the exogenous rate is not identifiable along an explicit μ–n ridge for k > 1, extend the threshold to sums of exponentials, and give a streaming estimator whose regret, for K=1, grows like (d/2) log T with a burn-in governed by the same near-critical time scale." Use only after steps 1–2 of §7 are done. Do **not** add a claim that the regret constant equals ½ log det I.

## 5. Evidence from my own experiments

**E1 — exact log(1/ε) coefficient of ½ log det I (K=1, stationary, β = 1, μ = k₀/2; `code/reposition_tests.py`).** Local slope of ½ log det I in log(1/ε) over ε = 0.3 … 0.001:
- k₀ = 0.2: 0.347, 0.397, 0.404, 0.403, 0.402 (→ (1−k₀)/2 = 0.40, essentially reached);
- k₀ = 0.6: 0.267, 0.274, 0.255, 0.236, 0.223 (→ 0.20, slowly);
- k₀ = 1: 0.207, 0.188, 0.151, 0.120, 0.096 (→ 0 like log log);
- k₀ = 2, 4: slopes → 0 and A_T = E[1/λ]E[λ] → 1.98, 1.33 (limits k₀/(k₀−1) = 2, 4/3).

**E3 — sums of exponentials** (`reposition_e3_se.py`, K=2, β = (1,4), w = (½,½), mean lifetime m = 0.625; 10 reps ± SE): A_T at ε = 0.3, 0.1, 0.03: k = 2: 1.29±0.02, 1.53±0.04, 1.79±0.05 (K=1 at k=2: 1.26, 1.55, 1.78); k = 0.5: 1.78±0.04, 3.56±0.05, 7.07±0.32 (slope 0.60±0.02 vs asymptotic 0.50). Supports the threshold 2μ·Σw_k/β_k = 1 and dependence on the dictionary only through m. One dictionary, small ε range; not a test of uniform integrability.

**KL lemma** (`code/lit_check_kl_lemma.py`): see §2; 8/8 pass.

**E2 — does the regret constant follow ½ log det I?** (K=1, exact F-norm projection, start at θ*, Tγ = 150; `code/reposition_e2.py`.)
- 16 reps, k₀ = 0.2: E[R − log T] = −3.47±0.19, −3.61±0.35, −2.94±0.75, −2.74±1.06 at ε = 0.3, 0.1, 0.03, 0.01, against ½ log det I = 0.51, 0.89, 1.37, 1.81. The difference is flat within error, but "no log(1/ε) term" is equally consistent: **low power, no conclusion.**
- k₀ = 2 control: −2.57±0.38, −2.22±0.67, −3.55±0.67; no discrimination.
- High-power run (100 reps, k₀ = 0.2, ε = 0.3 and 0.01): E[R − log T] = −3.44 ± 0.10 and −3.19 ± 0.32; difference +0.25 ± 0.34 against the +1.30 that tracking ½ log det I would require (≈ 3.1σ below), and consistent with a flat constant (0.7σ). **The conjecture is disfavoured for this algorithm and design.**
- Why this does not settle the lower-bound question: the estimator starts at θ*, so the test is superefficient at that point; a Bayes/redundancy bound is an average over a prior, and the ½ log det I term in such bounds is the prior-averaged cost. A fair test would average over random starts, but then the prior-cost term ½‖θ₁ − θ*‖²_{F₀} (which scales like 1/ε for a Fisher-scaled F₀, Problem 5) enters and the comparison depends on the prior design. I did not run that.

## 6. Correction to my Problem 2 note

In `problem02_identifiability.md` I reported, for K ≥ 2, that "λ_min(I) does **not** shrink with ε in the range tested" and that "cond ∼ 1/ε", and called the factor ε/μ in Lemma 3 "pessimistic for K ≥ 2". The ε range I used stopped at 0.05. Problem 12 shows why that was misleading: for β = (1,2) the middle-tier (shape-mixture) eigenvalue is about 0.067 and the ridge is about ε/(μ(k−1)); the ridge becomes the smallest only below ε* ≈ 0.067·μ(k−1) (≈ 0.03 in my setting). My data were therefore measuring the middle tier. **For ε → 0, λ_min(I) ∝ ε for every K** (P12; confirmed there by Monte Carlo down to ε = 0.01: 0.0053 at ε = 0.01 for K=2), so Lemma 3's ε/μ scaling is sharp in ε as ε → 0 and my "conservative in ε" remark applies only at moderate ε. (Erratum appended to `problem02_identifiability.md`.)

## 7. What would settle this (ordered by risk to the plan)

1. ~~Verify the load-bearing proofs~~ **Done 25 Oct for Theorems A, E, F (`verification_thmA_E_F.md`): correct, with two repairable justifications in A Step 2 and a gap in the *use* of E (initial-state law, Corollary E′).** Still unverified: the P12 ridge bound's use, P11 §3–4c, Problem 12's middle tier.
2. **E2 at high power** (pending) and then, if it supports the conjecture, a proof sketch of regret = (d/2) log T + ½ log det I + O(1).
3. **Targeted literature search** (not yet done): (a) Ben Alaya–Kebaier–Tran and other CIR-drift estimation work, to see whether the dichotomy in the *limit experiment* is classical and which of C4′'s content is therefore only the prelimit transfer; (b) the discrete analogues (Wei–Winnicki; Ispány et al.); (c) regret / redundancy lower bounds for point processes; (d) Solo (#6) and the online-EM rows for streaming prior art.
4. **Open mathematical pieces**: K ≥ 2 borderline k_K = 1, estimator-side matching rate for k < 1 (P11 §4c), middle-tier upper bounds (C3″).
5. **Overlap risk is live**: two unrefereed single-author preprints posted the same day; a revision could add the light-tailed converse. The parts neither preprint touches (online regret, K-dictionary, the μ–n ridge in the light-tailed regime) are where to move first.

## 8. E2 high-power result (final)

k₀ = 0.2 (μ = 0.1, β = 1), 100 replicates, Tγ = 150, start at θ*, ρ₀ = 1, exact F-norm projection:

| ε | events per run | E[R − log T] | ½ log det I | difference |
|---|---|---|---|---|
| 0.30 | 166 | −3.44 ± 0.10 | 0.51 | −3.95 |
| 0.01 | 150 000 | −3.19 ± 0.32 | 1.81 | −5.01 |

The regret constant moved by +0.25 ± 0.34 while ½ log det I moved by +1.30. Output file: `$TEMP/e2_highpower.out` (not kept in the project).

## Literature search update (25 Oct 2026)

Full details and access levels: `literature_search_report.md`.

- The K=1 information structure (bounded / log / divergent) is the Hawkes prelimit of classical CIR theory (Ben Alaya-Kebaier-Tran, fixed regimes). My exact Fisher matrix converges to their Thm 2.1 matrix as eps -> 0 (`code/lit_check_cir_fisher.py`). No local-to-critical CIR result was found.
- Discrete analogues found are all CLSE/OLS based (Li 2011; Anne-Lu-Yu-Zhou 2026, overlap Medium). None is Hawkes, none gives a likelihood information bound.
- No Hawkes-specific regret lower bound found. Rissanen / Merhav-Feder read at abstract level only (fixed-parameter bounds).
- "Solo" is Piggott and Solo (ICASSP 2018), abstract only: simulations, no regret bound.
- Unverified: Herrera-Marin's Wei-Winnicki attribution. Not searched: Zhang 2023, Wang 2023, online EM, Bacry/Cartea rows.
