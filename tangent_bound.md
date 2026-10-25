# A tangent-line minorant: sharp rigorous lower bounds for the whole spectrum of I, and an information level for the regret chain

Strengthens Problems 2, 12 and 5. Code: `code/p2_tangent_bound.py` (double precision), `code/p2_tangent_mp.py` (60-digit multiprecision), `code/p2_tangent_check.py` (independent simulation check), `code/p12_tangent_spectrum.py`. Tags: **[P]** proved on paper, **[P\*]** modulo listed inputs, **[N]** numerical, **[O]** open.

## 1. The minorant  **[P]**

1/λ is convex on (0, ∞), so for every λ₀ > 0: 1/λ ≥ 2/λ₀ − λ/λ₀² for all λ > 0. Since aaᵀ ⪰ 0 this holds as a matrix inequality, pathwise and in expectation:
 **I = E[aaᵀ/λ] ⪰ N(λ₀) := (2/λ₀)·M − (1/λ₀²)·M₃,  M = E[aaᵀ],  M₃ = E[λ aaᵀ] = μM + Σ_kα_k E[S_k aaᵀ].**
M and M₃ involve only moments of S up to order 3, which are **exact** solutions of linear systems (from E[𝒜f] = 0):
- order 1: E S_k = μ/(εβ_k); order 2: (β_i+β_j)Q_ij = E[λS_i] + E[λS_j] + Eλ, with E[λS_i] = μ m_i + Σ_kα_kQ_ki  (K² unknowns);
- order 3: 𝒜(S_iS_jS_l) = −(β_i+β_j+β_l)S_iS_jS_l + λ[S_iS_j + S_iS_l + S_jS_l + S_i + S_j + S_l + 1], so (β_i+β_j+β_l)T_ijl = Σ_{pairs} (μQ_pq + Σ_kα_kT_kpq) + Σ_p(μm_p + Σ_kα_kQ_kp) + Eλ  (K³ unknowns).
Integrability: all polynomial moments are finite (Lemma `lem:moments`/`lem:expmoment`), so E[𝒜f] = 0 is justified. Check (`p2_tangent_check.py`): second and third moments agree with a 4·10⁵-event simulation to 0.4–0.9 %.
By Weyl monotonicity, **λ_j(I) ≥ λ_j(N(λ₀)) for every ordered eigenvalue j and every λ₀ > 0**; one may take, for each j, the best λ₀.

## 2. Problem 2: the bound on λ_min(I) is now within a factor 1.3–3.5 of the truth  **[P]**

Multiprecision values (60 digits, `p2_tangent_mp.py`) and double-precision values agree; rigorous lower bound on λ_min(I) against the Monte Carlo value of the earlier notes (μ = 1, β_k = r^{k−1}, α_k = (1−ε)β_k/K):

| K, r | ε | MC λ_min(I) | tangent bound | ratio | old closed-form chain | old Theorem T (truncation) |
|---|---|---|---|---|---|---|
| 1, 2 | 0.3 | 0.139 | 0.0791 | 0.57 | 4.8·10⁻⁵ | 6.9·10⁻³ |
| 1, 2 | 0.1 | 0.0654 | 0.0249 | 0.38 | 3.7·10⁻⁶ | 1.4·10⁻³ |
| 2, 2 | 0.3 | 0.0337 | 0.0264 | 0.78 | 1.0·10⁻⁷ | 6.4·10⁻⁴ |
| 2, 2 | 0.1 | 0.0455 | 0.0240 | 0.53 | 7.8·10⁻⁹ | 4.6·10⁻⁴ |
| 2, 4 | 0.1 | 0.0677 | 0.0270 | 0.40 | 6.9·10⁻⁹ | 1.6·10⁻⁴ |
| 3, 2 | 0.3 | 0.0033 | 0.0024 | 0.73 | 1.1·10⁻⁹ | 1.7·10⁻⁵ |
| 3, 2 | 0.1 | 0.0045 | 0.00285 | 0.63 | 8.7·10⁻¹¹ | 1.7·10⁻⁵ |

The tangent bound also has the **correct scaling in ε**: for ε → 0 its ratio to ε settles at 0.24–0.27 (K = 1: 0.264, 0.249, 0.240, 0.237 at ε = 0.3, 0.1, 0.03, 0.01; K = 2: 0.260 at ε = 0.03 and 0.01; K = 3 reaches 0.268 at ε = 0.01 after the middle tier stops being the minimum), against the true λ_min/ε → 1/(μ(k−1))-type constants. The best λ₀ is ≈ 1.4–2.2·Eλ. The dependence on K and the grid ratio enters through the exact moment systems (no further estimate is lost), which meets the criterion of Problem 2 with constants that are no longer 10³–10⁸ off. Not claimed: a closed-form formula in (K, r) (the bound is computed from three explicit linear systems and one scalar optimisation).

## 3. Problem 12: every tier of the spectrum is bounded below rigorously  **[P] (bounds), [N] (comparison)**

Weyl's inequality gives a lower bound for each eigenvalue. Against the Monte Carlo spectra of Problem 12 (μ = 1.5; values: MC | rigorous lower bound | ratio):

| case | smallest | middle | top |
|---|---|---|---|
| K=2, β=(1,2), ε=0.1 | 0.032 \| 0.0143 \| 0.45 | 0.065 \| 0.0423 \| 0.65 | 18.2 \| 15.1 \| 0.83 |
| K=2, ε=0.03 | 0.015 \| 0.0045 \| 0.30 | 0.064 (middle) \| 0.0445 \| 0.70 | 57.2 \| 48.4 \| 0.85 |
| K=2, ε=0.01 | 0.0053 \| 0.0015 \| 0.28 | 0.066 \| 0.0455 \| 0.69 | 171 \| 144 \| 0.84 |
| K=3, β=(1,2,4), ε=0.1 | 0.0047 \| 0.0033 \| 0.71 | 0.0446 \| 0.0162 \| 0.36; 0.1495 \| 0.0943 \| 0.63 | 19.1 \| 15.1 \| 0.79 |
| K=3, ε=0.03 | 0.0055 \| 0.0035 \| 0.64 | 0.0225 \| 0.0051 \| 0.23; 0.1652 \| 0.1051 \| 0.64 | 60.9 \| 48.8 \| 0.80 |

Consequences, now proved (up to the validity of the exact moment systems and 60-digit evaluation): (i) the top eigenvalue is ≥ 0.8× its Monte Carlo value and ∝ 1/ε; (ii) the K−1 middle eigenvalues are bounded below by ε-independent constants set by the grid (e.g. 0.0455 for K = 2, β = (1,2), the same at ε = 0.1, 0.03, 0.01, against the heuristic 0.066 of the spectral formula); (iii) the smallest eigenvalue is bounded below by ≈ 0.25 ε for small ε. So the three-tier picture of Problem 12 is a theorem as regards **lower** bounds. Upper bounds remain: ridge Rayleigh quotient (proved, uses E[1/λ] = O(ε)) and the heuristic spectral formula for the middle tier (not proved). The statement that the smallest eigenvector converges to the ridge direction is still only numerical.

## 4. Problem 5: an information level without Assumption M  **[P\*]**

Replace, in Theorem FK (`K2_regret_chain.md`), the truncation device by the tangent minorant. Pathwise, J_T = ∫aaᵀ/λ\* dt ⪰ N_T(λ₀) := ∫₀^T aaᵀ(2/λ₀ − λ\*/λ₀²)dt, whose entries are time integrals of the cubic polynomials f_ij(S) = a_ia_j(2/λ₀ − λ(S)/λ₀²) with E f_ij = N(λ₀)_ij. Let λ_\* := sup_{λ₀}λ_min(N(λ₀)) > 0 (§2). Concentration of ∫f_ij: f_ij is not Lipschitz globally, so clip the state at s_max (componentwise), f̃_ij = f_ij(min(S, s_max)), Lipschitz in d_ℓ with a constant Λ_f(s_max) = O(s_max²) polynomial in s_max; on the event {max_{t≤T}S ≤ s_max} (probability ≥ 1 − (μ/ε)T·C e^{−θ_\* s_max/…} by Lemma `lem:expmoment`, an event of probability ≥ 1 − T^{−κ} for s_max = O((ln T)/ε)) f̃ = f, and |E f̃ − E f| ≤ E[|f|; S > s_max] is exponentially small. The Poisson-equation/Freedman argument of Theorem FK then gives, with probability ≥ 1 − 3(K+1)(K+2)e^{−c_fT} − T^{−κ},
 **J_T ⪰ T(λ_\* − τ) I_{K+1}**, τ = λ_\*/2,
i.e. information level λ_\*/2 instead of λ_L/(4L). For the K = 2 example of `K2_regret_chain.md` (β = (1,2), μ = 1.5): λ_\* = 0.0252 (ε = 0.4), 0.0264 (ε = 0.3) against the old level λ_L/(4L) ≈ 10⁻⁶ — a gain of about **10⁴** in the constant T₁ ∝ 1/level of the bootstrap and in the burn-in. For K = 1 this also makes Assumption M unnecessary (the chain can use the tangent level instead of the two-interval device), at the price of a smaller exponential rate c_f because Λ_f(s_max) is polynomial in s_max.
Limits: the rate c_f ≍ τ²γ²ε/(Λ_f(s_max)²μ) is small (Λ_f grows like s_max² ~ ε⁻²ln²T); this affects T₁, not the leading constant.

## 5. Status
| item | before | now |
|---|---|---|
| P2: bound in K and grid ratio | closed-form chain 10³–10⁸ loose | rigorous, within 1.3–3.5× of truth, correct ε-scaling, all K (computed from exact moment systems) |
| P12: middle tier | heuristic | lower bounds for **all** eigenvalues proved (0.23–0.85 of MC); upper bounds / eigenvector still open |
| P5: information level, K ≥ 2 | 10⁻⁶ (truncation) | 10⁻² (tangent), Assumption M no longer needed |
