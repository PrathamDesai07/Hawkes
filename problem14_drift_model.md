# Problem 14: choice of drift model for the tracking bound (C5)

## Decision

**Bounded variation, instantiated as a piecewise-Lipschitz class P(Δ, J, D; c_μ, c_ε):** the true parameter path θ*_t, t ∈ [0, T], satisfies
- (box) μ_t ∈ [μ, c_μ μ] and ε_t = 1 − n_t ∈ [ε, c_ε ε] for all t, for fixed constants c_μ, c_ε (so the path stays in Θ_ε and the Fisher metrics along it are uniformly equivalent, see "Metric" below);
- (variation) θ*_t is Lipschitz with ‖θ̇*_t‖ ≤ Δ between at most J jump times, and each jump has size at most D. The path length is V_T ≤ ΔT + JD.

Drift is measured in a fixed reference Fisher norm ‖·‖_{I(θ_ref)} with θ_ref in the box.

Target form for Problem 13 (average tracking MSE of a constant-gain estimator with gain η):
**η σ²(ε) + Δ²/η² + J D²/(η T)**, with σ²(ε) = r/(2·eff) in event units (r = μ/ε; eff = information about μ, equal to m = E[1/λ] when α is known), numerically ≈ 0.27–0.5·μ²/ε² for K = 1 depending on k = 2μ/β (see the follow-up section); the lag term is the response of a time-varying-gain filter (gain η·r_t), which equals Δ²/η² for slow drift; the data-driven α-unknown recursion adds an amplitude-dependent excess, about 1.3–1.75× at relative amplitude 0.2–0.5 and about 1 as the amplitude goes to 0 (see 'Resolution' at the end; the earlier c_lag ≈ 1.8 is retracted).
The first two terms are the C5 target; the third is the re-convergence cost of regime jumps. In path-length form: avg bias² ≤ min(Δ_max V_T/(η²T), V_T²/(2ηT)).

## Why not the alternatives

**Smooth only (|θ̇*| ≤ Δ).** It gives exactly the problem-list target Δ²/η², but (i) it cannot represent regime changes, the main source of apparent criticality in market data (Filimonov–Sornette; the change-point papers on the Literature sheet), and (ii) its η-exponent misleads step-size tuning when jumps are present (jump bias scales like 1/η, smooth like 1/η²), which matters for Problem 15: minimising aη + Δ²/η² gives η* ~ (Δ²/a)^{1/3}, minimising aη + ρD²/(2rη) gives η* ~ (ρD²/(2ra))^{1/2}.

**Piecewise-constant only.** (i) It excludes smooth intraday seasonality in μ_t; (ii) a clean analysis needs every epoch to be longer than the estimator's relaxation time 1/(ηr) and than the mixing time 1/(εβ) (Problem 3), and the jump cost saturates at about D²/4 when relaxation is slower than the epoch (formula below, confirmed); (iii) it needs J small or known. It is a special case of the chosen class.

**Bare BV budget.** Standard comparator class, but one budget cannot separate the exponents, and the Cauchy–Schwarz bound is loose by 1/η for pure jumps.

**Random-walk drift** (E‖w‖² = q per step). Needs a distributional assumption; gives q/η (same exponent as jumps). A remark, not the primary model.

## Where the three terms come from (linear recursion, Fisher-preconditioned coordinates)

e_{t+1} = (1 − η)e_t + ηξ_t − d_t with ξ_t a martingale-difference noise and d_t = θ*_{t+1} − θ*_t.
- Noise: stationary variance ≈ η σ²/2.
- Smooth drift (sinusoid, amplitude A, per-event frequency ω): bias² = (A²/2) ω²/(η² + ω²), which is ≈ Δ²/η² for η ≫ ω.
- Jumps: a jump of size D relaxes like D e^{−κt}, κ = ηr per unit time. For alternating levels with fixed spacing L, solving the periodic recursion gives avg bias² = **D² tanh(κL/2)/(2κL)**, which is ρD²/(2ηr) (ρ = 1/L) when κL ≫ 1 and saturates at D²/4 when κL ≪ 1.
- Path-length form: Young's inequality on Σ_s (1−η)^s d_{t−s} gives avg bias² ≤ (1/η²)·mean(d²) ≤ Δ_max V_T/(η²T) and ≤ V_T²/(2ηT).
The score is a martingale difference under the true time-varying intensity, so drift enters only through the parameter lag, not through a mismatch of the state law. The curvature is the information of the current (random) state, the near-critical difficulty for Problem 13.

## Evidence (K = 1, α known, μ_0 = β = 1; two constant-gain estimators, one update per event)

- **naive:** μ̂ += η (r/m)(1/λ̂ − w) with m = E[1/λ], r = μ/ε (preconditioner = stationary information).
- **newton:** F ← (1−η)F + ΔF_i, μ̂ −= g_i/F, with ΔF_i = ∫_gap du/λ̂ (Fisher information with exponential forgetting, window 1/η events; the constant-gain version of the Problem 5 recursion; not the Problem 5 code).
`code/p14_v2.py`; ensemble-mean bias with estimator noise variance/R subtracted (replicates share the drift path); `code/p14_filter_check.py` and `code/p14_bias_scaling.py` for the first runs.

**(A) Noise coefficient σ²(ε)** (no drift, per-replicate MSE ≈ a(ε)η at small η; T = 2500, 10 replicates):

| ε | a(ε) = MSE/η (naive, η = 0.003) | prediction r/(2m) | a·ε² |
|---|---|---|---|
| 0.4 | 2.7–2.9 | 2.66 | 0.46 |
| 0.2 | 10.4–10.6 | 9.13 | 0.42 |
| 0.1 | 43–45 | 32.3 | 0.44 |

So **σ²(ε) ≈ 0.44/ε²** (event units, μ = β = 1), within 10–35 % of r/(2m); with m ≈ 2ε/(2μ−β) this is ~ μ(2μ−β)/(4ε²). The two estimators agree at small η (e.g. ε = 0.2: 10.4 vs 10.6).

**(B) Smooth drift bias vs the filter formula** (P = 1200, T = 6000, 10–12 replicates; entries: prediction | naive | newton):

| ε | A | η | η/ε² | prediction | naive | newton |
|---|---|---|---|---|---|---|
| 0.4 | 0.25 | 0.001 | 0.006 | 0.02545 | 0.02496 | 0.02551 |
| 0.4 | 0.25 | 0.003 | 0.019 | 0.01024 | 0.00954 | 0.01040 |
| 0.4 | 0.5 | 0.001 | 0.006 | 0.10179 | 0.10261 | 0.11055 |
| 0.4 | 0.5 | 0.003 | 0.019 | 0.04096 | 0.03920 | 0.05017 |
| 0.2 | 0.25 | 0.001 | 0.025 | 0.01635 | 0.01566 | 0.01668 |
| 0.2 | 0.25 | 0.003 | 0.075 | 0.00339 | 0.00539 | 0.00393 |
| 0.2 | 0.5 | 0.001 | 0.025 | 0.06538 | 0.06021 | 0.07790 |
| 0.2 | 0.5 | 0.003 | 0.075 | 0.01358 | 0.01476 | 0.01837 |
| 0.1 | 0.25 | 0.001 | 0.10 | 0.00672 | 0.01114 | 0.01070 |
| 0.1 | 0.25 | 0.003 | 0.30 | 0.00092 | 0.03247 | 0.00532 |
| 0.1 | 0.5 | 0.001 | 0.10 | 0.02690 | 0.02983 | 0.04410 |
| 0.1 | 0.5 | 0.003 | 0.30 | 0.00370 | 0.04306 | 0.01140 |

- **Filter formula and A² scaling confirmed across ε ∈ {0.4, 0.2, 0.1} and A ∈ {0.25, 0.5}** as long as η/ε² ≲ 0.03: naive within about 5 %, newton within 0–20 %; doubling A multiplies the bias by about 4. Above η/ε² ≈ 0.05 both estimators leave the linear regime (naive more than newton: at ε = 0.1, η = 0.003 naive 0.032 vs newton 0.0053 vs prediction 0.0009).

**(C) Jumps** (levels 1 ↔ 1.5, D = 0.5, mean spacing 600, earlier run, ε = 0.2, naive, prediction = realised-rate asymptote ρD²/(2ηr), ρ = 0.00133):

| η | κL ≈ ηr·600 | measured | asymptote | fixed-spacing formula D² tanh(κL/2)/(2κL) at L = 600 |
|---|---|---|---|---|
| 0.0005 | 1.5 | 0.0222 | 0.0667 | 0.0529 |
| 0.001 | 3 | 0.0168 | 0.0333 | 0.0415 |
| 0.002 | 6 | 0.0142 | 0.0167 | 0.0208 |
| 0.004 | 12 | 0.0077 | 0.0083 | 0.0104 |

The asymptote holds (within about 15 %) once κL ≳ 20 and the measured value saturates well below it for κL ≲ 3; the random spacing (many short regimes) lowers the measured values below the fixed-spacing formula, which assumes L = 600 for every regime. Slope between η = 0.002 and 0.004: −0.89 (jump exponent −1), against −1.36 for smooth drift on the same range (filter slope −2η²/(η²+ω²) → −2 only for η ≫ ω).

## What the earlier "noise-induced bias" caveat turned out to be

At larger η the ensemble-mean bias was not drift: it is the nonlinear regime, where the noise standard deviation is comparable to μ. Since MSE ≈ 0.44η/ε², the linearisation needs 0.44η/ε² ≲ (μ/2)², and the bias test shows the filter formula holds only for η/ε² ≲ 0.03. The forgetting-Newton recursion degrades much more gently than the naive one beyond that (ε = 0.1: bias 0.005 vs 0.032), so the Problem 5 style adaptive preconditioning is worth keeping; but neither removes the window.

**Window rule (K = 1, μ = β = 1, empirical):** the tracking bound's linear regime needs ω ≲ η ≲ 0.03 ε², where ω is the per-event drift frequency (ω = 2π/(P r) for a sinusoid of period P). The window is non-empty only if drift is slow: ω = 2πε/(Pμ) ≲ 0.03ε², i.e. P ≳ 210/(εμ) (β = 1 units), a few hundred times the mixing time 1/(εβ). **Drift faster than this cannot be tracked in the linear regime at all**, and that threshold is part of what the model class has to say; *Superseded by the follow-up section: the constants depend on k = 2μ/β, and the dimensionless criterion is aη/μ² ≲ 0.013–0.03.*

## Metric: fixing a reference norm (resolves the "I(θ*_t) varies with t" caveat)

With the box condition, all Fisher matrices along the path are uniformly equivalent: κ_I := max over θ, θ' in the box of λ_max/λ_min of I(θ')^{−1/2} I(θ) I(θ')^{−1/2} (K = 1, exact I; `code/p14_metric_equiv.py`):

| ε_0 | μ range ×1.5 | ε range ×1.5 | both ×1.5 | both ×2 | cond(I(θ_0)) |
|---|---|---|---|---|---|
| 0.2 | 7.3 | 3.3 | 19.5 | 149 | 48 |
| 0.1 | 8.2 | 3.4 | 22.5 | 183 | 152 |
| 0.05 | 9.0 | 3.3 | 25.6 | 221 | 522 |

κ_I is a constant of the box and is essentially independent of ε (3.3 for ε varying by ×1.5 at all three ε_0), whereas cond(I) grows like ε^{−2}. So converting drift to a fixed reference norm costs κ_I^{1/2} ≈ 3 to 5 for a ×1.5 box, **not** the ε^{−2} I feared earlier. It grows fast with the box size (×2 in both: 150–220), so c_μ, c_ε should be small constants. K = 2 is computed in the follow-up (κ_I ≈ 10–11 for μ ×1.5, 3.6–4.1 for ε ×1.5, 27–30 for both).

## Literature check (the earlier caveat)

- **Hall–Willett (2016), verified from the arXiv full text (1409.0031):** they do *not* assume a drift model for the parameters. Algorithm 1 (dynamic mirror descent with a known dynamical model Φ_t) is analysed as adversarial online learning with an oblivious environment: Theorem 1 gives regret ≤ C(1 + Σ‖λ_{t+1} − Φ_t(λ_t)‖²)√T against *any* comparator sequence of rate vectors, with steps η_t ∝ 1/√t or 1/√(T/δ), constants depending on δ, p, x_max, λ_max, λ_min, data in [0, x_max]^p; Theorem 2 gives the analogue with the network matrix W unknown, using min over W. The complexity measure is the squared deviation of the comparator from the known dynamics, not a bounded-variation budget on θ*. They remark that previous tracking regret bounds are sublinear only for comparators that change slowly or at a few instants. So (a) my earlier claim that they use a path-length assumption was wrong; (b) their measure Σ‖λ_{t+1} − Φ_t(λ_t)‖² is a squared-increment sum, and for increments d_t it counts Σd_t²: for smooth drift ≈ Δ²δT, for jumps ≈ JD², the same two ingredients that appear in my class; (c) their result is an adversarial regret in the rate vectors, with no dependence on ε, no parameter-error bound and a 1/√t step, while C5 is a parameter MSE for a constant step, so the comparison must be written as "different objects", not "a stronger or weaker version".
- **Solo / Piggott–Solo (ICASSP 2018), not verified.** The IEEE page is not retrievable here, and the one arXiv identifier that a search snippet suggested (1809.03878) is an unrelated Solo co-authored paper (topological brain-network distances). The only confirmed content (search-snippet abstract paraphrase): a fixed-gain adaptive (online) distributed estimator for Hawkes network parameters on a causal Laguerre basis, a two-time-scale algorithm based on exponentially weighted least squares that keeps the non-negativity constraints satisfied. Its drift or tracking assumptions are unknown to me. **Action needed from you (or via Problem 22):** obtain the PDF and check whether it states a tracking-error bound and under what parameter-variation model.

## Which caveats are resolved

| caveat | status |
|---|---|
| K = 1, α known, one ε, one amplitude | **mostly done** (see the follow-up: K = 1 α known and unknown, K = 2, k ∈ {0.6, 1, 2, 4, 8}); other drift shapes and the jump class at other ε, k, K not done |
| noise-induced bias at larger η | **resolved:** it is the nonlinear regime; linear window η ≲ θ_k μ²/a with θ_k ≈ 0.013–0.03 (k-dependent; η ≲ 0.03ε² only at k = 2); forgetting-Newton degrades less |
| earlier uninformative run | **resolved:** replaced by the ensemble-bias design |
| Hall–Willett assumptions | **resolved** (verified; assumption differs from what I had assumed) |
| Solo assumptions | **not resolved** (paper not accessible) |
| Fisher norm varies with t | **resolved for K = 1 and K = 2** (κ_I tables) |
| estimator is not the Problem 5 estimator | **resolved:** an exact copy of Problem 5's loop (fidelity checked) was run under drift; see the follow-up |

## Follow-up (2026-10-25): K ≥ 2 and α unknown, Problem 5's own estimator, and the μ/β dependence

### F1. Problem 5's estimator under drift (α unknown, K = 1)

`code/p14_p5_drift.py` is a copy of the loop in `code/p5_online_newton_k1.py` (same Fisher increment, same Newton step, same clipping, same prior), with two additions: the truth μ_t = μ_s(1 + A sin(2πt/P)) and an optional forgetting factor F ← (1−η)F + ΔF; η = 0 is Problem 5's algorithm verbatim. **Fidelity check:** with no drift and η = 0, the normalised error e_n of the copy equals Problem 5's `run()` output at the same seed to every printed digit (difference 0.0 at ε = 0.3, 0.1, 0.03; e_n = 0.5990, 0.6053, 9.5809).

Results (ε = 0.2, μ_s = β = 1, A = 0.5, P = 1200, T = 6000, 12 replicates, both μ and α estimated, start at the truth; `code/p14_p5_experiment.py`):

- **Verbatim (η = 0): it does not track.** The ensemble-mean μ̂ stays in [0.70, 0.80] over t > 0.2T while the truth swings over [0.5, 1.5]; the bias² is 0.185, 0.193, 0.199, 0.200 in the four successive time windows. Because F accumulates, the gain decays like 1/n, so the recursion converges toward a (biased) time-average and cannot follow drift. Problem 5's regret-style analysis is for a fixed parameter; **to track drift the recursion needs forgetting (a constant gain)**, which is exactly what the C5 model assumes.
- **With forgetting (η > 0) it tracks,** but its bias exceeds the α-known filter formula (F2) and the noise coefficient is that of the α-unknown model (F2).

### F2. α unknown: noise coefficient confirmed; drift-bias ratios (the c_lag part is RETRACTED, see R0 and R1 at the end)

**Noise.** For α unknown the per-event noise coefficient should be a = r/(2·eff), with eff = (m − ε/μ)/(1 − 2ε + εμm) the information about μ with α unknown (Problem 11), instead of r/(2m). At ε = 0.2, μ = β = 1: r = 5, m = 0.2739, eff = 0.1128, so a = 22.2 (α unknown) against 9.1 (α known). Measured: **21.8 (η = 0.001, P = 1200), 17.6–18.0 (η = 0.0002, 0.0005, P = 2400)**; α-known recursion on the same data: 10.4. So the extra noise from estimating α is real and quantitatively predicted (within 20 %).

**Smooth-drift bias** (ratio of measured ensemble bias² to the filter formula (A²/2)ω²/(η²+ω²)):

| α | K | ε | η | ratio |
|---|---|---|---|---|
| known | 1 | 0.4 | 0.002, 0.004 | 1.17, 1.28 |
| unknown | 1 | 0.4 | 0.002, 0.004 | 1.79, 2.00 |
| unknown | 1 | 0.3 | 0.0015, 0.003 | 1.75, 1.98 |
| unknown | 1 | 0.2 (P = 2400) | 0.0002, 0.0005, 0.001 | 1.69, 1.83, 1.87 |
| unknown | 2 | 0.4 | 0.002, 0.004 | 1.64, 1.87 |

- The excess for α unknown is **about 1.6–2.0, the same at K = 1 and K = 2, and independent of ε, η and of the amplitude** (ε = 0.4, η = 0.004: ratios 1.76 at A = 0.5 and 1.82 at A = 0.2 with 64 replicates, split-half 0.0075 / 0.0079). So it is not nonlinearity in A.
- Hypotheses I tested and rejected: (i) worse conditioning of I — the factor does not shrink as ε grows (cond(I) = 18.5 at ε = 0.4 vs 48 at 0.2, same factor); (ii) a wrong mean gain — the time-average of the realised gain matrix F⁻¹ΔF over a stationary run has eigenvalues 0.975–1.0 times η (ε = 0.4: 0.975, 0.997; ε = 0.2: 0.979, 1.000; `code/p14_gain_matrix.py`). A phase/amplitude decomposition of the ensemble-mean response at ε = 0.4 gives amplitude ratios 0.63 and 0.91 (filter: 0.69 and 0.89), a phase lag larger than the filter's by about 0.17–0.19 rad (a constant extra delay of roughly 85 events), and a mean offset of −0.13 at η = 0.002 (−0.01 at η = 0.004). **I did not identify the mechanism**; empirically the lag term should carry a constant c_lag ≈ 1.8 ± 0.2 when α is estimated (c_lag ≈ 1.2 when α is known): Δ²/η² → c_lag·Δ²/η² in the C5 bound, until the cause is understood.

### F3. K = 2 (α unknown): tracking and metric equivalence

Truth β = (1, 2), α_k = (1−ε)β_k/2, μ_s = 1.5, sinusoidal μ_t with A = 0.3, P = 1200, T = 6000, 12 replicates, ε = 0.4; three-parameter recursion with forgetting (`code/p14_k2.py`). The information with α unknown from stationary MC is (I⁻¹)_μμ = 12.45, so the predicted noise coefficient is r/(2eff) = 23.3; measured 21.9 (η = 0.002) and 27.8 (η = 0.004, already past the linear window). Bias ratios to the filter formula: 1.64 and 1.87 (table above). So the K = 1 conclusions (noise coefficient r/(2 eff), bias excess ≈ 1.6–2) carry over to K = 2 at the single (β, ε, μ) tested.

**Metric equivalence for K = 2** (stationary MC, 250 000–450 000 events per point, I good to about 5 %, so κ_I to about 10–20 %; `code/p14_metric_k2.py`, β = (1, 2), μ_0 = 1.5):

| ε_0 | μ range ×1.5 | ε range ×1.5 | both ×1.5 | cond(I(θ_0)) |
|---|---|---|---|---|
| 0.2 | 9.7 | 4.1 | 27 | 268 |
| 0.1 | 10.8 | 3.6 | 30 | 599 |

As for K = 1 (7–9, 3.3, 20–26), κ_I is a modest box constant that does not grow as ε falls, while cond(I) does (268 → 599). K = 2 values are about 30 % above K = 1. So fixing a reference Fisher norm costs κ_I^{1/2} ≈ 2 to 5 for a ×1.5 box, for K = 1 and 2, independent of ε.

### F4. μ/β dependence of the noise constant and of the linear window

**Exact rescaling.** If t → t/c, μ → cμ, β → cβ, the point process is mapped exactly onto itself and the per-event recursion is unchanged. So only k = 2μ/β matters; varying μ at β = 1 covers everything (not run separately).

**Noise coefficient** (K = 1, α known, ε = 0.2, forgetting-Newton, no drift, η = 0.003 and 0.006; `code/p14_mu_beta.py`):

| μ (k) | measured a = MSE/η | prediction μ/(2εm) | a ε²/μ² (measured) |
|---|---|---|---|
| 0.3 (0.6) | 0.59 – 0.64 | 0.52 | 0.27 |
| 0.5 (1) | 1.99 – 2.08 | 1.79 | 0.32 |
| 1 (2) | 10.4 – 11.0 | 9.13 | 0.42 – 0.44 |
| 2 (4) | 50.1 – 53.5 | 42.5 | 0.50 – 0.54 |
| 4 (8) | 208 – 207 | 184 | 0.52 |

The formula a = μ/(2εm), with the exact m from Problem 10, predicts the measured value within 9–18 % (measured is above) across k = 0.6 to 8. So the constant "0.44/ε²" is specific to k = 2; the relative noise a/μ² · ε² rises from 0.27 (k = 0.6) to about 0.5 (k ≥ 4). For k < 1, m ~ Cε^k (Problem 10), so a ~ ε^{−(1+k)}, not ε^{−2}.

**Linear-window onset** (smooth drift, per-event frequency ω = 0.001, A = 0.3, ε = 0.2, T = 4P, 8 replicates; ratio of bias² to the filter formula):

| μ (k) | η = 0.001 | 0.002 | 0.004 | 0.008 | 0.016 |
|---|---|---|---|---|---|
| 0.3 (0.6): ratio | 0.90 | 0.87 | 1.03 | 5.2 | 77 |
| 2 (4): ratio | 1.26 | 1.57 | 2.19 | 4.5 | 38 |
| aη/μ² (k = 0.6) | 0.006 | 0.012 | 0.023 | 0.046 | 0.092 |
| aη/μ² (k = 4) | 0.011 | 0.021 | 0.043 | 0.085 | 0.170 |

(and, from the ε scans at μ = 1, k = 2: ratio ≈ 1.1–1.3 up to aη/μ² ≈ 0.013, growing beyond). The ratio leaves ~1 at **a η/μ² ≈ 0.013 (k = 2), ≈ 0.015 (k = 4), ≈ 0.03 (k = 0.6)**: the dimensionless criterion is the relative noise standard deviation √(aη)/μ ≲ 0.11 to 0.17, and the earlier "0.03ε²" is its k = 2 instance. Revised window rule:

**ω ≲ η ≲ η_max = θ_k μ²/a = 2θ_k ε m μ, with θ_k ≈ 0.013–0.03 (θ_k ≈ 0.013 at k = 2, ≈ 0.015 at k = 4, ≈ 0.03 at k = 0.6),** where m = E[1/λ] (exact from Problem 10; for k > 1, m ≈ 2ε/(2μ−β) up to O(1) corrections at moderate ε). Non-empty iff drift is slow enough: ω ≲ η_max. Check: μ = β = 1, ε = 0.2: η_max = 2·0.013·0.2·0.274 = 0.0014 = 0.036ε², matching the earlier 0.03ε².

### Which caveats are now resolved

| caveat | status |
|---|---|
| K ≥ 2, α unknown not run | **done at one setting** (K = 2, β = (1,2), μ = 1.5, ε = 0.4): noise coefficient and the bias excess agree with K = 1; K ≥ 3, other grids, other ε not done |
| metric equivalence for K ≥ 2 | **done for K = 2** (ε_0 = 0.2, 0.1): κ_I ≈ 10–11 / 3.6–4.1 / 27–30, ε-independent |
| Problem 5's own code | **done:** exact copy (fidelity 0.0 difference), shows it does not track as written (η = 0) and tracks with forgetting; the α-unknown noise coefficient is confirmed |
| μ/β dependence | **done for k ∈ {0.6, 1, 2, 4, 8} (noise) and {0.6, 2, 4} (window onset)**: constants depend on k; replaced by the dimensionless criterion a η/μ² ≲ 0.013–0.03 |
| (new) extra lag when α is estimated | **superseded:** see R1 at the end (amplitude effect of the data-driven preconditioner; c_lag retracted) |
| still open | mechanism of c_lag; K ≥ 3 and other grids; drift shapes other than a sinusoid, jump class at other ε, k, K; Solo (paper not accessible) |

## Resolution of the remaining open items (2026-10-25, later): mechanism of the extra lag, K ≥ 3 and other grids, other drift shapes, jump class, Solo

All simulations in this section use **causal sampling** (μ̂ is read at the last event at or before each grid time) and the general engine `code/p14_gen.py` (any K; estimators `adaptive` = α unknown, F ← (1−η)F + ΔF with quadrature Fisher increments, and `known` = α known).

### R0. Two corrections to the earlier part of this note

1. **Sampling artifact.** The first version of `code/p14_mech.py` and `code/p14_k2.py` recorded μ̂ at the *first event after* each grid time, which peeks at the long gap containing that time and makes the estimator look more accurate. `p14_mech.py` is fixed. **The K = 2 numbers in F3 (bias ratios 1.64 and 1.87) came from `p14_k2.py` and are superseded by R3**; the F3 noise coefficient (21.9η) and the metric tables are unaffected. The other tracking scripts (`p14_v2`, `p14_p5_drift`, `p14_filter_check`) sample causally.
2. **The "c_lag ≈ 1.8" claim in F2 is retracted.** The mean-square bias ratio, as a statistic, is too noisy and mixes offsets and harmonics: for the same algorithm at the same parameters (ε = 0.4, η = 0.004, P = 1200) different seed sets gave ratios of 1.2 to 1.9. The decomposition and the transfer function below replace it.

### R1. Mechanism of the extra lag

**Demodulated transfer function** (`code/p14_transfer.py`, ε = 0.4, μ = β = 1, A = 0.3, P = 1200, 40 replicates; H = response amplitude at the drift frequency relative to A, compared with the first-order filter η/(η+iω)):

| η | estimator | H (measured) | H (filter) |
|---|---|---|---|
| 0.004 | α known | 0.768 − 0.428i | 0.785 − 0.411i |
| 0.004 | α unknown, fixed preconditioner | 0.799 − 0.411i | 0.785 − 0.411i |
| 0.004 | α unknown, adaptive | 0.814 − 0.518i (modulus +9 %, extra lag 0.085 rad) | 0.785 − 0.411i |
| 0.002 | α known | 0.454 − 0.503i | 0.477 − 0.499i |
| 0.002 | α unknown, adaptive | 0.434 − 0.540i | 0.477 − 0.499i |
| 0.002 | α unknown, fixed | 0.505 − 0.486i | 0.477 − 0.499i |

(standard errors 0.016–0.032). Inside the noise window all three estimators follow the first-order filter within about 1–2 standard errors; the only distinguishable deviation is the adaptive α-unknown estimator at η = 0.004.

**Decomposition of the ensemble-mean error** into offset², fundamental-frequency error², second-harmonic² and residual (`code/p14_decomp.py`, ε = 0.4, η = 0.004; fundamental error² as a ratio to the constant-gain filter formula (A²/2)ω²/(η²+ω²)):

| A | α known | α unknown, fixed preconditioner | α unknown, adaptive | adaptive / fixed |
|---|---|---|---|---|
| 0.2 | 1.07 | 1.07 | 1.41 | 1.32 |
| 0.35 | 1.15 | 0.98 | 1.48 | 1.48 |
| 0.5 | 1.24 | 0.92 | 1.61 | 1.75 |

Offsets² and harmonics are small except at A = 0.2 (offset² 0.003–0.005 for α unknown, i.e. a mean offset of about +0.05 to +0.07, about 5–7 % of μ), which is what had inflated the earlier total ratios there.

**Findings.**
- The mean gain matrix of the adaptive recursion is ηI (eigenvalues 0.975–1.0 times η; `code/p14_gain_matrix.py`), and the lag excess does not shrink when conditioning improves (cond(I) = 18.5 at ε = 0.4 against 48 at ε = 0.2, same excess). Those two causes are ruled out.
- With a *fixed* (deterministic) preconditioner the response follows the filter at every amplitude (ratios 1.07, 0.98, 0.92). With the *data-driven* preconditioner there is an additional excess that **grows roughly linearly with the relative drift amplitude** (adaptive/fixed = 1.32, 1.48, 1.75 at A = 0.2, 0.35, 0.5, extrapolating to about 1 as A → 0). This identifies the mechanism as an amplitude effect of normalising a time-varying information matrix with a lagged window, not a constant lag.
- A second, smaller part of the amplitude dependence is not due to the preconditioner at all (R2).
- I did not derive the amplitude dependence analytically.

### R2. The right deterministic baseline: a time-varying-gain filter

A per-event gain η is a gain η·r_t per unit time, where r_t = μ_t/ε is the local event rate; for a drifting μ_t the filter is therefore dy/dt = η r_t (μ_t − y), not the constant-gain filter (`code/p14_tv_filter.py`). For the sinusoid (ε = 0.4, η = 0.004, P = 1200) it predicts a bias² of 1.01, 1.04, 1.09 times the constant-gain value at A = 0.2, 0.35, 0.5 (α-known measurement: 1.07, 1.15, 1.24), so it explains about half of the amplitude dependence for α known. For the triangle wave (A = 0.3) it predicts 0.0161 at η = 0.002 and 0.0068 at η = 0.004; measured (α known) 0.0154 and 0.0066; adaptive 0.0170 at η = 0.002. **The Δ²/η² term should be understood as the response of this time-varying-gain filter** (equal to Δ²/η² for slow drift and small relative amplitude).

### R3. K ≥ 2 and other grids (causal sampling; ε = 0.4, A = 0.3, P = 1200; `code/p14_genK_transfer.py`)

| setting | η | H (measured) | H (filter) | noise coefficient: measured vs r/(2 eff) |
|---|---|---|---|---|
| K = 2, β = (1,2), μ = 1.5 | 0.002 | 0.656 − 0.574i (\|H\| 0.871) | 0.672 − 0.469i (0.820) | 23.7η vs 23.2η |
| K = 2, β = (1,2), μ = 1.5 | 0.004 | 1.044 − 0.388i (\|H\| 1.113) | 0.891 − 0.311i (0.944) | 30.4η vs 23.2η |
| K = 2, β = (1,4), μ = 1 | 0.004 | 0.750 − 0.452i (0.876) | 0.785 − 0.411i (0.886) | 5.9η vs 6.1η |
| K = 3, β = (1,2,4), μ = 1.5 | 0.002 | 0.864 − 0.485i (0.991) | 0.672 − 0.469i (0.820) | 25.6η vs 22.6η |
| K = 3, β = (1,2,4), μ = 1.5 | 0.004 | 1.112 − 0.269i (1.144) | 0.891 − 0.311i (0.944) | 32.2η vs 22.6η |
| K = 3, β = (1,3,9), μ = 1.5 | 0.002 | 0.650 − 0.469i (0.802) | 0.672 − 0.469i (0.820) | 19.5η vs 18.7η |

(standard errors 0.04–0.06; stationary-MC eff and cond(I): K2a 0.0806 and 178; K2b 0.2045 and 69; K3a 0.0830 and 1790; K3b 0.1000 and 755.) With η inside the window (aη/μ² ≲ 0.02) the noise coefficient r/(2 eff) holds to 3–13 % and the response follows the filter to 2 % (K = 3, (1,3,9)), 6 % (K = 2, (1,2)) and 21 % in modulus (K = 3, (1,2,4)); outside the window (η = 0.004 at μ = 1.5, aη/μ² ≈ 0.04) the modulus overshoots and the noise rises 30–40 %, as the window rule predicts. The deviation is not monotone in cond(I) (21 % at cond 1790, 2 % at cond 755), so the grid dependence is not captured by cond(I) alone. Mean offsets are about ±0.05 (3–5 % of μ). **Not run:** K ≥ 4, and ε < 0.4 for K ≥ 2 (the window forces η ≲ 0.001, which needs ~10⁵ time units per run).

### R4. Other drift shapes (K = 1, ε = 0.4; `code/p14_ramp2.py`, `code/p14_shapes.py`)

- **Linear ramp.** Steady lag with the offset removed (half the difference between an up-ramp and a down-ramp): α known −0.042 (η = 0.002) and −0.019 (0.004); adaptive −0.048 and −0.022; predicted −(slope/r)/η with the local event rate r = μ/ε = 3.0 at the ramp midpoint: −0.0444 and −0.0222. So the lag law Δ_event/η holds within about 15 % (α known 95 %, 86 %; adaptive 108 %, 101 %). With the nominal r = 2.5 the prediction would be off by 20 %, which again shows the event rate must be the local one.
- **Triangle wave** (A = 0.3, P = 1200): measured bias² 0.0154 and 0.0066 (α known, η = 0.002, 0.004) vs the time-varying-gain filter 0.0161 and 0.0068; adaptive 0.0170 at η = 0.002. The steady-state formula (4Aμ/(P r η))² overestimates by 2.3× because a half-period of 1500 events is only about three relaxation times; the corner-aware filter is the correct baseline.
- **Random walk** (variance q = 4·10⁻⁵ per time unit, stationary MSE ≈ aη + (q/r)/(2η), a = 2.8 for α known and 7.0 for α unknown from earlier measurements): α known 0.0117 ± 0.0011 vs predicted 0.0096 (η = 0.002) and 0.0161 ± 0.0021 vs 0.0132 (η = 0.004); adaptive 0.0228 ± 0.0032 vs 0.0180 and 0.0427 ± 0.0066 vs 0.0300. Measured values are 20–40 % above the prediction, larger at the larger η where the noise coefficient is inflated. (My first run used the wrong a for α known and showed measured below predicted; that was an error in the script, now fixed.)

### R5. Jump class at other ε, k, K (`code/p14_jump.py`; alternating levels μ_0 and 1.3μ_0, fixed spacing with η·(spacing in events) = 5; ratio measured/(time-varying-gain filter); 12–16 replicates)

| setting | η | α known | α unknown, adaptive |
|---|---|---|---|
| K = 1, k = 2, ε = 0.4 | 0.004 | 1.00 | 1.43 |
| K = 1, k = 2, ε = 0.2 | 0.004 (η/ε² = 0.1, outside the window) | 1.49 | 9.0 |
| K = 1, k = 2, ε = 0.2 | 0.001 | | 1.37 |
| K = 1, k = 0.6, ε = 0.4 | 0.004 | 1.01 | 1.09 |
| K = 1, k = 4, ε = 0.4 | 0.004 | 1.10 | 2.60 |
| K = 1, k = 4, ε = 0.4 | 0.001 | | 1.02 |
| K = 2, β = (1,2), μ = 1.5, ε = 0.4 | 0.002 | 1.04 | 1.21 |

The jump response (exponential relaxation after each jump, saturating at D²/4 when relaxation is slow, formula in the derivation section) holds to 0–10 % for α known and to 2–40 % for α unknown whenever aη/μ² is below the window threshold, and fails (ratio 2.6–9) exactly where the window rule says η is too large. Not run: K ≥ 3, ε < 0.2 with small η, jump sizes other than 0.3μ.

### R6. Solo (Piggott–Solo, ICASSP 2018): still not read

Crossref gives the authors as Marc Piggott and Victor Solo (ICASSP 2018, pp. 6333–6337, DOI 10.1109/ICASSP.2018.8462138). Unpaywall reports the paper as closed access with no open copy; the IEEE page does not return text. Crossref also lists a publisher PDF link labelled as a version of record for similarity checking; I did not use it to get around the paywall. The only confirmed content is the abstract paraphrase (fixed-gain adaptive distributed estimator, causal Laguerre basis, two-time-scale exponentially weighted least squares preserving non-negativity). **What needs to be checked in the full text, which has to come from you or your library:** (1) whether a tracking-error or MSE bound is stated at all; (2) the parameter-variation model, if any (random walk, slow drift, none); (3) the forgetting factor, which would play the role of our η; (4) how the two-time-scale scheme handles non-negativity compared with our projection; (5) the Laguerre basis versus our exponential dictionary.

### R7. Status of every item that was open

| item | status |
|---|---|
| mechanism of the extra lag | **resolved empirically:** not a constant factor; an amplitude effect of the data-driven preconditioner (absent with a fixed preconditioner, about 1.0 as A → 0, 1.3–1.75× at A = 0.2–0.5); plus a time-varying-gain effect of the per-event gain that is common to all estimators (R2); an analytic derivation is not done |
| K ≥ 3 and other grids | **partly:** K = 3 on (1,2,4) and (1,3,9) and K = 2 on (1,2) and (1,4) at ε = 0.4; K ≥ 4 and ε < 0.4 for K ≥ 2 not run |
| K = 2 at one setting | **done** (two grids, two η, jump class) |
| drift shapes other than a sinusoid | **done for K = 1:** ramp, triangle, random walk; jumps R5 |
| jump class at other ε, k, K | **done for K = 1 at ε ∈ {0.4, 0.2}, k ∈ {0.6, 2, 4}, and K = 2** |
| Solo | **not resolved:** paper not accessible to me |

Re-derive the linear-recursion computation and the tanh jump formula yourself; the numerics check them in the settings above.
