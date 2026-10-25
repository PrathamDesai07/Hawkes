# Problem 9: lower bounds for online estimation, and where they meet the upper bound (mathematics, then code)

Tags: **[P]** proved here, **[P\*]** proved modulo the inputs named, **[H]** heuristic (prediction), **[N]** numerical (section 8, filled after the runs), **[O]** open.
Code: `code/p9_bayes.cs`, `code/p9_run.sh`, `code/p9_analyze.py`. Date: 25 Oct 2026. Nothing refereed.
Regret means `Reg_T = KL(P*_T ‖ Q_T)` of `problem08_regret_notion.md`, in time units (events `n`, `T = t_n`); `d = K+1`.

## 0. Where the problem stands, and what is redefined

Tracker row 9: "adapt the kernel-agnostic two-point bound to exponential dictionaries, convert it into a regret or minimax statement, compare with the upper bound in eps". Status before this session:
- The two-point KL bound itself is prior art (Herrera-Marin, arXiv 2609.38268, Thm 2.2 + Cor 2.3). `novelty_reassessment.md` therefore redirects the problem to (a) a regret lower bound for sequential estimators `>= (d/2) ln T + (1/2) ln det I - O(1)` and (b) a matching *upper* bound on `lam_min(I)`.
- Proposition LB (`box_cost_and_lower_bound.md`): `sup Reg >= (d/2) ln(n/ln n) - O(1)`, **conditional on (U)**, uniform scalar MSE of an estimator; no eps-explicit constant; non-consistency at `eps = c/T` (Problem 11, Thm NC).
Missing for Done (tracker): the comparison table in eps, mu, T and an eps-explicit lower bound for `eps = eps(n) -> 0`.
This note (i) removes (U) as a hypothesis (it is a statement about our own estimator, proved from the Theorem R chain), (ii) replaces the scalar MSE by a *matrix* covariance, which produces the term `(1/2) ln det I` that the upper bound has, (iii) gives the eps-explicit comparison, (iv) treats the crossover `eps ~ 1/T` through the Gaussian-experiment benchmark, (v) records the exponential-dictionary two-point bound in Fisher form.

## 1. Two-point bound in Fisher form for exponential dictionaries [P]

**Proposition 9.0.** Fix `theta, theta'` in a ratio set (`c^- <= lam_theta'/lam_theta <= c^+` on the whole record) and let both laws start from the **same** initial state `S_0 = s_0`. Then
 `KL(P_theta^T ‖ P_theta'^T) = E_theta int_0^T lam_theta psi(lam_theta'/lam_theta) dt`  (Problem 8, Definition 1 with the constant predictor `theta'`)
 `<= (1/(2 c~^{-2})) (theta'-theta)^T E_theta[J_T(theta)] (theta'-theta)`, `c~^- = min(1,c^-)`,
and `>= (1/(2 c~^+)) (theta'-theta)^T E_theta[J_T] (theta'-theta)`. With a stationary start `E_theta J_T = T I(theta)` (and `T I + O(1)` from a fixed state).
*Proof.* Lemma `lem:regid` (the likelihood-ratio formula) and Lemma `lem:regsand` of Problem 8. ∎
Comparison with Herrera-Marin Thm 2.2 (kernel-agnostic, `KL <~ delta^2 T mu (1-n)`): ours is **direction-resolved**: `u^T I u` is `~ eps/(mu(k0-1))` along the ridge (matches `mu(1-n) = mu eps` up to constants) but `~ 1/eps` along the top eigenvector, which the isotropic bound cannot see. Numerical check of the bound in section 8 (T7).
**Caveat (final audit).** With a *stationary and unobserved* past, `S_0 ~ pi_theta` differs in law between `theta, theta'`, so the marginal term `TV(pi_theta, pi_theta')` or `KL(pi_theta ‖ pi_theta')` must be added (`theoremA_step2_and_stationary_tv.md`). All statements below are for the experiment with `S_0` observed or ancillary, or carry that extra term, which is bounded by the KL-diameter of `pi_theta` over `Theta_0` (O(delta^2) for a box of relative width delta).

## 2. A matrix-covariance Bayes lower bound (replaces scalar (U)) [P]

**Proposition 9.1.** Let `Theta_0` be a box of positive volume `|Theta_0|`, `Theta ~ Unif(Theta_0)`, `X` the record on `[0,T]` (with `S_0` ancillary). Suppose a statistic `theta_hat = theta_hat(X)` satisfies the matrix bound
 **(U\*)** `E_theta[(theta_hat - theta)(theta_hat - theta)^T] ⪯ (kappa_T/T) A^{-1}` for all `theta in Theta_0`, `A ≻ 0` fixed.
Then for every predictor `Q` (every sequential or batch predictor)
 **(9.1)** `sup_{theta in Theta_0} Reg_T(Q; theta) >= avg_{Theta_0} KL(P_theta ‖ Q) >= I(Theta; X) >= ln|Theta_0| + (d/2) ln( T/(2 pi e kappa_T) ) + (1/2) ln det A`.
*Proof.* (i) sup >= average; (ii) `E_Theta KL(P_Theta ‖ Q) = I(Theta;X) + KL(P_bar ‖ Q) >= I(Theta;X)` (compensation identity); (iii) data processing `I(Theta;X) >= I(Theta;theta_hat)`; (iv) `I(Theta;theta_hat) = h(Theta) - h(Theta|theta_hat) >= ln|Theta_0| - h(Theta - theta_hat)` and, by the Gaussian maximum-entropy bound, `h(Theta - theta_hat) <= (1/2) ln((2 pi e)^d det E[(Theta-theta_hat)(Theta-theta_hat)^T])`; (v) the determinant is monotone in the PSD order: `det E[...] <= det((kappa_T/T)A^{-1}) = (kappa_T/T)^d / det A`. ∎
With a scalar bound `E‖theta_hat - theta‖^2 <= C_1 ln n/n` (Prop LB) the last step uses AM-GM and loses the term `(1/2) ln det A`; the matrix form keeps it, and `det A` is where `eps` enters.

## 3. The online Newton estimator satisfies (U\*) with explicit constants [P\*]

**Proposition 9.2.** Let `theta_hat = theta_{n+1}` (Newton iterate after the last event before `T`), `F_n` its Fisher matrix. On the good event `G_T = {F_n ⪰ (T/(2 c^+)) N(theta)}` (information concentration: tangent minorant `N(theta) ⪯ I(theta)` of `tangent_bound.md` with `J_T ⪰ (T/2) N`, probability `>= 1 - C e^{-c_f T}`),
 `v v^T ⪯ ‖v‖^2_{F_n} F_n^{-1}`,  `F_n^{-1} ⪯ (2 c^+/T) N^{-1}`,
and on `G_T` the Fisher-norm error satisfies `E[D_n 1_{G_T}] <= c_D d ln T` (Theorem C-exp + Lemma C1: `D_n <= C(d ln T + x)` with probability `>= 1 - 5 e^{-x}`). Hence
 `E[v v^T] ⪯ (2 c^+ c_D d ln T / T) N^{-1} + diam^2 C e^{-c_f T} I`, i.e. (U\*) holds with `A = N(theta)` (uniform over `Theta_0` if `N(theta) ⪰ A` there; `N` is continuous in `theta`), `kappa_T = 2 c^+ c_D d ln T (1 + o(1))`.
*Proof.* Two displays above plus `E[D_n F_n^{-1} 1_G] ⪯ (2c^+/T)N^{-1} E[D_n 1_G]`. The inputs are exactly those of Theorem R (Theorem C-exp, concentration of `J_T`), so the status is the same: **[P\*]**. ∎
Consequence: (U) of Prop LB is *not an extra hypothesis*; it is a consequence of the upper-bound chain. The price is `kappa_T ∝ ln T` (hence the `(d/2) ln ln T` loss); removing it needs `E D_n = O(d)` (measured: E D_n = 1.9-2.0 for d = 2, `p8_out`), i.e. `E X_n = E G_n + O(1)`, [O].

## 4. The lower bound with its eps-dependence [P\*]

**Theorem 9.3.** For `Theta_0` in the interior of the ratio set, `T >= t_*(eps)` (mixing and leverage burn-in of Problem 6, `T eps >> ln`),
 **(9.2)** `sup_{Theta_0} Reg_T >= (d/2) ln T + (1/2) ln det N_* + ln|Theta_0| - (d/2) ln( 2 pi e kappa_T ) ,  N_* := min_{Theta_0} N(theta)` (in the PSD order, or at the centre up to continuity).
Matching upper bound (Problem 6, Theorem 6 with a prior `F_0` matched to `Theta_0`, `F_0 = c_0 diag(|width_j|^{-2})`, so `-(1/2) ln det F_0 = ln|Theta_0| + O(d)` and `P_0 = O(d)`):
 **(9.3)** `sup_{Theta_0} Reg_T <= (1+O(delta+eta)) { (d/2) ln T + (1/2) ln det I_* } + ln|Theta_0| + O(d c^+ ln t_*) + O(1)`,  `I_* = max I`.
**Reading in eps** (K = 1, `Theta_0 = [mu(1-w), mu(1+w)] x [beta(1-3 eps/2), beta(1 - eps/2)]`, the natural stability band, `|Theta_0| = 2 w mu beta eps`):
- the common term `ln|Theta_0| = ln eps + ln(2 w mu beta)`: **the minimax regret over the stability band carries `- ln(1/eps)` in time units**, which is the prior-volume cost of a parameter set that shrinks with eps (its Fisher radius stays O(1) because the Fisher information in alpha is `~ 1/eps`);
- `(1/2) ln det I(eps)` and `(1/2) ln det N(eps)`: for `k0 > 1`, `det I -> 1/(beta^2(k0-1))` and `det N` is bounded above and below by eps-free constants (K = 1: `ln det N` = -1.39, -1.62, -1.70, -1.72 for eps = 0.3, 0.1, 0.03, 0.01), so the gap between (9.2) and (9.3) in eps is only `(1/2) ln(det I/det N) = O(1)` plus `(d/2) ln ln T`: **the lower and upper bounds match in the eps-dependence of the constant, for k0 > 1, up to eps-free terms**;
- for `k0 < 1`: `det I ≈ (C eps^{k0-1} - 1)/beta^2` (Problem 6 section 17) gives `+((1-k0)/2) ln(1/eps)` in (9.3); the tangent minorant `N` (moments of order <= 3) cannot see the heavy tail of `1/lam` near `mu`, so (9.2) lacks this term: the bounds differ by `((1-k0)/2) ln(1/eps)` unless `N` is replaced by a minorant that uses `P_pi(lam <= 2 mu) ≍ eps^{k0}` [O];
- mu and beta enter only through `k0` (and `ln mu, ln beta` in the common volume term): `(1/2) ln det I = -ln beta - (1/2) ln(k0-1)`.

## 5. The regime eps ~ 1/T: crossover and the Gaussian-experiment benchmark [H], [P] for the benchmark

When `T lam_min(I) = O(1)` the ridge direction is not learnable and (9.2) does not apply (Theorem NC: no consistent estimator of mu). The experiment is then not LAN in the ridge, but the *Bayes information* of the Gaussian experiment with the same Fisher matrix is explicit:
 **(9.4)** `I_G(T) := (1/2) ln det( I_d + T Sigma^{1/2} I Sigma^{1/2} )`,  `Sigma = diag(width_j^2/12)` (prior covariance of the box).
It interpolates exactly: `I_G = (1/2) sum_j ln(1 + T lam_j sigma_j^2)` in the eigenbasis of `Sigma^{1/2} I Sigma^{1/2}`; as `T lam_min sigma^2 -> 0` the ridge term vanishes and `I_G -> ((d-1)/2) ln T`, as `T lam_min sigma^2 >> 1` it is `(d/2) ln T + (1/2) ln det I + ln|Theta_0| + const`.
- *Upper side* [P\*]: the tier-resolved Newton bound `E R_T <= (1/2) sum_j ln(1 + T lam_j / f_j)` of Problem 6 section 4 is (9.4) with `f_j = 1/sigma_j^2` for `F_0 = Sigma^{-1}` and the tangent/true eigenvalues.
- *Lower side*: for `T lam_min sigma^2 >> 1` Theorem 9.3; for `T lam_min sigma^2 = O(1)` the lower bound is `((d-1)/2) ln T - O(1)` from (9.1) applied to the learnable directions (the covariance bound (U\*) holds for the projection onto the span of the d-1 non-ridge eigenvectors), plus the ridge contribution `>= 0`. The remaining gap is the ridge term `(1/2) ln(1 + T lam_min sigma_r^2)` itself, which is `O(1)` in this regime: [O] to certify below by a better-than-0 term (needs non-asymptotic LAN or a two-point argument in the ridge direction).

## 6. The comparison table to be filled (K = 1, then K >= 2)

| regime (`eps`, `T`) | upper (9.3) | lower (9.2) | gap | status |
|---|---|---|---|---|
| fixed `eps`, `k0 > 1`, `T >> t_*` | `(d/2) ln T + (1/2) ln det I + ln|Theta_0| + O(1)` | `(d/2) ln T + (1/2) ln det N + ln|Theta_0| - (d/2) ln(c ln T)` | `(d/2) ln ln T + (1/2) ln(det I/det N) + O(1)`; no eps-dependence | [P\*] |
| `eps -> 0`, `T eps >> ln T`, `k0 > 1` | same, constants uniform in eps | same | same | [P\*], uniformity in eps from Problem 6 |
| `k0 < 1` | extra `((1-k0)/2) ln(1/eps)` | none (tangent minorant) | `((1-k0)/2) ln(1/eps)` | [O] |
| `T eps = O(1)` | `(1/2) sum_j ln(1+T lam_j f_j^{-1})` | `((d-1)/2) ln T - O(1)` | `(1/2) ln(1 + T lam_min sigma_r^2) = O(1)` | [H], benchmark (9.4) |
| two-point (Prop 9.0) | `KL <= T u^T I u/(2 c~^{-2})` | `KL >= T u^T I u/(2 c~^+)` | ratio `c~^+ c~^{-2}` (-> 1 as delta -> 0) | [P] |

## 7. Predictions tested in the code [H]

- **T1.** The measured covariance of the Newton iterate over a prior draw satisfies `T Cov ⪯ kappa_eff N^{-1}` with `kappa_eff = O(1)`, in fact `≈ I^{-1}` (efficiency), uniformly in eps once `T eps >> 1`.
- **T2.** The bound of Prop 9.1 evaluated with the *measured* covariance, `L_emp = ln|Theta_0| - (1/2) ln((2 pi e)^d det Cov)`, is a valid lower bound on the Bayes regret of our own predictor: `L_emp <= E_Theta Reg_T(Newton)`, for all eps and T.
- **T3.** The gap `E_Theta Reg_T(Newton) - L_emp` is `O(1)` and flat in eps (k0 > 1) and in T (it should converge to a constant of order d/2).
- **T4.** `E_Theta Reg_T(Newton) ≈ I_G(T)` (Gaussian-experiment benchmark (9.4)) + `O(1)`, including `T eps ~ 1`.
- **T5.** The ln|Theta_0| term: Bayes regret over the stability band decreases by `ln 3` per factor 3 decrease of eps at fixed `T eps >> 1`.
- **T6.** Crossover: at fixed `T`, as `T eps` goes from `>> 1` to `< 1`, the measured Bayes regret loses `~ (1/2) ln(T eps)`-type terms: slope in `ln T` drops from `d/2` to `(d-1)/2`.
- **T7.** Two-point (Prop 9.0): `E_theta R_T(theta') / (T u^T I u/2)` lies in `[1/c~^+, 1/c~^{-2}]` and tends to 1 as `delta -> 0`.

## 8. Numerical results (25 Oct 2026; `code/p9_bayes.cs`, `code/p9_run.sh`, `code/p9_run2.sh`, `code/p9_analyze.py`; raw tables `code/p9_out/`)

Setup (K = 1): `Theta_0 = [mu(1-w), mu(1+w)] x [beta(1-3 eps/2), beta(1-eps/2)]`, `w = 0.25`, `Theta ~ Unif(Theta_0)` redrawn per replicate, stationary start under the drawn theta (state observed), online Newton with `F_0 = 12 diag(1/width_mu^2, 1/width_alpha^2)` (precision of the uniform), start at the box centre, clip; 60-300 replicates (200-300 up to T = 3000, 60-120 to T = 1e4 for eps <= 0.01). `ln|Theta_0| = ln(2 w mu beta eps)`. Columns: E R = Bayes regret of the Newton predictor (which is >= I(Theta;X)), `L_emp` = Prop 9.1 with the *measured* covariance (valid, being the bound for this very statistic), `I_G` = (9.4) with the measured `Ibar`.

**T2 — `L_emp <= E R`: holds in all 31 rows tabulated** (`p9_analyze.py` (A)), eps in {0.3, 0.1, 0.03, 0.01, 0.003}, `T eps` from 0.28 to 285.

**T3 — the gap `E R - L_emp` is `O(1)` and flat:** 0.36-0.46 (eps = 0.3), 0.27-0.38 (0.1), 0.05-0.37 (0.03), 0.24-0.33 (0.01), 0.24-0.40 (0.003). So the lower bound (9.1) with the measured covariance is within about 0.4 of the regret of the Newton predictor over the whole range, including `T eps < 1`. (The 0.05 at eps = 0.03, T = 9487 is the best case: the Newton predictor is within 0.05 of the information-theoretic lower bound.)

**T1 — efficiency of the covariance:** `kappa_eff = lambda_max((T Ibar)^{1/2} Cov (T Ibar)^{1/2})` rises from 0.04 (`T eps = 0.3`, prior-limited) through 0.54 (T eps ≈ 10) and 0.6-0.86 (T eps ≈ 30) to 0.8-1.26 (T eps ≈ 100-300): the estimator is efficient (`kappa_eff ≈ 1`), not `kappa_T ∝ ln T`. So the `(d/2) ln ln T` loss of Theorem 9.3 is an artefact of the proof (Prop 9.2 uses the Fisher-norm tail `D_n <= C d ln T`); empirically `kappa ≈ 1`.

**T4 — benchmark `I_G`:** `E R - I_G` = -0.03...+0.00 for `T eps <= 3` (all eps), 0 to -0.14 at `T eps = 10-30`, then -0.09 (eps = 0.3), -0.19 (0.1), -0.38 (0.03), -0.12 (0.01) at `T eps ≈ 100-300`. The asymptotic offset is explained: `I_G` is the information of a *Gaussian* prior; for the uniform prior of the same variance `h` is smaller by `(d/2) ln(2 pi e/12) = 0.35` (d = 2), so `E R ≈ I_G - 0.35 + (Newton excess over I)`, consistent with -0.38 at eps = 0.03 (s.e. 0.09). The shape of the whole curve is reproduced by (9.4) from `T eps = 0.3` to `300`.

**Scaling collapse (new, T5/T6 revised).** With the band `Theta_0(eps)` of width `∝ eps` in alpha, both Fisher directions have `x_j = T lam_j sigma_j^2 ∝ T eps`, and the Bayes regret is a function of `x = T eps` alone:

| eps \ x = T eps | 0.3 | 1 | 3 | 10 | 30 | 100 | 300 |
|---|---|---|---|---|---|---|---|
| 0.3 | - | 0.075 | 0.170 | 0.394 | 0.849 | 1.532 | 2.344 |
| 0.1 | 0.024 | 0.067 | 0.165 | 0.394 | 0.830 | 1.500 | 2.400 |
| 0.03 | 0.021 | 0.060 | 0.156 | 0.380 | 0.787 | 1.385 | 2.247 |
| 0.01 | 0.025 | 0.069 | 0.165 | 0.396 | 0.897 | 1.625 | - |
| 0.003 | 0.026 | 0.070 | 0.181 | 0.438 | 0.843 | - | - |

At `x = 3` the five values span 0.156-0.181 for eps ranging over a factor 100; at `x = 30` 0.79-0.90. So **the whole eps-dependence of the minimax regret over the stability band is the dependence through `T eps` (in events, `n eps^2/mu`)**, and for large x it is `ln x + const` (slope `d/2 = 1` in `ln T`, with `ln|Theta_0| = ln eps + const` supplying the `ln eps`). The earlier prediction "slope drops from d/2 to (d-1)/2 as T eps falls below 1" (T6) is **wrong for this band**: both directions switch off together because the prior width in alpha shrinks with eps; the local slope of `E R` in `ln T` at `x = 1, 3, 10, 30, 100` is 0.03, 0.06, 0.14, 0.29, 0.48 (all eps), approaching `d/2` only for `x >> 100`. A `(d-1)/2` plateau appears only for a prior that is wide in the learnable direction; the benchmark (9.4) covers both cases.

**Additive constant** `c = E R - (d/2) ln T - (1/2) ln det Ibar - ln|Theta_0|` (`p9_analyze.py` (B)): at `T eps = 30`: -1.58, -1.78 (mu = 0.3, k0 = 0.6, eps = 0.1, 0.03), -1.68, -1.83 (mu = 0.6), -1.62, -1.83, -1.92, -1.92, -1.99 (mu = 1, eps = 0.3...0.003), -1.86, -1.90 (mu = 2); at `T eps = 100`: -2.07 to -2.47. Flat in eps and in mu (including k0 < 1, where `(1/2) ln det Ibar` measured is 0.59 and 0.94 and **is needed to make c flat**: the term `((1-k0)/2) ln(1/eps)` is real). The expected asymptote is `-ln(2 pi e) = -2.84` plus the Newton excess (0.05-0.4): the drift from -1.6...-2.0 (x = 30) to -2.1...-2.5 (x = 100) is the approach to it.

**Rigorous-form bound `L_N`** (Prop 9.1 with `A = N`, the tangent minorant, `kappa_N` measured; `p9_analyze.py` (C)) at the largest T: `E R`, `L_emp`, `L_N` = 3.364, 2.935, 2.579 (eps = 0.3); 2.400, 2.082, 1.572 (0.1); 2.247, 2.195, 1.674 (0.03); 1.625, 1.383, 0.525 (0.01); 0.843, 0.606, -0.154 (0.003). Ordering `L_N <= L_emp <= E R` in all 5; `E R - L_N` = 0.78, 0.83, 0.57, 1.10, 1.00 (the part `(1/2) ln(det I/det N)` = 0.50, 0.75, 0.90, 0.97, 1.05 saturates in eps; `ln det N` = -1.40, -1.63, -1.71, -1.73, -1.74 against `ln det I` = -0.40, -0.13, 0.08, 0.21, 0.36). The rigorous lower bound is therefore within about 1 of the regret uniformly in eps over a 100-fold range.

**T7 — two-point (Prop 9.0), valid configurations** (`tpv_*`, alpha-perturbation `delta <= 0.4 eps`, so theta' stays stable): ratio `E R_T/(T u^T I u/2)` is 0.965-1.000 for the mu-direction and the ridge, 0.92-0.998 for the alpha-direction, tending to 1 as delta shrinks (`delta = 0.1 eps`: 0.993-1.000), and **always inside the proved bracket** `[1/c~^+, 1/c~^{-2}]` (e.g. [0.893, 1.0] at eps = 0.3, delta = 0.4 eps). (Configurations with `delta > eps` in the alpha-direction put theta' outside the stable set, the code clamps it, and the ratio drops to 0.04-0.34; those rows are outside the hypothesis of Prop 9.0 and were discarded.)

## 9. The comparison table, filled (K = 1, mu = beta = 1, k0 = 2)

Regret over the stability band, `x = T eps`. "Upper" = Newton with matched prior (measured); "L_emp" = Prop 9.1 with the measured covariance; "L_N" = Prop 9.1 with the tangent minorant N; formulas in sections 4 and 6.

| eps | T | x = T eps | upper (E R, Newton) | L_emp | L_N | upper - L_N | formula (9.2)/(9.3) constant `c` |
|---|---|---|---|---|---|---|---|
| 0.3 | 3000 | 900 | 3.364 | 2.935 | 2.579 | 0.78 | `c = -2.84 + (0 to 0.4)` |
| 0.1 | 3000 | 300 | 2.400 | 2.082 | 1.572 | 0.83 | " |
| 0.03 | 9487 | 285 | 2.247 | 2.195 | 1.674 | 0.57 | " |
| 0.01 | 9487 | 95 | 1.625 | 1.383 | 0.525 | 1.10 | " |
| 0.003 | 9487 | 28 | 0.843 | 0.606 | -0.154 | 1.00 | " |

mu-dependence (eps = 0.1 / 0.03, `T eps = 30`, `c` as above): mu = 0.3 (k0 = 0.6): -1.58 / -1.78; 0.6 (1.2): -1.68 / -1.83; 1 (2): -1.83 / -1.92; 2 (4): -1.86 / -1.90.
T-dependence: slope of `E R` in `ln T` is 0.48 at `x = 100` and rises toward `d/2 = 1` (I_G slope); the regret at fixed `x` does not depend on `T` separately.

## 10. Net result for Problem 9

1. **Statements now proved or reduced to the Theorem R inputs:** Prop 9.0 (two-point in Fisher form, exponential dictionaries, direction-resolved; sandwich verified), Prop 9.1 (matrix-covariance lower bound), Prop 9.2/Theorem 9.3 (lower bound without the hypothesis (U), `(d/2) ln T + (1/2) ln det N + ln|Theta_0| - (d/2) ln(2 pi e kappa_T)`).
2. **The eps-dependence:** upper and lower bound have the same eps-structure `ln(T eps) + const` for `k0 > 1` (the stability band contributes `ln eps` through `ln|Theta_0|`, `ln det I` and `ln det N` are eps-flat); for `k0 < 1` the upper bound has the extra `((1-k0)/2) ln(1/eps)` (needed for flat `c`), the tangent minorant does not see it [O]. **For `eps = eps(n) -> 0` with `T eps_T -> infinity`: `sup Reg_T = ln(T eps_T) + O(1)` up to the lower-side `ln ln T`**; for `T eps = O(1)` the Bayes regret is `O(1)` (0.17 at x = 3), the non-learnable regime of Theorem NC. Empirically the lower bound with the measured covariance is within 0.4 of the regret.
3. **Corrections to earlier statements:** the "(d-1)/2 plateau" prediction (section 5, T6) does not apply to the natural band (collapse in `T eps`); Prop LB's (U) is not an extra hypothesis (it is delivered by the Theorem R chain); the `(d/2) ln ln T` loss is a proof artefact (measured `kappa_eff ≈ 1`).
4. **Open:** the `ln ln T` removal (needs `E D_n = O(d)`, i.e. `E X_n = E G_n + O(1)`); a minorant of `1/lam` that captures `P_pi(lam <= 2 mu) ≍ eps^{k0}` to make the lower bound see the `((1-k0)/2) ln(1/eps)` term; K >= 2 numerics (the propositions are K-general, only K = 1 was simulated here); a certification of the ridge term `(1/2) ln(1 + T lam_min sigma_r^2)` in the lower bound for fixed (not shrinking) prior width at `T eps = O(1)`; uniformity of the constants of Theorem 9.3 over `Theta_0` was stated, not written out.

| "Done when" (tracker) | status |
|---|---|
| a lower-bound theorem | **Theorem 9.3 / Prop 9.1 [P\*]**, unconditional given the Theorem R inputs, with `ln det N` in the constant; `latex/part2b_extensions.tex` updated |
| a comparison table showing where upper and lower bounds match or differ in eps, mu and T | **section 9** (numerical) and section 6 (formulas): match up to `(1/2) ln ln T + O(1)` for `k0 > 1`; differ by `((1-k0)/2) ln(1/eps)` for `k0 < 1` (rigorous side) |
