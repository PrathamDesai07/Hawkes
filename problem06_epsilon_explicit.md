# Problem 6: making the epsilon-dependence of the regret bound explicit

Tags: **[P]** derived here from earlier proved statements, **[H]** heuristic derivation (a prediction to test), **[O]** open. Numerical results are in the last section (filled after the runs). Code: `code/p6_scaling.cs`, `code/p6_run.ps1`.
Dates: 25 Oct 2026. Nothing here is refereed.

## 0. What the tracker asks, and what is already there

Target: `E R_T <= (d/2) ln T + C ln(1/eps) + burn-in(mu eps)`, d = K+1, every constant a function of (eps, mu, K, beta-grid).
Already there (Theorem R, `theoremR_proof_complete.md`): `E R_T <= (1/2+O(delta+eta)) d ln T + C`, C = O(mu d ln(1/(eps delta)) / (beta eps^2 delta^2))`, an order only, plus the remark that the eps^-2 is probably a proof artefact.
This note (a) separates the bound into pieces and tracks the eps-power of each piece, (b) identifies which pieces are real and which are artefacts, (c) states predictions that the simulation tests.

## 1. Units: events versus time [P]

The stationary intensity has mean `E lam = mu/eps`, so n events take about `T = eps n/mu` time. Hence
 `(d/2) ln n = (d/2) ln T + (d/2) ln(mu/eps)`.
A statement `E R_n = (d/2) ln n + c(eps)` in events is the statement `E R_T = (d/2) ln T + c'(eps)` in time with `c(eps) = c'(eps) - (d/2) ln(mu/eps)`. So if the time-unit constant is flat in eps (it is in the earlier runs), the event-unit constant must **decrease like (d/2) ln eps**. The "C ln(1/eps)" of the target form is therefore not a cost in time units; it is the conversion between the two clocks. All statements below are in time units unless said otherwise.

## 2. Master inequality with every constant kept [P]

From `theoremR_proof_complete.md` section 2 (K arbitrary), with `A := (1/2 + C_2 delta)/(1 - delta)`:

 `E R_T <= A [ P_0 + (1+eta)(1+delta)^2 E ln det(F_T/F_0) + (1+eta) C_rho + C_B / lam_min(F_0) ] + C_B`,

valid once `B := (1+1/eta) lbar(t_*) / mu_lo <= delta`. The pieces:

| piece | meaning | eps-dependence (derived below) |
|---|---|---|
| P_0 = ‖theta_1 - theta*‖^2_{F_0} | start/prior cost | depends on how F_0 is scaled, section 3 |
| E ln det(F_T/F_0) | information term, carries the leading constant | section 4 |
| t_* | burn-in (leverage small) | section 5 |
| C_rho | truncation remainder | section 6 |
| C_B | bad gaps (probability e^{-c1 gamma delta^2 t}) | section 6 |

## 3. The prior term [P]

Take `F_0 = rho_0 diag(1/lam_bar, lam_bar/beta_1^2, ..., lam_bar/beta_K^2)` with `lam_bar = mu/eps` (the scaling used in the code). Then
 `det F_0 = rho_0^{K+1} / prod_k beta_k^2` (K = 1: `rho_0^2/beta^2`), **free of eps**,
 `P_0 = rho_0 [ eps (Dmu)^2/mu + (mu/eps) sum_k (D alpha_k / beta_k)^2 ]`.
So `-1/2 ln det F_0` costs nothing in eps, but a start with an O(1) relative error in alpha_k costs `P_0 ~ rho_0 mu/eps` in the **bound** (Newton pseudo-regret bound `R <= X`). The measured regret from a far start is much smaller (the bound `psi(r) <= (r-1)^2/r` is loose for ratios far from 1; log-loss grows linearly there). [H] The true start cost is `O(ln(1/eps))`, not `1/eps`; the simulation checks this.

## 4. The information term with explicit eps [P] for K = 1, [H] for K >= 2

Concavity of ln det: `E ln det(F_T) <= ln det(F_0 + E[F_T - F_0]) <= ln det(F_0 + T Ibar_T / c^-)`, `Ibar_T = E J_T / T -> I(theta*)`.

**K = 1, closed form** (`K1_main_results.md` A1, A2): `m = E[1/lam] = 2eps/(2mu-beta) (1+o(1))` for `k0 = 2mu/beta > 1`, and
 I_mumu = m, I_mualpha = (1 - mu m)/alpha, I_alphaalpha = (E lam - 2 mu + mu^2 m)/alpha^2, E lam = mu/eps, alpha = (1-eps) beta.
Leading orders as eps -> 0:
 `I_mumu ~ 2 eps/(2mu-beta)`, `I_mualpha ~ 1/beta`, `I_alphaalpha ~ mu/(eps beta^2)`,
 `det I ~ 1/(beta(2mu-beta)) = 1/(beta^2 (k0-1))`  (bounded in eps),
 `lam_max(I) ~ mu/(eps beta^2)`, `lam_min(I) ~ det/tr = eps/(mu (k0-1))`  (the ridge, units: time).
Therefore
 `1/2 ln det(F_0 + T I) - 1/2 ln det F_0 = (d/2) ln T + 1/2 ln det I - 1/2 ln det F_0` only for `T >> f_ridge/lam_min(I)`; before that the ridge contributes `1/2 ln(1 + T lam_min/f_ridge)` instead of `1/2 ln T`. In tier-resolved form, with `f_j` the eigenvalue of `F_0` in the j-th eigendirection of I:

 **`E R_T ~ 1/2 sum_j ln(1 + T lam_j(I)/f_j)`**   (prediction [H], used as the comparison curve below).

For the F_0 above, `f_ridge ~ rho_0 eps/mu` along mu, so the ridge saturates when `T eps/(mu(k0-1)) ~ rho_0 eps/mu`, i.e. at `T ~ rho_0 (k0-1)`, **independent of eps**: with an eps-scaled prior the log-det curve reaches slope d/2 at O(1) time. With an eps-independent prior `F_0 = rho_0 diag(1,1/beta^2)` the saturation is at `T ~ rho_0 mu (k0-1)/eps` (time), i.e. `n ~ rho_0 mu^2 (k0-1)/eps^2` events. This is the origin of the eps^-2 events: it is the price of a prior that does not know the scale `mu/eps`.

**K >= 2** (`problem12_upper_bounds.md`, Thm 4-6 + tangent bound): eigenvalues `Theta(1/eps)` (top), `[lam_j(N), C_mid]` eps-independent (K-1 middle), ridge `Theta(eps)`; so `ln det I` is bounded in eps and the saturation times are `~ eps`, `~ 1`, `~ 1/eps` for the three tiers. The K-1 middle directions saturate at O(1) time, the top at `O(eps)`; the ridge at `1/lam_ridge ~ 1/eps` (fixed prior).

## 5. The burn-in t_* [P] K = 1, [H] K >= 2

Needed: `lbar(t_*) <= delta mu_lo / (1+1/eta)`. For K = 1 and `F ~ t I / c^+`, `lbar(t) = sup_u phi^T F^{-1} phi / lam*` with `phi = (1, y)`, `y in [0,S]`:
 `phi^T I^{-1} phi / lam* = (I_aa - 2 y I_ma + y^2 I_mm) / (det I (mu + alpha y))`.
At y = 0: `I_aa/(mu det I) = (mu/(eps beta^2)) beta (2mu-beta)/mu = (k0-1)/eps`. At y = S = mu/(eps beta): ratio `~ 1`. So the sup is at y = 0 (the rare state S ~ 0) and
 **`lbar(t) = c^+ (k0-1)/(eps t) (1+o(1))`  [P]**, hence
 **`t_* = (k0-1)(1+1/eta) c^+ / (eps delta mu_lo)`  (time)**, i.e. `~ (k0-1)/(mu eps delta)`.
This is the "1/(mu eps)" of the tracker's guess, with the coefficient `(k0-1)/delta`: the burn-in in **time** is `(k0-1)/(mu eps delta)`, in **events** `n_* = (mu/eps) t_* = (k0-1)/(eps^2 delta)` (mu cancels in events). The earlier note's "1/(beta eps)" agrees when mu ~ beta. For K >= 2 the proved leverage has an extra `(ln t)^2/eps` (K2 chain section 2.3): `t_* ~ (ln)^2/(eps^2 ...)`; the K >= 2 leverage from the tangent minorant is [O].
Mixing also needs `gamma delta^2 t >> 1`: `t_mix = ln(1/q)/(c_1 beta eps delta^2)`; the larger of `t_*` and `t_mix` governs.

## 6. Where the eps^-2 in C comes from, and the improved bookkeeping [P/H]

(a) *C_rho.* The remainder per gap is `rho_i = 3 lbar_i^3 ln(e/lbar_i)/mu*` with `m_i = 2 ln(1/lbar_i)`. Summed over events after `t_*` (rate `mu/eps`): `C_rho ~ (mu/eps) a^3/(2 t_*^2)` with `a = c^+(k0-1)/eps` -- order `mu delta^2/eps^2`. **Fix [P]:** take `m_i = p ln(1/lbar_i)`; the remainder becomes `lbar^{p+1}` and `C_rho ~ (mu/eps^2) (delta mu_lo)^p`. With `p = 2 ln(1/eps)/ln(1/(delta mu_lo))` this is `O(mu)` at the price `x̄ = lbar m = O(lbar ln(1/eps))`, i.e. `t_*` larger by a factor `ln(1/eps)`. So C_rho = O(1) and `t_* = O((k0-1) ln(1/eps) / (mu eps delta))`.
(b) *C_B.* `sum_i P(g_i = 0)^{1/2} ~ (mu/eps) int_{t_*}^inf e^{-c_1 gamma delta^2 t/2} dt = 2 mu e^{-c_1 gamma delta^2 t_*/2}/(c_1 beta eps^2 delta^2)`. This is `O(1)` iff `c_1 gamma delta^2 t_* >= 2 ln(1/eps)+...`, i.e. iff `t_* >= t_mix` (above). With `t_* = max(t_lev, t_mix) = O(ln(1/eps)/(beta eps delta^2))` the sum is O(mu/(c_1 beta delta^2)) -- **free of eps**.
(c) *Early gaps.* The proof counts the `(mu/eps) t_*` early gaps at O(1) each: `(mu/eps) t_* = O(mu ln(1/eps)/(beta eps^2 delta^2))`. This is the remaining eps^-2 and it **is** an artefact of counting events instead of using the potential: during `[0, t_*]` the regret is bounded by the Newton bound `P_0 + E G` and `E G <= (1+eta) E sum Q_i + ...` with `Q_i` dominated by the log-det increments *if* the typical (not sup) leverage is used. [H] With typical leverage `phi^T F^-1 phi/lam ~ 1/t` (checked above at y = S) the early phase costs `O(ln(det ratio)) = O(d ln(1/eps))`. Making this rigorous needs the "leverage from occupation" bound (`code/p5_leverage_from_occupation.py`) in the potential argument: [O].

**Conjectured sharp form [H]:**
 `E R_T <= (d/2 + o(1)) ln T + 1/2 ln det I(theta*) + 1/2 P_0 - 1/2 ln det F_0 + C_1 d ln(1/eps) + C_2`   for `T >= t_* = O((k0-1) ln(1/eps)/(mu eps delta))`,
with `ln det I` bounded in eps, `det F_0` free of eps for the scaled prior, and the only unbounded term the `ln(1/eps)` from the burn-in potential. Compare with the earlier `O(mu d ln(1/(eps delta))/(beta eps^2 delta^2))`: the factor `eps^-2 delta^-2` disappears.

## 7. Predictions to test [H]

- **P1.** In time units `E R_T - (d/2) ln T` is flat in eps for `T >= T_sat`; in event units `E R_n - (d/2) ln n` decreases like `(d/2) ln eps` (K=1: slope +1 against ln eps).
- **P2.** `E R_T ~ 1/2 ln det(F_T/F_0)` and `~ 1/2 sum_j ln(1 + T lam_j(I)/f_j)` (tier-resolved curve) at all T, K = 1, 2.
- **P3.** Saturation time (slope reaches d/2): `O(1)` for the eps-scaled prior, `~ rho_0 mu (k0-1)/eps` for the fixed prior. So the burn-in is a statement about the prior scale, not about the algorithm.
- **P4.** Regret accrued on `[0, t_*]` grows like `ln(1/eps)` (not eps^-2) for a far start.
- **P5.** `lam_min(I) = eps/(mu(k0-1))`: measured from the accumulated information `J_T/T`.
- **P6.** mu-dependence: with fixed eps, the fixed-prior saturation time is proportional to `mu (k0-1) = mu(2mu-beta)/beta`.
- **P7.** K = 2: slopes approach 1.5 and ln det I stays bounded as eps decreases.

## 8. Status of the problem after this note

| item | status |
|---|---|
| explicit t_*, K = 1 | [P] `(k0-1)(1+1/eta)c^+/(eps delta mu_lo)`; with remainder fix `x ln(1/eps)` |
| explicit information term, K = 1 | [P] up to `o(1)` in eps (Theorem A, `m` asymptotics) |
| removing eps^-2 from C | [P] for C_rho, C_B; [O] for the early-gap count (needs typical-leverage potential) |
| K >= 2 explicit constants | [O] (three tiers known up to constants; leverage `(ln t)^2/eps` worse) |
| simulation of P1-P7 | see section 9 |

## 9. Numerical results (25 Oct 2026; `code/p6_scaling.cs`, `code/p6_run.sh`, `code/p6_analyze.py`; raw tables in `code/p6_out/`, `code/p6_out_long/`)

Setup: K = 1, mu = beta = 1 (k0 = 2) unless stated, alpha* = (1-eps) beta, stationary start (warm-up 20/eps^2 events), Newton recursion with clip (not the F-projection), 32 replicates to T = 3000, 16 replicates to T = 3e4 (long runs). K = 2: beta = (1,2), mu = 1.5, near start (1.3 mu*, 0.8 alpha*), 24 replicates. "scaled" prior F_0 = rho_0 diag(1/lam_bar, lam_bar/beta^2), "fixed" prior F_0 = rho_0 diag(1, 1/beta^2), rho_0 = 1. Caveat: the simulator is C# (no Python at first; Python exists at `C:/msys64/ucrt64/bin/python.exe` and is used for the analysis).

**P1 (flat constant in time units) — supported.** Truth start, scaled prior, T = 1e4 (16 replicates, s.e. 0.5-1.1): `E R_T - ln T` = -2.04, -1.71, -0.84, -1.33 for eps = 0.3, 0.1, 0.03, 0.01. Over a 30-fold range of eps the constant stays within [-2.2, -0.8]; an eps^-2 law would span a factor 1000. Local slopes at T = 1e4 are 0.87, 1.06, 1.31, 0.83 (noisy), at the target d/2 = 1. In event units the constant is `c_time - ln(mu/eps)`, i.e. falls by `ln 30 = 3.4` over the same range (check: table column c_events). So **the "C ln(1/eps)" of the target form is the clock conversion of section 1, not a regret cost.**

**Transient overshoot.** For small eps the regret is below the asymptote early and then rises faster: local slope at T = 1000 is 1.58 (eps = 0.03) and 1.95 (eps = 0.01) against 1.0 (eps = 0.3). This is the ridge direction switching on; it does not change the final constant but is the visible burn-in. [N]

**P2 (regret vs log det) — partly.** `E R_T` stays below `(1/2) E ln det(F_T/F_0)` and the tier-resolved curve `(1/2) sum_j ln(1 + T lam_j/f_j)` (upper envelope) and its slope matches for eps >= 0.1 (slopes 1.0, 1.0). For eps <= 0.03 with the fixed prior and truth start the slope of `E R_T` over T in [100, 3000] is 0.48-0.70 against 0.85-0.9 for the curve: the regret lags the information. So `E R_T ≈ (1/2) ln det` is an upper bound, not a sharp description at small eps before the lag closes. For K = 2 the curve and regret slopes are 1.4-1.5 vs 0.9-1.8 (noisy at 24 replicates).

**P3 (saturation time of the prior) — refuted as stated.** Predicted slope d/2 from `T ~ rho_0 mu (k0-1)/eps` (fixed prior), i.e. T ~ 100 at eps = 0.01. Measured slope over [100, 3000] is 0.48-0.75 there. The saturation is later than the information crossover because the estimator error needs `T >> 1/(gamma delta^2)` (mixing), not only `T lam_min >> f`. The mixing time `1/(beta eps)` = 100 at eps = 0.01 times a factor of order 10 reproduces the observed time; the section 5 statement `t_* = O(ln(1/eps)/(beta eps delta^2))` (max of leverage and mixing times) is consistent with this, the information-only crossover is not enough.

**P4 (start cost is O(ln 1/eps)) — refuted; the bound P_0 is sharp for the scaled prior.** Excess of the far start (mu_0 = 2 mu*, alpha_0 = 0.5 beta) over the truth start at T = 3000, scaled prior, against `(1/2) P_0 = (rho_0/2)[eps (Dmu)^2/mu + (mu/eps)(Dalpha/beta)^2]`:

| eps | measured excess | (1/2) P_0 |
|---|---|---|
| 0.3 | 0.9 | 0.22 |
| 0.1 | 1.5 | 0.85 |
| 0.03 | 4.1 | 3.7 |
| 0.01 | 12.5 | 12.0 |

For small eps the start cost is exactly `(1/2) P_0 ∝ 1/eps` (agreement to 4%), so a far start with a prior scaled to the natural information level costs `1/(2 eps)`. With the fixed prior the excess is 0.8 at eps = 0.01 against `(1/2) P_0 = 0.6`. Conclusion: the eps-dependence of the additive constant is **entirely through P_0 = ‖theta_1 - theta*‖^2_{F_0}**; the P_0 term in the master inequality is the real source, with `‖.‖_{F_0}` weighting the alpha error by `lam_bar/beta^2 = mu/(eps beta^2)`. It is not removable by the method, but it is also not eps^-2: it is `rho_0 mu (Dalpha/beta)^2/(2 eps)` and vanishes with a prior of weight `rho_0 ∝ eps` (then the far-start cost is O(1) in eps but the early wandering grows, see `box_cost_and_lower_bound.md` A2).

**P5 (lam_min(I) = eps/(mu(k0-1))) — right order, biased by 35-55 %.** Estimated from accumulated oracle information `J_T/T` (tail average): 0.14, 0.066, 0.025, 0.0064 for eps = 0.3, 0.1, 0.03, 0.01 against 0.30, 0.10, 0.030, 0.010. The ratio is 0.46, 0.66, 0.82, 0.64: linear in eps (log-log slope 1.0 of the data) but the leading-order formula overestimates and the estimate at eps = 0.01 is noisy (mean of 1/lam over 30 mixing times). `det I` measured 0.51, 0.69, 0.94, 0.72 against the leading-order `1/(beta^2 (k0-1)) = 1`: **det I is bounded and O(1) in eps** (also noisy at eps = 0.01), as claimed. Top eigenvalue times eps = 1.10, 1.06, 1.14, 1.12 (predicted 1.0 + O(eps)): confirms `lam_max = Theta(1/eps)`.

**P6 (mu-dependence) — inconclusive.** mu = 0.6 (k0 = 1.2), 2 (k0 = 4), fixed prior, truth start: measured lam_min 0.20 (mu = 0.6, eps = 0.1), 0.0145 (mu = 2, eps = 0.1) against leading-order 0.83 and 0.0167. Agreement for k0 = 4 (87 %), poor for k0 = 1.2 where `m = E[1/lam] ≈ 2 eps/(2mu - beta)` is not accurate (corrections `eps^{min(1,k0-1/2)}`). Slopes at T = 3000 are not at saturation (0.58-1.24).

**P7 (K = 2) — consistent.** Top eigenvalue times eps = 1.94-2.0 (= K-fold, since two exponentials share the mean intensity); lam_min(J/T) = 0.032, 0.030, 0.0125 for eps = 0.3, 0.1, 0.03: not yet in the O(eps) regime at these T. Slopes of `E R_T` over the last 1.5 decades: 1.79, 1.44, 1.01 (fixed prior) against d/2 = 1.5: below target at eps = 0.03 (not saturated at T = 3000).

## 10. Net result for Problem 6

1. **Resolved (explicit):** the dictionary events <-> time (section 1); explicit `lbar(t)` and `t_*` for K = 1 (section 5); the mechanism of every eps-power in the proof's C (section 6): `C_rho` and `C_B` can be made eps-free by taking truncation level `p ln(1/lbar)` with `p ∝ ln(1/eps)` and `t_* >= t_mix`; the start term is exactly `(1/2) P_0` (measured, section 9 P4).
2. **Numerically supported:** `E R_T - (d/2) ln T` is O(1) uniformly in eps in time units (section 9 P1), so the true additive constant is eps-flat, not eps^-2 ln(1/eps); `det I` bounded, `lam_max ∝ 1/eps`, `lam_min ∝ eps`.
3. **Not proved:** the early-gap count in C (section 6c) still gives eps^-2 ln(1/eps) in the written proof; the conjectured sharp form of section 6 needs the typical-leverage potential argument. Not covered: K >= 2 explicit constants, k0 <= 1.
4. **Corrections to earlier notes:** burn-in in time is `(k0-1)/(mu eps delta)` for leverage (not `1/(mu eps)` without the factor, not `1/(beta eps)` in general) plus `ln(1/eps)/(beta eps delta^2)` for mixing; P3 (information-only saturation) and P4 (O(ln 1/eps) start cost) were wrong as first stated.

---
# Completion pass (25 Oct 2026, second session on Problem 6)

Code: `code/p6_scaling.cs` (extended: potential defect, count of gaps with x>1, exit time `T_exit`), `code/p6_run2.sh`, `code/p6_analyze2.py`; raw tables `code/p6_out2/`.

## 11. Removing the event count from C: the potential-defect lemma [P], and its size [N]

**Lemma 6.4 (leverage-free conversion).** For a PSD increment `DF_i` and `F = F_{i-1} ≻ 0`, put `X_i = F^{-1/2} DF_i F^{-1/2}` with eigenvalues `x_j ≥ 0`. Then, for every gap, with no condition on leverage,
 `tr X_i = ln det(F_i/F_{i-1}) + def_i`,  `def_i := sum_j [x_j - ln(1+x_j)]  <=  (tr X_i)^2`.
*Proof.* `ln det(F_i/F_{i-1}) = sum_j ln(1+x_j)` exactly. For `x <= 1`: `x - ln(1+x) <= x^2/2`; for `x > 1`: `x - ln(1+x) <= x <= x^2`. Sum and use `sum x_j^2 <= (sum x_j)^2`. ∎

*Consequence.* The score part of the proof (F4) needs `q_i = tr(F_{i-1}^{-1} E DJ_i) <= c^+ E tr X_i`. By Lemma 6.4 this is `<= c^+ ( E ln det(F_i/F_{i-1}) + E def_i )`. So the sharp conversion `q_i <= (1+delta) r(xbar) E ln det`, with its leverage condition `xbar -> 0`, is not needed in the early phase: one pays `c^+` instead of `(1+delta)` there, and the price is `D_early := E sum_i def_i <= E sum_i (tr X_i)^2`. This replaces the event count `(mu/eps) t_*` of the written proof (section 6c) by a second moment of the gap trace.

**Size of D_early [H].** `tr X_i = int_0^tau phi^T F_{i-1}^{-1} phi / lam_theta du`. Once `F_{i-1} ⪰ F_0 + J_t/c^+` dominates the prior (time `t_J`), a typical gap has `tr X_i ~ d tau_i / t_{i-1}` and `E sum_i (tr X_i)^2 ~ (2 d^2/lam_bar) int_{t_J}^inf t^-2 dt = 2 d^2 eps/(mu t_J)`, which is tiny. Before `t_J`, crudely `tr X_i <= O(1)` per gap, giving `(mu/eps) t_J` events. For the **scaled prior** `t_J = O(rho_0)` (J dominates F_0 in every direction after O(1) time), so the crude count is `O(mu rho_0/eps)`, one power of eps better than the written proof. For the **fixed prior** `t_J = rho_0 mu (k0-1)/eps` and the crude count is `eps^-2`: the written proof's count is exactly the fixed-prior, scale-unaware case.

**Measured [N]** (`code/p6_out2`, K = 1, mu = beta = 1, truth start, T = 6000; `ln det` = ln det(F_T/F_0)):

| eps | prior | events | ln det | E sum def_i | E #(gaps with tr X > 1) |
|---|---|---|---|---|---|
| 0.3 | scaled | 1.9e4 | 16.6 | 1.07 | 0.7 |
| 0.1 | scaled | 5.6e4 | 16.8 | 0.61 | 0.3 |
| 0.03 | scaled | 1.8e5 | 17.0 | 0.17 | 0.0 |
| 0.01 | scaled | 5.7e5 | 17.1 | 0.05 | 0.0 |
| 0.3 | fixed | 1.9e4 | 16.6 | 1.38 | 1.0 |
| 0.1 | fixed | 5.6e4 | 16.8 | 1.10 | 0.7 |
| 0.03 | fixed | 1.8e5 | 17.0 | 1.20 | 0.8 |
| 0.01 | fixed | 5.7e5 | 17.1 | 1.19 | 0.8 |

The defect is O(1). For the scaled prior it decreases with eps (1.07 -> 0.05); for the fixed prior it is flat at about 1.2. The event count grows 30-fold. So the early-phase cost is `c^+ ln det(F_{t_*}/F_0) + O(1)`, not `(mu/eps) t_*`.

**Early term, explicit.** `ln det(F_{t_*}/F_0) <= d ln t_* + ln det(I/c^-) - ln det F_0 + o(1)` with `t_* = O(ln(1/eps)/(beta eps delta^2))` (section 12), so the early term is `c^+ (d ln(1/eps) + ...)`. This is the **C ln(1/eps)** of the tracker's target form, here with coefficient `c^+ d (1+eta)`, in time units, from the early phase (rather than from the clock conversion of section 1).

**What is still not proved.**
(i) The moment bound on `E (tr X_i)^2` along the real state process. It needs the occupation-measure leverage bound; the measured defect supports it.
(ii) The regret of **bad** early gaps: `R_i <= C_3 X_i` with `X_i <= D_{i-1} lam_max(X_i) <= D_{i-1} tr X_i`, so `E sum_{bad,early} R_i <= C_3 (sup D)(ln det + D_early)`. This is `O(d^2 ln^2(1/eps))` if `sup D = O(d ln(1/eps))` on the early window, but Theorem C-exp covers only t >= T_1.
These two items are the remaining gap between the written proof and the numerically observed `O(1)`.

## 12. Burn-in corrected by measurement [N]

`T_exit(delta)` is the last time the iterate is outside `R(delta)` (delta = 0.25). 32 replicates, scaled prior, K = 1, mu = beta = 1.

| eps | 0.3 | 0.1 | 0.03 | 0.01 |
|---|---|---|---|---|
| T_exit median | 123 | 214 | 1207 | 3215 |
| p90 | 496 | 855 | 2968 | 5539 |
| `1/(beta eps delta^2)` | 53 | 160 | 533 | 1600 |
| ratio | 2.3 | 1.3 | 2.3 | 2.0 |

So **T_exit ≈ 2/(beta eps delta^2)**: proportional to `1/eps`, with constant between 1.3 and 2.3 (mixing-time scaling `1/gamma`, `gamma = beta eps`). The information-based prediction `(k0-1)/(mu eps delta^2)` is *not* what limits the burn-in. At eps = 0.1 and `mu = 0.3, 0.5, 0.6, 1, 1.5, 2, 4` (k0 = 0.6 ... 8) the medians are 296, 320, 294, 214, 306, 326, 440: nearly flat in mu and not divergent at k0 <= 1, where the information-based formula is undefined.

**Corrected burn-in:** `t_* = max(t_lev, t_mix)` with `t_mix = c_mix/(beta eps delta^2)` (c_mix ≈ 2 empirically) and `t_lev = (k0-1)(1+1/eta)c^+/(eps delta mu_lo)` from section 5. `t_lev` has the same 1/eps order with a coefficient (k0-1)/(mu delta) that is not resolved by the exit-time data. The section 5 statement is a correct bound for the *leverage* but not the binding constraint.

## 13. Other regimes [N]

- **k0 <= 1 (mu <= beta/2).** mu = 0.3, 0.5 (k0 = 0.6, 1.0): at T = 6000, `E R_T - ln T` = +0.96 / +0.15 (mu = 0.3, eps = 0.1 / 0.03) and -0.91 / -0.80 (mu = 0.5); local slope 0.98-1.17. The log slope is still `d/2 = 1` and the constant is bounded in eps, but it increases as k0 drops (k0 = 3: -2.7; k0 = 1.2: -1.1; k0 = 0.6: +1.0), as expected from `E[1/lam]` no longer being O(eps). Not proved; outside Theorem A.
- **K = 2** (beta = (1,2), mu = 1.5, T = 3e4, 32 replicates): `E R_T - 1.5 ln T` = -5.97, -6.08, -5.74 for eps = 0.3, 0.1, 0.03, **flat in eps** (spread 0.3, s.e. 0.6). Local slope 1.18, 1.27, 1.33 against 1.5: still below target at T = 3e4 (the ridge is slow), so the constant is meaningful but the slope is not yet asymptotic. Defect 1.03, 0.60, 0.16: O(1) and decreasing.
- **K = 3** (beta = (1,2,4), mu = 2, 24 replicates, T = 3e4): constants -0.6 (s.e. 2.9), -7.3, -10.1 with slope 1.65 against 2. **Not asymptotic** at this horizon (the extra middle tier has a smaller eigenvalue; `lam_min(M)` collapses geometrically in K). No conclusion about eps-flatness beyond K = 2.

## 14. The eps-explicit statement (what can now be written) and status

**Theorem 6 (K = 1, k0 > 1; constants explicit; one hypothesis).** Let `A = (1/2 + C_2 delta)/(1-delta)` and `t_* = max(t_lev, t_mix)` as above. For `T >= t_*`:

 `E R_T <= A { P_0 + (1+eta)(1+delta)^2 [ 2 ln T + ln det(I/c^-) - ln det F_0 ] + (1+eta) c^+ [ 2 ln t_* + ln det(I/c^-) - ln det F_0 ] + (1+eta)(D_early + C_rho) + C_B/lam_min(F_0) } + C_B'`,

with, for the scaled prior and `d = 2`:
- `ln det F_0 = 2 ln rho_0 - 2 ln beta`;
- `ln det I = -2 ln beta - ln(k0-1) + O(eps^min(1,k0-1/2))` (Theorem A; measured `beta^2 det I` = 0.51-0.94 against 1.0, noisy);
- `P_0 = rho_0 [eps (D mu)^2/mu + (mu/eps)(D alpha/beta)^2]` (exact; measured to 4 % for small eps);
- `t_* = O(ln(1/eps)/(beta eps delta^2))`; `C_rho = O(1)` and `C_B = O(mu/(c_1 beta delta^2))` for truncation level `p ∝ ln(1/eps)` and `t_* >= t_mix` (section 6);
- hypothesis **(D)**: `D_early = O(d^2 ln^2(1/eps))` (the two items of section 11, "what is still not proved").

Leading terms: `E R_T <= (1 + O(delta+eta)) { ln T - ln rho_0 - (1/2) ln(k0-1) } + (P_0/2)(1+...) + O(c^+ ln(1/eps))`. Against the measurement (truth start, scaled prior, rho_0 = 1, k0 = 2, where the bracket's eps-free part is 0), the measured `E R_T - ln T` is -2.9 to -1.2 (eps = 0.3, 0.1, 0.03, 0.01: -2.2, -1.5, -1.2, -2.9). The bound is valid and its eps-free part is loose by about 1-3 (the Gaussian-mixture constant `-d/2 = -1` accounts for one unit).

| item (Problem 6 "Done when") | status |
|---|---|
| every constant a function of (eps, mu, K, grid): K = 1, k0 > 1 | **done up to hypothesis (D)** and the conservative universal constants c_1, C_1 of (H3) (explicit in `thmT_explicit_constants.md` and the K2 chain) |
| burn-in: what is really needed | **measured**: `t_* ≈ 2/(beta eps delta^2)` (mixing-limited), not `1/(mu eps)` |
| `C log(1/eps)` form | two sources identified: clock conversion `(d/2) ln(mu/eps)` (events) and early-phase potential `c^+ d ln t_*` (time) |
| removal of eps^-2 | **numerically** yes (flat constant, O(1) defect); **rigorously** only for C_rho, C_B; (D) open |
| K >= 2 | same statement with d = K+1 and the three-tier `ln det I`; constants not explicit; K = 2 constant flat in eps (numerics); K = 3 not asymptotic |
| k0 <= 1 | slope d/2 holds numerically, constant larger; no proof |

---
# Resolution pass (25 Oct 2026, third session on Problem 6): hypothesis (D), K >= 2 constants, k0 <= 1

Code: `code/p6_moments.cs` (stationary negative moments and I(theta*)), `code/p6_scaling.cs` (now also E G and E X), `code/p6_chat.py` (leverage-scale constant), `code/p6_tangent_det.py` (rigorous lower bounds on every eigenvalue of I), `code/p6_run3.sh`; raw tables `code/p6_out3/`.
Tags as before. This pass changes the status of (D): it is now **proved for K = 1, k0 > 1** (given the information event of Theorem R); for K >= 2 two of its three ingredients are proved and the third is reduced to a mode-adapted leverage bound; for k0 <= 1 it is open.

## 15. Proof of (D)

Notation: `lam_bar = mu/eps`; `n = 1-eps`; `b = min beta_k`, `B = max beta_k`, `r = B/b`; `Z = sum_k (alpha_k/beta_k) S_k` (so `dZ = -(lam - mu) dt + n dN`); `w := mu + b Z <= lam <= mu + B Z <= r w`. Box constants `c^- <= lam_theta/lam* <= c^+`, `kappa_1 = max(1, 1/c^- - 1)`. `DF_i = int_0^{tau_i} phi phi^T / lam_theta du`, `DJ_i` the same with lam*. `X_i^{tr} := tr(F_{i-1}^{-1} DF_i)` (written tr X_i in section 11).

### 15.1 Lemma D1 (negative moments of the intensity), every K [P]

For `p >= 0` with `(p+1) B n^2 < 2 mu`,
 `E_pi w^{-(p+1)} <= c_p E_pi w^{-p}`,   `c_p := eps r / (mu - (p+1) B n^2 / 2)`.
Hence `E_pi lam^{-q} <= E_pi w^{-q} <= prod_{j} c_{q-1-j}` for integer q, and for `q = j + s`, `s in (0,1]`, `E w^{-q} <= c_{q-1} ... c_s (c_0)^s`; all finite exactly when `q < k_eff := 2 mu/(B n^2)`.
*Proof.* Let `f = w^{-p}`. The generator of `(Z, N)` gives `L f = f'(Z)(-(lam - mu)) + lam [f(Z+n) - f(Z)]` with `f' = -p b w^{-p-1}`. By the alternating Taylor bound `(1+x)^{-p} - 1 <= -p x + p(p+1)x^2/2` (x >= 0),
 `L f <= p b w^{-p-1}[(lam - mu) - n lam] + (p(p+1)/2) b^2 n^2 lam w^{-p-2}`
 `= p b w^{-p-1}(eps lam - mu) + (p(p+1)/2) b^2 n^2 lam w^{-p-2}`   (since 1 - n = eps).
Use `lam <= r w` in both terms: `L f <= p b [ eps r w^{-p} - (mu - (p+1) b r n^2/2) w^{-p-1} ]`, and `b r = B`. Stationarity `E_pi L f = 0` (f is bounded by mu^{-p}; L f is dominated by polynomial moments, finite by Lemma `lem:moments`) gives `(mu - (p+1)Bn^2/2) E w^{-(p+1)} <= eps r E w^{-p}`. ∎
For K = 1 this is `E lam^{-(p+1)} <= [gamma/(beta mu - (p+1)alpha^2/2)] E lam^{-p}`, `gamma = beta eps`; at p = 0 it is Theorem A's `m <= 2 eps/(2mu - n^2 beta)`.
**Numerical check** (`code/p6_moments.exe`, 8 replicates of 3e5 time units; bound / Monte Carlo):

| case (k_eff) | q = 1 | 1.5 | 2 | 2.5 | 3 |
|---|---|---|---|---|---|
| K=1, mu=1, eps=0.3 (4.08) | 1.05 | 1.19 | 1.33 | 1.80 | 2.69 |
| K=1, mu=1, eps=0.1 (2.47) | 1.09 | 1.43 | 2.25 | n/a | n/a |
| K=1, mu=1, eps=0.03 (2.13) | 1.06 | 1.46 | 4.42 | n/a | n/a |
| K=1, mu=2, eps=0.1 (4.94) | 1.01 | 1.07 | 1.07 | 1.19 | 1.31 |
| K=2, beta=(1,2), mu=1.5, eps=0.1 (1.85) | 2.94 | 10.3 | n/a | n/a | n/a |
| K=2, mu=3, eps=0.1 (3.70) | 2.26 | 3.86 | 6.54 | 14.0 | 35.9 |

The bound holds in all 23 cases where it is finite, and "n/a" is exactly `q >= k_eff`. It is within 1-9 % of the truth for K = 1, q = 1.

### 15.2 Lemma D2 (the gap quantities, any K, any horizon) [P]

(a) Pathwise: `DJ_i <= c^+ DF_i`, `DF_i <= DJ_i/c^-`, and `F_{i-1} >= F_0 + J_{t_{i-1}}/c^+` (since DF >= DJ/c^+).
(b) For `M = F_{i-1}^{-1}`: `E[ g_i^T M g_i | G_{i-1} ] <= 2(1 + kappa_1^2) tr(M E[DJ_i | G_{i-1}]) <= 2(1+kappa_1^2) c^+ E[tr X_i | G_{i-1}]`.
*Proof.* `g = s + h`, `g^T M g <= 2 s^T M s + 2 h^T M h`. The score identity `E[s s^T | past] = E DJ_i` (Lemma `lem:score`) gives `E s^T M s = tr(M E DJ_i)`. The drift part is rank one, `h = a (1/lam* - 1/lam_theta)(tau)`, so `h^T M h = (a^T M a/lam*(tau)^2)(1 - lam*/lam_theta)^2 <= kappa_1^2 a^T M a/lam*^2`, and since tau has density `lam* e^{-Lambda*}`, `E[a a^T/lam*(tau)^2 | past] = int e^{-Lambda*} phi phi^T/lam* du = E DJ_i`. ∎
No leverage condition appears in (b): this replaces the lemma F3/F4 pair in the early phase. Measured (`code/p6_out3`, K = 1, truth start, scaled prior, T = 6000): `E G / sum_i tr X_i` = 0.99, 1.01, 1.11, 1.09 for eps = 0.01, 0.03, 0.1, 0.3 (G = sum g^T F_{i-1}^{-1} g), against the proved constant 2(1+kappa_1^2)c^+; and `E X <= P_0 + E G` holds in all 8 runs (ratio 0.60-0.98).

### 15.3 Proposition D3 (global decomposition without a bad-gap count) [P]

Let `t_*` be deterministic, `N_e = #{i : t_{i-1} < t_*}` (a stopping time; `1{t_{i-1} < t_*}` is predictable), `L_e := E ln det(F_{N_e}/F_0)`, `D_e := E sum_{i <= N_e} def_i`, and `Gam := 2(1+kappa_1^2) c^+ (L_e + D_e)`. Then
 (i) **early window** (Proposition `prop:reg` stopped at N_e; optional stopping for the square-integrable martingale M, E M = 0): `E sum_{i<=N_e} X_i <= P_0 + Gam`;
 (ii) **whole horizon**: `E X_n <= P_0 + Gam + E G_late`, `G_late` the sum over gaps with `t_{i-1} >= t_*`;
 (iii) **regret**: `R_i <= X_i` always and `R_i <= (1/2 + C_2 delta) X_i` on good gaps, so
 `E R_n <= (1/2 + C_2 delta) E X_n + (1/2 - C_2 delta)[ E sum_{i<=N_e} X_i + E sum_{late, bad} X_i ]`
 `<= (1/2 + C_2 delta)(P_0 + Gam + E G_late) + (1/2 - C_2 delta)(P_0 + Gam) + C_B`.
Here `E G_late <= (1-delta)^{-1}-type bound of Theorem R` on `ln det(F_n/F_{N_e})` with the sharp prefactor `(1+eta)(1+delta)^2`, plus `C_rho + C_B/lam_min(F_0)` exactly as in `theoremR_proof_complete.md` section 2 (that argument needs `lbar(t) <= delta mu_lo/(1+1/eta)` only for `t >= t_*`).
**What this removes:** the term `C_B^{early}` = (number of early gaps) x (bound per gap) = `(mu/eps) t_* x O(1)`, and the need for a high-probability bound on sup D in the early window. The early window now costs `P_0 + Gam` with `Gam = O(c^+ d ln t_*) + O(D_e)`.

### 15.4 Proposition D4 (D_e = O(eps^a)), K = 1, k0 > 1 [P]

Statement. Let the prior be scaled, `F_0 = rho_0 diag(1/lam_bar, lam_bar/beta^2)`, assume stationary start, the information event `E_t := {J_s >= (s/2) N for all s >= t_J}` of probability `>= 1 - C e^{-c t_J}` (Theorem FK/tangent minorant; N = N(lam_0), the explicit tangent matrix), and `0 < a < (k_eff - 1)/2`. Then there are explicit constants with
 `D_e <= C_D eps^a`,  `C_D = C_D(a, k0, mu, beta, rho_0, c^+, c^-, t_J, C_hat)`.
Proof (steps; constants kept symbolic).
1. *Pathwise reduction* (Lemma 6.4 + D2a): `def_i <= (tr X_i)^{1+a}` (`min(1,x) <= x^a`), and `tr X_i <= (1/c^-) tr(F_{i-1}^{-1} DJ_i) <= (1/c^-) min{ tr(F_0^{-1}DJ_i), tr((F_0 + J_{t_{i-1}}/c^+)^{-1} DJ_i) }` -- a functional of the **true** process only.
2. *Leverage along a gap.* For any weight matrix `W >= 0`: `tr(W^{-1}DJ_i) = int_0^{tau_i} phi^T W^{-1} phi / lam* du`. Suppose `phi^T W^{-1} phi <= (C/omega) lam_bar (1 + (beta y/lam_bar)^2)` with `y = S e^{-beta u}`. Then, using `y^2 <= (lam - mu)^2/alpha^2 <= lam^2/alpha^2` and `int_0^tau lam du = E_i` (E_i = Lambda*(tau_i) ~ Exp(1), independent of the past),
 `tr(W^{-1} DJ_i) <= (C/omega)[ lam_bar T_i + (beta/alpha)^2 E_i/lam_bar ]`,  `T_i := int_0^{tau_i} du/lam*(u) <= E_i/lam_end^2`, `lam_end = lam*(tau_i^-)`.
 - *F_0 window:* `W = F_0`, `omega = rho_0`, `C = 1`.
 - *N window:* `W = (t/(2c^+)) N` on E_t, `omega = t/(2c^+)`, `C = C_hat := 1/(lam_bar lambda_min(S N S))`, `S = diag(1, beta/lam_bar)`; this is exactly the quantity tabulated by `code/p6_chat.py`: **C_hat = 6.6, 8.2, 9.2, 9.6 (f = 1.5) or 6.6, 7.5, 7.8, 7.9 (best f) for eps = 0.3, 0.1, 0.03, 0.01 (mu = beta = 1)**, i.e. O(1) in eps. (The cancellation behind this: `phi^T I^{-1} phi = lam_bar[(k0-1) - 2(k0-1)s + k0 s^2]`, `s = beta y/lam_bar`, so the bad direction y = 0 has leverage (k0-1)/eps only at atypically small states.)
3. *One-gap moment.* Conditional on the start intensity `lam_0`, for E < (lam_0 - mu)/beta one has `lam_end >= lam_0 - beta E` (compare `Lambda(u) >= (lam_0 - mu)(1 - e^{-beta u})/beta`), so for `E <= lam_0/(2 beta)`, `lam_end >= lam_0/2`. Hence
 `E[ T_i^{1+a} | lam_0 ] <= (2/lam_0)^{2(1+a)} Gamma(2+a) + mu^{-2(1+a)} Gamma(2+a+m) (2 beta/lam_0)^m`  (Markov of order m on {E > lam_0/2beta}), any m > 0.
4. *Sampling by events.* The sum over events of a function of the pre-event state is a time integral against `lam`: `E sum_i h(lam_{t_i^-}) = int E[lam_t h(lam_t)] dt`, and the post-event intensity is `lam_0 = lam_{t^-} + alpha >= lam_{t^-}`. So
 `E_palm[ lam_0^{-q} ] <= E_pi[ lam^{1-q} ]/lam_bar <= B_{q-1}/lam_bar`,
 with `B_q` the bound of Lemma D1 (finite for `q - 1 < k_eff`). Taking `q = 2(1+a) - ...`: the first term needs `B_{1+2a}`, finite since `a < (k_eff - 1)/2`; the second with `m = 1 + 2a`.
5. *Assemble.* Per gap, with `omega` the weight of step 2,
 `E[(tr X_i)^{1+a}] <= 2^a (C/(c^- omega))^{1+a} [ lam_bar^{1+a} E_palm T_i^{1+a} + ((beta/alpha)^2 lam_bar^{-1})^{1+a} Gamma(2+a) ]`,
 and by steps 3-4 `E_palm T_i^{1+a} <= C_3 B_{1+2a}/lam_bar`, while Lemma D1 gives `B_{1+2a} = O(lam_bar^{-(1+2a)})`. The first term is therefore `O(lam_bar^{1+a} lam_bar^{-1} lam_bar^{-(1+2a)}) = O(lam_bar^{-(1+a)})`, the second is `O(lam_bar^{-(1+a)})` too. The event rate is `lam_bar`, so per unit time the cost is `lam_bar^{-a} omega^{-(1+a)}` times a constant:
 - F_0 window `[0, t_J]` (`omega = rho_0`): `sum <= C' t_J (eps/mu)^a rho_0^{-(1+a)}`;
 - N window `[t_J, infinity)`, `omega = t/(2c^+)`: `sum <= C' (2c^+ C_hat/c^-)^{1+a} (eps/(mu t_J))^a / a` (the `int t^{-(1+a)} dt` converges);
 - complement of E_t: use the F_0 bound with weight `P(E_t^c)^{1/2} <= C e^{-c t/2}`; sum against `lam_bar dt`, per-gap `(E (tr X)^{2(1+a)})^{1/2} <= C'' (eps/(mu rho_0))^{1+a}`: `<= C'' (2C/c) rho_0^{-(1+a)} (eps/mu)^a e^{-c t_J/2}`, still `O(eps^a)` for any c > 0.
Sum: `D_e <= C_D eps^a`. ∎ With a = 1 allowed when `k_eff > 3` and the Gamma factors; for the main example (k0 = 2, k_eff = 2.1-2.5) any `a < 0.5`. For the fixed prior `t_J = rho_0 mu (k0-1)/eps` is forced and the same computation gives `O(eps^{a-1})`, a vacuous bound for `a < 1`: **the proof of (D) needs a prior scaled to the information level** (the measured defect is still O(1) for the fixed prior, 1.1-1.4).
**Measured exponent** (`code/p6_out2`, scaled prior): `E sum def_i` = 1.07, 0.61, 0.17, 0.05 for eps = 0.3, 0.1, 0.03, 0.01; log-log slope `1.09` against the proved `a < 0.5` (conservative: the proof uses `min(1,x) <= x^a`).

### 15.5 K >= 2: what carries over, what does not [P] / [O]

- Lemma D1 (every K, with `k_eff = 2 mu/(B n^2)`), Lemma D2, Proposition D3 (steps 1-3) and steps 1, 3-5 of D4 carry over verbatim.
- **Step 2 does not**: `C_hat_K = 1/(lam_bar lambda_min(S N S))` (`code/p6_chat.py`, tangent N, best f):

| K, grid, mu | eps = 0.3 | 0.1 | 0.03 | 0.01 | trend |
|---|---|---|---|---|---|
| 1, 1, 1 | 6.6 | 7.5 | 7.8 | 7.9 | bounded |
| 1, 1, 2 | 10.7 | 11.5 | 11.9 | 12.0 | bounded |
| 2, (1,2), 1.5 | 72 | 157 | 461 | 1329 | `~ 1/eps` |
| 3, (1,2,4), 2 | 469 | 1043 | 3071 | 8864 | `~ 1/eps` |

With a single scale for all state coordinates, the K >= 2 leverage constant grows like `1/eps` (the fast modes `y_j` of Problem 12 are O(1) at typical states but the bound allows O(lam_bar/beta)). Then `D_N <= C eps^{a-1}` is vacuous. A proof for K >= 2 needs the mode-adapted estimate `E[(y_j^2/lam)^{1+a}] = O(1)` (Lemma 2 of `problem12_upper_bounds.md` gives the first moment `a = 0`) -- **open**; measured defect for K = 2, 3 is O(1) and decreasing (K = 2: 1.03, 0.60, 0.16 for eps = 0.3, 0.1, 0.03; K = 3: 1.36, 0.42, 0.14).

## 16. Explicit constants for K >= 2

What the K-general theorem needs: `P_0`, `ln det F_0`, `ln det I`, `t_*`, `c^+, c^-, kappa_1`, `C_rho, C_B` (K-general by `K2_regret_chain.md`), `D_e`. Formulas and numbers for beta = (1,2), mu = 1.5 (K = 2) and beta = (1,2,4), mu = 2 (K = 3), alpha_k* = (1-eps) beta_k/K:

- `P_0 = rho_0 [ eps (D mu)^2/mu + (mu/eps) sum_k (D alpha_k/beta_k)^2 ]` (scaled prior).
- `ln det F_0 = (K+1) ln rho_0 + (K-1) ln(mu/eps) - 2 sum_k ln beta_k`. For K = 1 this is eps-free; **for K >= 2 it contains `-(K-1) ln(1/eps)`**, so `-(1/2) ln det F_0 = ((K-1)/2) ln eps + ...` *decreases* the bound by `((K-1)/2) ln(1/eps)`. An eps-free alternative is `F_0 = rho_0 diag(1/lam_bar, lam_bar/beta_1^2, 1/beta_2^2, ...)` (ridge and top weights scaled, middle weights fixed); not run.
- `ln det I`: **bracket**. Rigorous lower bound = sum of the best tangent-minorant eigenvalue lower bounds (Weyl); upper = stationary Monte Carlo (8 replicates, `p6_moments`):

| K | eps | rigorous lower bound on ln det I | Monte Carlo ln det I | eigenvalues of I (MC) |
|---|---|---|---|---|
| 2 | 0.3 | -4.75 | -3.89 | 0.032, 0.095, 6.6 |
| 2 | 0.1 | -4.70 | -3.28 | 0.031, 0.065, 18.9 |
| 2 | 0.03 | -4.65 | -2.94 | 0.013, 0.064, 62.8 |
| 2 | 0.01 | -4.63 | -2.85 | 0.0046, 0.066, 191 |
| 2 | 0.003 | (n/a, double precision) | -2.74 | 0.0015, 0.066, 641 |
| 3 | 0.3 | -9.73 | -8.67 | 0.0036, 0.042, 0.126, 9.2 |
| 3 | 0.1 | -9.37 | -7.68 | 0.0048, 0.024, 0.151, 26.6 |
| 3 | 0.03 | -9.18 | -7.18 | 0.0054, 0.0096, 0.168, 88 |
| 3 | 0.01 | -9.11 | -7.03 | 0.0032, 0.0060, 0.173, 266 |
| 3 | 0.003 | (n/a) | -6.87 | 0.0012, 0.0059, 0.174, 813 |

`ln det I` is bounded and converges (K = 2: -2.74 at eps = 0.003, K = 3: -6.87), while `lam_max ≈ (1.9 to 2.7)/eps` and `lam_min ≈ (0.4 to 0.5) eps`. This confirms the three-tier prediction: one eigenvalue `~ 1/eps`, K-1 of order one, one `~ eps`.
- `k_eff = 2 mu/(B n^2)`: K = 2: 1.85 (eps = 0.1), 3.70 for mu = 3; the condition for (D) is `a < (k_eff - 1)/2`.
- Burn-in: `t_* = max(t_lev, t_mix)`, `t_mix ≈ 2/(beta_1 eps delta^2)` is the K = 1 measurement; for K = 2 the iterate exit time `T_exit` (delta = 0.25, scaled prior, 32 replicates) is 1013, 545, 777 for eps = 0.3, 0.1, 0.03 (the ridge sets `1/(beta_1 eps)` = 3.3, 10, 33 times a factor): the 1/eps scaling is not resolved there (non-monotone, p90 up to 2500), so **t_* for K >= 2 is only an order**.
So for K = 2, 3 the statement of Theorem 6 holds with `d = K+1`, `P_0`, `ln det F_0` as above, `ln det I` from the bracket, `k_eff` as above; its only unproved input is the mode-adapted version of step 2 of D4.

## 17. k0 <= 1

Stationary information vs eps, K = 1, beta = 1 (`code/p6_out3/I_k1_*`, 8 replicates; det I):

| mu (k0) | eps = 0.3 | 0.1 | 0.03 | 0.01 | 0.003 | fit |
|---|---|---|---|---|---|---|
| 0.3 (0.6) | 1.50 | 2.70 | 5.24 | 9.20 | 15.3 | `det I ≈ 1.5 eps^{-0.4} - 1` (the implied C at the five eps, 0.3 to 0.003: 1.54, 1.47, 1.54, 1.62, 1.57) |
| 0.5 (1.0) | 1.00 | 1.57 | 2.45 | 3.40 | 4.41 | `det I ≈ ln(1/eps) - 1` (predicted from m = 2 eps ln(1/eps)/beta; values 0.20, 1.30, 2.51, 3.61, 4.81) |
| 0.8 (1.6) | 0.65 | 0.89 | 1.16 | 1.31 | 1.53 | limit 1/(k0-1) = 1.67 |
| 1 (2) | 0.52 | 0.68 | 0.82 | 0.89 | 0.93 | limit 1 |
| 2 (4) | 0.26 | 0.30 | 0.32 | 0.32 | 0.32 | limit 1/3 |

Exact algebra: `det I = (m E lam - 1)/alpha^2` with `E lam = mu/eps`; Theorem A's dichotomy for `m = E[1/lam]` (`2 eps/(2mu - beta)`, `(2/beta) eps ln(1/eps)`, `C eps^{k0}`) gives
 `det I ≈ 1/(beta^2 (k0 - 1))` (k0 > 1),  `≈ (ln(1/eps) - 1)/beta^2` (k0 = 1),  `≈ (C eps^{k0-1} - 1)/beta^2` (k0 < 1).
**Consequence for the regret constant:** `(1/2) ln det I` grows like `((1-k0)/2) ln(1/eps)` for k0 < 1 and `(1/2) ln ln(1/eps)` for k0 = 1: **for k0 < 1 there is a genuine `C ln(1/eps)` term in time units with `C = (1 - k0)/2`** (for k0 > 1, none). The ridge eigenvalue is also larger (`lam_min ≈ eps^{0.5 to 0.7}` for k0 = 0.6, 1.0 against eps), so the ridge saturates earlier, as seen in the slopes. The three ingredients of the theorem that use k0 > 1 are: Lemma D1 (needs `k_eff > 1`; for k0 <= 1 only the crude `E lam^{-q} <= mu^{-q}` is available), the bounded-det-I statement, and Theorem A's `m`. The remaining chain (Theorem R/D3, tangent information level, mixing) does not use it. So for k0 <= 1: **(D) is not proved** (D1 fails), the log-det term is explicit as above, and the numerics (section 13) show slope d/2 and defect O(1) (1.9 and 2.8 for k0 = 0.6, eps = 0.03 and 0.1).
Prediction vs data: the constant `E R_T - ln T` at T = 6000, mu = 0.3, should rise by `(1/2) ln(det I(0.003)/det I(0.3)) = 1.16` from eps = 0.3 to 0.003; measured at eps = 0.1, 0.03: 0.96, 0.15 (s.e. 0.6-0.8): not resolved by the data.

## 18. Final status of Problem 6

| item | status |
|---|---|
| (D)(ii) bad early gaps | **proved** (D3: no sup D, no gap count), any K |
| (D)(i) `D_e` | **proved for K = 1, k0 > 1**, scaled prior: `O(eps^a)`, `a < (k_eff-1)/2` (measured exponent 1.09); open for K >= 2 (mode-adapted leverage), k0 <= 1 (D1 fails), fixed prior (bound vacuous) |
| negative moments of lam (Lemma D1) | **proved, every K**; checked against simulation in 23/23 cases |
| explicit constants K >= 2 | formulas and tables for K = 2, 3: `P_0`, `ln det F_0`, bracket on `ln det I`, `k_eff`; burn-in only an order |
| k0 <= 1 | explicit log-det law (`(1-k0)/2 · ln(1/eps)`), checked against 5 values of mu; the D-chain is not covered |
| Theorem 6 (K = 1, k0 > 1) | **no hypothesis (D) any more**; inputs: info event (Theorem FK/tangent), conservative universal constants of (H3) |
