# Problem 8: which regret, and how it relates to parameter error (mathematics, then code)

Tags: **[P]** proved here, **[N]** numerical (section 6, filled after the runs), **[O]** open. Code: `code/p8_link.cs`, `code/p8_run.sh`, `code/p8_analyze.py`.
Dates: 25 Oct 2026. Nothing refereed. The task (tracker row 8): "A precise definition used consistently in all theorems, and a lemma connecting the two notions if both are kept."

## 0. The notions currently in use, and the confusion between them

The notes use, without ever relating them in one place:

| symbol | where | what |
|---|---|---|
| `r_n` | `box_cost_and_lower_bound.md` B1 | realised log-loss regret `sum_i [l_i(theta_i) - l_i(theta*)]` (random) |
| `E r_n = KL(P_theta*^{(n)} ‖ Q^{(n)})` | B1, Prop LB | the headline of the lower bound |
| `R_n` | Prop 1 (`problem05_online_newton.md`), Theorem R | pathwise pseudo-regret `sum_i int_gap lam* psi(lam_theta/lam*) du`, `psi(r) = r - 1 - ln r` |
| `X_n = sum v_i^T DF_i v_i`, `D_n = ‖theta_{n+1} - theta*‖^2_{F_n}`, `e_n = (theta - theta*)^T J_n (theta - theta*)` | Lemma C1, Theorem R, `problem14` | parameter error in various metrics |
| "best fixed theta*" | tracker wording | ambiguous: the **truth** theta* or the **hindsight** maximiser `theta_hat_n` |

These are four different objects. Section 1 fixes one definition; sections 2-4 prove the lemmas that connect them; section 5 says what is lost in each direction.

## 1. The definition to use throughout [P]

**Setting.** Observation of the point process on `[0, T]`. A *sequential predictor* is a predictable family of intensities `lam_hat_t = phi(t)^T theta_{N_{t-}+1}` (parameters updated at events; `theta_i` is `G_{i-1}`-measurable). `Q` denotes the law of a point process with this predictive intensity, `P*` the law with the true intensity `lam*`.

**Definition (headline).** The *regret* of the predictor against the truth is
 `Reg_T := E_{P*} [ log (dP*/dQ)|_{[0,T]} ] = KL(P*_T ‖ Q_T)`.
By the likelihood-ratio formula for point processes,
 `log (dP*/dQ)|_{[0,T]} = int_0^T (lam_hat - lam*) dt - sum_{t_i <= T} log(lam_hat(t_i)/lam*(t_i)) =: r_T`
(the realised log-loss regret against theta*), and
 **(1.1)** `Reg_T = E r_T = E R_T`,  `R_T := int_0^T lam*(t) psi(lam_hat_t/lam*_t) dt`.
Events-indexed version: `Reg_n = E sum_{i<=n} r_i = E sum_{i<=n} R_i` with the same identity gap by gap.

*Proof of (1.1).* `r_T - R_T = int_0^T (lam_hat - lam*) dt - sum log(lam_hat/lam*)(t_i) - int_0^T [lam_hat - lam* - lam* log(lam_hat/lam*)] dt = int_0^T log(lam_hat/lam*)(lam* dt - dN_t)`, the integral of a predictable bounded-on-compacts integrand (`lam_hat`, `lam*` are bounded below by `mu_lo > 0` and have all moments) against the compensated martingale `dN - lam* dt`; it has mean zero. ∎ (Gap by gap this is the computation of Proposition 1(a): given `G_{i-1}`, `tau_i` has density `lam* e^{-Lambda*}`, so `E int_0^tau f du = int e^{-Lambda*} f du`.)

**Consequences.**
- `Reg_T` is the one-line object of the minimax lower bound (Prop LB), `inf_Q sup_theta KL(P_theta ‖ Q)`; the regret of *every* predictor, not only online Newton, has this form (B1).
- `R_T` is nonnegative pathwise (`psi >= 0`); `r_T` is not. The proofs use `R_T` (pathwise, and a Bregman divergence), the statements use `Reg_T = E R_T`. This removes the ambiguity in "regret" between Prop 1, Theorem R and Prop LB: they are all `Reg_T`.
- **Not** the headline: regret against the hindsight maximiser `theta_hat_n` (section 4: it differs by the Wilks constant `d/2`).
- Horizon: time `T` is the primary index; `n` events is the same statement with `T = t_n` (`Reg_n` above), and `ln T_n = ln n - ln(mu/eps) + O(1)` (Problem 6, section 1). Every theorem should say which clock it uses.

## 2. Two-sided comparison of the Bregman divergence with the quadratic form [P]

**Lemma 8.1 (psi).** For `r > 0`, with `r_- := min(1,r)`, `r_+ := max(1,r)`:
 `(r-1)^2/(2 r_+^2) <= psi(r) <= (r-1)^2/(2 r_-^2)`, and also `(r-1)^2/(2 r_+) <= psi(r)` and `psi(r) <= (r-1)^2/r`.
*Proof.* `psi(1) = psi'(1) = 0`, `psi'' = 1/r^2`; Taylor with Lagrange remainder gives `psi(r) = (r-1)^2/(2 xi^2)`, `xi` between 1 and r, hence the first pair. The bound `psi(r) >= (r-1)^2/(2 r_+)` and `psi <= (r-1)^2/r` are `log r >= 1 - 1/r` and its integrated form (checked on a grid in section 6). ∎

**Lemma 8.2 (gap sandwich).** On a ratio set (`c^- <= lam_theta/lam* <= c^+` on the gap) put `c~^+ = max(1,c^+)`, `c~^- = min(1,c^-)`, `v = theta - theta*`, `Q_i := v^T DJ_i v = int_0^{tau_i} (phi^T v)^2/lam* du`. Then, pathwise,
 **(2.1)** `Q_i/(2 c~^+) <= R_i <= Q_i/(2 c~^{-2})` (lower bound from `psi(r) >= (r-1)^2/(2 max(1,r))`, upper from the Lagrange form; the weaker lower bound from the Lagrange form is `Q_i/(2 c~^{+2})`).
On a localisation region `|lam_theta/lam* - 1| <= delta` (delta <= 1/4): `R_i = (1/2 + O(delta)) Q_i`; precisely `(1+delta)^{-2} Q_i/2 <= R_i <= (1-delta)^{-2} Q_i/2`.
*Proof.* `lam* (r-1)^2 = (phi^T v)^2/lam*` with `r = lam_theta/lam*`; insert Lemma 8.1 and integrate over the gap. ∎
Relation to `X_i = v^T DF_i v`: `DF_i` and `DJ_i` differ by `lam_theta/lam* in [c^-, c^+]`, so `Q_i/c^+ <= X_i <= Q_i/c^-`. (This is the `R_i <= X_i` of Theorem R, F1, with the constant `1/(2 c~^-)`.)

## 3. Regret controls parameter error: the information-weighted average [P]

The problem statement asks for a lemma that goes in both directions. The direction *regret => error* is not available for the last iterate (the last iterate can be bad while the regret is small), but it holds exactly for an information-weighted average.

**Lemma 8.3 (regret => error of the averaged iterate).** Let `J_n = sum_i DJ_i` (assume `J_n > 0`) and the information-weighted average
 `theta_bar_n := J_n^{-1} sum_{i<=n} DJ_i theta_i`.
Then pathwise (any predictor, any realisation)
 **(3.1)** `sum_{i<=n} Q_i = (theta_bar_n - theta*)^T J_n (theta_bar_n - theta*) + sum_i (theta_i - theta_bar_n)^T DJ_i (theta_i - theta_bar_n) >= ‖theta_bar_n - theta*‖^2_{J_n}`,
and on a ratio set, with Lemma 8.2,
 **(3.2)** `‖theta_bar_n - theta*‖^2_{J_n} <= 2 c~^{+} R_n`,  so  `E ‖theta_bar_n - theta*‖^2_{J_n} <= 2 c~^{+} Reg_n`.
*Proof.* (3.1) is the exact bias-variance decomposition for the matrix-weighted mean: `sum DJ_i (theta_i - theta*) = J_n (theta_bar - theta*)`, so cross terms vanish. (3.2) combines (3.1) with the left inequality of (2.1). ∎
Interpretation: `Reg_n ~ (d/2) ln T` therefore bounds the Fisher-normalised error of the averaged iterate, `‖theta_bar - theta*‖^2_{J_n} <= (d + o(1)) c~^{+} ln T`, which is weaker than the optimal `chi^2_d` but is a *pathwise* consequence of the regret with no further estimate.

**Remark (why not the last iterate).** `D_n = ‖theta_{n+1} - theta*‖^2_{F_n}` stays `O(d)` for the Newton iterate (it is a `chi^2_d`-like quantity, E `~ d`, Problem 5 numerics 1.3-2.7), while `Reg_n` grows like `(d/2) ln T`: the two are on different scales and there is no inequality of the form `D_n <= C Reg_n`. What holds is that both are controlled by the same quantity `G`:

**Lemma 8.4 (common bound; restatement of Lemma C1 / Prop 1).** For the Newton recursion, pathwise `D_n + X_n = D_0 + sum_i g_i^T F_i^{-1} g_i - M_n`, and `E X_n <= P_0 + E G_n`, `E D_n <= P_0 + E G_n - E X_n`; hence `E R_n <= E X_n/(...)` by Lemma 8.2 and `E D_n <= P_0 + E G_n`, with `E G_n <= (1+o(1)) ln det(F_n/F_0)` (Theorem R). So the *rate* of the regret, `d/2 · ln T`, and the *size* of the last-iterate error, `O(d)`, are two projections of one information-accumulation statement.

## 4. Regret versus the hindsight maximiser: Wilks constant [P] heuristic, [N] check

Let `L_n(theta) = sum_i l_i(theta)` and `theta_hat_n = argmin L_n` (batch MLE). The *hindsight regret* is `Reg^hind_n := E[ sum_i l_i(theta_i) - L_n(theta_hat_n) ]`. Then
 **(4.1)** `Reg^hind_n = Reg_n + E[ L_n(theta*) - L_n(theta_hat_n) ]` (the comparator `theta_hat_n` is at least as good as the truth on the sample, so the hindsight regret is the *larger* one) and `L_n(theta*) - L_n(theta_hat_n) = (1/2) (theta_hat_n - theta*)^T H (theta_hat_n - theta*) (1 + o(1)) -> (1/2) chi^2_d` under local asymptotic normality, so `Reg^hind_n = Reg_n + d/2 + o(1)` when LAN holds (ridge saturated, `T lam_min(I) >> 1`; Problem 11 shows LAN fails at `eps = c/T`).
So the two definitions differ by the additive constant `d/2`; the leading term `(d/2) ln T` and all `eps`-dependence of the additive constant are unaffected except through this `d/2` and the ridge-saturation correction. **Headline stays against the truth theta*:** it is the one equal to a KL divergence and the one for which lower bounds exist; the hindsight version equals `Reg_n + d/2` and should be quoted as such.

## 5. The converse (error => regret) and what is lost [P]

- *Error => regret.* Lemma 8.2: `R_n <= (1/(2 c~^{-2})) sum Q_i <= (c^+/(2 c~^{-2})) X_n`; and `E X_n <= P_0 + E G_n` (Prop 1). So any bound on the cumulative Fisher-norm error `X_n` gives a regret bound, with the sharp factor `1/2 + O(delta)` on a localisation region.
- *Regret => error.* Only for the information-weighted average (Lemma 8.3). A last-iterate error bound needs the Newton structure (Lemma 8.4) or an additional stability property of the algorithm.
- *Minimax link (already proved).* Prop LB: an estimator with `E‖theta_hat - theta‖^2 <= C ln n/n` forces `sup Reg_n >= (d/2) ln(n/ln n) - O(1)` (Bayes mixture, mutual information). Together with Lemma 8.3 this closes the loop: regret lower bound <= error upper bound (Prop LB), error of the averaged iterate <= regret upper bound (Lemma 8.3).

## 6. Predictions tested in the code [H]

- **T1.** `E r_n = E R_n` (identity 1.1): difference zero within standard error, for every start and every eps.
- **T2.** Lemma 8.1 on a fine grid of r, and Lemma 8.2 gap by gap: no violation; empirical constants `min R_i/(Q_i/2)`, `max R_i/(Q_i/2)` against `1/c~^+`, `c~^{-2}`-bounds from the visited ratio range.
- **T3.** `sum R_i/(Q-sum/2) = 1 + O(delta)` once the iterates are in R(delta).
- **T4.** (3.1) holds as an identity (residual machine precision), and `‖theta_bar - theta*‖^2_{J_n} <= 2 c~^{+} R_n`.
- **T5.** `E D_n ≈ d` and `E ‖theta_bar - theta*‖^2_{J_n}` is `O(d)`, while `E R_n ≈ (d/2) ln T`: different scales (Remark in section 3).
- **T6.** Hindsight regret `= Reg_n + d/2` when LAN holds (eps not tiny, T large).

## 7. Numerical results (25 Oct 2026; `code/p8_link.cs`, `code/p8_onestep.py`, `code/p8_run.sh`, `code/p8_analyze.py`; raw tables `code/p8_out/`)

Setup: mu = beta = 1 (K = 1) and mu = 1.5, beta = (1,2) (K = 2), alpha*_k = (1-eps) beta_k/K, stationary start, Newton recursion with clip and scaled prior rho_0 = 1, "near" start (1.3 mu*, 0.8 alpha*) or at the truth, n = 20000 events, 64 replicates (48 for K = 2). Regret columns are in events.

**T1, identity `E r = E R` — supported.**
- *Single gap, exact sampling* (`p8_onestep.py`, 4e5 gaps per cell, 28 cells over three parameter sets, start states s = 0, 1, 8, four theta's including far ones): `E(r - R)` within 2 standard errors of zero in every cell (e.g. +0.0012 +- 0.0009, -0.0010 +- 0.0008; the cells where r and R are 0.58 and 0.58 agree to 3 digits). One parameter set (alpha* = 1.2 > beta) is not stable; the identity does not need stability.
- *Whole run, high power* (4000 replicates, n = 100, 300): `E r - E R` = -0.013 +- 0.018, -0.018 +- 0.026 (K=1, eps=0.1), +0.016 +- 0.028, -0.021 +- 0.038 (eps=0.3), -0.005 +- 0.018, -0.024 +- 0.025 (K=2). The few 2-3 sigma excursions seen at 16-64 replicates (e.g. +0.32 +- 0.10) vanish at 4000 replicates, so they were sampling noise of the heavier-tailed realised regret r.
- At n = 20000 (64 replicates): `E r - E R` = -0.12 (0.37), -0.06 (0.29), +0.34 (0.49), +0.20 (0.48), +0.11 (0.48), +0.14 (0.47) for K = 1; -0.16 (0.49), -0.21 (0.63) for K = 2. All within 1 standard error.

**T2, Lemma 8.1 and 8.2 — supported, no violation.** Lemma 8.1 holds at all 4e5 grid points r in [e^-8, e^8] for all four forms (squared Lagrange lower/upper, linear lower `(r-1)^2/(2r_+)`, and `(r-1)^2/r`). Lemma 8.2 gap by gap: 0 violations in 28 single-gap cells x 4e5 gaps, and 0 violations in the full runs (about 1.9 million gaps each) once gaps with `Q_i < 1e-9` are excluded; the earlier "violations" (about 2e-5 of gaps) all had `R, Q ~ 1e-13` at `theta ≈ theta*` where `R = (b^T theta - b^T theta*) - quadrature` loses its last digits. Empirical `R_i/(Q_i/2)` over all gaps lies in [0.55, 1.87] (K = 1), [0.67, 2.13] (K = 2), inside the proved `[1/c~^+, 1/c~^{-2}]` for the visited ratio ranges ([0.30, 2.9] and [0.38, 2.1]).

**T3, sharp factor — supported.** `sum R_i / (sum Q_i/2)` at n = 20000: 1.008, 0.996 (eps = 0.03), 0.951, 0.957 (0.1), 0.950, 0.962 (0.3), K = 2: 0.981, 0.956. So `R_n = (1/2)(1 - 0.05) sum Q_i` or better: the factor 1/2 is attained to 1-5 %, in line with `R_i = (1/2 + O(delta)) Q_i` on R(delta).

**T4, averaged iterate — supported.** The decomposition (3.1) is an algebraic identity (the code reports `sum Q - ‖theta_bar - theta*‖^2_J >= 0` with minimum over runs from about 0.05 (n = 100) to about 3 (n = 20000), always positive). Inequality (3.2): `‖theta_bar - theta*‖^2_{J_n} / (2 c~^+ R_n)` over all runs: at most 0.62-0.72 (eps = 0.03), 0.54-0.60 (eps = 0.1), 0.39-0.40 (eps = 0.3), K = 2: 0.61, 0.34. So (3.2) holds with margin at least 1.4 (the constant 2 c~^+ with `c~^+ = max ratio` along the run is the pessimistic part; typical `‖theta_bar - theta*‖^2_J` is 3.2-4.3 for K = 1 against `2 R_n` = 8.2-15).

**T5, scales — supported.** `E D_n` (last iterate, Fisher norm of F_n) at n = 20000: 2.03, 1.95, 1.89, 1.88, 1.98, 1.97 (K = 1; d = 2) and 2.64, 2.75 (K = 2; d = 3), flat in n from n = 316 on (K = 1, eps = 0.1: 2.06, 2.09, 2.18, 2.01, 1.89 for n = 316 ... 20000), while `E R_n` rises 1.2 -> 2.3 -> 3.8 -> 5.3 -> 6.0 over the same range; the information-weighted average has `E ‖theta_bar - theta*‖^2_J` 0.5 -> 3.9 (saturating at roughly 2 d, because it carries the early iterates). So no inequality `D_n <= C Reg_n` can exist (Remark, section 3), and the regret rate and the last-iterate error are two views of the same accumulation `P_0 + G`.

**T6, Wilks constant — supported.** `E[L_n(theta*) - L_n(theta_hat_n)]` = hindsight regret minus truth regret = 1.08, 1.08 (eps = 0.03), 0.95, 0.95 (0.1), 0.99, 0.99 (0.3) against `d/2 = 1` for K = 1; 1.37, 1.40 against `d/2 = 1.5` for K = 2 (the 8-9 % deficit is the unsaturated ridge: `T lam_min(I)` is not yet >> 1 at n = 20000 for these eps, cf. Problem 6). Hindsight regret is therefore larger than `Reg_n` by about `d/2`, as in (4.1) after the sign correction made in section 4.

## 8. Consistency audit and the wording to use in the paper

| where | quantity called "regret" | is it `Reg` of Definition 1? | action |
|---|---|---|---|
| `problem05_online_newton.md` Prop 1 / `latex/part2_regret.tex` Prop `prop:reg` | pseudo-regret `R_n`; "E R_n <= X <= P_0 + G" | yes: `E R_n = Reg_n` by Lemma 8.0 (identity) | refer to Definition `def:regret` and Lemma `lem:regid` |
| Theorem R, `thm:R`, `thm:RK`, Problem 6 | `E R_n`, `E R_T` | yes | state "regret in the sense of Definition `def:regret`"; say which clock (events vs time), cf. Problem 6 section 1 |
| `box_cost_and_lower_bound.md` B1, Prop LB `prop:LB` | `KL(P_theta ‖ Q)` | yes, it is the definition | unchanged |
| `problem14_tracking_theorem.md` | tracking error in the Fisher norm; no regret | not applicable | state explicitly that it is a *parameter-error* result; the dynamic-regret reading would need a comparator path (Hall-Willett) |
| Hall-Willett (2016) regret (Problem 22) | adversarial dynamic regret of mirror descent against a comparator sequence, path-length variation | **no**: different object, no probabilistic model | the "difference from" paragraph must say: ours is the expected log-loss regret (KL) against the truth under the model, theirs is adversarial |
| Tracker row 8, "best fixed theta*" | ambiguous | against the truth: yes; against the hindsight maximiser: differs by `+d/2 + o(1)` (Wilks) | quote both, headline is the truth |
| `novelty_reassessment.md` E2 (conjecture `Reg = (d/2) ln T + (1/2) ln det I + O(1)`) | `E R_T` | yes | consistent with Problem 6 results |

**Paper-ready text** (to sit at the top of the model/regret section; implemented as `latex/part2a_regret_notion.tex`, included before `part2_regret` in `hawkes_master.tex`):
> *Regret.* A predictor issues a predictable intensity `lam_hat_t`. Its regret against the truth up to time `T` is `Reg_T = KL(P*_T ‖ Q_T) = E r_T`, with `r_T` the realised log-loss excess of `lam_hat` over `lam*`; equivalently `E R_T` with the nonnegative pathwise pseudo-regret `R_T = int_0^T lam* psi(lam_hat/lam*) dt`. All upper bounds in this paper bound `E R_T`; the lower bound is for `KL(P*_T ‖ Q_T)`. Regret against the hindsight maximiser exceeds `Reg_T` by `d/2 + o(1)` under LAN. Parameter-error statements are made in the Fisher norm and are related to the regret by Lemmas `lem:regsand`, `lem:avgerr` (regret controls the error of the information-weighted average) and by the common bound `P_0 + E G` (last iterate).

## 9. Status

| "Done when" (tracker) | status |
|---|---|
| a precise definition used consistently in all theorems | **done**: Definition 1 (KL = expected log-loss regret against the truth, with the pseudo-regret as its pathwise version); audit table in section 8; `latex/part2a_regret_notion.tex` added and included. The theorems themselves already bound `E R`; the only edit left in the existing text is replacing the parenthetical remark after `R_n` in `part2_regret.tex` by a reference to the definition (not done, to avoid touching the proofs). |
| a lemma connecting the two notions if both are kept | **done**: Lemma `lem:regid` (identity), `lem:regsand` (two-sided comparison with the Fisher form), `lem:avgerr` (regret => error of the information-weighted average, pathwise), Remark (last iterate is not comparable), (4.1) hindsight comparator. All checked numerically (section 7). |
| not covered | regret against a *time-varying* comparator (Problem 14, Hall-Willett type): a different notion, not defined here; high-probability regret (all statements are in expectation); a converse "last-iterate error => regret" beyond Lemma 8.2 for the averaged iterate |
