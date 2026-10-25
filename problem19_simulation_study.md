# Problem 19: simulation study (protocol, predictions, results)

Tags: **[P]** from proved statements, **[H]** prediction to test, **[N]** measured here. Code: `code/p19_sim.cs`, `code/p19_run.sh`, `code/p19_analyze.py`, `code/p19_thinning_check.py`; figures `figures/p19_*.png`; numbers `code/p19_summary.txt`. Date: 25 Oct 2026.
Tracker row 19: "Simulate with Ogata thinning using known (mu, alpha, beta). Sweep eps and T and check that regret and error follow the predicted scalings, including the burn-in. Done when: plots of error and regret versus T for several eps, with fitted slopes compared with theory."

## 0. Starting point

Before this session the simulation scripts existed (`code/p5_*.py`, `p14_*.py`, `p11_obs_info.py`, `p13*.py`) and the numbers were in tables inside the notes; **no figure had been produced**, no slope had a confidence interval, and the runs used different starts, priors and horizons. This study runs one protocol for all of them and produces the figures.

## 1. Protocol

**Process.** Exponential-dictionary Hawkes process `lam_t = mu + sum_k alpha_k S_k(t)`, `dS_k = -beta_k S_k dt + dN`, known `(mu, beta)`, truth `alpha*_k = (1-eps) beta_k/K` (branching ratio `1-eps`), simulated by **Ogata thinning** (the intensity is non-increasing between events, so `lam` at the gap start dominates). Stationary start (warm-up of `20/eps^2` events), state observed.
**Validity of the sampler** (`p19_thinning_check.py`, 3e5 events each, five parameter sets): compensator increments `Lambda(t_{i-1}, t_i)` satisfy the time-rescaling theorem, i.e. are i.i.d. Exp(1): KS distance 0.0008-0.0012 (p = 0.79-0.99), mean 1.001-1.002, variance 1.005-1.011, lag-1 autocorrelation within +-0.004; time-average intensity 3.32, 9.92, 31.3, 19.7, 4.94 against `mu/eps` = 3.33, 10, 33.3, 20, 5 (the 6 % deficit at `eps = 0.03` is the slow mixing, `1/gamma = 33`); `E[1/lam]` = 0.378, 0.157, 0.059 for `eps = 0.3, 0.1, 0.03` (`mu = beta = 1`), in line with `Theorem A` and the independent estimates of `p6_moments`.
**Estimator.** The Fisher-preconditioned online Newton recursion of Problem 5 (`F += int phi phi^T/lam_theta`, `theta <- theta - F^{-1} g`) with a clip to a box and the stability cap `sum alpha_k/beta_k <= 1 - eps/2` (the clip is not the exact `F`-projection of `problem04_projection.md`; it is active in a small fraction of steps). Prior `F_0 = rho_0 diag(1/lam_bar, lam_bar/beta_k^2)`, `lam_bar = mu/eps`, `rho_0 = 1`. Start "near" `(1.3 mu*, 0.8 alpha*)` or "far" `(2 mu*, 0.5 beta/K)`.
**Recorded** at 8 geometric checkpoints per decade of `T` (time units, `T` from 1): pseudo-regret `R_T = sum_i int lam* psi(lam_theta/lam*) du` (`E R_T` is the expected regret, Problem 8), `D = (theta - theta*)^T F (theta - theta*)`, `e = (theta - theta*)^T J (theta - theta*)`, the error vector `v`, the oracle information `J`, and per replicate the last exit time `T_exit(delta)` from the region `R(delta)` (delta = 0.25).
**Sweep.** K = 1 (mu = beta = 1): eps in {0.3, 0.1, 0.03, 0.01}, near start; far start for eps in {0.1, 0.03, 0.01}; T up to `3e4` (`2e4` for eps = 0.01); 96-200 replicates. K = 2 (beta = (1,2), mu = 1.5): eps in {0.3, 0.1, 0.03}, T up to `3e4`, 64-96 replicates. Roughly `1e9` simulated events.
**Fits.** Slopes of replicate-mean `E R_T` against `ln T` on the last decade `[T_max/10, T_max]` and on the previous decade; log-log slope of the mean squared error; 95 % intervals by bootstrap over replicates (B = 400).

## 2. Predictions [P]/[H]

- **P1 (regret).** `E R_T = (d/2) ln T + c(eps) + o(1)`, `d = K+1`; slope `d/2` against `ln T` once `T >> t_*`, and the constant `c` is flat in eps in time units (Problems 5, 6; events: `c - (d/2) ln(mu/eps)`). Before saturation the slope is below `d/2` (ridge not yet learned); the tier-resolved curve `(1/2) sum_j ln(1 + T lam_j(I)/f_j)` (Problem 6, section 4) is an upper envelope.
- **P2 (raw error).** `E (theta_hat_j - theta*_j)^2 = (I^{-1})_jj/T (1 + o(1))`: log-log slope `-1`, and `MSE_j * T -> (I^{-1})_jj` with `I` measured from `J_T/T`. For K = 1: `(I^{-1})_{mu mu} ≈ mu(k0-1)/eps` (ridge), `(I^{-1})_{alpha alpha} ≈ 2 eps beta` (Problem 6, section 4).
- **P3 (normalised error).** `E D_n` and `E e_n` converge to `d` (a `chi^2_d` mean), flat in `T` after the burn-in (Problem 8, T5; Problem 5 numerics).
- **P4 (burn-in).** The iterate stays in `R(0.25)` after `T_exit ≈ c_exit/(beta eps delta^2)` with `c_exit ≈ 2` (Problem 6, section 12), i.e. `32/eps`; the extra regret of a far start is `P_0/2` with `P_0 = rho_0[eps (D mu)^2/mu + (mu/eps) sum_k (D alpha_k/beta_k)^2]` (Problem 6, section 9).

## 3. Results (25 Oct 2026; every interval is a 95 % bootstrap over replicates)

Figures: `figures/p19_k1_regret.png` (regret vs T, fitted slopes; constant `E R - ln T`), `figures/p19_k1_error.png` (MSE of mu and alpha vs T, log-log, with `(I^-1)_jj/T` overlaid), `figures/p19_k1_normalised.png` (`E D_n`, `E e_n`), `figures/p19_burnin.png` (exit time vs `1/eps`, start cost), `figures/p19_k2.png` (K = 2: regret, MSE, `E D_n`). Full numbers: `code/p19_summary.txt`.

### P1 regret slope [N]: agrees with `d/2` in all 10 runs

| run | late-window slope of `E R_T` vs `ln T` [95 % CI] | target | mid-window slope |
|---|---|---|---|
| K=1, eps = 0.3, near | 0.95 [0.86, 1.05] | 1 | 1.01 |
| K=1, eps = 0.1, near | 0.96 [0.85, 1.08] | 1 | 1.01 |
| K=1, eps = 0.03, near | 1.10 [0.98, 1.23] | 1 | 1.26 |
| K=1, eps = 0.01, near | 1.08 [0.93, 1.22] | 1 | 1.41 [1.22, 1.62] |
| K=1, eps = 0.1, far | 1.08 [0.95, 1.21] | 1 | 1.08 |
| K=1, eps = 0.03, far | 0.92 [0.80, 1.08] | 1 | 1.49 [1.25, 1.75] |
| K=1, eps = 0.01, far | 1.05 [0.86, 1.22] | 1 | 2.38 [1.92, 2.89] |
| K=2, eps = 0.3 | 1.44 [1.26, 1.64] | 1.5 | 1.64 |
| K=2, eps = 0.1 | 1.57 [1.41, 1.75] | 1.5 | 1.50 |
| K=2, eps = 0.03 | 1.41 [1.19, 1.62] | 1.5 | 1.26 |

(late window = last decade of T, 1.8e3-2.4e4 for the small eps, up to 2.4e4.) All ten intervals contain the theoretical slope `d/2`. **The burn-in is visible as a mid-window slope above `d/2`** that grows as eps shrinks and the start gets worse (far start, eps = 0.01: 2.38): the regret first lags and then catches up with the start cost `P_0/2` (Problem 6, sections 9 and 12).
*Constant* `c = E R - ln T` at the end of the run (K = 1, near start): -1.32, -1.65, -1.11, -0.55 for eps = 0.3, 0.1, 0.03, 0.01 (s.e. about 0.3-0.5; right panel of `p19_k1_regret.png`). So c is **flat only up to a drift of about +0.8 at eps = 0.01**, where `T eps = 180` is still in the saturation climb (mid slope 1.41). With a start at the truth or a weak prior c is flatter (Problem 6, section 9: -2.2 to -1.2 at the same T); the claim "flat in eps in time units" is therefore supported as an *asymptotic* statement, with a finite-T correction that grows as eps shrinks and shrinks as `T eps` grows.

### P2 raw error [N]: slope -1 and MSE*T = (I^-1)_jj

Log-log slope of `E (theta_hat_j - theta*_j)^2` in the late window, K = 1: mu: -1.08 [-1.19, -0.96], -0.98 [-1.11, -0.85], -0.95 [-1.08, -0.83], -1.62 [-1.82, -1.37] for eps = 0.3, 0.1, 0.03, 0.01; alpha: -1.04, -1.02, -0.92, -1.29. The first three eps agree with -1; eps = 0.01 is steeper (-1.6 / -1.3) because the ridge has not saturated (`T eps = 180`, the MSE is still falling faster than 1/T while the information accumulates). K = 2: mu -0.98, -1.02, -1.10; alpha_1 -1.01, -0.96, -0.88; alpha_2 -1.07, -0.98, -0.85 (intervals about +-0.17 contain -1).
`MSE_j * T / (I^-1)_jj` (pointwise over the late window, `I` from the measured `J_T/T`): K = 1, mu: 0.88, 1.01, 1.15, 1.90; alpha: 0.92, 1.14, 1.24, 1.36 (eps = 0.3...0.01, near start); far start: mu 1.07-1.32, alpha 1.14-1.69; K = 2: 0.82-0.95 (eps = 0.3), 1.02-1.09 (0.1), 0.83-1.08 (0.03). So **the estimator attains the Cramer-Rao scale, `MSE*T = (1 +- 0.2) (I^-1)_jj`, for `T eps >~ 700`**; at `T eps = 180` the ratio is 1.4-1.9.
*Efficiency check at the last checkpoint* (K = 1, T = 23714, 160-200 replicates): the whitened error covariance `(T I)^{1/2} E[v v^T] (T I)^{1/2}` has eigenvalues (0.94, 1.00), (0.77, 1.03), (0.97, 1.26) for eps = 0.3, 0.1, 0.03; the per-replicate statistic `v^T (T I) v` has mean 1.94, 1.80, 2.24 and `P(. > 6) = 0.050, 0.050, 0.062` against the `chi^2_2` values 2 and 0.050; the mu-error correlation with the alpha-error is -0.77, -0.80, -0.77 against -0.78, -0.75, -0.73 from `I^{-1}`. So the late-time error law is the `chi^2_2`/Gaussian law with covariance `(T I)^{-1}`, the efficiency (LAN) statement behind Problem 9.
*Note on a statistic.* A first version of the analysis used `<MSE>_window * <T>_window`, which is biased by 1.52 for a log-uniform decade (Jensen) and gave a spurious ratio of about 1.6; the table above uses `<MSE_t T_t>` pointwise.

### P3 normalised error [N]

`E D_n` (late window) = 1.92, 1.89, 2.17, 2.11 for K = 1, eps = 0.3...0.01, near start; 2.13, 1.85, 2.11 far start (target 2); `E e_n` identical to within 0.01; K = 2: 2.89, 3.18, 2.95 (target 3). Flat in `T` after a burn-in of the order of `T_exit` (`p19_k1_normalised.png`, `p19_k2.png`; K = 2 right panel). The last-iterate Fisher-norm error stays at `d` while the regret grows like `(d/2) ln T`: the two-scales statement of Problem 8.

### P4 burn-in [N]

*Exit time* `T_exit(0.25)` (last time outside `R(0.25)`), median [p90], K = 1, near start: 129 [418], 331 [906], 932 [2525], 2373 [7651] for eps = 0.3, 0.1, 0.03, 0.01 (far start: 372, 966, 2052 for eps = 0.1, 0.03, 0.01). Proportional to `1/eps` (log-log slope of the median against `1/eps` = 0.86): `eps T_exit` = 39, 33, 28, 24. The earlier "mixing form" `2/(beta eps delta^2) = 32/eps` is right to within 25-40 % over this range, but the data fit an **information form** better:
 `T_exit ≈ z^2 max_j (I^{-1})_jj / (delta theta*_j)^2`, `z^2 = 1.4`,
the time for the standard error of the worst coordinate to fall to `delta/z` of its value. Observed/predicted over **all 10 runs**: 0.85, 0.98, 1.06, 1.02, 1.12, 0.99, 0.86 (K = 1) and 1.23, 1.11, 0.89 (K = 2); the binding coordinate is mu (the ridge) for K = 1 and for K = 2 at eps = 0.03, and alpha_1 for K = 2 at eps = 0.1 and 0.3. This explains the K = 2 behaviour noted in Problem 6 (burn-in 500-1100, not proportional to `1/eps`): there the *middle-tier* directions bind (`(I^-1)_{alpha_1 alpha_1}` is eps-independent) and the burn-in is set by the information in the middle tier, `~ 1/(delta^2 lambda_mid)`. [This is a refinement of Problem 6, section 12, which read the K = 1 data as mixing-limited; with `(I^{-1})_{mu mu} = mu(k0-1)/eps` the two forms coincide up to the constant for `mu = beta = 1`.]
*Start cost.* Far-start excess over the near start at the end of the run: 1.77 +- 0.88 (eps = 0.1), 3.41 +- 1.11 (0.03), 10.66 +- 2.06 (0.01) against `(P_0^{far} - P_0^{near})/2` = 0.68, 3.07, 10.05 (`p19_burnin.png`, right): agreement within the intervals for eps = 0.03 and 0.01 and 1.2 standard errors for eps = 0.1. The start cost is `∝ 1/eps` and is repaid only logarithmically slowly (the curves flatten, they do not close).

## 4. Net result and limits

| "Done when" | status |
|---|---|
| plots of error and regret versus T for several eps | **done**: 5 figures (K = 1: eps = 0.3, 0.1, 0.03, 0.01; K = 2: eps = 0.3, 0.1, 0.03; near and far starts) |
| fitted slopes compared with theory | **done**: regret slope vs `d/2` (10/10 intervals contain it), MSE slope vs -1 (inside for all runs with `T eps >~ 700`), `MSE*T/(I^-1)_jj` = 1 +- 0.2, `E D_n` vs d, exit time and start cost vs prediction |
| including the burn-in | **done**: mid-window slope above `d/2`, `T_exit`, far-start excess |

Limits (stated, not hidden): (i) the estimator uses the clip, not the exact `F`-projection of `problem04_projection.md`; (ii) `eps = 0.01` has only `T eps = 180` at the horizon (`T = 1.8e4`, 64-96 replicates), so its error slopes (-1.6, -1.3) and constant (+0.8 drift) are not asymptotic; larger `T` would need `~ 1e9` events per replicate; (iii) the stationary start, the known `beta`, and the prior scale `rho_0 = 1` follow the earlier notes; other priors change the finite-T constant and the start cost (Problem 6); (iv) K = 3 and `k0 <= 1` were not run in this study (`k0 <= 1` is in Problem 6, section 13, K = 3 not asymptotic at `T = 3e4`); (v) 64-200 replicates give slope intervals of about +-0.13 (K = 1) and +-0.2 (K = 2): statements of the form "slope = d/2" are consistency checks at that resolution, not sharp estimates; (vi) the regret is the pseudo-regret `R_T` whose expectation is the KL regret (Problem 8; identity checked there to within 1 s.e.).
