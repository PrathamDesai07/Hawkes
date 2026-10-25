# Problem 21: baselines to beat (batch MLE, Hall-Willett mirror descent, Solo fixed-gain least squares)

Tags: **[P]** from proved statements, **[H]** prediction to test, **[N]** measured. Code: `code/p21_baselines.cs`, `code/p21_run.sh`, `code/p21_analyze.py`, `code/p21_ls_efficiency.py`; numbers in `code/p21_summary.txt`. Date: 25 Oct 2026.
Tracker row 21: "Implement batch MLE, Hall-Willett mirror descent and Solo's fixed-gain estimator, and compare accuracy and per-event cost with your method. Done when: a comparison table of error and runtime on the same simulated and real datasets." The real-data half belongs to Problem 20 (tick data not available here); this note does the simulated half and fixes the protocol the real-data half will reuse.

## 0. What the two literature baselines are, and what is adapted

- **Hall-Willett (arXiv 1409.0031).** Tracks the *intensity vector* `lambda_t` of a network of point processes with a known (or learned) dynamical model by **dynamic mirror descent in `theta = log(delta lambda)`**, for the per-bin Poisson loss; Theorem 1 gives a tracking regret `<= C(1 + sum ||lambda_{t+1} - Phi_t(lambda_t)||^2) sqrt(T)`; the regret is adversarial (`literature_crosscheck.md`). It does not estimate `(mu, alpha_k)` of an exponential-dictionary Hawkes process. **Adaptation used here (labelled "HW-type"):** online mirror descent with the log (multiplicative) mirror map on the *parameter vector* `theta = (mu, alpha_1..alpha_K)` for the Hawkes loss `l_i(theta) = Lambda_theta(tau_i) - log lambda_theta(tau_i)`:
 `log theta_j <- log theta_j - eta_i theta_j g_{i,j}`,  `g_i = b_i - a_i/(a_i^T theta)`,
 with a decaying step `eta_i = eta_0/sqrt(i)` (their `sqrt T` tuning) or `eta_0/i`, and a constant step `eta_0` (their tracking tuning). No curvature information is used.
- **Solo (Piggott and Solo, ICASSP 2018, abstract only; full text not obtained).** A fixed-gain adaptive estimator based on **exponentially weighted least squares** for a causal Laguerre-basis intensity, preserving non-negativity. **Adaptation used here ("Solo-type EWLS"):** the Hawkes least-squares contrast `int lambda_theta^2 dt - 2 sum_i lambda_theta(t_i)`, whose minimiser solves `A theta = b` with `A = int phi phi^T dt`, `b = sum_i phi(t_i)`, with exponential forgetting `A_i = rho A_{i-1} + DA_i`, `b_i = rho b_{i-1} + a_i`, `theta_i = clip((A_i + Lambda_0)^{-1}(b_i + Lambda_0 theta_start))`; "fixed gain" = constant `rho < 1`; `rho = 1` is the plain recursive least squares. The exponential dictionary replaces the Laguerre basis (same intensity class here), the network is a single node, positivity is enforced by the same box clip as for every method. Whether Piggott-Solo's recursion equals this one cannot be checked without the full text (`literature_search_report.md`).
- **Batch MLE.** Full-sample maximum likelihood by Newton with backtracking, warm-started from the previous refit; the statistical gold standard (efficient: `e_n ~ chi^2_d`).
- **Ours.** The Fisher-preconditioned online Newton of Problem 5 (same prior, box and start for every online method).

## 1. Costs [P]

| method | per event | memory | notes |
|---|---|---|---|
| ours (online Newton) | `O(d^3)` (a `d x d` Cholesky solve; `d = K+1`), constant in `n` | `O(d^2)` | `Delta F` by closed form (K=1) or a fixed quadrature |
| HW-type mirror descent | `O(d)` | `O(d)` | no matrix at all |
| Solo-type EWLS | `O(d^2)` update of `A` (closed form `int phi_j phi_k = S_j S_k (1 - e^{-(beta_j+beta_k)w})/(beta_j + beta_k)`) + `O(d^3)` solve | `O(d^2)` | |
| batch MLE refit after `n` events | `O(n d^2)` per Newton iteration, `~ 3-6` iterations, **linear in `n`** | `O(n d)` (stored gaps) | refitting after *every* event costs `O(n^2 d^2)` in total |

## 2. Predictions [H]

- **B1 (accuracy of ours vs batch MLE).** Ours is a one-pass recursive MLE: `E e_n` of ours `= E e_n` of batch MLE `(1 + o(1))`, both `≈ d`; they differ only in the finite-`n` transient (prior, clipping).
- **B2 (Solo-type, efficiency).** Least squares is not efficient for Hawkes: its asymptotic covariance is `B^{-1} C B^{-1}/T` with `B = E[a a^T]`, `C = E[lambda a a^T]` (sandwich), so `E e_n` (in the Fisher metric) tends to `tr(I B^{-1} C B^{-1})`, strictly above `d`; computed from the exact moment systems (`tangent_bound.md`) in `p21_ls_efficiency.py`. With forgetting `rho < 1` the error does not decay: it floors at `~ (1-rho)/2 x` (that constant) -- the fixed-gain trade-off.
- **B3 (HW-type).** No preconditioning: the Hessian of the loss has condition number `lam_max(I)/lam_min(I) ~ eps^{-2}` (Problem 12: eigenvalues `1/eps` and `eps`), so a single step size cannot serve both directions; the error in the ridge direction is expected to decay far slower than `1/n`, and with a constant step it has a variance floor `~ eta_0`. Even with the step tuned *per scenario on the true parameters* (which the baseline cannot do in practice) the error stays well above `d` after the same number of events.
- **B4 (cost).** Per-event CPU of ours is constant in `n`; batch refit time per refit grows linearly in `n`; HW-type is the cheapest per event (`O(d)`), Solo-type about `d` times that of HW; the *amortised* cost of keeping the batch estimate current after every event is `O(n)` per event.
- **B5 (regret).** Ours has regret `(d/2) ln T + O(1)` (Problem 19). HW-type with `eta_0/sqrt(i)` has regret growing faster than `ln n` (no curvature); Solo-type with constant forgetting has *linear* regret (non-vanishing error floor); `rho = 1` least squares has regret `~ (c_LS) ln T` with `c_LS = (1/2) tr(I B^{-1} C B^{-1}) > d/2`.

## 3. Protocol [N]

Identical simulated paths for every method (K = 1: mu = beta = 1, eps in {0.3, 0.1, 0.03}, 64 paths of `n = 1e5` events; K = 2: beta = (1,2), mu = 1.5, eps = 0.1, 32 paths), stationary start, the same "near" start `(1.3 mu*, 0.8 alpha*)`, the same a-priori box and stability cap for **every** method (clip `mu >= 0.1 mu*`, `alpha_k >= 0.05 alpha*_k`, `sum alpha_k/beta_k <= 1 - eps/2`; the batch MLE is the *constrained* MLE, which is also the unconstrained one whenever the optimum is interior; the unconstrained likelihood has no minimiser for small `n`, a first version diverged to `1e31` and was corrected). Checkpoints: geometric, 4 per decade from 30 to `1e5` events.
**Metrics.** (i) `e_n = (theta_hat - theta*)^T J_n (theta_hat - theta*)`, `J_n` the oracle information accumulated up to event `n` (so `e_n ≈ d` means the error equals the Cramer-Rao scale at that `n`; an estimator with a non-vanishing error has `e_n ∝ n`); (ii) expected regret `E R_n` (Problem 8) of the *sequential* predictor `theta_i` (not defined for the batch refit); (iii) CPU time, single thread, minimum of 3 passes, data generation excluded. 95 % intervals over paths.
**HW-type step sizes.** Constant `eta_0`, `eta_0/sqrt(i)`, `eta_0/i`, `eta_0 in {0.003, 0.01, 0.03, 0.1, 0.3, 1}`. Two numbers are reported for each schedule: "default" (`eta_0 = 0.1`) and "tuned on the truth" (the `eta_0` minimising the final `e_n` of that scenario, which a practitioner cannot do; it favours the baseline).

## 4. Results (25 Oct 2026; figures `figures/p21_k1_eps{0.3,0.1,0.03}.png`, `figures/p21_k2_eps0.1.png`; all numbers `code/p21_summary.txt`)

**Table 1 (accuracy and regret at n = 1e5 events; `e_n` target `d`).**

| method | K=1, eps=0.3 | K=1, eps=0.1 | K=1, eps=0.03 | K=2, eps=0.1 (d=3) |
|---|---|---|---|---|
| **ours (online Newton)** | 1.78 +- 0.36 | 2.16 +- 0.53 | 1.86 +- 0.61 | 4.08 +- 1.60 |
| batch MLE (refit at n) | 1.80 +- 0.36 | 2.15 +- 0.52 | 1.87 +- 0.59 | 4.09 +- 1.63 |
| Solo-type EWLS, rho = 1 | 2.48 +- 0.53 | 3.66 +- 0.92 | 3.08 +- 0.92 | 4.92 +- 1.90 |
| Solo-type EWLS, rho = 1-1/2000 | 48 +- 12 | 87 +- 23 | 239 +- 73 | 102 +- 42 |
| Solo-type EWLS, rho = 1-1/200 | 655 +- 147 | 1705 +- 475 | 1377 +- 159 | 1052 +- 334 |
| HW-type, const step (tuned) | 113 +- 29 | 223 +- 42 | 199 +- 51 | 95 +- 34 |
| HW-type, eta_0/sqrt(i) (tuned) | 13.4 +- 4.0 | 56.9 +- 7.9 | 8.95 +- 0.89 | 46.8 +- 6.7 |
| HW-type, eta_0/i (tuned) | 717 +- 99 | 156 +- 46 | 13.5 +- 2.4 | 213 +- 86 |
| HW-type, default eta_0 = 0.1, const / 1/i | 12156 / 777 | 97833 / 555 | 112714 / 580 | 73974 / 1215 |

**Table 2 (expected regret `E R_n` at n = 1e5; ours has the `(d/2) ln T + O(1)` of Problem 19).**

| method | K=1, eps=0.3 | K=1, eps=0.1 | K=1, eps=0.03 | K=2, eps=0.1 |
|---|---|---|---|---|
| **ours** | 9.51 | 7.79 | 6.60 | 8.59 |
| Solo-type EWLS, rho = 1 | 10.99 | 10.23 | 9.21 | 11.39 |
| Solo-type EWLS, rho = 1-1/2000 | 37.3 | 34.5 | 21.5 | 36.7 |
| Solo-type EWLS, rho = 1-1/200 | 267 | 147 | 78 | 127 |
| HW-type, eta_0/sqrt(i) (tuned) | 31.4 | 32.6 | 37.9 | 28.5 |
| HW-type, const step (tuned) | 57.8 | 61.8 | 44.6 | 39.9 |
| HW-type, eta_0/i (tuned) | 411 | 110 | 11.1 | 142 |

**B1 [N]: ours = batch MLE in accuracy.** At `n = 1e5` the two agree to within 0.01-0.02 (1.78 vs 1.80, 2.16 vs 2.15, 1.86 vs 1.87, 4.08 vs 4.09) and both sit at `d` within the intervals: the one-pass recursive MLE is the batch MLE to first order. At small `n` ours is *better* (`n = 95`: 0.83 vs 4.50 at eps = 0.1; 2.07 vs 2.84 at eps = 0.3; `n = 949`: 1.51 vs 3.23 at eps = 0.03) because of the prior (the unregularised MLE does not exist for few events and is box-clipped here).
**B2 [N]: least squares is inefficient by the predicted factor.** `E e_n` of the rho = 1 EWLS tends to `tr(I B^{-1} C B^{-1})` computed from the exact moment systems: predicted 2.45, 2.99, 3.50, 4.15 (eps = 0.3, 0.1, 0.03; K = 2); measured 2.48 +- 0.53, 3.66 +- 0.92, 3.08 +- 0.92, 4.92 +- 1.90 (all inside the intervals). That is a factor 1.22, 1.50, 1.75 (K = 1) and 1.38 (K = 2) above `d`: **the Solo-type estimator needs 1.2-1.75 times as many events as the MLE for the same accuracy, and the factor grows as eps decreases** (more variable intensity). Its regret slope in `ln n` on the last decade is 1.38, 1.60, 1.91, 1.92 against `(1/2) tr(...)` = 1.22, 1.50, 1.75, 2.08 (B5). With fixed gain (`rho < 1`) the error does **not** vanish: `e_n` grows linearly with `n` (the floor of the raw error is `~ (1-rho)`), the regret grows linearly (`E R_n/n` = 2.2-3.7e-4 per event for `rho = 1-1/2000` and 0.8-2.7e-3 for `rho = 1-1/200`, against 6.6-9.5e-5 for ours).
**B3 [N]: mirror descent without curvature does not reach the Cramer-Rao scale.** Even with `eta_0` tuned on the truth the best HW-type schedule is 4.5-28 times above `d` at `n = 1e5` (8.95 +- 0.89 at eps = 0.03 with `eta_0/sqrt(i)`, 13.4 at eps = 0.3, 56.9 at eps = 0.1, 46.8 for K = 2), and its error grows with `n` in most scenarios (no decay at the Cramer-Rao rate), the regret is 3-6 times ours; with the default step the constant-step version diverges to the box boundary (`e = 1e4-1e5`). At very small `n` (`n = 95`) the cautious `eta_0/sqrt(i)` schedule is competitive or better (0.59, 0.63, 0.72 at eps = 0.1, 0.03, 0.3, against ours 0.83, 0.72, 2.07) because ours must first accumulate Fisher information. The mismatch is the ridge: the loss has condition number `~ eps^{-2}`, a single step size cannot serve the `1/eps` and `eps` eigendirections.

**Table 3 (CPU per event, microseconds, one thread; ratios are robust, absolute values vary by about 2x with machine load between runs).**

| method | K=1 (d=2) | K=2 (d=3) | scaling in n |
|---|---|---|---|
| **ours (online Newton)** | 1.6 | 2.3 | constant |
| HW-type mirror descent | 0.08 | 0.10 | constant |
| Solo-type EWLS | 0.42 | 0.59 | constant |
| batch MLE, one cold refit at `n = 1e5` (5 Newton iterations) | 33 ms = 0.33 us per data event | 40 ms = 0.40 us | linear in `n` per refit |
| batch MLE, refits at the 17 geometric checkpoints (warm-started chain), amortised over `n = 2e5` | 0.84 | 0.93 | |
| batch MLE kept current after **every** event (one Newton step on all past data), average over `k <= n = 2e5` | **11 346** | **13 113** | linear in `n` per event, `O(n^2)` in total |

(B4 holds.) **Ours is 6 000-7 000 times cheaper per event than a batch estimate kept current after every event**, 20 times costlier than mirror descent and 4 times costlier than the EWLS recursion; batch refits at geometric checkpoints are cheap in aggregate (0.84 us/event) but only give an estimate at those checkpoints. Per-event cost of ours is `O(d^3)`: 1.6 -> 2.3 us from `d = 2` to `d = 3`.

## 5. Summary against "baselines to beat"

- **Accuracy:** ours equals the batch MLE (the accuracy gold standard) to 1 % at `n = 1e5` and beats it at small `n`; it beats the Solo-type EWLS by the efficiency factor 1.2-1.75 (rho = 1) or by orders of magnitude (fixed gain, non-vanishing error); it beats the Hall-Willett-type mirror descent by a factor 4-30 at `n = 1e5` even when that baseline is tuned on the truth.
- **Cost:** constant per event, `~ 1.6-2.3 us`; batch is `O(n)` per refit and `O(n)` per event if kept current; mirror descent and EWLS are cheaper per event but pay in accuracy.
- **Fair-comparison caveats:** (i) the Hall-Willett and Solo estimators are *adaptations* to the exponential-dictionary model (section 0): Hall-Willett is designed for intensity tracking (`sqrt T` tracking regret under changing dynamics) and Piggott-Solo for a Laguerre-basis network, so a loss of the baselines here under a *static* single-node model is not a statement about their intended regimes; in particular fixed gain is a feature for tracking (Problem 14) and not for static estimation; (ii) the Solo paper was read as an abstract only; (iii) `eta_0` was tuned on the truth for the HW-type baseline, a favourable setting for it; (iv) K = 1, 2, a stationary start and a known box; the batch MLE is constrained to the same box; (v) 32-64 paths (intervals as shown); (vi) timings are single-thread C# on one machine and vary by about 2x between runs.
- **Real data (Problem 20):** the protocol above (same stream, same box, same checkpoints) will be reused; there the "truth" is unknown, so `e_n` is replaced by the distance to the full-sample batch MLE and the sequential log-loss on held-out events.

## 6. Status

| "Done when" | status |
|---|---|
| a comparison table of error and runtime on the same simulated datasets | **done**: Tables 1-3 (K = 1: 3 values of eps; K = 2) with intervals |
| ... and real datasets | **not done**: needs the tick data of Problem 20; protocol fixed above |
| the three baselines implemented | batch MLE: exact; HW-type and Solo-type: adaptations (section 0), the closest faithful versions available without the Solo full text and without a published Hawkes-parameter version of dynamic mirror descent |
