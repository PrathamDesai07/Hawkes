# Problem 20: real-data check (Deribit BTC option trades; NSE / Bank Nifty tick data is not available here)

Tags: **[N]** measured here, **[H]** interpretation. Code: `code/p20_prep.py`, `code/p20_fit.cs`, `code/p20_run.sh`, `code/p20_analyze.py`; figures `figures/p20_*.png`; numbers `code/p20_summary.txt`; raw fits `code/p20_out/`. Date: 25 Oct 2026.
Tracker row 20: "Test on NSE or Bank Nifty tick or order-flow data. Fit the dictionary, report the estimated n, show the effect of misspecification, and show tracking behaviour across regime changes. Needs access to suitable tick data. Done when: a figure of n_hat(t) over a trading period, with a comment on regime shifts and on bias relative to a batch fit."

## 0. What data exist, and what was used

- **NSE / Bank Nifty tick data: not available.** `Desktop/BankNifty Options/icharts_option_banknifty_monthly_atm_up_to_5strikes.csv` (866 MB, 2019-2025) is **1-minute option-chain snapshots** (spot, future, ATM-strip prices, Greeks), not trade or order-flow events: event arrival times cannot be recovered from it. A previous attempt (24 Oct 2026, `BankNifty Options/hawkes/`) built *directional-change events from 1-minute future bars* and fitted a two-type marked Hawkes model; that is a different object (events are manufactured from minute bars, with tied events spread evenly inside the minute) and is not used here.
- **Tick data that do exist:** `Desktop/market-data/btc-1y/options_trades/`: **5,411,992 Deribit BTC option trades, 2024-11-01 .. 2025-10-31, millisecond timestamps**, all strikes and expiries (and 12 sample days of 5-level order-book updates at microsecond resolution, not used). This is genuine event-time data of the kind the tracker asks for, from a different market (24/7 crypto options, not NSE).
So this study is the real-data check **on BTC option trade arrivals**; the Bank Nifty part of the tracker row remains open (it needs tick or trade data from NSE).

## 1. The event stream

**Event** = a distinct millisecond timestamp over all instruments (trades sharing a millisecond are almost always fills of one aggressing order; they are one event). Total 5,411,992 trades -> **4,770,908 events** (tie fraction 11.8 %; full-day counts: 5th-95th percentile 5,243-24,469, minimum 3,337, median 12,058, the busiest day is 2025-10-10 with 36,596). 90 % of gaps are shorter than 18 s, the median is 1.5 s, the 99th percentile is 75 s. Block trades 4,700-15,000 per month, flagged liquidations 430-2,500 per month (not modelled). Intraday rate varies by a factor 3 (306-948 events per UTC hour in March 2025), weekends are quiet. **Window used for the main analysis:** 13:00-20:00 UTC (US trading hours) every day, which is the flattest part of the intraday profile; the state `S` is initialised from the events of the preceding 10 minutes and the likelihood includes the final partial gap. 365 daily windows, median 4,435 events (range 742-15,544).

## 2. Model, estimators

Exponential-dictionary Hawkes `lam(t) = mu + sum_k alpha_k S_k(t)`, `dS_k = -beta_k S_k dt + dN`, `beta` fixed (per second), `n_hat = sum_k alpha_k/beta_k`. **Batch:** constrained Newton MLE per window with backtracking (`mu >= 1e-6`, `alpha_k >= 1e-8`, `n <= 0.999`); `n_hat` standard error by the delta method from the observed information (valid only if the model is correct). **Online:** the forgetting Fisher-Newton of Problem 14, `F_i = (1-eta) F_{i-1} + DeltaF_i`, `theta_{i+1} = clip(theta_i - F_i^{-1} g_i)`, `eta = 1/(m N_prev)` with `N_prev` the events of the previous window (`m = 5`, `m = 1`), or `eta = 0` (growing window); one pass through the year, prior `F_0 = diag(1/lam_bar, lam_bar/beta_k^2)`.

## 3. Dictionary selection (train: days 0-29, validation: days 30-59, pooled fit, log-likelihood per event)

| dictionary beta (per s) | d | train nll/ev | validation ll/ev | n_hat (pooled) |
|---|---|---|---|---|
| 0.03 | 2 | 2.147 | -2.309 | 0.911 |
| 0.3 | 2 | 2.026 | -2.196 | 0.676 |
| 3 | 2 | 1.992 | -2.199 | 0.395 |
| 30 | 2 | 2.067 | -2.300 | 0.213 |
| (0.1, 1) | 3 | 1.955 | -2.121 | 0.713 |
| (1, 10) | 3 | 1.922 | -2.116 | 0.482 |
| (0.3, 3) | 3 | 1.927 | -2.103 | 0.608 |
| (3, 30) | 3 | 1.944 | -2.154 | 0.373 |
| **(0.1, 1, 10)** | 4 | 1.858 | **-2.017** | 0.730 |
| (0.03, 0.3, 3) | 4 | 1.885 | -2.041 | 0.823 |
| (0.3, 3, 30) | 4 | 1.853 | -2.026 | 0.623 |
| (0.1, 1, 10, 100) | 5 | 1.819 | -1.982 | 0.728 |

Held-out log-likelihood improves monotonically with the number of components and with the span of time scales (the trade clustering is multi-scale, not a single exponential), so **no finite dictionary is correct**; the estimate `n_hat` for K = 1 is **not even stable in the choice of beta** (0.21 to 0.91 for a single exponential). I use `beta = (0.1, 1, 10)/s` (K = 3, d = 4) as the main dictionary (decade-spaced, covers 0.1-10 s) and compare K = 1, 2, 4 for misspecification.

## 4. Results

**Figure `figures/p20_nhat_year.png` (the requested n_hat(t)):** daily batch `n_hat` (points, 95 % delta-method bands) with the online forgetting estimate (red) and the growing-window estimate (blue) over 365 days, with events per window below.

**Level and variation [N].** `n_hat` per day (K = 3): median **0.668**, 10-90 % range 0.561-0.768, minimum 0.440 (2025-09-27), maximum 0.930 (2025-03-02); median s.e. 0.017. The day-to-day dispersion is **5.0 times the s.e.** (sd of `(n_hat - median)/s.e.`): the variation is far larger than estimation noise, i.e. a real regime variation *or* model-misspecification variation (they cannot be separated from this data alone). Correlations: `corr(log events/window, n_hat) = +0.77`, `corr(log events, log mu_hat) = +0.92`: busy days have both more baseline activity and a larger branching ratio. Weekend median `n_hat` 0.596 vs weekday 0.684 (events per window 1,930 vs 5,138). Monthly medians 0.61-0.70 with no trend (Nov 2024 0.61, Apr 2025 0.70, Sep 2025 0.66).
**Regime shifts [N, interpretation H].** The six busiest days: 2025-03-02 (15,544 events, `n_hat` 0.930, `mu_hat` 0.043), 2025-10-10 (14,591, 0.797, 0.118), 2024-11-13 (14,548, 0.770, 0.133), 2025-10-03 (13,989, 0.781, 0.121), 2025-04-22 (13,341, 0.830, 0.090), 2025-05-21 (13,317, 0.811, 0.100). The two events with the largest rate responses are not alike: **on 2025-10-10 the activity spike enters mainly through the baseline** (`mu_hat` 0.118 vs a typical 0.05, `n_hat` 0.80, near the busy-day level) whereas **on 2025-03-02 the clustering itself jumps** (`n_hat` 0.93, the maximum of the year, with a *low* `mu_hat`): endogenous burst vs exogenous surge. (That these dates coincide with well-known crypto-market events is from my background knowledge, not from the data; I have not checked news sources.)

**Misspecification [N].** The same daily windows fitted with four dictionaries:

| dictionary | median n_hat | 10-90 % | in-sample nll/ev (median) | median KS | days with KS > critical | median lag-1 corr of rescaled gaps |
|---|---|---|---|---|---|---|
| K=1, beta = 1/s | 0.501 | 0.361-0.603 | 2.176 | 0.117 | 100 % | 0.113 |
| K=2, (0.3, 3) | 0.554 | 0.433-0.660 | 2.098 | 0.092 | 100 % | 0.094 |
| K=3, (0.1, 1, 10) | **0.668** | 0.561-0.768 | 1.998 | 0.057 | 98 % | 0.067 |
| K=4, (0.1, 1, 10, 100) | 0.664 | 0.560-0.766 | 1.935 | **0.016** | **15 %** | 0.072 |

(a) `n_hat` **increases with the dictionary and saturates**: the K = 1 value is 0.17 below K = 3 on the median day (IQR 0.14-0.20), the K = 3 and K = 4 values agree (0.668 vs 0.664). This is exactly the bias-versus-K picture of Problem 18 (the KL pseudo-truth `n(theta_dagger)` of a too-small dictionary is biased downward for a heavy-tailed kernel). (b) The time-rescaling test (compensator increments should be i.i.d. Exp(1)) **rejects every dictionary with K <= 3 on 98-100 % of days; K = 4 (with the 100/s component, the millisecond scale) passes on 85 % of days**: the missing ingredient in K = 3 is the sub-100-ms time scale, not the long memory. (c) The delta-method standard errors (0.017 daily) assume the model is correct; for K <= 3 they are therefore optimistic.
**Seasonality as misspecification [N].** Fitting K = 3 to the full 24 h day instead of the 13-20 UTC window changes the median `n_hat` only from 0.668 to 0.682 and the KS from 0.057 to 0.053 (`mu_hat` 0.043 vs 0.056): at the time scales of the dictionary (0.1-10 s) the 24-hour cycle is absorbed by the baseline; the misspecification that matters is the dictionary (the time-scale content), not the intraday cycle.

**Tracking and bias relative to a batch fit [N].**

| online estimator | mean(online - daily batch) | sd of difference | corr with daily batch | sd of online `n_hat(t)` (batch: 0.079) |
|---|---|---|---|---|
| forgetting `eta = 1/(5 N_prev)` | +0.077 | 0.075 | 0.34 | 0.034 |
| forgetting `eta = 1/N_prev` | +0.039 | 0.069 | 0.50 | 0.049 |
| growing window (`eta = 0`) | +0.090 | 0.079 | 0.14 | 0.015 |

(i) **The growing-window online estimate equals the pooled batch MLE**: after 30 / 120 / 365 days, pooled batch `n_hat` = 0.7300, 0.7600, 0.7590 vs online 0.7287, 0.7576, 0.7565 (difference **-0.0013, -0.0024, -0.0025**): the one-pass recursion reproduces the batch fit on real data, as on simulated data (Problem 21). (ii) The +0.04 to +0.09 offset to the *daily* batch median is **not estimator error**: a pooled (or long-memory) fit weights the busy days, and busy days have larger `n_hat` (event-weighted mean of the daily values 0.697, pooled MLE 0.759, median daily 0.668): the difference is regime heterogeneity. (iii) The forgetting estimators are smoothed versions of the daily series (sd 0.034-0.049 against 0.079); with a window of about one day (`eta = 1/N_prev`) they correlate 0.50 with the daily batch estimate, with a five-day window 0.34: they follow the slow variation (the peak at 2025-03-02 and the clusters of busy days are visible in the red curve of the figure) but not the day-to-day scatter, which is the expected bias-variance trade-off of a forgetting window (Problem 14: error `≈ eta d/2 + (drift/eta)^2`).
**Intra-day tracking [N].** `figures/p20_rolling.png`: hourly batch fits (K = 3) around 2025-10-10, around 2025-03-02 and in a quiet week (2025-09-15..25): median hourly `n_hat` 0.619, 0.628 and 0.563, median s.e. 0.045-0.055, sd of `n_hat`/s.e. = 2.0-2.4; max events per hour 5,813, 4,487, 2,244. At one-hour resolution the estimates are too noisy (s.e. about 0.045 on an effect of 0.1-0.2) to resolve the regime change *within* a day; the daily window is the resolution at which the regime structure is visible.

## 5. Reading against the theory [H]

- `n_hat` in this data is 0.5-0.8: **far from criticality** (`eps = 1 - n_hat = 0.2-0.5`), the regime where Problems 5-6, 19 give a fast burn-in (`T_exit ≈ 1.4 max_j (I^{-1})_jj/(delta theta_j)^2`) and where the near-critical results (the `eps -> 0` dichotomy of Problem 11) are not tested. The nearly-critical regime would need either a different asset/time scale or another dictionary (`n_hat` approaches 1 for power-law-like kernels as `K` grows with a slow component; here the slowest component is 0.1/s, i.e. 10 s).
- The finding that `n_hat` is stable between K = 3 and K = 4 (0.668 / 0.664) while K = 1, 2 are lower is the real-data counterpart of the Problem 18 conjecture (bias of `n(theta_dagger)` decays with `K`), here already saturated at K = 3 for this stream.
- The online estimator behaves on real data as on simulated data on the points that can be tested without the truth: it reproduces the batch MLE (-0.0013..-0.0025), it is cheap (one pass over the 1.6M window events, including the Fisher quadrature, took 18 s in single-threaded C#), and its smoothing is controlled by `eta`.

## 6. Limits (stated, not hidden)

(i) **Not the NSE/Bank Nifty data**: BTC option trades, a different market microstructure (24/7, single venue, trade events not order-flow with marks). (ii) Event definition = distinct millisecond, 11.8 % of trades collapsed; block trades and liquidations are ordinary events here. (iii) The state is reset/initialised from 10 minutes before each window; windows start at 13:00 UTC with the intraday rate not exactly flat. (iv) Standard errors are model-based and optimistic for K <= 3 (KS rejects on 98-100 % of days); the 5.0x dispersion statement is therefore an upper bound on "real" regime variation. (v) `beta` is fixed by a coarse grid chosen on the first 60 days; the K = 3 vs K = 4 comparison is the only sensitivity check. (vi) No independent ground truth exists, so "bias relative to a batch fit" is the only available bias statement; the baselines of Problem 21 (Hall-Willett-type, Solo-type) were *not* run on this stream here. (vii) Regime-shift attributions to market events are from background knowledge and were not verified.

## 7. Status

| "Done when" (tracker row 20) | status |
|---|---|
| a figure of n_hat(t) over a trading period | **done on BTC option trades**: `figures/p20_nhat_year.png` (365 daily windows), `p20_rolling.png` (hourly), `p20_misspecification.png`, `p20_seasonality_activity.png` |
| a comment on regime shifts | **done** (section 4): variation 5x the standard error, `corr(n_hat, log activity) = +0.77`, two distinct burst types (2025-03-02 clustering jump, 2025-10-10 baseline surge) |
| a comment on bias relative to a batch fit | **done**: growing-window online = pooled batch MLE to 0.0025; offset to the daily batch explained by heterogeneity and window length |
| NSE / Bank Nifty data | **open**: no tick data available; the minute-bar CSV cannot supply event times |
| baselines of Problem 21 on the real stream | **open** (protocol in `problem21_baselines.md`; the event file `code/p20_data/events.bin` makes it a short run) |
