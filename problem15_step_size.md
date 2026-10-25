# Problem 15: optimal constant step size and its ε-dependence

**Status: draft. Conditional on Problem 13.** Problem 15 depends on the tracking bound of Problem 13, which is not proved. Everything here builds on the *measured* structure of the tracking error from Problem 14 (noise term a(ε)η, plus the response of a time-varying-gain filter to the drift). The formulas are derived from that structure and checked by simulation; they are not consequences of a theorem. Evidence is K = 1 (α known and α unknown), μ = β = 1, at ε between 0.2 and 0.4.

## 1. Setup and the two terms

Per-event gain η, event rate r = μ/ε, drift measured in the parameter μ. Average tracking MSE ≈ **a η + B(η)**:
- noise: a = r/(2·eff) = μ/(2ε·eff), where eff is the information about μ (eff = m = E[1/λ] if α is known; eff = (m − ε/μ)/(1 − 2ε + εμm) if α is estimated). Measured a_meas ≈ 1.07–1.15 × this. For k = 2μ/β > 1, as ε → 0: m ≈ 2ε/(2μ−β), eff ≈ εβ/(μ(2μ−β)), so **eff ≈ m/k** (estimating α costs a factor k = 2μ/β in information about μ) and a ∝ ε^{−2}. For k < 1, eff ~ ε^k and a ~ ε^{−(1+k)} (Problem 10).
- drift: B(η) = the squared lag of the filter dy/dt = η r_t(μ_t − y), r_t = μ_t/ε (Problem 14, R2). For a sinusoid of relative amplitude A and period P, per-event frequency ω = 2π/(P r) = 2πε/(Pμ): B ≈ (A²μ²/2) ω²/(η² + ω²). For a jump of size Dμ every L events: B ≈ (Dμ)² tanh(ηL/2)/(2ηL).

## 2. Closed-form optimal steps (derived)

**Smooth drift, fast-tracking regime (η* ≫ ω).** B ≈ Δ_e²/(2η²), with Δ_e = Aμω the per-event slope amplitude (Δ_e = Δ·ε/μ for a drift rate Δ per unit time). Minimising aη + Δ_e²/(2η²):

η* = (Δ_e²/a)^{1/3},  MSE* = (3/2) a η*.

With a = μ/(2ε eff) and eff = κ ε (k > 1; κ = 2/(2μ−β) for α known, β/(μ(2μ−β)) for α unknown):
**η* = 2 (π²A²κ/(P²μ))^{1/3} ε^{4/3}** for a sinusoid (equivalently η* ∝ Δ^{2/3} ε^{4/3}). So the optimal per-event gain falls like ε^{4/3}; in time units the tracking bandwidth η* r falls only like ε^{1/3} (the event rate grows like 1/ε), and the averaging window 1/η* grows like ε^{−4/3} events. For k < 1 (eff ~ ε^k): η* ∝ ε^{1+k/3}.

**Jump drift.** Relaxation faster than the regime (η*L ≫ 1): minimise aη + (Dμ)²/(2ηL):
η* = Dμ/√(2aL_e) = D·√(ρ_t eff)·ε, with L_e = L_t r events per regime and ρ_t = 1/L_t jumps per unit time. So η* ∝ ε^{3/2} for k > 1 (ε^{1+k/2} for k < 1), ∝ D, ∝ L_t^{−1/2}.

**Window cap.** The linear (small-noise) regime needs aη/μ² ≲ θ with θ ≈ 0.013–0.03 (Problem 14, F4): η ≤ η_max = θμ²/a = 2θ ε eff μ ∝ ε². Since η* ∝ ε^{4/3} (or ε^{3/2}) falls more slowly than η_max ∝ ε², the cap becomes binding as ε → 0 unless the drift is slow enough: for the sinusoid, η* ≤ η_max iff P·eff ≥ πA/(θ^{3/2}μ²) (≈ 1100 A/μ² at θ = 0.02). Near criticality, drift on timescales shorter than a few hundred times 1/eff cannot be tracked in the linear regime at all.

## 3. When is tracking worth doing at all? (derived, then checked)

A static (η → 0) estimator has MSE ≈ A²μ²/2 for a sinusoid. In the fast-tracking regime the tracker's best MSE is (3/2)aη* = (3/2) a^{2/3}(A²μ²ω²)^{1/3}. It beats the static estimator iff a ω ≲ (Aμ)²/3^{3/2}, and aω = π/(P·eff), so

**tracking helps iff P · eff · (Aμ)² ≳ 3^{3/2} π ≈ 16**:

(information accumulated over one drift period) × (squared drift amplitude) must exceed about 16, i.e. the drift amplitude must exceed roughly 4 standard errors of the estimate over one period. Below that, the best constant gain is η → 0 (give up tracking). This threshold is a consequence of the asymptotic formula and is only as good as it; I checked it at three points below.

## 4. Simulation checks (K = 1, μ = β = 1, forgetting-Newton recursion, causal sampling, per-replicate MSE; prediction = c_a r/(2 eff)·η + time-varying-filter bias with c_a = 1.14 for α known and 1.15 for α unknown)

`code/p15_eta_opt.py`, `code/p15_jump_opt.py`, `code/p15_alpha_unknown.py`.

**(a) Smooth drift, α known** (A = 0.2, P = 1200, T = 6000, 12–16 replicates): the full MSE(η) curve agrees with the prediction within 3–7 % at ε = 0.4 (e.g. at η = 0.002, 0.004, 0.008: simulated 0.0159, 0.0154, 0.0245 against 0.0166, 0.0165, 0.0256). Optimal step, simulated (parabola in log η) vs predicted argmin vs the clean asymptotic formula:

| ε | simulated η_opt | predicted argmin | asymptotic (Δ_e²/a)^{1/3} | η_opt/ω_e |
|---|---|---|---|---|
| 0.4 | 0.00292 | 0.00295 | 0.00387 | 1.4 |
| 0.3 | 0.00179 | 0.00193 | 0.00269 | 1.1 |
| 0.2 | 0.00089 | 0.00095 | 0.00162 | 0.85 |

The simulated optimum follows the prediction within 7 %. It falls between ε = 0.4 and 0.2 by a factor 3.3 (log-slope 1.7), steeper than the asymptotic ε^{4/3} (factor 2.5) because at fixed drift period the optimum moves from the tracking regime (η/ω_e = 1.4) toward the saturated regime (η/ω_e < 1) as ε falls, where the clean power law does not apply. I did not run the regime η* ≫ ω_e for several ε (it needs P of order 10⁴ time units).

**(b) Jump drift, α known** (D = 0.3, regime length L_t = 800 time units, alternating levels):

| ε | L_e | simulated η_opt | predicted argmin | asymptotic Dμ/√(2aL_e) |
|---|---|---|---|---|
| 0.4 | 2000 | 0.00195 | 0.00248 | 0.00272 |
| 0.2 | 4000 | 0.00077 | 0.00090 | 0.00104 |

Simulated η_opt is 15–27 % below the prediction (the simulated MSE is above the prediction once η is beyond the linear window, which pulls the optimum lower). The simulated ratio between the two ε is 2.5 (log-slope 1.34) against the predicted ε^{3/2} (2.8). Regime length (ε = 0.4): L_t = 1600 gives 0.00163 against 0.00195 at L_t = 800 (ratio 0.84; the asymptotic L^{−1/2} gives 0.71; predicted argmin 0.77). A smaller jump (D = 0.15, L_t = 800) put the simulated optimum at the edge of the grid (0.0005; predicted 0.00025): there ηL_e ≈ 1 to 3, the saturated regime, as the theory says.

**(c) α unknown (adaptive recursion), ε = 0.4, P = 1200:**
- **A = 0.2** (P·eff·(Aμ)² = 8.7, below 16): there is no interior optimum. Simulated MSE is 0.0219, 0.0233, 0.0243, 0.0247, 0.0266 at η = 0.0002, 0.0004, 0.0007, 0.001, 0.002, rising monotonically from the smallest step tried, while the static value is A²μ²/2 = 0.02. The prediction (0.0218 … 0.0265) agrees within 3 %. For α known at the same A the criterion value is 22.6 (above 16) and tracking does help: best MSE 0.0154 against 0.02.
- **A = 0.4** (criterion 34.7, above 16): an interior optimum appears at simulated η = 0.0033 (predicted 0.0038), best MSE 0.059 against the static 0.08 (26 % lower); the α-known prediction for the same drift is 0.0058.
So the criterion separates the two cases correctly at all three points (8.7: no gain; 22.6: modest gain; 34.7: 26 % gain), and estimating α does shift the optimum down (the asymptotic factor is k^{−1/3}; predicted 0.0038 vs 0.0058 at finite ε).

## 5. Practical tuning rule

1. **Do not fix η across ε.** The noise MSE aη is ∝ η/ε²: η = 0.01 is about 16 times noisier at ε = 0.1 than at ε = 0.4. Scale η like ε^{4/3} (slow drift) or ε^{3/2} (regime shifts), or better compute it from the formula below with the current estimates.
2. Estimate ε̂ = 1 − n̂ and the information eff from the recursion's own information matrix F (eff ≈ the μ-μ entry of F⁻¹ inverted, per event times r), so a ≈ r/(2 eff).
3. Estimate the drift scale (Δ_e for smooth drift, or D and the regime length for jumps).
4. **Check tracking is worthwhile**: P · eff · (Aμ)² ≳ 16 (for a sinusoid; for generic drift, drift variance over a period × information accumulated over that period). If not, use η → 0 (a decaying gain).
5. Otherwise set η = min{ (Δ_e²/a)^{1/3} [smooth] or Dμ/√(2aL_e) [jumps], θ μ²/a }, θ ≈ 0.013–0.03, and use the forgetting version of the Newton recursion (Problem 14: Problem 5's recursion with F accumulating forever cannot track).
6. If α is estimated, expect the step to be smaller by about k^{−1/3} and the tracking threshold harder to meet (information is divided by k).

## 6. What is not established

- The formulas rest on the measured structure (a(ε), the time-varying-filter response), not on a proved bound (Problem 13 pending, and Problem 5's analysis is itself in progress).
- Only K = 1 was simulated for the step-size optimum; K ≥ 2 follows through eff and the window, but η_opt was not run for K ≥ 2.
- The clean asymptotic laws η* ∝ ε^{4/3} and ε^{3/2} were not reached in the smooth-drift runs (optimum near or below ω_e); the scaling checks are consistent with the prediction including finite-ε effects, not with the pure power law. The jump runs cover two values of ε; no run at ε < 0.2.
- Constants (c_a ≈ 1.14, θ ≈ 0.013–0.03, the criterion 16) are empirical or follow from the asymptotic formula; the amplitude effect of the adaptive preconditioner (Problem 14, R1) makes α-unknown tracking up to 1.3–1.75 times worse in bias than the filter formula at amplitudes 0.2–0.5, which shifts the optimum slightly lower than predicted (A = 0.4: simulated 0.0033 against 0.0038).
- The k < 1 exponents (ε^{1+k/3}, ε^{1+k/2}) follow from eff ~ ε^k and were not simulated.
- Re-derive the minimisations in §2 and §3 yourself; they are two-line calculus but the criterion in §3 depends on the asymptotic regime.
