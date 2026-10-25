# Problem 10: Fisher information for K = 1 as n → 1

Model: λ = μ + αS, dS = −βS dt + dN, n = α/β = 1 − ε, a = (1, S), θ = (μ, α). I(θ) = E[a aᵀ/λ] (stationary, per unit time).

## 1. Exact reduction to one number

Stationarity gives E λ = μ/ε and E S = (E λ − μ)/α. With S = (λ − μ)/α and m := E[1/λ]:

- I_μμ = m
- I_μα = (1 − μm)/α
- I_αα = (μ/ε − 2μ + μ²m)/α²
- det I = (m·μ/ε − 1)/α²   (≥ 0 by Jensen: m ≥ 1/Eλ = ε/μ)

Everything about I(θ) is therefore determined by m = E[1/λ] (and I_μμ ≥ ε/μ always).

## 2. Exact formula for m (quadrature, not closed form)

The stationary Laplace transform φ(u) = E e^{−uλ} solves φ' = −uβμ φ/(uβ + e^{−uα} − 1), so

φ(u) = exp(−βμ ∫₀^u G(s) ds), G(s) = s/(εβs + h(αs)), h(x) = e^{−x} − 1 + x, **m = ∫₀^∞ φ(u) du.**

I could not find an elementary antiderivative of G, so the exact answer is this single integral, which `code/p10_fisher_k1_exact.py` evaluates to ~1e-10. The near-critical behaviour, which is what the problem needs, is in closed form (next section).

## 3. Near-critical asymptotics (ε → 0, μ and β fixed)

Let k = 2μ/β (the shape of the CIR/Gamma limit at criticality). **All three asymptotics below are now proved, not just matched: see Theorem A in [problem10_rigorous.md](problem10_rigorous.md).** Heuristic behind them: for u ~ ε, G ≈ 1/(εβ + α²s/2), so φ(u) ≈ (1 + α²u/(2βε))^{−k}. This Gamma-type power law has an integrable tail iff k > 1; otherwise the integral is dominated by u = O(1), where the true φ matters.

| regime | m = E[1/λ] | I_μμ | I_αα | det I | efficient info for μ (α unknown) = det/I_αα |
|---|---|---|---|---|---|
| k > 1 (2μ > β) | 2ε/(2μ−β) | ≈ 2ε/(2μ−β) | ≈ μ/(εβ²) | → 1/(β(2μ−β)) | ≈ εβ/(μ(2μ−β)) |
| k = 1 (β = 2μ) | (2/β) ε ln(1/ε) | ε ln(1/ε)·2/β | ≈ μ/(εβ²) | ~ ln(1/ε)·(const) | ≈ m |
| k < 1 (2μ < β) | C(μ,β) ε^k | ≈ C ε^k | ≈ μ/(εβ²) | ≈ Cμ ε^{k−1}/β² | ≈ C ε^k |

with I_μα → 1/β in the first row, and for k < 1:

C(μ,β) = (2/β)^k ∫₀^∞ u^{−k} F(u) du, F(u) = exp(−βμ ∫₀^u [g₀(s) − 2/(β²s)] ds), g₀(s) = s/(βs + e^{−βs} − 1).

In every regime I_αα ≈ μ/(εβ²) ~ 1/ε and the smallest eigenvalue of I is I_μμ-like (order ε for k > 1, ε^k for k < 1, ε ln(1/ε) at k = 1). For k > 1 this gives λ_min(I) ≈ εβ/(μ(2μ−β)); for μ = β = 1 that is ε, matching the K=1 table in the Problem 1 note.

## 4. Numerical check (exact m from the ODE, β = 1)

| μ | ε | m | m/ε | m/ε^k | local slope d ln m/d ln ε |
|---|---|---|---|---|---|
| 1.0 (k=2) | 1e-2 | 1.899e-2 | 1.899 | | 0.94 |
| | 1e-4 | 1.998e-4 | 1.998 | | 0.998 |
| 0.5 (k=1) | 1e-2 | 8.74e-2 | 8.74 | | 0.75 |
| | 1e-4 | 1.79e-3 | 17.9 | | 0.88 |
| 0.3 (k=0.6) | 1e-2 | 0.331 | | 5.25 | 0.525 |
| | 1e-4 | 0.0236 | | 5.93 | 0.589 |

- k = 2: m/ε → 2.000 (prediction 2/(2μ−β) = 2), and det·β²(1−ε)² → 1.000 (prediction 1).
- k = 1: m/(ε ln(1/ε)) at ε = 1e-4 is 1.94, against the predicted 2/β = 2; the slope creeps toward 1 as expected for a logarithm.
- k = 0.6: predicted C = 6.052 from the integral; m/ε^k = 5.93 at ε = 1e-4 and still rising toward it, and the slope approaches k = 0.6 (0.589).
- Every row also checks the λ_min formula: for μ = β = 1, λ_min(I)/ε = 0.998 at ε = 1e-4.

## 5. What this says for Problems 9 and 11 (corrected 2026-10-25)

Two different observation models give different answers, and my earlier version of this section conflated them.

**(a) Stationary observation** (the process is observed after a long burn-in, so λ_0 ~ π). Then the information over [0, T] is exactly T·I(θ*) (additivity under stationarity). With ε = c/T and k = 2μ/β:
- information about μ with α known: T·m → 2c/(2μ−β) (bounded) if k > 1; ~ 2c·ln(T/c)/β (log-divergent) if k = 1; ~ C c^k T^{1−k} (divergent) if k < 1.
- information about α (μ unknown) ≈ T/(2βε) ~ T²/c for k > 1.

**(b) Cold start** (S_0 = 0, λ_0 = μ, the natural model for an observed event stream from time 0). The nonstationary early phase matters: information about μ **always diverges, even for k > 1**. Rigorous: for any K, E λ_t ≤ μ(1 + β_max t), hence I_T(μ,μ) ≥ ln(1 + β_max T)/(μ β_max). Numerically (K=1, exact formula, T = 2/ε, β = 1): d M_T/d ln T → 2.03 for k = 2 (conjectured limit 2/(β(2μ−β)) = 2), has a slope that itself grows with ln T (roughly M_T ~ ln² T) for k = 1, and M_T ~ T^{1−k} for k < 1. Details and the exact transient formula are in [problem10_rigorous.md](problem10_rigorous.md) §4.
So the "bounded vs divergent" dichotomy for baseline information is really **logarithmic vs polynomial divergence** under cold start, with the same threshold k = 1, and it is **bounded vs divergent** only under stationary observation. Which of the two models contribution C4 intends needs to be decided.

For Problem 9, the information scale per unit time is ε (k > 1) or ε^k (k < 1), so a single ε-exponent will not match the lower bound across all μ.

## 6. Status and caveats (updated 2026-10-25)

- Proved: §1, §2, the three asymptotics in §3 (Theorem A), the exact cold-start formula, the rigorous log lower bound for general K, and the stationary-start identity. All in [problem10_rigorous.md](problem10_rigorous.md).
- Still conjectural (numerics only): the exact cold-start coefficients (2/(β(2μ−β)) for k > 1, ln² T for k = 1, T^{1−k} for k < 1).
- Not done: stationary asymptotics for K ≥ 2 (no closed-form Laplace transform; the Gamma-majorant argument has no direct analogue); a rigorous lower bound on the efficient information for K ≥ 2 (K = 2 Monte Carlo is consistent with K = 1); the regime where μ varies with ε.
- Re-derive §1 and the φ ODE yourself, as the problem list asks; the full derivation is written out step by step in problem10_rigorous.md §1–2, and both identities were also checked by independent simulation.
- Errata: my Monte Carlo thinning loop had a bug in the first version of some scripts (see problem10_rigorous.md §6). Nothing in this note's K=1 tables depends on a simulation: they come from the exact ODE.
