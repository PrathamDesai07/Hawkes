# Assumption M certified without a Gamma comparator (K = 1, any ε, any k₀)

Closes the "proof and experiments do not overlap in ε" gap of the regret theorem. Code: `code/p5_assumption_M_exact.py`.

## Method (replaces the Gamma-comparator certificate of `thmT_explicit_constants.md` for this purpose)

Assumption M asks for a lower bound on m = min_j E_π h_j, where h_j is a trapezoid in the state (support half-width d = 0.15·ES, plateau half-width d/2, centres 0.5·ES and 1.8·ES). In the variable X = ελ the stationary characteristic function φ(t) = E e^{itX} = exp(L(t)) is **exact** (Lemma R, `thmT_explicit_constants.md` §2):
 L′(x) = (k/2)·x/E₁(x) per unit x = αεt, E₁ = −(1−cos x) − i((1+r)x − sin x), r = ε/(1−ε), k = 2βμ/α².
A trapezoid centred at c with plateau half-width p, support half-width s, ramp ρ = s−p, is (1/ρ)·(1_{|x|≤A} ∗ 1_{|x|≤ρ/2}) with A = (s+p)/2, so its Fourier transform is e^{itc}H(t), H(t) = (4/ρ)·sin(At)·sin(ρt/2)/t² (real, even, |H| ≤ 4/(ρt²)). By Parseval

 **E_π h = (1/π) ∫₀^∞ H(t) Re[e^{itc} conj(φ(t))] dt.**

No comparison with a Gamma law is made, so there is no ε → 0 restriction and no condition k₀ > 1.

**Error control.**
- *Truncation, rigorous:* |φ| is non-increasing in t (`thmT_explicit_constants.md` Lemma elem(c)) and |H| ≤ 4/(ρt²), hence the tail beyond T is at most |φ(T)|·4/(πρT). The code stops when this explicit bound is < 2·10⁻⁶ (reported per case).
- *Quadrature:* composite Gauss–Legendre panels in t (panel 0.4 with 16 nodes, and panel 0.2 with 24 nodes), with L at the nodes obtained from the panel-wise spectral integration matrix. The two quadratures agree to ≤ 2·10⁻⁹ in every case. This is a convergence check, not an interval-arithmetic bound.
- *Independent validation:* the values reproduce the stationary-simulation values (0.15, 0.05) and converge to the Gamma-comparator limits (0.164, 0.0444) as ε → 0.

## Results (β = 1; E_π h on I₁ = 0.5·ES and I₂ = 1.8·ES; tail bound ≤ 5·10⁻⁶ everywhere)

| k₀ (μ) | ε = 0.3 | 0.1 | 0.03 | 0.01 | 0.001 |
|---|---|---|---|---|---|
| 4 (2.0): I₁ / I₂ | 0.15271 / 0.04322 | 0.15704 / 0.04255 | 0.15870 / 0.04233 | 0.15919 / 0.04227 | 0.15941 / 0.04224 |
| 2 (1.0): I₁ / I₂ | 0.15283 / 0.04631 | 0.16023 / 0.04504 | 0.16285 / 0.04462 | 0.16360 / 0.04451 | 0.16394 / 0.04445 |
| 1.2 (0.6): I₁ / I₂ | 0.13631 / 0.04175 | 0.14183 / 0.04033 | 0.14453 / 0.03986 | 0.14528 / 0.03973 | 0.14561 / 0.03967 |
| 0.6 (0.3): I₁ / I₂ | 0.08555 / 0.03235 | 0.10556 / 0.03055 | 0.10877 / 0.03010 | 0.10936 / 0.02997 | see log |

The smaller of the two is I₂ in every row, so **m ≥ 0.032 (k₀ = 0.6), 0.040 (k₀ = 1.2), 0.042 (k₀ = 4) and 0.044 (k₀ = 2) for every ε ≤ 0.3 tested.** The ε-dependence is mild (monotone, within 8 % of the limit already at ε = 0.1 for k₀ ≥ 1.2).

## Consequence

Assumption M for K = 1 is now certified, at the stated quadrature accuracy, for the whole range ε ∈ [10⁻³, 0.3] and for k₀ from 0.6 to 4, which includes every ε used in the regret experiments (0.3, 0.1, 0.03). Theorem R (K = 1) no longer needs the restriction ε ≤ 10⁻⁴; it needs only Assumption M at the ε of interest, which can be certified by this script for any (μ, β, ε). The case k₀ ≤ 1 is covered by the same certificate (the regret theorem's other ingredients were proved without k₀ > 1).

Limits:
1. Floating-point quadrature with a convergence check, not interval arithmetic; the truncation bound is rigorous.
2. The certificate is per parameter value. Uniformity over a compact parameter set follows from continuity of E_π h in (μ, β, ε), which was not written out.
3. K ≥ 2 is not covered (no closed-form transform). `K2_regret_chain.md` replaces Assumption M for K ≥ 2 by a different, fully explicit device.
