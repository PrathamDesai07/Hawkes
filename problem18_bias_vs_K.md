# Problem 18: bias of n(θ†) versus K for a power-law truth on a geometric β-grid (contribution C6)

Date: 25 Oct 2026. Tags: **[P]** proved on paper, **[N]** numerical, **[C]** conjecture with numerical support, **[O]** open. Builds on `problem16_kl_pseudotruth.md` (Prop. 16.4 rate matching, Cor. 2(iv), Theorem 16, section 6). Not refereed.

## 0. What the tracker asks

"For a power-law truth on a geometric β-grid, give an explicit rate for the bias of n(θ†) as K grows. This is a rational-approximation problem for a completely monotone kernel with a branch point. Done when: an explicit rate (or a numerically supported conjecture) for the bias as a function of K and grid ratio."

Result in one line: **b_K := n* − n(θ†_K) ≈ κ(η) · n* (τβ_1)^η / Γ(1+η) = κ(η) · (mass of the truth below the slowest pole), with β_1 = β_top r^{−(K−1)}, so b_K ∝ r^{−η(K−1)} decays geometrically in K at rate η ln r; κ(η) ∈ (0,1) is a universal number (≈ 0.43 at η = ½) that does not depend on n*, μ*, or the distance to criticality.** It is a conjecture [C] supported by a semi-analytic computation that agrees with the simulations of Problem 16; the reductions leading to it are proved.

## 1. Setting

Truth: stationary Hawkes with μ*, kernel φ*(t) = ∫ e^{−βt} ν(dβ), n* = ∫φ* = ∫ν(dβ)/β < 1, Λ* = μ*/(1 − n*). Power law: φ*(t) = (n*η/τ)(1 + t/τ)^{−(1+η)}, i.e. (Laplace/gamma representation) ν(dβ) = c β^η e^{−τβ} dβ with c = n* τ^η/Γ(η) (the completely monotone kernel with a branch point at 0: φ̂*(s) = ∫ν/(β+s) ≈ n* − C s^η as s → 0).
Model: K exponentials, β_k = β_1 r^{k−1}, θ† the KL pseudo-truth of Problem 16 (also the surrogate minimiser θ_Q). Bias b_K = n* − n(θ†).

## 2. Proved reductions

**Proposition 18.1 (bias of n equals bias of the baseline) [P].** By Prop. 16.4 / Cor. 2 (valid for KL and for the quadratic surrogate), if μ† > μ_min and the stability constraint is slack, μ† + n†Λ* = Λ*. Since μ* + n*Λ* = Λ*:
 **n* − n† = (μ† − μ*)/Λ*.**
So the bias in the branching ratio is exactly the excess baseline divided by the mean rate; the fitted model absorbs the missing slow memory into extra immigration. (Verified: μ† = 1 + 5·b in the numerical tables of Problem 16: e.g. K = 2, grid {0.3, 1.5}: b = 0.1788, μ† = 1.894.)

**Proposition 18.2 (frequency-domain form of the surrogate pseudo-truth) [P; KNOWN: it is Rong–Nair (arXiv 2602.20795) Lemma 4(b) (Bartlett spectrum Λ/|1 − φ̂₀|²), Lemma 13(b) (Parseval form of R*, R*^(1,0)) and Theorem 15 (pseudo-true α* = R*⁻¹R*^(1,0), c* = Λ(1 − Γ*)), checked in their PDF on 25 Oct 2026; it is not a contribution here].** Let G(s) = φ̂*(s) = ∫ν(dβ)/(β + s) and R_α(s) = Σ α_k/(β_k + s). For the stationary Hawkes truth the Bartlett spectrum of N is S(ω) = Λ*/|1 − G(iω)|² (Hawkes 1971). The surrogate objective E(λ*_t − λ_θ(t))² equals
 (μ* − μ + Λ*(n* − n))² + (1/π) ∫₀^∞ |G(iω) − R_α(iω)|² · Λ*/|1 − G(iω)|² dω,
because λ*_t − λ_θ(t) = (μ* − μ) + ∫(φ* − φ_α)(t − s) dN_s, whose mean is (μ* − μ) + Λ*(n* − n) and whose fluctuation is a linear filter of the centred N. The first term is zero at the optimum (choose μ), so **α_Q is the weighted L² approximant of G on the imaginary axis by the pole sums R_α, with weight w(ω) = Λ*/|1 − G(iω)|²**, and b_K = G(0) − R_α(0) is the DC value of the approximation error. The weight is huge at low frequency (Λ*/ε*² at ω = 0) and tends to Λ* at high frequency: the near-critical regime pushes accuracy to low frequencies, exactly where the branch point is.
*Consequence.* Computing θ_Q needs only deterministic integrals (no simulation), for any K, in double precision with a log-frequency grid. This is the tool used below.

**Lemma 18.3 (relation to the KL pseudo-truth) [N].** KL and Q pseudo-truths differ little: paired on the same path, n_KL − n_Q = −0.017 (K=2), +0.001 (K=3), 0.000 (K=4), −0.001 (K=5), −0.001 (K=6), −0.002 (K=7), statistical error ≈ 0.001 (grid β_k = 3·3^{−(K−k)}, n* = 0.8). So b_K^{KL} = b_K^{Q} within about 0.002 for K ≥ 3 (about 10–40 % of b when b < 0.02).

## 3. The rate (numerically supported)

**Truncation mass.** The best a pole set can do is limited by the kernel mass at rates below β_1: n_low(β_1) = ∫₀^{β_1} ν(dβ)/β = cβ_1^η/η = n* (τβ_1)^η/Γ(1+η).

**Conjecture 18.4 [C].** For a geometric grid with ratio r in a moderate range and K large enough that the truncation term dominates the discretization floor,
 **b_K = κ(η, r) · n* (τβ_1)^η / Γ(1+η)·(1 + o(1)),  β_1 = β_top r^{−(K−1)}.**
Equivalently b_K ∝ r^{−η(K−1)}: geometric decay in K with exponent η ln r per added exponential; κ(η, r) ∈ (0, 1) is the fraction of the truncated mass that cannot be recovered. Observed properties:
1. *Exponent η.* For η = ½, r = 2: b = 1.306·10⁻², 4.158·10⁻³, 1.352·10⁻³, 4.261·10⁻⁴ at β_1 = 10⁻³, 10⁻⁴, 10⁻⁵, 10⁻⁶: ratios 0.318, 0.325, 0.315 per decade, against 10^{−η} = 0.316.
2. *Geometric in K.* r = 2, η = ½: b_K = 0.0592, 0.0307, 0.0158, 0.0079, 0.0020 at K = 8, 10, 12, 14, 18; per-step ratio ≈ 0.71 = 2^{−1/2}. For r = 3: ratio per step 0.58 = 3^{−1/2}; r = 5: 0.43 vs 0.45.
3. *Fixed range, K → ∞.* For the grid fixed on [10⁻³, 10], b_K saturates at 0.0129 = κ n* (τ·10⁻³)^η/Γ(1+η) with κ ≈ 0.45: refining the grid cannot remove the truncation at β_1 (finer r at fixed β_1 gains nothing).
4. *κ is independent of n*, i.e. of the distance to criticality.* η = ½, r = 1.5: κ = 0.462, 0.461, 0.461, 0.460, 0.458 for n* = 0.3, 0.5, 0.8, 0.9, 0.95. The weight w only enters through its value on the frequency band ω ≈ β_1, which is flat when β_1 is far below the near-critical relaxation rate. This is why a self-similar local problem governs the constant.
5. *κ(η).* Signed least squares (no positivity needed for r ≥ 3, K ≤ 12, see 7): η = 0.15: 0.69–0.73; 0.25: 0.61–0.66 (r = 2…4); 0.35: 0.53–0.60; 0.5: 0.40–0.50 (r = 2…4, the spread is the discretization floor); dense-grid positive fit (r = 1.3): 0.46 at η = ½. So κ decreases from about 0.7 to 0.45 as η goes from 0.15 to 0.5 and depends weakly on r. A one-pole heuristic (the lowest pole absorbs the truncated mass optimally in flat-weight L²) gives the closed form κ_1(η) = 1 − 2η∫₀¹ x^η/(1+x) dx = 0.82, 0.73, 0.66, 0.57 for η = 0.15, 0.25, 0.35, 0.5, an upper bound on all the numerical κ (the extra poles shrink the unrecovered fraction by roughly a fifth).
6. *Discretization floor and sign change.* For coarse grids the bias goes negative (n† > n*): the discrete poles overshoot. At fixed top β = 3 and K = 6 the bias is 0.103, 0.045, 0.023, 0.013, 0.007, 0.0016, −0.0018, −0.0057 for r = 2, 3, 4, 5, 6, 8, 10, 13. The zero crossing ratio decreases with K (r ≈ 17 at K = 4, ≈ 8 at K = 6, ≈ 6 at K = 8); beyond it |b| grows again. Heuristically the floor is the trapezoid error of the Stieltjes integral on the log scale, exp(−2π d / ln r) with half-width d = π/2 of the analytic strip at imaginary frequency, i.e. ≈ exp(−π²/ln r); it is only seen when β_1^η-term is below it.

**Design consequence [N].** For a fixed K the best grid ratio is not the smallest: with the top pole fixed at 3, the minimum of |b| is at r ≈ 40 (K=2), 30 (K=3), 17 (K=4), 13 (K=5), 8 (K=6), 6 (K=8). To get a bias ≤ δ with ratio r the number of exponentials needed is K ≈ 1 + ln(κ n*(τβ_top)^η/(Γ(1+η)δ)) / (η ln r). For η = ½, n* = 0.8, β_top = 3, r = 3, δ = 0.01: K ≈ 9.

## 4. Why the exponent is η (argument, not a proof)

G(s) = n* − C s^η + regular terms near s = 0 (branch point), while every R_α is analytic at 0 with poles at −β_k ≤ −β_1. In the variable x = β/β_1 the part of G coming from rates below β_1 is f(σ) = ∫₀¹ x^η/(x + σ) dx times cβ_1^η, a function of σ = s/β_1 with no dependence on any other scale; the weight is flat there when β_1 ≪ the criticality scale (point 4 above). The projection of f onto span{1/(σ + x_j)}, x_j = r^j ≥ 1, gives a remainder whose DC value is κ(η, r)·f(0) = κ/η. Scaling s → β_1 s then yields b_K = κ(η, r) cβ_1^η/η. What is not proved is that κ(η, r) is positive and r-independent in the limit, and its value: this is a Wiener–Hopf-type problem for a geometric ladder of poles (Newman/Stahl-type rational approximation of functions with a branch point), and I did not solve it.

## 5. What is and is not established

Established [P]: Props. 18.1 and 18.2 (bias of n = excess baseline over Λ*; the pseudo-truth is a weighted L² pole approximation with the Hawkes spectrum as weight); the semi-analytic code reproduces the simulations of Problem 16.
Supported numerically [N,C]: Conjecture 18.4 and its five observed properties; KL ≈ Q.
Not established: a proof of the exponent η and of κ(η, r) > 0 (a rigorous lower bound would follow from a quantitative Müntz-type argument for non-negative exponent sums); the sign change/floor law ≈ exp(−π²/ln r); the dependence on the top of the grid; general ν beyond the Gamma mixing measure (the argument only uses ν ~ cβ^η at small β, tested only for the gamma form); and everything for the unconstrained-sign vs α ≥ 0 distinction beyond K ≤ 12: for r ≤ 2 and K ≥ 5 the positivity constraint binds (1 zero α at r = 2) and the signed fit slightly under-estimates the bias (r = 2, K = 8: 0.0592 nonneg vs 0.0543 signed).

## 6. Numerical record

Code: `code/p18/fd18.py` (frequency-domain pseudo-truth; log-frequency trapezoid, 760 frequencies from e^{−24} to e^{14}; NNLS/lstsq on the √w-weighted real and imaginary parts), `e2.py`…`e8.py`. Truth for the scaling study: the gamma mixing measure ν ∝ β^η e^{−β} (τ = 1, n* = 0.8, μ* = 1, Λ* = 5) on a fine log grid (step 0.02 from e^{−30} to e^{8}, truncation mass ≈ 3·10⁻⁷).
- **Validation vs simulation** (24-exponential truncated truth, as in Problem 16): n_Q = 0.459, 0.6212, 0.7005, 0.7361 for β = {1}, {0.3, 1.5}, {0.1, 0.5, 2.5}, {0.05, 0.2, 0.8, 3.2} against simulated KL n† = 0.463, 0.621, 0.697, 0.733.
- **Fixed top β = 3, bias n* − n†** (K = 1, 2, 3, 4, 6, 8, 10, 14, 18): r = 2: 0.556, 0.436, 0.307, 0.217, 0.113, 0.059, 0.031, 0.0079, 0.0020; r = 3: …0.0445, 0.0155, 0.0053, 0.0005, 0.0001; r = 5: …0.0128, 0.0021, 0.0004, 0.0000; r = 10: …−0.0018 for K ≥ 6 (floor).
- **KL vs Q at r = 3** (above): bias_KL = 0.338, 0.181, 0.100, 0.052, 0.023, 0.0085 at K = 2…7, Q: 0.322, 0.182, 0.100, 0.051, 0.022, 0.006.
- Scripts reproduce every number above (seeded).
