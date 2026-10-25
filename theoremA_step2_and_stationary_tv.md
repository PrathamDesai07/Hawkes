# (1) Repaired Step 2 of Theorem A, and (2) the initial-state gap in the impossibility result (K = 1)

Written 25 Oct 2026. Part 1 repairs two inadequate justifications found in `verification_thmA_E_F.md`. Part 2 attempts to close the gap in Corollary E′ / Theorem B of `problem11_dichotomy.md`. **These proofs are mine, written in one pass; re-derive them before relying on them.** Code: `code/verify_step2_repair.py`, `code/tv_stationary_laws.py`.

---

## Part 1. Repaired Step 2 of Theorem A

Notation (Problem 10): α = (1−ε)β, h(x) = e^{−x} − 1 + x, x = αs, r := εβ/α = ε/(1−ε), and
D_ε(s) = (α²s²/2 − h(αs)) / ((εβs + h(αs))(εβ + α²s/2)).
In terms of x: εβs = r x, εβ + α²s/2 = α(r + x/2). Elementary facts used (all for x ≥ 0): h is convex, h(0) = 0, so h(x)/x is nondecreasing; (x²/2)(1 − x/3) ≤ h(x) ≤ x²/2; x²/2 − h(x) ≤ x³/6; h(x) ≤ x − c with c := 1 − e⁻¹ for x ≥ 1.

**Lemma 2.** For every ε ∈ (0,1):
- (a) D_ε(s) ≥ 0 for all s > 0.
- (b) D_ε(s) ≤ 2.08/α for all s > 0.
- (c) If ε ≤ 1/4, then D_ε(s) ≥ 0.166/α ≥ 0.166/β for all s ≥ 1/α.

*Proof.* (a) h ≤ x²/2 makes the numerator ≥ 0; both denominator factors are positive.

(b) The denominator is ≥ h(x)·αx/2 (drop the non-negative terms rx and r). *Case x ≤ 1.5:* h(x) ≥ (x²/2)(1 − x/3) ≥ x²/4, and the numerator is ≤ x³/6, so D_ε ≤ (x³/6)/((x²/4)(αx/2)) = 4/(3α). *Case x ≥ 1.5:* by convexity h(x) ≥ x·h(1.5)/1.5, so the denominator is ≥ αh(1.5)x²/3, while the numerator is ≤ x²/2; hence D_ε ≤ 3/(2h(1.5)α) = 2.074/α. (The original text bounded the denominator by h(1.5)αx/2, which only gives D_ε ≤ x/(αh(1.5)).)

(c) *Numerator.* (x²/2 − h(x))/x² = ½ − h(x)/x² is nondecreasing: d/dx[h/x²] = q(x)/x³ with q(x) = 2 − x − (x+2)e^{−x}, q(0) = q′(0) = 0 and q″(x) = −xe^{−x} < 0, so q′ < 0 and q < 0 for x > 0. Hence for x ≥ 1 the numerator is ≥ (½ − e⁻¹)x² = 0.1321x². *Denominator.* For ε ≤ 1/4 we have r ≤ 1/3, and both factors of the denominator are increasing in r, so for x ≥ 1
Den ≤ α((4/3)x − c)(1/3 + x/2) = α[(2/3)x² + (4/9 − c/2)x − c/3] ≤ αx²[2/3 + (4/9 − c/2)/x] ≤ 0.7951·αx²
(using h(x) ≤ x − c so that rx + h(x) ≤ (4/3)x − c, dropping −c/3 ≤ 0, and 4/9 − c/2 = 0.1283 with x ≥ 1). (The original text replaced the exact constant 1/3 in the second factor by x/3, which overshoots: that gives 10/9·αx² > αx².) Therefore D_ε ≥ 0.1321x²/(0.7951αx²) = 0.1662/α. ∎

**Consequences used in Steps 3 of Theorem A (replacing the corresponding lines).**
- F_ε(u) = exp(−βμ∫₀^u D_ε) satisfies, for **all** u > 0 and all ε, F_ε(u) ≥ exp(−2.08βμu/α); this one inequality replaces the two-case lower bound and is what the dominated-convergence arguments (k₀ > 1, k₀ = 1) use.
- For ε ≤ 1/4: F_ε(u) ≤ exp(−0.166μ(u − 1/α)), which supplies the exponential decay needed for k₀ < 1 (the original constant 0.13 is also correct, with room).
- D_ε → D₀ pointwise and is uniformly bounded (b), so F_ε → F₀ locally uniformly.

**Numerical check (`verify_step2_repair.py`, β = 1, fine grids up to x = 10⁷):** q(x) < 0 on the grid (max 0); ratio Den/(αx²) ≤ the closed-form bound with max(ratio − bound) = 0 and max ratio 0.667; sup D_εα = 1.000 (≤ 2.08); min D_εα over s ≥ 1/α and ε ≤ 1/4 = 0.226 (≥ 0.166). With the repairs, Step 2, and therefore Theorem A, is fully justified.

---

## Part 2. The initial-state gap (K = 1, k₀ = 2μ/β > 1)

**What was missing.** Corollary E′ uses the KL lemma KL ≤ ρ·uᵀI_Tu, which compares the two path laws with the *same* initial state. For a stationary record with unobserved past the initial state is drawn from π_θ under one hypothesis and π_θ′ under the other. The note did not bound this discrepancy.

**Setting.** θ = (μ, α), α = (1−ε)β; ridge perturbation θ′ = (μ+δ, α(1 − εδ/μ)), so u = θ′ − θ = (δ, −αεδ/μ). Q_θ, Q_θ′ are the laws of the event record on [0,T] of the stationary processes (S₀ unobserved).

### Lemma R (reduction to total variation of the stationary laws)
TV(Q_θ, Q_θ′) ≤ √(ρ·uᵀI_Tu / 2) + TV(π_θ, π_θ′), where I_T = T·I(θ) and ρ = max(μ/μ′, α/α′).
*Proof.* Let P̃ be the law of the record when S₀ ~ π_θ and the dynamics are those of θ′. (i) Q_θ and P̃ are mixtures over S₀ ~ π_θ of P^θ_s and P^θ′_s; therefore KL(Q_θ‖P̃) ≤ E_{π_θ}KL(P^θ_s‖P^θ′_s) ≤ ρ·uᵀI_Tu by the KL lemma (a stationary-start expectation, so no initial-law term), and Pinsker gives TV(Q_θ, P̃) ≤ √(ρuᵀI_Tu/2). (ii) P̃ and Q_θ′ are mixtures of the same kernels P^θ′_s over S₀ ~ π_θ versus π_θ′, so TV(P̃, Q_θ′) ≤ TV(π_θ, π_θ′) by data processing. Triangle inequality. ∎
With ε = c/T and Theorem A (E[1/λ] ∼ 2ε/(2μ−β) for k₀ > 1), ρ·uᵀI_Tu = ρδ²T[(1+ε)²m − (ε/μ)(1+2ε)] stays bounded: the first term of Lemma R is O(δ), independent of T. **Only TV(π_θ, π_θ′) = O(δ), not KL, is needed.**

### Lemma F (total variation from characteristic functions)
For probability densities f, g with characteristic functions ϕ_f, ϕ_g, Δ := ϕ_f − ϕ_g, if Δ ∈ H¹(R):
TV(f, g) ≤ (1/(2√2)) · (∫|Δ|²dt + ∫|Δ′|²dt)^{1/2}.
*Proof.* With h = f − g, Cauchy–Schwarz gives ‖h‖₁ ≤ (∫(1+x²)⁻¹dx)^{1/2}(∫(1+x²)h²dx)^{1/2} = √π(‖h‖₂² + ‖xh‖₂²)^{1/2}; by Plancherel ‖h‖₂² = ‖Δ‖₂²/(2π) and ‖xh‖₂² = ‖Δ′‖₂²/(2π); TV = ½‖h‖₁. ∎

### Theorem T (K = 1, k₀ > 1)
Fix β and μ with k₀ > 1. There are ε₀, δ₀, C depending only on (β, μ) such that for ε ≤ ε₀, |δ| ≤ δ₀:
TV(π_θ, π_θ′) ≤ C|δ| + C ε^η, η = min{1, k₀ − ½} > 0.
**Corollary (impossibility for the stationary record with unobserved past).** For K = 1 and k₀ > 1, with ε = c/T, there is δ₁ > 0 such that for δ ≤ δ₁ and T large, TV(Q_θ, Q_θ′) ≤ 1/2, hence by Le Cam max_{θ,θ′} P(|μ̂ − μ| ≥ δ/2) ≥ 1/4 for every estimator: μ is not consistently estimable (α known or not), also without access to S₀.

*Proof of Theorem T.*

**Step 1 (characteristic function).** Since c(u) = βu − (1 − e^{−αu}) and |1 − e^{−αu}| ≤ α|u| < β|u| for Re u ≥ 0, c has no zero in {Re u ≥ 0}∖{0}; φ(u) = E e^{−uλ} = exp(−βμ∫₀^u s/c(s) ds) on u > 0 extends analytically to Re u > 0 and by continuity to the imaginary axis. For X = ελ, ϕ(t) := E e^{itX} = exp(L(t)), L(t) = βμ ε²∫₀^t η dη / c(−iεη).

**Step 2 (Gamma factorisation on the imaginary axis).** Using s/c = 1/(εβ + α²s/2) + D_ε(s): ϕ(t) = G(t)·F̃(t), G(t) = (1 − iat)^{−k}, a = α²/(2β), k = 2βμ/α² (the characteristic function of Gamma(k, scale a)), and F̃(t) = exp(iβμε∫₀^t D_ε(−iεη)dη).

**Step 3 (imaginary-axis bound for D_ε).** For all τ ∈ R and ε ∈ (0,1): |D_ε(−iτ)| ≤ 31.5/α. With x = α|τ|: numerator N = (1 − cos x − x²/2) + i(x − sin x); E₁ = εβs + h(αs) = −(1 − cos x) − i((1+r)x − sin x); E₂ = α(r − ix/2). *x ≤ 1:* |N| ≤ x³/6 + x⁴/24 ≤ 0.209x³, |E₁| ≥ 1 − cos x ≥ 0.458x², |E₂| ≥ αx/2, so |D| ≤ 0.91/α. *x ≥ 1:* |N| ≤ x²/2 + 2x ≤ 2.5x², |E₁| ≥ x − sin x ≥ 0.158x, |E₂| ≥ αx/2, so |D| ≤ 31.5/α. (Numerically the true sup of |D_ε(−iτ)|α over τ and ε ∈ {0.9, …, 10⁻⁵} is 1.143, at x ≈ 7.4; a first numerical attempt gave 52.7, which was catastrophic cancellation in the complex evaluation of h at |τ| < 10⁻⁴, fixed by using the series.) Consequently |F̃(t)| ≤ e^{M₁ε|t|}, |F̃′(t)| ≤ M₁ε|F̃(t)|, with M₁ = 31.5βμ/α.

**Step 4 (monotone decay).** Re L′(t) = −βμ ε² t (1 − cos x)/|c(−iεt)|² ≤ 0 with x = αεt, |c(−iεt)|² = (1 − cos x)² + ((1+r)x − sin x)². So |ϕ| is non-increasing in |t|.

**Step 5 (tail beyond t₁ = 1/(M₁ε)).** For t ≥ t₁: |ϕ(t)| ≤ |ϕ(t₁)| ≤ e|G(t₁)| ≤ e(at₁)^{−k} = O(ε^k). Beyond x ≥ X₀ (t ≥ t_* := X₀/(αε), X₀ large) use |c|² ≤ (1+r)²x²(1 + 3/x) to get Re L(t) − Re L(t_*) ≤ −q[ln(t/t_*) − 2/x_*] with q = k/(2(1+r)²(1 + 3/X₀)) (the integral ∫(1 − cos y)/y dy ≥ ln(x/x_*) − 2/x_*). Since k₀ > 1, choosing ε₀ small and X₀ large gives q > ½. Then ∫_{t₁}^{t_*}|ϕ|² ≤ e²(at₁)^{−2k}t_* = O(ε^{2k−1}) and ∫_{t_*}^∞|ϕ|² ≤ O(ε^{2k−1})/(2q − 1). For the derivative, ϕ′ = ϕL′ with |L′(t)| = βμε²|t|/|c| ≤ 1.1k/t for x ≤ 1 (|c| ≥ 0.458x²) and ≤ 8.4με for x ≥ 1 (|c| ≥ |Im c| ≥ 0.158x); so ∫_{|t|>t₁}|ϕ′|² = O(ε^{2}·ε^{2k−1}) + O(ε^{2k+1}). The Gamma terms G, G′ on |t| > t₁ are O(t₁^{1−2k}) = O(ε^{2k−1}).

**Step 6 (head |t| ≤ t₁).** There M₁ε|t| ≤ 1, so |F̃ − 1| ≤ (e−1)M₁ε|t| and |Δ| = |G||F̃ − 1| ≤ (e−1)M₁ε|t||G|; |Δ′| ≤ |G′||F̃ − 1| + |G||F̃′| ≤ M₁ε((e−1)k + e)|G|, using |G′(t)| = ka|1 − iat|^{−k−1} = (ka/√(1 + a²t²))|G| ≤ (k/|t|)|G|. Hence ∫_{|t|≤t₁}|Δ|² ≤ Cε²(1 + t₁^{3−2k}) = O(ε^{min(2,2k−1)}) and ∫|Δ′|² = O(ε²).

**Step 7 (conclusion).** ‖Δ‖²_{H¹} = O(ε^{min(2, 2k−1)}); Lemma F gives TV(π_θ, Γ_θ) ≤ Cε^{min(1,k−1/2)} for the Gamma comparator Γ_θ = Gamma(k, scale a/ε) in λ-space. The same holds for θ′ (its own scale ε′ ≍ ε, k′ near k, uniformly on the compact family). Finally Γ_θ and Γ_θ′ have shapes k = 2βμ/α², k′ = 2βμ′/α′² and scales a/ε, a′/ε′ whose logarithms differ by O(δ) uniformly in ε; their KL is C² in these parameters with minimum 0 at equality, hence ≤ Cδ², and Pinsker gives TV(Γ_θ, Γ_θ′) ≤ Cδ. The triangle inequality completes the proof. ∎

### Numerical check (exact characteristic functions inverted by FFT; β = 1; `tv_stationary_laws.py`)

μ = 1 (k₀ = 2), ridge pair θ′ = (μ+δ, α(1 − εδ/μ)):

| ε | δ | TV(π_θ, π_θ′) | TV(π_θ, Γ_θ) | TV(Γ_θ, Γ_θ′) | ‖ϕ − G‖_{H¹} |
|---|---|---|---|---|---|
| 0.1 | 0.1 / 0.3 | 0.0334 / 0.0875 | 0.0201 | 0.0301 / 0.0836 | 0.111 |
| 0.03 | 0.1 / 0.3 | 0.0267 / 0.0731 | 0.0060 | 0.0263 / 0.0723 | 0.0368 |
| 0.01 | 0.1 / 0.3 | 0.0256 / 0.0702 | 0.0020 | 0.0253 / 0.0693 | 0.0128 |
| 0.003 | 0.1 / 0.3 | 0.0257 / 0.0704 | 0.0006 | 0.0250 / 0.0683 | 0.0039 |

TV(π_θ, π_θ′) is bounded in ε (converging to the Gamma–Gamma value) and roughly linear in δ; TV(π_θ, Γ_θ) and the H¹ norm both decay ∝ ε, as predicted for k = 2.4 > 3/2 (rate ε^{min(1, k−1/2)} = ε). Lemma F's inequality TV ≤ H¹/(2√2) holds in every row (e.g. 0.0393 ≥ 0.0201 at ε = 0.1). Mass check Σf·dx = 1.0000 in all cases. Near the threshold, μ = 0.6 (k₀ = 1.2): TV(π_θ, π_θ′) = 0.103 / 0.176, 0.055 / 0.124, 0.045 / 0.113, 0.043 / 0.111 at ε = 0.1, 0.03, 0.01, 0.003 (δ = 0.1 / 0.3); H¹ norm 0.534, 0.275, 0.139, 0.063 (decay ≈ ε^{0.6–0.7}, matching the worst-case rate ε^{k−1/2} = ε^{0.7}); TV(π_θ, Γ_θ) decays faster, 0.075 → 0.0018.

### Status and limits
- **Closed (K = 1, k₀ > 1):** the gap in Corollary E′ / Theorem B for the stationary record with unobserved past, given Lemma R, Lemma F and Theorem T. Constants (ε₀, δ₀, C) are not explicit; the argument is asymptotic in ε → 0 at fixed (β, μ).
- **Not closed:** K ≥ 2 (no closed-form Laplace transform; Theorem E only gives E[1/λ] = O(ε), and no analogue of the Gamma comparator is available), the threshold k₀ = 1 and k₀ < 1 (Step 5 needs q > ½, i.e. k₀ > 1; for k₀ < 1 there is no impossibility to prove, information diverges), and explicit constants.
- **Checks:** Steps 3–7 were hand-derived and each analytic bound was tested numerically (Lemma 2, the imaginary-axis bound, the rates, the TV values), but the chain of constants in Steps 5–6 was not formalised line by line; treat it as a proof outline with all key estimates verified numerically.

---

## Update (25 Oct 2026): K ≥ 2

The direct extension of Theorem T to K ≥ 2 was not obtained (no closed-form transform; no scalar comparator). The impossibility itself, however, is proved for every K by a route that never compares π_θ with π_θ′: `general_K_impossibility.md` (Theorem V: TV(Q_θ, Q_θ′) ≤ Cδ(ln(1/δ) + 1)^{1/2} for k_K > 1). Theorem T remains the only place where TV(π_θ, π_θ′) = O(δ) is established (K = 1).

---
## Update (25 Oct 2026): explicit constants
Steps 3–7 of Theorem T are now written out with explicit constants and checked against the exact characteristic function in `thmT_explicit_constants.md` (certified TV(π,Γ) ≤ Φ/(2√2), Φ computed from hand-proved majorants). The earlier remark "the chain of constants was not formalised" no longer applies; what remains informal is listed in its §9.
