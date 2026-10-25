# Theorem T (K = 1): the constant chain written out

Replaces the "outline" status of Steps 3–7 in `theoremA_step2_and_stationary_tv.md`. Code: `code/thmT_constants.py`, `code/thmT_final_checks.py`.
Convention: β = 1 (time scaling; the only parameters are k₀ = 2μ/β and ε). α = 1−ε, r = ε/(1−ε) (so 1+r = 1/(1−ε)), k = 2βμ/α² = k₀/(1−ε)², a = α²/(2β), X := ελ. Then k/(1+r)² = k₀ **exactly**.
Standing: k₀ > 1, ε ≤ ½. Dimensionless variable x = αεt.
Tags: [H] proved by hand below, [C] computed by floating-point quadrature of an explicit majorant (not interval arithmetic), [V] verified numerically against the exact characteristic function.

## Statement (Theorem T′)

Let Φ(k₀, ε) be the explicit number defined in §6 (a finite sum of closed-form terms and two one-dimensional integrals of explicit functions). Then for every k₀ > 1 and ε ≤ ½
 TV(π_θ, Γ_θ) ≤ Φ(k₀, ε)/(2√2),
where Γ_θ = Gamma(k, scale a/ε) in λ-space. Φ(k₀, ε) = O(ε^{min(1,k−½)}) as ε → 0 (for k = 3/2 an extra log). Values (k₀ = 2) of the TV bound Φ/(2√2): 1.96·10⁻² (ε = 10⁻⁴), 3.9·10⁻³ (3·10⁻⁵), 1.0·10⁻³ (10⁻⁵), 8.6·10⁻⁵ (10⁻⁶), 8.4·10⁻⁷ (10⁻⁸). Other k₀ (ε = 10⁻⁶ / 10⁻⁸ / 10⁻¹⁰): k₀ = 1.2: 0.10 / 4.0·10⁻³ / 1.6·10⁻⁴; k₀ = 1.5: 1.4·10⁻² / 1.4·10⁻³ / 1.4·10⁻⁴; k₀ = 3: 1.3·10⁻⁴ / 1.3·10⁻⁶ / 1.3·10⁻⁸.
The true values (exact characteristic function, [V]) are 10²–10³ times smaller (k₀ = 2: H¹ norm 1.3·10⁻² at ε = 10⁻², 1.3·10⁻⁶ at 10⁻⁶), i.e. the bound is valid and loose but has the right ε-scaling (linear in ε for k₀ = 2).

## 1. Lemma F (TV from characteristic functions) [H]

For probability densities p, q with characteristic functions ϕ_p, ϕ_q (ϕ(t) = ∫e^{itx}p), Δ = ϕ_p − ϕ_q, ρ = p − q:
∫|ρ| = ∫(1+x²)^{−1/2}(1+x²)^{1/2}|ρ| ≤ (∫dx/(1+x²))^{1/2}(‖ρ‖₂² + ‖xρ‖₂²)^{1/2} = √π · ((‖Δ‖₂² + ‖Δ′‖₂²)/(2π))^{1/2} (Plancherel, ‖ϕ‖₂² = 2π‖p‖₂², xρ ↔ −iΔ′).
So TV = ½∫|ρ| ≤ ‖Δ‖_{H¹}/(2√2), ‖Δ‖²_{H¹} := ∫(|Δ|² + |Δ′|²)dt. (Used in the variable X = ελ, whose scale is O(1); the t-variable is the right one for the weight (1+x²)⁻¹.)

## 2. The exact characteristic function and its Gamma factorisation [H]

Lemma R (Laplace transform of the stationary intensity, `problem11_dichotomy.md`/Theorem T Step 1): E e^{−uλ} = exp(−βμ∫₀^u s ds/c(s)), c(s) = βs − (1−e^{−αs}) = εβs + h(αs) (h(y) = e^{−y}−1+y), analytic with no zero in Re s ≥ 0, s ≠ 0 (since |1−e^{−αs}| ≤ α|s| < β|s|). On the imaginary axis s = −iτ, τ = εt: ϕ(t) = E e^{itX} = exp(L(t)), L(t) = −βμ∫₀^{−iεt}s ds/c(s).
Write s/c(s) = 1/(εβ + α²s/2) + D_ε(s) (so the first term integrates to −k ln(1−iat): ϕ = G·F̃, G(t) = (1−iat)^{−k} the characteristic function of Gamma(k, a); F̃ = exp(iβμε∫₀^tD_ε(−iεη)dη)). Algebra: D_ε = [α²s²/2 − h(αs)]/(c(s)(εβ+α²s/2)).
With s = −iτ, y = ατ ≥ 0 (negative τ by conjugation):
 numerator N(y) = (1−cos y − y²/2) + i(y − sin y), c = E₁(y) = −(1−cos y) − i((1+r)y − sin y), εβ+α²s/2 = α(r − iy/2).
Define d(y) := αD_ε(−iy/α) = N/(E₁(r − iy/2)). Then F̃ = exp(z), z(x) = i(k/2)∫₀^x d(y)dy (x = αεt), and L′(x) = (k/2)x/E₁(x).

## 3. Elementary bounds [H]

(a) **|d(y)| ≤ 19 for all y ≥ 0, ε ∈ (0,1).**
 y ≤ 1: |Re N| = y⁴/24 − … ≤ y⁴/24, |Im N| ≤ y³/6, so |N| ≤ y³(1/6 + 1/24) ≤ 0.2084y³; |E₁| ≥ 1−cos y ≥ y²/2 − y⁴/24 ≥ 0.4583y²; |r − iy/2| ≥ y/2. So |d| ≤ 0.2084/(0.4583·0.5) ≤ 0.91.
 y ≥ 1: since 1−cos y ≤ y²/2, Re N = (1−cos y) − y²/2 ∈ [−y²/2, 0], so |Re N| ≤ y²/2; |Im N| = y − sin y ≤ y. Hence |N| ≤ y²/2 + y ≤ 1.5y². |E₁| ≥ |Im E₁| = (1+r)y − sin y ≥ y − sin y ≥ (1 − sin 1)y ≥ 0.158y (because sin y/y ≤ sin 1 for y ≥ 1). So |d| ≤ 1.5y²/(0.158y·y/2) = 18.99 < 19. ∎
 [V] The true sup over y, ε is 1.1429 (ε → 0; attained near y ≈ 7.4); 0.0999 at ε = 0.9.
 Consequences: |z(x)| ≤ 9.5k x, i.e. with M₁ := 19kα/2 (units of t): |F̃(t)| ≤ e^{M₁ε|t|}, |F̃′(t)| ≤ M₁ε|F̃(t)|.

(b) **|ϕ| is non-increasing on t ≥ 0.** Re L′(x) = (k/2)x Re(1/E₁) = −(k/2)x(1−cos x)/|E₁|² ≤ 0. [H]

(c) **|L′| majorants** (L′ per unit t is αε(k/2)x/E₁): for x ≤ 1: |E₁| ≥ 0.4583x² gives |L′(t)| ≤ k/(0.9166t) ≤ 1.1k/t. For x ≥ 1: |E₁| ≥ 0.158x gives |L′(t)| ≤ (k/2)αε/0.158 ≤ 3.2kαε. [H]. [V]: max of |L′|t/k on x ≤ 1 is 1.028 (< 1.1), and both majorants hold on the whole grid for ε = 10⁻², 10⁻³, 10⁻⁴.

(d) **Upper bound for |E₁|² when x ≥ 5**: |E₁|² ≤ 4 + ((1+r)x+1)² = (1+r)²x²(1 + 2/((1+r)x)… ) ≤ (1+r)²x²(1+1/x)² + 4 ≤ (1+r)²x²(1 + 2/x + 5/x²) ≤ (1+r)²x²(1+3/x) (x ≥ 5). [H]

## 4. Head region |t| ≤ t₁ := 1/(M₁ε) [H]

There |z| ≤ M₁ε|t| ≤ 1 so |F̃| ≤ e and |F̃−1| ≤ e^{|z|}−1 ≤ (e−1)|z| ≤ (e−1)M₁ε|t|.
 |Δ| = |G||F̃−1| ≤ (e−1)M₁ε|t||G|.
 Δ′ = G′(F̃−1) + GF̃′, |G′| = ka(1+a²t²)^{−(k+1)/2} ≤ (k/|t|)|G|, so |Δ′| ≤ M₁ε(k(e−1)+e)|G|.
Hence ∫_{|t|≤t₁}|Δ|² ≤ (e−1)²M₁²ε²J, J = ∫_{|t|≤t₁}t²(1+a²t²)^{−k}dt (= a⁻³B(3/2,k−3/2) if bounded by the full line, k > 3/2; for k ≤ 3/2, J ≤ 2a^{−2k}t₁^{3−2k}/(3−2k)); ∫|Δ′|² ≤ M₁²ε²(k(e−1)+e)²·a⁻¹B(½,k−½). [V]: both pointwise inequalities hold on the exact functions for ε = 10⁻², 10⁻³, 10⁻⁴ (note the Δ′ inequality is in t, i.e. αε × the x-derivative).
Order: M₁ ∝ k, t₁ ∝ 1/ε, so head² = O(ε²·(1 + t₁^{3−2k})) = O(ε^{min(2,2k−1)}).

## 5. Tail region |t| > t₁ [H]

At t₁: M₁εt₁ = 1, so by (b) and the head bound, for all t ≥ t₁: |ϕ(t)| ≤ |ϕ(t₁)| ≤ e|G(t₁)| ≤ e(at₁)^{−k} =: A, with at₁ = (1−ε)/(19kε) (large).
For x ≥ x_* := X₀ := max(5, 6/(k₀−1)) (t_* = X₀/(αε)): from Re L′ in (b) and (d),
 Re L(t) − Re L(t_*) ≤ −(k/(2(1+r)²(1+3/X₀)))∫_{X₀}^x (1−cos y)/y dy ≤ −q[ln(x/X₀) − 2/X₀], **q = k₀/(2(1+3/X₀)) ≥ k₀/(k₀+1) > ½** (this is where k₀ > 1 enters; recall k/(1+r)² = k₀ exactly), using ∫_{X₀}^x cos y/y dy ≤ 2/X₀ by parts.
So |ϕ(t)| ≤ A e^{k₀/X₀}(t_*/t)^q for t ≥ t_*, and ≤ A on [t₁, t_*].
Integrals (two sides, t₁ < t_* since 2/(19k) < 5):
 ∫_{|t|>t₁}|ϕ|² ≤ 2A²[(t_* − t₁) + e^{2k₀/X₀}t_*/(2q−1)] =: I_ϕ,
 ∫_{|t|>t₁}|ϕ′|² ≤ 2[(1.1k)²A²(1/t₁ − αε… ) + (3.2kαε)²A²((t_* − t_b) + e^{2k₀/X₀}t_*/(2q−1))] =: I_ϕ′ (t_b = 1/(αε); piece x ≤ 1 uses |ϕ′| = |ϕL′| ≤ A·1.1k/t, piece x ≥ 1 uses 3.2kαε),
 ∫_{|t|>t₁}|G|² ≤ 2a^{−2k}t₁^{1−2k}/(2k−1) =: I_G, ∫_{|t|>t₁}|G′|² ≤ 2k²a^{−2k}t₁^{−2k−1}/(2k+1) =: I_G′.
Orders: A² = O(ε^{2k}), t_* = O(1/ε): I_ϕ = O(ε^{2k−1}), I_ϕ′ = O(ε^{2k+1}), I_G = O(ε^{2k−1}), I_G′ = O(ε^{2k+1}). [V]: |ϕ| ≤ majorant on the exact function wherever t ≥ t₁, ε = 0.1…10⁻⁸.

## 6. Assembly

Δ = ϕ − G, so on the tail ‖Δ‖² ≤ 2(I_ϕ + I_G) and ‖Δ′‖² ≤ 2(I_ϕ′ + I_G′). Define
 **Φ(k₀, ε)² := (e−1)²M₁²ε²J + M₁²ε²(k(e−1)+e)²a⁻¹B(½,k−½) + 2(I_ϕ + I_G) + 2(I_ϕ′ + I_G′).**
Then ‖ϕ−G‖_{H¹} ≤ Φ and by Lemma F, TV(π_θ, Γ_θ) ≤ Φ/(2√2). The leading behaviour as ε → 0 is the head term, ∝ ε (k₀ = 2: Φ = 2.4·10⁻² at ε = 10⁻⁴ is 99.8 % head). Φ ≥ true ‖ϕ−G‖_{H¹} in every case tested (ε from 0.1 to 10⁻⁸), by factors 75 (10⁻⁸) to 7·10⁴ (0.1).
Caveats: Φ is evaluated by floating-point quadrature of closed-form majorants (and one closed-form beta function), not by interval arithmetic; "proved" means the chain of inequalities in §§3–5 is proved by hand and the final number is computed. The tail is not tight at moderate ε (k₀ = 2, ε = 0.1: Φ = 8.5·10³, vacuous; the bound is informative for ε ≲ 10⁻⁴).

## 7. Step 7 (two Gamma laws) and the full TV bound [C]

For the ridge pair θ′ = (μ+δ, α(1−εδ/μ)) (β = 1): k′ = 2μ′/α′², scales a/ε, a′/ε′, ε′ = 1−α′ = ε + (1−ε)εδ/μ. KL(Γ(k,s)‖Γ(k′,s′)) = (k−k′)ψ(k) − ln Γ(k) + ln Γ(k′) + k′ ln(s′/s) + k(s−s′)/s′ (closed form). Computed for μ = 1 (k₀ = 2), ε = 10⁻² and 10⁻⁴: KL/δ² = 0.299, 0.291, 0.282, 0.265, 0.251 (δ = 0.01, 0.05, 0.1, 0.2, 0.3; ε = 10⁻²) and 0.288 … 0.240 (ε = 10⁻⁴): ≤ 0.30. Pinsker: TV(Γ_θ, Γ_θ′) ≤ √(KL/2) ≤ 0.39δ for δ ≤ 0.3. (The ridge direction keeps the mean μ/ε almost fixed, ratio 1+O(εδ), so the KL is dominated by the shape change k′/k ≈ 1+δ/μ.)
The KL bound is a numerical evaluation of an exact closed form on a δ-grid, not an analytic bound; a proof that sup_{δ≤0.3}KL/δ² ≤ 0.31 would use the Bregman-divergence representation (log-partition Hessian bounded on the segment), unwritten.

**Theorem T″ (explicit TV for the ridge pair, k₀ = 2, μ = β = 1).** For ε ≤ 10⁻⁵ and |δ| ≤ 0.3:
 TV(π_θ, π_θ′) ≤ [Φ(2,ε) + Φ(k₀′,ε′)]/(2√2) + 0.39|δ| ≤ 2.4·10⁻³ + 0.39|δ|,
where the first term is [Φ(2,10⁻⁵) + Φ(k₀′, ε′)]/(2√2) evaluated at ε = 10⁻⁵ for δ = 0, 0.1, 0.3 (k₀′ = 2(1+δ), ε′ = ε(1+(1−ε)δ)): 1.01·10⁻³ + {1.01, 1.00, 1.38}·10⁻³. E.g. TV ≤ 0.041 at δ = 0.1. (For δ between the tested values the same check was not run.)

## 8. Consequences

**Assumption M (K = 1, k₀ = 2, μ = β = 1) is certified for ε ≤ 10⁻⁴.** E_π h_j ≥ E_Γ h_j − TV(π,Γ) ≥ E_Γ h_j − 0.0196, with E_Γ h₁ = 0.164 and E_Γ h₂ = 0.0444 (`code/p5_assumption_M_gamma.py`), so m := min_j E_π h_j ≥ 0.0248 for ε ≤ 10⁻⁴, and ≥ 0.0405 for ε ≤ 3·10⁻⁵ (TV bound 3.9·10⁻³ there). For other k₀ the same subtraction gives a certified ε-range from the table in the Statement.
**Impossibility for the stationary record, K = 1.** The corollary of `theoremA_step2_and_stationary_tv.md` needs TV(Q_θ,Q_θ′) ≤ ½ with ε = c/T. By Theorem T″, for ε ≤ 10⁻⁵ the stationary-law term TV(π_θ,π_θ′) is ≤ 2.4·10⁻³ + 0.39δ, which is below 0.05 for δ ≤ 0.12. The other term of the reduction, √(ρ uᵀI_Tu/2), comes from the KL lemma and the Fisher computation of Problem 10 and was not redone here. So the stationary-law part of the corollary now has explicit constants in this range (k₀ = 2, μ = β = 1, ε ≤ 10⁻⁵); the KL part keeps the constants of the earlier notes.

## 9. What is and is not "fully formal" now

Fully written, hand-proved: Lemma F; the factorisation and N, E₁ identities; bound |d| ≤ 19; monotonicity of |ϕ|; the L′ majorants; the head estimates; the tail decay with q = k₀/(2(1+3/X₀)) > ½; the assembly into Φ; every pointwise inequality was also tested on the exact function.
Not formal: (i) Φ's numerical value comes from floating-point quadrature, not interval arithmetic; (ii) Theorem T″'s KL step is a closed-form evaluation on a grid (analytic Hessian bound unwritten); (iii) the statement Φ(k₀′,ε′) ≤ Φ(2,ε) for the perturbed pair was checked on the displayed values only; (iv) the full corollary also needs the KL lemma and Fisher bound (Problem 10 notes) — unchanged; (v) k₀ ≤ 1 and K ≥ 2 are out of scope (Step 5 needs q > ½).
