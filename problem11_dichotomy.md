# Problem 11: information dichotomy for the baseline μ (K = 1 done, K ≥ 2 reduced)

Builds on [problem10_rigorous.md](problem10_rigorous.md) (Theorem A, exact transient formula). Setup: K = 1, λ = μ + αS, dS = −βS dt + dN, ε = 1 − α/β = c/T (nearly unstable), **k = 2μ/β** (Gamma shape of the CIR limit; equivalently 2μ times the mean delay 1/β). Two observation models: **(S)** stationary start, with S_0 revealed to the statistician (a stronger statistician than one who does not see S_0, so lower bounds transfer); **(C)** cold start S_0 = 0.

## 1. Fisher information: exact statements

"Efficient information for μ" = 1/(I⁻¹)_μμ, the information when α is unknown.

**Proposition 1 (stationary).** eff = (m − ε/μ)/(1 − 2ε + εμm) exactly, with m = E_π[1/λ]. (Proof: eff = det/I_αα with det = (mEλ − 1)/α², I_αα = (Eλ − 2μ + μ²m)/α², Eλ = μ/ε.) With Theorem A and ε = c/T:

| k | T·I_μμ (α known) | T·eff (α unknown) | verdict |
|---|---|---|---|
| > 1 | → 2c/(2μ−β) | → cβ/(μ(2μ−β)) | bounded |
| = 1 | ~ (2c/β) ln(T/c) | same leading order | log-divergent |
| < 1 | ~ C c^k T^{1−k} | same leading order | polynomially divergent |

Numerics (c = 2, β = 1, exact m): k = 2, T·eff = 1.74, 1.95, 1.99 at T = 10², 10³, 10⁴ (limit 2); k = 1: 11.2, 19.9, 29.1 (+9.2 per decade, predicted 4·ln 10 = 9.2); k = 0.6: 42.8, 129, 349 (ratio per decade ≈ 10^{0.4} = 2.5, observed 3.0, 2.7).

**Proposition 2 (cold start).** M_T = I_μμ diverges for every k: M_T ≥ ln(1 + βT)/(μβ) (rigorous, Problem 10). Numerically M_T ≈ A ln T with A → 2/(β(2μ−β)) for k > 1, ≈ ln² T for k = 1, ≈ T^{1−k} for k < 1 (conjectural coefficients; Problem 10 §4).

## 2. Impossibility and lower bounds (rigorous for K = 1)

**Lemma (KL bound).** For θ = (μ, α) and θ' = (μ', α') with the same initial state,
KL(P_θ^T ‖ P_θ'^T) = E_θ ∫₀^T [λ log(λ/λ') − (λ − λ')] dt ≤ E_θ ∫ (λ − λ')²/λ' dt ≤ ρ · uᵀ I_T(θ) u,
where u = θ' − θ and ρ = max(μ/μ', α/α') ≥ sup λ/λ'. (First inequality: x log(x/y) − x + y ≤ (x−y)²/y. Second: 1/λ' ≤ ρ/λ, and I_T = E∫ a aᵀ/λ dt.)

**Theorem B (K = 1).** Let θ' = (μ + δ, α − wδ) with w ≥ 0, so uᵀI_Tu = δ²(I_μμ − 2wI_μα + w²I_αα). If ρ·uᵀI_Tu ≤ 1/2 then by Le Cam and Pinsker, for every estimator μ̂, max_{θ,θ'} P(|μ̂ − μ| ≥ δ/2) ≥ 1/4.
1. **(S), k > 1, ε = c/T.** Take w = I_μα/I_αα (≈ εβ/μ, so ρ = 1 + O(ε)): uᵀI_Tu = δ²·T·eff → δ² c β/(μ(2μ−β)), bounded in T. So for δ below an explicit constant (δ₀² = μ(2μ−β)/(2ρcβ)), **no estimator is consistent for μ**, even when α is unknown, and even more so when α is known (δ₀² = (2μ−β)/(4ρc) with w = 0).
2. **(S), k < 1:** uᵀI_Tu ≲ δ² C c^k T^{1−k}, so error ≳ T^{−(1−k)/2} (lower bound only); k = 1: error ≳ (ln T)^{−1/2}.
3. **(C), k > 1:** uᵀI_Tu = δ²(M_T + O(1)); numerically 6.75, 11.4, 16.0 per δ² at T = 10², 10³, 10⁴ (≈ +4.6 per decade = 2 ln 10). **Conditional on M_T ≲ A ln T** (not proved), the error is ≳ (ln T)^{−1/2}: μ is estimable only at a logarithmic rate.

For comparison, `code/p11_twopoint.py` prints the largest δ with ρ·uᵀI_Tu ≤ 1/2 under the cold start (C) for k = 2: 0.27, 0.21, 0.18 at T = 10², 10³, 10⁴ (ρ taken as 1). It shrinks with T like (ln T)^{−1/2}, whereas under (S) the threshold δ₀ above does not depend on T at all.

## 3. Attainability (numerics only)

`code/p11_mle_mc.py`: MLE of μ (cold start, ε = 2/T, β = 1, thinning simulation, pure Python), standard deviation of μ̂ against 1/√M_T:

| μ (k) | T | R | α known: sd (robust sd) | 1/√M_T | joint (μ, α): sd (robust) | √(I⁻¹)_μμ |
|---|---|---|---|---|---|---|
| 1 (2) | 100 | 150 | 0.332 (0.335) | 0.327 | 0.669 (0.512) | 0.399 |
| 1 (2) | 400 | 75 | 0.314 (0.298) | 0.285 | 0.472 (0.306) | 0.334 |
| 0.3 (0.6) | 100 | 150 | 0.179 (0.123) | 0.114 | 0.247 (0.175) | 0.120 |
| 0.3 (0.6) | 400 | 75 | 0.115 (0.106) | 0.081 | 0.159 (0.130) | 0.083 |

For α known the MLE sits at the information scale within sampling noise at k = 2 (and shows the predicted very slow decay with T), and within about 10–30 % at k = 0.6. The joint MLE has finite-sample bias (mean μ̂ = 1.38 at T = 100 for μ = 1) and heavy tails, so efficiency is not demonstrated for α unknown; this is likely a finite-T boundary effect (n̂ near 1) but I did not test that. R is small (75–150), so these are ±10 % numbers.

## 4. K ≥ 2

**Proposition 3 (general K, rigorous).** Let a = (1, S_1..S_K), m_K = E_π[1/λ]. For the direction v = (1, −εα/μ): vᵀa = 1 + ε − ελ/μ, so
vᵀ I v = (1+ε)² m_K − (ε/μ)(1 + 2ε),
and hence **eff(μ) ≤ (1+ε)² m_K − (ε/μ)(1+2ε)** for every K (checked on K = 1: bound 0.0676 vs eff 0.0675 at ε = 0.1, μ = β = 1; on K = 2 by MC the identity holds to 2 %). With the KL lemma (ρ = 1/(1 − εδ/μ) for the perturbation α ↦ α(1 − εδ/μ)), the K = 1 impossibility argument goes through for any K as soon as **m_K ≲ ε**. So the K ≥ 2 dichotomy reduces to the behaviour of E_π[1/λ].

**Conjecture.** With ν₁ = Σα_k/β_k² (mean delay of the kernel at n = 1) and k_K = 2μν₁: m_K ≈ 2ε/(2μ − 1/ν₁) for k_K > 1 (this reduces to Theorem A for K = 1), ~ ε^{k_K} for k_K < 1. Evidence (K = 2, β = (1, 2), α_k = (1−ε)β_k/2, μ = 1, so ν₁ = 0.75, k_K = 1.5, predicted m/ε → 3): MC gives m/ε = 1.81 (ε = 0.1) and 2.16 (ε = 0.03). It is rising, but the K = 1 sequence also approaches its limit slowly (0.78 → 0.89 of the limit over the same range, against 0.60 → 0.72 here), so this is consistent with the conjecture, not a confirmation (250 000-event runs, ±10 %).

The Laplace transform of λ depends on the whole K-dimensional state and has no closed-form ODE, so the Gamma factorisation of Theorem A does not carry over. A Lyapunov argument does, in a sub-region:

**Theorem C (any K, stationary).** Let b = β_min, B = β_max, n = 1 − ε. If 2μ > n²B then
**E_π[1/λ] ≤ 2εB / (b(2μ − n²B)).**

*Proof.* Let c_k = α_k/β_k, Z = Σc_kS_k (so A Z = μ − ελ, Problem 3), z₀ = μ/B, W = Z + z₀ ≥ z₀.
1. The drift of Z is −Σα_kS_k = −(λ − μ) and Z jumps by n at rate λ, so A log W = −(λ−μ)/W + λ log(1 + n/W).
2. log(1+x) ≥ x − x²/2 for x ≥ 0 gives A log W ≥ (μ − ελ)/W − n²λ/(2W²).
3. λ = μ + Σβ_kc_kS_k ≤ μ + BZ = BW, so λ/W ≤ B. Hence (μ−ελ)/W ≥ μ/W − εB and n²λ/(2W²) ≤ n²B/(2W).
4. Stationarity: E_π[A log W] = 0 (justified below). Therefore 0 ≥ (μ − n²B/2) E[1/W] − εB, i.e. E[1/W] ≤ 2εB/(2μ − n²B).
5. λ ≥ μ + bZ ≥ bW (because μ ≥ bμ/B), so 1/λ ≤ 1/(bW) and E[1/λ] ≤ E[1/W]/b. ∎

*Justification of E_π[A log W] = 0.* log W is integrable (all moments of Z are finite by Problem 3, and W ≥ z₀ > 0). A log W is bounded: |(λ−μ)/W| ≤ B and 0 ≤ λ log(1+n/W) ≤ nλ/W ≤ nB. Jumps of log W are at most n/z₀. So log W(S_t) − log W(S_0) − ∫₀^t A log W ds is a true martingale (E∫λ ds < ∞), and stationarity gives E_π[A log W] = 0.

For K = 1 (B = b = β) this reads m ≤ 2ε/(2μ − n²β), valid for every ε with 2μ > n²β, and its ε → 0 limit is exactly Theorem A's constant 2/(2μ−β). It reproduces the k > 1 upper bound of Theorem A without the Gamma factorisation. Numerically (exact m, μ = β = 1): m/bound = 0.95, 0.92, 0.97 at ε = 0.3, 0.1, 0.01. For K = 2, 3 the bound holds but is loose by a factor of 2.5 to 9 (K = 2, β = (1,2), μ = 1.5: MC m = 0.103, 0.039 vs bound 0.290, 0.107 at ε = 0.1, 0.03; K = 3, β = (1,2,4), μ = 2.5: 0.057, 0.021 vs 0.455, 0.194; `code/p11_lyapunov_check.py`).

**Corollary D (impossibility for general K, rigorous in the region 2μ > β_max).** Under stationary observation with ε = c/T, perturb θ ↦ θ' = (μ + δ, α(1 − εδ/μ)). Then ρ = 1/(1 − εδ/μ), and by the KL lemma, Proposition 3 and Theorem C,
KL ≤ ρ T δ² [(1+ε)² m_K − (ε/μ)(1+2ε)] ≤ ρ(1+ε)² c δ² C_K, with C_K = 2B/(b(2μ − n²B)),
which does not grow with T. Le Cam gives max P(|μ̂ − μ| ≥ δ/2) ≥ 1/4 once ρ(1+ε)² c δ² C_K ≤ 1/2, so **μ is not consistently estimable** for any K in the region 2μ > β_max(1+o(1)), even when α is unknown. For K = 1 this is the whole range k > 1.

*Theorem C and Corollary D cover only the region 2μ > β_max. The gap to the conjectured threshold 2μν₁ > 1 is closed in §4b below, which supersedes the "conjecture" paragraph above.*

## 4b. Closing the K ≥ 2 gap: sharp threshold, both sides (Theorems E and F)

**Setup.** α_k > 0, β_1 < … < β_K, n = Σα_k/β_k = 1 − ε, b = β_1, B = β_K. Let γ ∈ (0, b) be the root of Σα_k/(β_k − γ) = 1 (it exists and is unique: the left side increases from n < 1 at γ = 0 to +∞ at γ = b). Put
- ω_k = α_k/(β_k − γ), so Σω_k = 1 and β_kω_k − α_k = γω_k; Z = ωᵀS;
- ν̃ = Σω_k/β_k, c₀ = n/(2ν̃), and **k_K := μ/c₀ = 2μν̃/n**. As ε → 0: γ → 0, ν̃ → ν₁ = Σα_k/β_k², k_K → 2μν₁.
- Exact relation: 1 = n + γν̃, so **γ = ε/ν̃**.
- κ_k = (c₀ω_k − α_k/2)/β_k (so Σκ_k = 0), ψ = κᵀS, κ̄ = max|κ_k|/ω_k.

**Lemma G.** A Z = μ − γZ (since Σ(α_k − β_kω_k)S_k = −γZ), hence E_π Z = μ/γ; and (b−γ)Z ≤ λ − μ = Σ(β_k − γ)ω_kS_k ≤ (B−γ)Z.

**Lemma H (corrected Lyapunov function).** For z > 0 let W = Z + z, f_z = log W + ψ/W², B_z = max(μ/z, B−γ), C₂ = 2κ̄μ + 3B_zκ̄, C₄ = B_z/3 + C₂. Then, pathwise,

(μ − c₀ − 2κ̄γ)/W − γ − (μ/2 + C₂)/W²  ≤  A f_z  ≤  (μ − c₀ + γz + 2κ̄γ)/W − γ + (c₀z + C₄)/W².

*Proof.* (i) A log W = (μ − γZ)/W + λ[log(1 + 1/W) − 1/W], because W jumps by Σω_k = 1 and has drift −Σβ_kω_kS_k. (ii) With q = ψ/W² and Σκ_k = 0 (ψ does not jump): A q = −Σκ_kβ_kS_k/W² + 2ψ(Σβ_kω_kS_k − λ)/W³ + λψ(3W+2)/(W³(W+1)²), using 1/W² − 1/(W+1)² = 2/W³ − (3W+2)/(W³(W+1)²). (iii) Σβ_kω_kS_k − λ = γZ − μ, |ψ| ≤ κ̄Z ≤ κ̄W, λ ≤ B_zW (since λ ≤ μ + (B−γ)Z), and (3W+2)/(W+1)² ≤ 3/W, so the last two terms of A q are bounded in absolute value by 2κ̄γ/W + 2κ̄μ/W² + 3B_zκ̄/W². (iv) The choice of κ gives λ/2 + Σκ_kβ_kS_k = μ/2 + c₀Z, hence −λ/(2W²) − Σκ_kβ_kS_k/W² = −μ/(2W²) − c₀/W + c₀z/W². (v) −x²/2 ≤ log(1+x) − x ≤ −x²/2 + x³/3 for x > 0, and (μ−γZ)/W = μ/W − γ + γz/W. Collect the terms. ∎
Checked pathwise on 200 000 random states for each of 5 parameter sets (K = 1, 2, 3, covering both k_K < 1 and k_K > 1; `code/p11_corrector_check.py`): the exact generator always lay between the two bounds (all minimum slacks ≥ 0).

**Why the corrector matters.** The leading noise term of A log W is λ/(2W²), and λ/W varies between b and B over the state space. Theorem C bounded it by its worst case B. The corrector ψ/W² cancels the fluctuation of λ/W around its typical value and leaves the coefficient c₀ = n/(2ν̃) ≈ 1/(2ν₁) (the CIR value), which is where the sharp threshold comes from.

**Justification of E_π A f_z = 0.** f_z is integrable (log W with W ≥ z, and ψ/W² is bounded), A f_z is bounded, and the jumps of f_z are bounded, so f_z(S_t) − f_z(S_0) − ∫A f_z is a true martingale (E∫λ < ∞); stationarity gives E_π A f_z = 0. The same argument allows differentiating E[1/(Z+z)] in z.

**Theorem E (any K, stationary, upper bound).** If k_K > 1, put δ₀ = μ − c₀ = c₀(k_K − 1). For ε small enough that 2κ̄γ ≤ δ₀/4, and z ≥ 4(μ/2 + C₂)/δ₀,
E_π[1/W] ≤ 2γ/δ₀, hence **E_π[1/λ] ≤ max(z/μ, 1/(b−γ)) · 2ε/(ν̃δ₀) = O(ε).**
*Proof.* The lower bound of Lemma H, 1/W² ≤ 1/(zW), and E A f_z = 0 give 0 ≥ (δ₀ − δ₀/4 − δ₀/4)E[1/W] − γ. Finally λ ≥ μ + (b−γ)Z, so W/λ ≤ max(z/μ, 1/(b−γ)), and γ = ε/ν̃. ∎

**Corollary E′ (impossibility, any K, the whole region k_K > 1).** With the perturbation of Corollary D, KL ≤ ρ(1+ε)²δ²·T·E_π[1/λ] ≤ Cρ(1+ε)²cδ², which does not depend on T. By Le Cam, **μ is not consistently estimable whenever k_K > 1**, with or without α known, for every K. This replaces the region 2μ > β_max of Corollary D by the conjectured one, 2μν₁ > 1 (as ε → 0). For K = 1 it is the whole range k > 1.

**Theorem F (any K, stationary, lower bound).** If k_K < 1 then for every η > 0 there are c_η > 0 and ε_η > 0 such that for ε ≤ ε_η
**E_π[1/λ] ≥ c_η ε^{k_K + η}**, so in particular E_π[1/λ]/ε → ∞.
*Proof.* Let a(z) = E_π[1/(Z+z)], so −a′(z) = E_π[1/W²]. The upper bound of Lemma H, E A f_z = 0, and −γ < 0 give (c₀z + C₄)(−a′(z)) ≥ γ + (c₀ − μ − γz − 2κ̄γ)a(z). For 1 ≤ z ≤ η/γ the bracket is at least θ := c₀ − μ − η − 2κ̄γ, which is positive for small η and γ because c₀ > μ. Integrating −a′/a ≥ θ/(c₀z + C₄) from 1 to η/γ:
a(1) ≥ a(η/γ)·[(c₀η/γ + C₄)/(c₀ + C₄)]^{θ/c₀}.
By Jensen and E Z = μ/γ (Lemma G), a(η/γ) ≥ 1/(E Z + η/γ) = γ/(μ + η). So a(1) ≥ c·γ^{1 − θ/c₀}, and 1 − θ/c₀ = k_K + (η + 2κ̄γ)/c₀. Finally 1/λ ≥ 1/(max(μ, B)(Z + 1)), so E[1/λ] ≥ a(1)/max(μ, B), and γ = ε/ν̃. ∎
Consequence: the α-known information satisfies T·I_μμ = T·E[1/λ] ≥ c·T^{1 − k_K − η} → ∞ at ε = c/T. The exponent k_K matches the K = 1 constant of Theorem A and the conjecture. Theorem F gives a lower bound only; it does not prove E[1/λ] ≲ ε^{k_K − η}.

**Numerical checks** (`code/p11_theoremEF_check.py`, `code/p11_efficient_K2.py`; stationary Monte Carlo with the corrected simulator):
- Theorem E vs MC m (explicit constants, z optimised): K = 2, β = (1,2), μ = 0.8 (the old gap region, k_K = 1.21): MC m = 0.040 = 4.0ε at ε = 0.01, bound = 596ε (valid but very loose because the margin μ − c₀ is small). μ = 1.5: bound 11ε vs MC 1.1ε. K = 3, β = (1,2,4), μ = 2.5: bound 7.2ε vs MC 0.60ε.
- Theorem F: K = 2, β = (1,2), μ = 0.4 (k_K ≈ 0.61): m = 0.80, 0.45, 0.29 at ε = 0.1, 0.03, 0.01, local slopes of ln m against ln ε of 0.48 and 0.41, below k_K as a lower bound allows. K = 3, β = (1,3,9), μ = 1 (k_K ≈ 1.0, borderline): slopes 0.71, 0.52. MC noise is about ±10 % per point (±0.1 on a slope), so this is consistent with the theorem, not a test of its sharpness.
- Efficient information for K ≥ 2 (α unknown): eff/m = 0.74, 0.82, 0.87 (K = 2, μ = 0.4) and 0.72, 0.81, 0.88 (K = 3, β = (1,3,9), μ = 0.5) at ε = 0.1, 0.03, 0.01, rising toward 1, with eff/ε = 5.9, 11.2, 17.7 for K = 2. This supports eff ≈ m − O(ε) and hence divergence of the α-unknown information for k_K < 1, but **I do not have a proof of that step for K ≥ 2** (for K = 1 it is exact: eff = (m − ε/μ)/(1 − 2ε + εμm)).

## 4c. The matching upper bound: what is and is not shown

For k < 1 (K = 1, α known) the Fisher information diverges like T^{1−k}, and Theorem B gives only the lower bound T^{−(1−k)/2} on the error. A matching rate for an estimator is **not proved here**, and the obstruction is structural:
- The score ∫(1/λ)(dN − λ dt) is a martingale with predictable quadratic variation 𝓘_T = ∫₀^T dt/λ_t. With ε = c/T the state has a single time scale 1/(βε) = T/(βc), so ∫dt/λ_t does not average: it stays random at O(1) relative size. The expected information E 𝓘_T is therefore the wrong normaliser; the right one is the observed 𝓘_T (a mixed-normal, LAMN-type structure, as for the CIR drift).
- Numerics (`code/p11_lamn_check.py`, K = 1, cold start, T = 100, ε = 0.02, 400 paths): the coefficient of variation of 𝓘_T is 0.45 (μ = 1, k = 2) and 0.52 (μ = 0.3, k = 0.6). Normalising the MLE error by the observed information gives sd 1.07 and 1.04 with excess kurtosis +1.5 and −0.02; normalising by the expected information gives sd 1.27 and 1.63 with excess kurtosis about +7.7. So the MLE is close to mixed normal under the observed-information scaling and heavy-tailed under the expected one. This also explains why the earlier MLE standard deviations were noisy and above 1/√(E 𝓘). (These runs are cold start, K = 1; the mixed-normal structure is the point, not the stationary rate.)
- A proof would need (i) a stable martingale CLT for the score under the time change 𝓘_T, (ii) 𝓘_T → ∞ in probability (not just in expectation) when k < 1, and (iii) control of the observed second derivative ∫dN/λ² near μ̂. These steps are standard in spirit (Jaisson–Rosenbaum scaling, CIR-drift LAMN theory) but I have not carried them out, so the lower bounds for k ≤ 1 are not shown to be tight.

## 5. What is and is not established for Problem 11 (updated)

- **Established, rigorous, stationary observation.**
  - K = 1: exact information formulas, the full dichotomy (T·eff bounded iff k > 1), impossibility for k > 1 (Theorem B.1), and divergence of the information for k ≤ 1 (Theorem A).
  - Any K, with k_K := 2μν̃/n → 2μΣα_k/β_k²: E_π[1/λ] = O(ε) when k_K > 1 (Theorem E), and E_π[1/λ] ≥ c ε^{k_K+η} when k_K < 1 (Theorem F). Hence for every K: **no consistent estimator of μ when k_K > 1** (Corollary E′), and **divergence of the α-known information when k_K < 1**. The role of the shape parameter is explicit: k_K = 2μ × (mean delay of the kernel at n = 1).
- **Not established:**
  - Divergence of the α-unknown (efficient) information for K ≥ 2 and k_K < 1 (numerically eff/m → 1).
  - A matching upper bound: an estimator with a theorem for k_K < 1 (§4c).
  - The borderline k_K = 1 for K ≥ 2.
  - The regime where μ or β varies with ε.
  - Cold-start coefficients (conditional, Problem 10).
- **Decision (2026-10-25, by you): C4 uses stationary observation (model S).** The dichotomy to state is *bounded vs divergent*: with ε_T = c/T, the information about μ stays bounded iff k_K > 1, and for k_K > 1 no estimator of μ is consistent. Cold-start results (§1 Prop 2, Theorem B.3) are remarks showing the statement is sensitive to the observation model. In the paper, state the model explicitly (the process is observed in its stationary regime, S_0 known for the lower bound).
- As with the other notes, re-derive the KL lemma, Proposition 1 and the generator computation of Lemma H yourself before building on them; the pathwise numerical check supports the algebra but is not a proof.


---
## Audit correction (25 Oct 2026, final audit): what the impossibility results do and do not prove

The sentences in §2 (Theorem B.1), §4 (Corollary D), §4b (Corollary E′) and §5 that say "no estimator of μ is consistent ... given S₀ / under stationary observation" are **stronger than what the KL lemma proves**. The KL lemma compares path laws that start from the same initial state; the stationary experiment starts from S₀ ~ π_θ versus π_θ′, whose laws differ. The correct two-point bound is
 TV(Q_θ, Q_θ′) ≤ TV(π_θ, π_θ′) + √(ρ T uᵀI(θ)u / 2)  (chain rule for TV; stationarity gives E_{π_θ}KL_cond = T·E_π[λ log(λ/λ′) − λ + λ′] ≤ ρT uᵀIu).
The path term is what Theorems B/E bound (it stays bounded when k_K > 1); the marginal term TV(π_θ, π_θ′) is needed for an impossibility statement, and treating S₀ as ancillary changes the information (from S₀ = 0 it diverges like ln T). Rigorous impossibility statements: K = 1, stationary record: `theoremA_step2_and_stationary_tv.md` (TV(π_θ, π_θ′) = O(δ), explicit constants in `thmT_explicit_constants.md`); any K, stationary record with unobserved past: `general_K_impossibility.md` (Theorem V). For K ≥ 2 with S₀ *observed*, inconsistency is open. Details: `final_audit.md` item 1.
