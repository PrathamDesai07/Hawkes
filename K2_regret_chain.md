# The regret chain for general K ≥ 2

Goal: carry Theorem R (expected regret (d/2 + o(1)) ln T + C, `problem05_E5_expected_regret.md`) from K = 1 to any K, d = K + 1.
Tags: **[P]** proved on paper, **[P\*]** proved modulo listed inputs, **[N]** numerical, **[O]** open. Code: `code/p5_k2_regret.py`.
Outcome in one line: every K = 1 ingredient that used the scalar structure has a K-general replacement that I can prove; the one K-specific obstacle of the K = 1 chain, Assumption M, is **removed** (also for K = 1) by a fully explicit device (Theorem FK below). What does not carry over is the quality of the constants: they are astronomically conservative for K ≥ 2.

## 1. Ingredients that carry over verbatim (checked line by line)

Setting: Θ_c = {μ ∈ [μ_lo, μ_hi], α_k ∈ [α_lo,k, α_hi,k], Σα_k/β_k ≤ 1 − ε/2} ∋ θ\*, a convex polytope; F-norm projection Π^F onto Θ_c (a convex QP; existence and uniqueness are all the proofs need).

- **Regret reduction** (Prop. `prop:reg`): uses only convexity of Θ_c, θ\* ∈ Θ_c and non-expansiveness of the metric projection in the F_i-norm. **[P]**
- **Ratio set** (Prop. `prop:ratio`): for every state s ≥ 0, (μ + Σα_k s_k)/(μ\* + Σα\*_k s_k) lies between the extreme coordinate ratios, so c⁻λ\* ≤ λ_θ ≤ c⁺λ\* with c⁺ = max(μ_hi/μ\*, max_k α_hi,k/α\*_k), c⁻ analogous. **[P]**
- **Score/drift/log-det lemmas and E5a–E5d, Proposition E5, Lemma E4, Lemma C1**: the proofs use only λ\* ≥ μ\*, the compensator structure of the first event in a gap (time change: E_i = Λ\*(τ_i) ~ Exp(1)), and linear algebra in dimension d. The point leverage is defined as ℓ̄_i ≥ sup_{u ≥ 0} φ(u)ᵀF_{i−1}⁻¹φ(u)/λ\*(u) (φ(u) = (1, S e^{−βu}), a path in R^{K+1}), which is the K = 1 definition after the substitution y = S e^{−βu}. **[P]**
- **Log-det growth**: Σ_i log det(F_i/F_{i−1}) = log det(F_n/F_0) ≤ d ln(λ_max(F_n)/λ_min(F_0)) and λ_max(F_n) ≲ n·E|a|²/μ\*, so log det(F_n/F_0) = d ln T_n + O(1) once J_n/T → I (see §2 for the concentration). **[P\*]**
- **Moments and exponential tails of the state**: Lemma `lem:moments`, Lemma `lem:expmoment` (E exp(εZ/e) ≤ 2e^{4μ/(eb)}, Z = Σ(α_k/β_k)S_k), so max_{i ≤ n} S_k(t_i) ≤ (e/ε)(ln C + ln n + x)/ℓ_min w.p. ≥ 1 − e^{−x}. This replaces Lemma E0 (which was K = 1 only). **[P]**
- **Contraction with the exact rate** γ in the ℓ-weighted L¹ metric (Lemma `lem:ellcontr`) is K-general. **[P]**

## 2. Replacement of the two K = 1 pieces

### 2.1 The identities of Lemma E1 for general K  **[P]**
With c_k = α_k/β_k, Z_c = Σc_kS_k (so Σc_k = n): dZ_c = −(λ − μ)dt + n dN. Writing N = Λ + M′ (M′ = N − Λ) gives the **exact** identity
 εΛ_T = μT + n M′_T − (Z_c(T) − Z_c(0)),
the integrated form of the key identity A Z_c = μ − ελ (Problem 3), valid for every K and identical to the K = 1 formula (a) of Lemma E1. Hence the lower tail of Λ_T (Lemma E2) holds verbatim:
 P(Λ_T ≤ (1 − η)μT/ε, |Z_c(T) − Z_c(0)| ≤ y) ≤ exp(−(ημT/n − y)²ε/(2μT)),
and the upper tail is Bernstein, with Z_c having the exponential moment of Lemma `lem:expmoment`. (An earlier draft of this note used Z_ω = ωᵀS and the inequalities (b − γ)∫Z ≤ Σα_k∫S_k ≤ (B − γ)∫Z, which lose a constant factor B·ν̃ ≥ 1 in the upper tail for K ≥ 2; the exact Z_c identity makes that unnecessary.)

### 2.2 Theorem FK: concentration of the accumulated information without Assumption M  **[P]**
Let ψ_L(λ) = min(1, max(0, 2 − λ/L)) (equal to 1 for λ ≤ L, 0 for λ ≥ 2L, Lipschitz). Since λ\* ≤ 2L on the support of ψ_L,
 J_T = ∫₀^T aaᵀ/λ\* dt ⪰ (1/(2L)) ∫₀^T aaᵀψ_L(λ_t)dt.
Let A_L = E_π[aaᵀψ_L(λ)] and f_ij(s) = a_ia_jψ_L(λ(s)) (entries of the integrand).
1. **Mean.** ψ_L ≥ 1{λ ≤ L}, so λ_min(A_L) ≥ λ_min(M) − (E|a|⁴)^{1/2}(Eλᵖ/Lᵖ)^{1/2} =: λ_L (Theorem `thm:trunc`; M, the moments and the optimisation over (p, L) are explicit and computable). λ_L > 0 for L large enough whenever λ_min(M) > 0.
2. **Lipschitz constant.** On {λ ≤ 2L}: S_k ≤ 2L/α_k, so |a_i| ≤ A_s := max(1, 2L/α_min). Hence |∂f_ij/∂s_k| ≤ 2A_s + A_s²α_k/L and f_ij is Lipschitz in the weighted metric d_ℓ(s, s′) = Σℓ_k|s_k − s′_k| with constant Λ_f = max_k(2A_s + A_s²α_k/L)/ℓ_k.
3. **Poisson equation.** û(s) = ∫₀^∞(P_tf(s) − E_πf)dt satisfies |û(s) − û(s′)| ≤ Λ_f d_ℓ(s, s′)/γ (Lemma `lem:ellcontr`), |û(s)| ≤ Λ_f(ℓᵀs + ℓᵀE S)/γ, and an event changes the state by 𝟙, d_ℓ(S + 𝟙, S) = ℓᵀ𝟙 = 1. So ∫(f − m)dt = û(S_0) − û(S_T) + 𝓜_T with jumps of 𝓜 bounded by b₀ = Λ_f/γ and ⟨𝓜⟩_T ≤ b₀²Λ_T (exactly as in Theorem `thm:E3`).
4. **Freedman.** With Λ_T ≤ (1 + η)μT/ε (Lemma E2, K-general by 2.1) and the boundary terms ≤ τT/4: P(|∫₀^T(f_ij − m_ij)dt| ≥ τT) ≤ 3 exp(−c_f T), c_f = τ²/(32 b₀²(1 + η)μ/ε + 8b₀τ/3) ≍ τ²γ²ε/(Λ_f²μ).
5. **Union.** Over the (K+1)(K+2)/2 entries with τ = λ_L/(2(K+1)): with probability ≥ 1 − (K+1)(K+2)·3e^{−c_f T}, ‖(1/T)∫aaᵀψ_L − A_L‖_F ≤ λ_L/2, hence
 **J_T ⪰ (T λ_L/(4L)) · I_{K+1}.**
Thus Assumption M (a statement about π-mass of two or K+1 balls, impossible to certify for K ≥ 2 with the K = 1 Fourier method) is replaced by the **computable positive number λ_L**. No separation of K + 1 well-placed regions in R^K is needed. The same argument gives, for K = 1, an alternative to `assumptionM_exact.md`.

### 2.3 Leverage from Theorem FK  **[P]**
On the event of Theorem FK and under (L1): F_{i−1} ⪰ J_{t_{i−1}}/c⁺ ⪰ (t_{i−1}λ_L/(4Lc⁺)) I, so the point leverage satisfies
 ℓ̄_i ≤ (4Lc⁺/(t_{i−1}λ_L μ\*))·(1 + K s_max²), s_max(t) ≤ (e/ε)(ln C + 2 ln t)/ℓ_min (Lemma `lem:expmoment`),
i.e. ℓ̄(t) = O((ln t)²/(ε² t)). This is a factor (ln t)²/ε worse than the K = 1 bound O(1/(εt)) (which used the structure of Q), but it tends to zero as t → ∞, which is all that Theorem R uses (the bad events are charged to their own time t, not to the horizon n).

## 3. Theorem R^(K)  **[P\*]**

Under the standing assumptions (known polytope Θ_c ∋ θ\*, exact F-projection, L ≥ L(K, ε, β, μ) with λ_L > 0), for every δ, η > 0 there is C_K = C_K(δ, η, ε, μ, β, Θ_c, F_0) independent of n with
 **E R_n ≤ (½ + O(δ + η))(K + 1) ln T_n + C_K**, hence limsup E R_n/ln T_n ≤ (K + 1)/2.
Proof: the proof of Theorem R with the replacements of §1–2: per-time failure events 𝓑_t = {θ_t ∉ R(δ)} ∪ {λ_min(F_t) < tλ_L/(4Lc⁺)} ∪ {ℓ̄(t) > ℓ̄_t} ∪ {max_k S_k(t) > s_max(t)}, each with probability ≤ C₁e^{−c₁γδ²t} for t ≥ t\* (Theorem FK + Theorem Cexp with q = e^{−c′γt}); predictable good-gap indicators; one global telescoping; Proposition E5 on the good part; Cauchy–Schwarz on the bad part. ∎ (same level of detail as for K = 1).

**What this does not give.**
1. *Constants.* For β = (1, 2), μ = 1.5, ε = 0.4 (and 0.3) the explicit bound on λ_min(I) from Theorem `thm:trunc` is 5.6·10⁻⁴ (5.5·10⁻⁴) against a true λ_min of order 6·10⁻² (i.e. 100× conservative, as before), and the factor 1/(4L) with L ≈ 50 μ/ε makes the information level λ_L/(4L) ≈ 10⁻⁶. Since T₁ ∝ 1/(that level), C_K is far too large to be informative; the theorem asserts the leading constant, not a usable burn-in. For tight grids (K ≥ 3, ratio 1.5–2) λ_min(M) itself collapses geometrically in K.
2. *Slow transients.* The ridge direction (information ≈ ε/(μ(k − 1)) per unit time) and the middle tier (O(1)) enter λ_L only through λ_min(M); the three-tier structure is not exploited, so the proof treats all K + 1 directions as if they were equally slow.
3. *K-general Assumption-free regime.* The theorem needs no condition k_K > 1 (the chain never used it once M is replaced by λ_L).

## 4. Numerics (K = 2) **[N]**

`code/p5_k2_regret.py`: β = (1, 2), μ\* = 1.5, α\*_k = (1 − ε)β_k/2, box [0.5μ\*, 2μ\*] × [0.4α\*_k, 2α\*_k] ∩ {Σα_k/β_k ≤ 1 − ε/2}, start (1.3μ\*, 0.8α\*), F_0 = diag, projection (SLSQP QP) active in 0.2 % of steps, 6 replicates, n = 40 000 events.

| ε | E R_n at n/16, n/8, n/4, n/2, n | slope dE R/d ln T between checkpoints | dE R/d log det F | target |
|---|---|---|---|---|
| 0.4 | 4.41, 5.29, 6.14, 7.28, 7.99 | 1.32, 1.25, 1.65, 1.01 (mean 1.3) | 0.44, 0.42, 0.55, 0.34 | (K+1)/2 = 1.5; 0.5 |
| 0.3 | 4.85, 5.86, 6.60, 7.64, 8.28 | 1.52, 1.10, 1.50, 0.93 (mean 1.26) | 0.52, 0.37, 0.51, 0.31 | same |

The growth is logarithmic with a slope of the order (K + 1)/2 (mean 1.3 against 1.5, with ±0.3 scatter between checkpoints at 6 replicates; the slope with respect to log det F, which removes the T-scaling, is 0.43 against 0.5). This supports the order and does not distinguish 1.5 from 1.25; the slow approach is expected because the ridge and shape directions have not saturated at these T. No run at ε < 0.3 was done for K = 2 (the transient scales like 1/ε² events).

## 5. Status

| item | status |
|---|---|
| projection onto a polytope, regret reduction, ratio set | [P] |
| E5, bootstrap, Theorem R structure for K ≥ 2 | [P\*] (as K = 1) |
| information concentration (Theorem FK, replaces Assumption M) | [P] |
| leverage bound | [P], weaker than K = 1 by (ln t)²/ε |
| explicit constants | uninformative (10⁶ conservative) |
| K = 2 regret slope | [N], consistent with (K+1)/2 within scatter |
| three-tier-aware preconditioning, sharp transients | [O] |
| K ≥ 3 and ε < 0.3 numerics | not run |

---
## Update (25 Oct 2026): better information level
The information level of Theorem FK can be taken from the tangent minorant of `tangent_bound.md` (lambda_*/2 = 0.0126 for the K = 2 example, against about 1e-6 for the truncation device), and Assumption M is not needed. The expected-regret proof is written in full in `theoremR_proof_complete.md`.
