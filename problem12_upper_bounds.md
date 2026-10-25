# Problem 12: rigorous upper bounds for the middle tier and convergence of the smallest eigenvector (mathematics)

Status of this file: mathematics (statements and proofs). The symbolic and numerical checks are in the next section of the same file once the code has been run.

## 0. The linear structure behind the three tiers

Let A_mat = diag(β) − 𝟙αᵀ (K × K), the matrix of the mean dynamics: E[S_t | S_0] − m = e^{−A_mat t}(S_0 − m). Because α_k > 0 and β_1 < … < β_K are distinct, det(A_mat − aI) = ∏(β_k − a)·(1 − Σ_kα_k/(β_k − a)), so the eigenvalues are the roots of the secular equation Σ_kα_k/(β_k − a) = 1. The left side increases from −∞ to +∞ on each interval (β_{j−1}, β_j), and on (−∞, β_1) it increases from 0 to +∞ with value n < 1 at a = 0. Hence

 **A_mat has K real, simple eigenvalues 0 < γ = a_1 < β_1 < a_2 < β_2 < … < a_K < β_K**, so a_j ∈ (β_{j−1}, β_j) ⊂ [b, B] for j ≥ 2 **independently of ε**, while a_1 = γ = ε/ν̃ → 0.

The left eigenvector for a_j is ℓ_{j,k} = α_k/(β_k − a_j) (normalised so that σ_j := ℓ_jᵀ𝟙 = Σ_kα_k/(β_k − a_j) = 1, which is the secular equation). Define the modes
 y_j := ℓ_jᵀS = Σ_kα_kS_k/(β_k − a_j), j = 1..K, with y_1 = Z (the Perron variable of Problem 11).

**Lemma 0 (exact generator of the modes).** For every j, with the generator 𝒜 of S: 𝒜y_j = −a_jy_j + μ. Moreover each event increases every y_j by exactly 1, and 𝒜(y_j²) = −2a_jy_j² + 2μy_j + λ.
*Proof.* The continuous drift of y_j is −ℓ_jᵀBS, B = diag β; the jumps add σ_j = 1 at rate λ = μ + αᵀS. Since ℓ_jᵀB = ℓ_jᵀ(A_mat + 𝟙αᵀ) = a_jℓ_jᵀ + (ℓ_jᵀ𝟙)αᵀ, the drift of y_j is −a_jy_j − αᵀS and the jump term adds λ = μ + αᵀS, giving −a_jy_j + μ. The second formula: 𝒜(y²) = 2y(−a_jy + μ) + λ((y+1)² − y² − 2y) = −2ay² + 2μy + λ. ∎
(Consequence: E_πy_j = μ/a_j; the slow mode has E Z = μ/γ ~ 1/ε, the fast modes have O(1) means. The fast modes are Ornstein–Uhlenbeck-type with no coupling to the slow level in the *mean* drift; the level enters only through the jump intensity λ.)

The matrix Lₓ with rows ℓ_j is invertible (distinct eigenvalues), S = Lₓ⁻¹y, so each S_k is a combination of Z and the fast modes.

## 1. A Lyapunov function for the fast modes

Let W := Z + z with a parameter z > 0 and set
 c_λ := min(μ/z, b − γ), B_z := max(μ/z, B − γ), c₁ := max(μ, B_z).
From λ = μ + Σ(β_k − γ)ω_kS_k·… (Lemma G of Problem 11): μ + (b−γ)Z ≤ λ ≤ μ + (B−γ)Z, hence **c_λ W ≤ λ ≤ B_z W**.

**Lemma 1 (generator of y_j²/W).** For j ≥ 2 and f_j := y_j²/W,
 𝒜f_j = −2a_j y_j²/W + y_j²[(γZ − μ)/W² + λ/(W²(W+1))] + λ/(W+1) + 2y_j[μ/W − λ/(W(W+1))].
*Proof.* Products of functions of the state: 𝒜(gh) = (𝒜g)h + g(𝒜h) + λ(Δg)(Δh) with Δ the one-jump increment. Take g = y_j², h = 1/W, so that 𝒜g = −2a_jy² + 2μy + λ (Lemma 0), 𝒜h = (λ − μ + γZ)/W² − λ/(W(W+1)) [the drift of W is −(λ − μ + γZ) and W jumps by 1], and Δg = 2y + 1, Δh = −1/(W(W+1)). Collecting the terms gives the displayed formula. ∎ (Verified symbolically in the next section.)

**Lemma 2 (bound).** Assume z ≥ 1 and γ ≤ a_j/4 and B_z/z ≤ a_j/4. Then
 𝒜f_j ≤ −a_j f_j + B_z + 2c₁²/(a_jz).
*Proof.* (i) (γZ − μ)/W² ≤ γ/W (since Z ≤ W); λ/(W²(W+1)) ≤ B_z/(W(W+1)) ≤ B_z/(zW). So the bracket multiplying y² is ≤ (γ + B_z/z)/W ≤ (a_j/2)/W. (ii) λ/(W+1) ≤ B_z. (iii) 0 ≤ μ/W − λ/(W(W+1)) is not guaranteed, but |μ/W − λ/(W(W+1))| ≤ max(μ/W, B_z/(W+1)) ≤ c₁/W, so |2y_j[…]| ≤ 2c₁|y_j|/W ≤ (a_j/2)y_j²/W + (2c₁²/a_j)/W ≤ (a_j/2)y_j²/W + 2c₁²/(a_jz). Sum: 𝒜f_j ≤ (−2a_j + a_j/2 + a_j/2)y²/W + B_z + 2c₁²/(a_jz). ∎

**Proposition 3 (uniform bound).** Under the assumptions of Lemma 2, E_π[y_j²/W] ≤ [B_z + 2c₁²/(a_jz)]/a_j and therefore
 **E_π[y_j²/λ] ≤ C_{Y,j} := [B_z + 2c₁²/(a_jz)] / (a_jc_λ), independent of ε** (for ε ≤ ε_0 such that γ ≤ min_ja_j/4).
*Proof.* Stationarity E_π𝒜f_j = 0 (justified by localisation: f_j is integrable and 𝒜f_j is bounded by polynomial moments, which are finite by Lemma `lem:moments`; the stopped martingale argument is that of Problem 11, "Justification of E_π𝒜f_z = 0"), combined with Lemma 2, gives a_jE f_j ≤ B_z + 2c₁²/(a_jz). Then λ ≥ c_λW. ∎

## 2. Consequences for the spectrum of I

Let Lₘ be the K × (K−1) matrix with columns ℓ_2,…,ℓ_K and s_min its smallest singular value (explicit, depends only on α/β; it stays bounded below as ε → 0 because ℓ_{j,k} → α_k/(β_k − a_j) with a_j → the K−1 roots at ε = 0).
**Theorem 4 (middle tier, upper bound).** For v = (0, w), w = Σ_{j≥2}c_jℓ_j,
 vᵀIv = E[(Σc_jy_j)²/λ] ≤ (K−1)Σc_j²E[y_j²/λ] ≤ (K−1)C_Y|c|² ≤ ((K−1)C_Y/s_min²)|v|², C_Y := max_{j≥2}C_{Y,j}.
By Courant–Fischer the (K−1)-th smallest eigenvalue satisfies **λ_{K−1}(I) ≤ C_mid := (K−1)C_Y/s_min²**, independent of ε. If the ridge Rayleigh quotient ρ_ε := e_rᵀIe_r (e_r the unit vector ∝ (1, −εα/μ)) is added: on span{e_r} ⊕ {(0,w)}, since the angle between e_r and the middle subspace is ≥ arccos(c_ε) with c_ε = O(ε), vᵀIv ≤ 2(s²ρ_ε + wᵀIw) ≤ 2max(ρ_ε, C_mid)(s²+|w|²) ≤ [2max(ρ_ε,C_mid)/(1 − c_ε)]|v|², hence λ_K(I) ≤ 2max(ρ_ε, C_mid)/(1−c_ε).
**Corollary 5 (top eigenvalue).** E[Z²/λ] ≤ E Z/c_λ = μ/(γc_λ) = O(1/ε) and E[y_j²/λ] ≤ C_{Y,j}; since S = Lₓ⁻¹y, tr I = E[1/λ] + ΣE[S_k²/λ] ≤ 1/μ + K‖Lₓ⁻¹‖²_F(μ/(γc_λ) + (K−1)C_Y) = O(1/ε), so **λ_{K+1}(I) ≤ tr I = O(1/ε)**; together with the tangent lower bound ~1/ε: λ_max(I) = Θ(1/ε).
Combined with the lower bounds (tangent bound): **top Θ(1/ε); K−1 middle eigenvalues in [λ_j(N), C_mid] (both ε-independent); the ridge ≥ 0.25ε-ish from below and ≤ ρ_ε = O(ε) from above** — the three-tier structure of Problem 12 is a theorem for k_K > 1 (the ridge upper bound uses E[1/λ] = O(ε), Theorem E).

## 3. Convergence of the smallest eigenvector to the ridge direction
**Theorem 6.** Let k_K > 1, let v_1 be the unit eigenvector of the smallest eigenvalue λ_1 of I, and suppose λ_2(I) ≥ g > 0 (for example g = λ_2(N(λ₀)) from the tangent bound, a computable ε-dependent number, or a fixed lower bound uniform in ε ≤ ε_0). Then
 **sin²∠(v_1, e_r) ≤ (ρ_ε − λ_1)/(g − λ_1) ≤ ρ_ε/g = O(ε/g).**
*Proof.* Expand e_r in the eigenbasis, e_r = Σc_jv_j: ρ_ε = Σλ_jc_j² ≥ λ_1c_1² + g(1 − c_1²), so 1 − c_1² ≤ (ρ_ε − λ_1)/(g − λ_1), and x ↦ (ρ_ε − x)/(g − x) is decreasing for ρ_ε < g, so the maximum over x ≥ 0 is at x = 0. ∎ (This is the Davis–Kahan sin θ bound for a rank-one subspace.) For k_K ≤ 1 only ρ_ε = O(E[1/λ]) is available, and E[1/λ] is bounded below, not above, by Theorem F, so the statement is not claimed there.

## 4. What the proof needs and what is left
- Needs: the eigenstructure of A_mat (proved above), Lemma 0–2 (algebra, checked symbolically below), E_π𝒜f_j = 0 (localisation), c_λ>0 (z fixed), γ ≤ a_j/4, B_z/z ≤ a_j/4 (choose z ≥ 4B/b and ε small).
- The bound on ρ_ε uses Theorem E (explicit constants, loose by 18×).
- g is computer-assisted for each ε (tangent bound); "g bounded away from 0 uniformly in ε → 0" is verified numerically (ε = 0.1, 0.03, 0.01) but not proved analytically.
- Not claimed: upper bounds for k_K ≤ 1; the exact value of the middle eigenvalues (the spectral integral); the case of coinciding β_k (excluded).

## 5. Symbolic and numerical checks (code/p12_upper_check.py, code/p12_eigvec_check.py)
- **Lemma 0 and Lemma 1 verified symbolically** (sympy, K = 2 with symbolic β₁, β₂, γ, a₂; α solved from the two secular equations): 𝒜y₁ + γy₁ − μ = 0, 𝒜y₂ + a₂y₂ − μ = 0, both modes jump by exactly 1, and 𝒜(y₂²/W) equals the formula of Lemma 1 (residual exactly 0 at 12 random rational points).
- **Lemma 2 pathwise inequality:** 5 721 random (state, mode) pairs (K = 2–4), 0 violations.
- **Constants (K = 2, β = (1,2), μ = 1.5):** a₂ = 1.5183, 1.5051, 1.5017 and γ = 0.1317, 0.0399, 0.0133 at ε = 0.1, 0.03, 0.01. Rigorous C_Y = 6.06, 6.78, 7.00 (Monte Carlo E[y₂²/λ] = 0.305, 0.320, 0.333; the bound is ~20× loose but **ε-independent**, as the theorem says); C_mid = 1.427, 1.424, 1.423 (constant in ε) against the Monte Carlo middle eigenvalue 0.064–0.066 (22× loose; the ratio is |ℓ₂|² ≈ 4.2 times the looseness of C_Y). The theorem therefore proves that the middle eigenvalue is O(1) uniformly in ε, with explicit constants, in the regime γ ≤ min a_j/4.
- **Eigenvector (Theorem 6), K = 2:** tangent gap g = λ₂(N) = 0.0423, 0.0445, 0.0456, 0.0460, 0.0461 at ε = 0.1, 0.03, 0.01, 0.001, 0.0003 (stable ⇒ uniform gap numerically). With the *rigorous* Theorem-E bound on E[1/λ] (which is ~18× loose) the bound sin² ≤ ρ_ε/g is 0.49 at ε = 10⁻³ and 0.146 at 3·10⁻⁴, decreasing linearly in ε as the theorem says; it is vacuous at ε ≥ 0.003 because of that looseness (measured sin² = 0.28, 0.020, 0.0016 at ε = 0.1, 0.03, 0.01). So the theorem is asymptotically correct and numerically informative only for ε ≲ 10⁻³ unless the E[1/λ] bound is sharpened.

## 6. Status after the check
| item | status |
|---|---|
| E[y_j²/λ] = O(1), middle-tier eigenvalues bounded above uniformly in ε (k_K arbitrary for the fast modes; needs γ ≤ min a_j/4) | **proved**, symbolic + pathwise checks |
| top eigenvalue Θ(1/ε) (upper bound tr I = O(1/ε); lower bound from the tangent bound) | **proved** |
| three-tier structure, k_K > 1 (ridge ≤ ρ_ε = O(ε) by Theorem E, ridge ≥ tangent bound) | **proved** (constants loose) |
| smallest eigenvector → ridge direction, rate O(ε/g) | **proved** given g ≥ g_0 > 0 (computer-assisted for each ε; uniform in ε numerically); informative for ε ≲ 10⁻³ with the rigorous E[1/λ] bound |
| k_K ≤ 1 | upper bounds on E[1/λ] unavailable (Theorem F is a lower bound); middle tier bounded as above, eigenvector statement not claimed |
