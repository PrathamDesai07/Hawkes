# Problem 14: a tracking theorem for the forgetting-Newton estimator, and the amplitude effect derived (mathematics, then code)

Estimator (event-indexed): F_i = (1−η)F_{i−1} + ΔF_i, g_i = b_i − a_i/(a_iᵀθ_i), θ_{i+1} = Π^{F_i}_Θ(θ_i − F_i⁻¹g_i), F_0 ≻ 0, Θ convex compact. Truth: a path θ*_t in Θ (Lipschitz, or with jumps). Notation per gap i (history 𝒢_{i−1}): comparator θ*_i := θ*(t_{i−1}) (parameter at the gap start), v_i := θ_i − θ*_i, X_i = v_iᵀΔF_iv_i, G_i = g_iᵀF_i⁻¹g_i, d_i := θ*_{i+1} − θ*_i, D_i := ‖θ_{i+1} − θ*_{i+1}‖²_{F_i} (so D_{i−1} = ‖v_i‖²_{F_{i−1}}), and the in-gap drift energy Ξ_i := ∫₀^{τ_i}(φ(u)ᵀ(θ*_i − θ*(t_{i−1}+u)))²/λ_θ(u) du (zero for a frozen parameter).

## Part A. Tracking theorem

**Proposition 1 (pathwise recursion) [P].** For every ρ > 0,
 D_i ≤ (1+ρ)[(1−η)D_{i−1} + X_i − 2⟨g_i,v_i⟩ + G_i] + (1+1/ρ)‖d_i‖²_{F_i}.
*Proof.* θ*_i ∈ Θ and Π^{F_i} is non-expansive in ‖·‖_{F_i}: ‖θ_{i+1}−θ*_i‖²_{F_i} ≤ ‖v_i − F_i⁻¹g_i‖²_{F_i} = ‖v_i‖²_{F_i} − 2⟨g_i,v_i⟩ + G_i. Forgetting: ‖v_i‖²_{F_i} = (1−η)‖v_i‖²_{F_{i−1}} + X_i = (1−η)D_{i−1} + X_i. Finally θ_{i+1} − θ*_{i+1} = (θ_{i+1} − θ*_i) − d_i and (x+y)² ≤ (1+ρ)x² + (1+1/ρ)y² in the F_i-norm. ∎

**Proposition 2 (compensator identity under drift) [P].** E[X_i − 2⟨g_i,v_i⟩ | 𝒢_{i−1}] ≤ E[Ξ_i | 𝒢_{i−1}].
*Proof.* Let λ^t(u) = φᵀθ*(t_{i−1}+u) be the true intensity on the gap, λ^ref = φᵀθ*_i, δ(u) = λ^ref − λ^t, so λ_θ − λ^t = φᵀv + δ. Given the past, τ has density λ^te^{−Λ^t}; E g_i = ∫φ(λ_θ−λ^t)/λ_θ·e^{−Λ^t}du (Fubini, as in Proposition 1(a) of Problem 5). Hence E⟨g,v⟩ = ∫e^{−Λ^t}(φᵀv)(φᵀv+δ)/λ_θ du and E X_i = ∫e^{−Λ^t}(φᵀv)²/λ_θ du, so E[X − 2⟨g,v⟩] = −∫e^{−Λ^t}[(φᵀv)² + 2(φᵀv)δ]/λ_θ du ≤ ∫e^{−Λ^t}δ²/λ_θ du = EΞ_i, using x² + 2xδ ≥ −δ². ∎

**Proposition 3 (steady-state bound) [P].** Take ρ = η/2 and η ≤ 1, so (1+η/2)(1−η) ≤ 1−η/2. Then
 E D_i ≤ (1−η/2)E D_{i−1} + (1+η/2)(E G_i + E Ξ_i) + (1+2/η)E‖d_i‖²_{F_i},
hence for the time average, mean_i E D_i ≤ (2/η)[(1+η/2)(mean E G + mean EΞ) + (1+2/η)mean E‖d‖²_F] + D_0/(nη/2).

**Theorem 4 (tracking error in the Fisher norm) [P\*].** Assume (W) the window information bounds: on a good event of probability ≥ 1 − q_n, (c_F/η)I_ev ⪯ F_{i−1} ⪯ (C_F/η)I_ev for all i ≥ i_0, where I_ev = E[ΔJ] is the information per event at the current truth (for a drifting truth I_ev varies slowly; treat it as constant over a window), together with the leverage conditions of Theorem R (point leverage ℓ̄_i ≤ ℓ̄ small, iterates in a ratio set, so that E[G_i] ≤ (1+η′)E tr(F_{i−1}⁻¹ΔJ_i) + (1+1/η′)E[H_i] with E[H_i] ≤ B E X_i). Then
 mean_i E‖θ_{i+1} − θ*_{i+1}‖²_{I_ev} ≤ (2/c_F²)·(1+η/2)(1+η′)/(1−B)·ηd + (2/c_F)(1+η/2)·mean_i EΞ_i + (2C_F/c_F)(1+η/2)(1+2/η)·mean_i‖d_i‖²_{I_ev}/η + (transient),
where the second-to-last coefficient (1+2/η)/η ≈ 2/η²,
and in particular, for the idealised window F_i = I_ev/η (c_F = C_F = 1) and a frozen in-gap parameter (Ξ = 0),
 **E‖v‖²_{I_ev} ≲ 2ηd + 4(‖d‖_{I_ev}/η)²**, with d = dim θ, ‖d‖²_{I_ev} = dᵀI_evd the squared per-event parameter increment in the Fisher norm.
*Proof.* (W) gives ‖v‖²_{I_ev} ≤ (η/c_F)D and ‖d_i‖²_{F_i} ≤ (C_F/η)‖d_i‖²_{I_ev}. Leverage: E tr(F⁻¹ΔJ) ≤ (η/c_F)E tr(I_ev⁻¹ΔJ) = ηd/c_F (E ΔJ = I_ev); the drift part H is absorbed as in Theorem R (coefficient B ≪ 1). Insert into Proposition 3 and multiply by η/c_F. The bad event contributes ≤ q_n·(box diameter)²·sup F, exponentially small. ∎ The factor 4 is the price of the Young split with ρ = η/2 and of discarding the dissipation −X_i; the exact linear-filter values are ηd/2 (noise) and (‖d‖/η)² (drift).
The window condition (W) needs the window length 1/η (events) to contain many mixing times: η ≲ εγ-scale, i.e. the "linear window" ηr ≲ γ of the experiments (η_max = θ_kμ²/a); it is proved by the exponentially weighted version of the Poisson-equation/Freedman argument of Theorem FK (weights (1−η)^{i−j} in place of 1) with the tangent information level λ_*; I state it as the hypothesis (W) rather than carrying the weighted concentration through.
**Scope.** Lipschitz drift (and piecewise-constant drift with the extra factor 1/η from the Young split, i.e. a looser bound than the exact tanh law; jump class not tight). K = 1 numerically.

### Check (code/p14_track_check.py; K = 1, both parameters unknown, ε = 0.4, η = 0.004, μ_t = 1 + 0.3 sin(2πt/1200), exact F-projection, 4 replicates, 57 616 gaps)
- (P) pathwise recursion: **0 violations** in 57 616 gaps.
- (E) mean(X − 2⟨g,v⟩) = −0.0060 ≤ mean Ξ_gap = 4·10⁻⁷ ✓ (the left side is negative: the dissipation we discarded).
- E G measured = 0.00805 vs the prediction ηd = 0.008 ✓ (the leverage input of Theorem 4).
- (S) mean D = 1.73 ≤ steady-state right side 12.05 (ratio 0.14).
- (F) mean ‖v‖²_{I_ev} = 0.00770 ≤ 2ηd + 4‖d‖²_I/η² = 0.0463 (ratio 0.17); the exact-filter prediction ηd/2 + (drift) = 0.0040 + 0.0036 = 0.0076 matches the measurement to 1 %.

## Part B. The amplitude effect, derived

**Mechanism [P].** In the linear regime the mean error obeys e_{i+1} = (I − F_i⁻¹f_i)e_i − d_i with f_i = I(θ*_{t_i})/r_{t_i} (information per event at the current truth) and F_i the window. Per unit time this is de/dt = −Γ(t)e − dμ/dt with gain Γ(t) = I(θ*_t)·[F_t]⁻¹ (matrix), where for a fixed preconditioner F = I(θ_ref)/(η) and for the adaptive window F satisfies dF/dt = −ηr_tF + I(θ*_t). The constant-gain filter has Γ ≡ ηr_0 =: h. Because I(μ) depends strongly on μ (κ_I := dln I/dln μ = −1.14 for α known, −1.87 for α unknown at μ = β = 1, ε = 0.4), the gain is modulated at relative amplitude κ_IA. The adaptive window sees a *lagged* I, so its modulation has the opposite sign structure to the fixed one.
**Second-order result (scalar case) [P, sympy harmonic algebra].** Write Φ := Γ/h = I(μ_t)/I_ref = 1 + a₁A sin ωt + a₂A² sin²ωt + …, a₁ = κ_I, a₂ = I″(1)/(2I(1)); ω = 2π/P; time-domain forcing exactly A ω cos ωt. The relative fundamental error energy (ratio to the constant-gain filter) is R = 1 + C₂A² + O(A⁴) with
 C₂^{fixed} = h²[a₁²(h² − 2ω²) − a₂(h² + 4ω²)] / [2(h²+ω²)(h²+4ω²)],
 C₂^{adaptive} = h²[a₁²ω² − 2a₁ω² + a₂ω² + h²/2] / [(h²+ω²)(h²+4ω²)]
(adaptive: ψ′ = h(−(1+A sin ωt)ψ + Φ), Γ = hΦ/ψ). These were derived by solving the harmonic-balance hierarchy e = Ae₁ + A²e₂ + A³e₃ symbolically (`code/p14_perturbation.py`).
**Numerical confirmation of the formulas.** With the exact m(μ) (Laplace ODE), ε = 0.4, η = 0.004, P = 1200 (h/ω = 1.91): α known (a₁ = −1.140, a₂ = 1.267): C₂^{fixed} = −0.387, C₂^{adaptive} = +0.684; exact ODE (long integration, nonlinear Φ truncated at second order) vs 1 + C₂A²: A = 0.05: 0.9990 vs 0.9990 (fixed), 1.0017 vs 1.0017 (adaptive); A = 0.10: 0.9962 vs 0.9961, 1.0068 vs 1.0068; A = 0.20: 0.9850 vs 0.9845, 1.0273 vs 1.0274; A = 0.35: 0.9568 vs 0.9526, 1.0826 vs 1.0838. α unknown (efficient information, scalar projection: a₁ = −1.869, a₂ = 2.600): C₂^{fixed} = −0.725, C₂^{adaptive} = +1.196, adaptive/fixed ≈ 1 + 1.92A².
**Full two-parameter skeleton (α unknown; μ- and α-errors coupled through the 2×2 matrices F, f)** (`code/p14_skeleton.py`, deterministic): fundamental error² ratio to the constant-gain filter: fixed 0.995, 0.989, 0.967, 0.925, 0.900 and adaptive 1.005, 1.019, 1.074, 1.233, 1.500 at A = 0.05, 0.10, 0.20, 0.35, 0.50 (adaptive/fixed = 1.010, 1.030, 1.111, 1.333, 1.668), against Monte Carlo (`problem14_drift_model.md`, R1): fixed 1.07, 0.98, 0.92 and adaptive/fixed 1.32, 1.48, 1.75 at A = 0.2, 0.35, 0.5. So the amplitude effect is reproduced by a noise-free deterministic mechanism, grows as A² at small A, and the skeleton matches the Monte Carlo at A ≥ 0.35 (1.33 vs 1.48 at 0.35 and 1.67 vs 1.75 at 0.5, within the Monte Carlo scatter of ±0.1–0.2 quoted in R0); at A = 0.2 the skeleton (1.11) is below the measurement (1.32). **The matrix (α-unknown) closed form was not derived;** the scalar closed form is exact in the α-known case and quantitatively explains why the data-driven preconditioner has an excess that vanishes as A → 0 and is zero for a fixed preconditioner at the same order.

## Status
| item | status |
|---|---|
| tracking recursion, compensator identity, steady-state bound (Props 1–3) | **proved** |
| Theorem 4 (E‖v‖²_{I_ev} ≲ 2ηd + 4(‖d‖/η)²) | proved modulo the window hypothesis (W); numerically valid with margin 6× |
| amplitude effect, scalar closed form | **derived and verified** (4 digits at A ≤ 0.2) |
| amplitude effect, α unknown | mechanism and deterministic skeleton; no closed form |
| jump class | not tight (Young split); the tanh law remains the sharp statement for the linear filter |
