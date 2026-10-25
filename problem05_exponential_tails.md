# Problem 5: exponential tails for the bootstrap and for Theorem F (K = 1)

Interpretation of the request: "exponential tails (Freedman/Bernstein)" was the first open item in `problem05_bootstrap.md`. This note replaces the Doob/Chebyshev/Markov steps by exponential-supermartingale (Freedman/Bernstein) steps. Rigorous: Lemmas E1, E2, E3, E4. Sketch only: the control of G (E5). Not done: the expected-regret statement with the sharp constant (see the end).
Code: `code/p5_exp_tails_check.py`, `code/p5_ito_check.py` (reuse `p5_*` modules).

## E1. Exact identities (K = 1, γ = βε, n = 1 − ε, ES = μ/γ)

With N the counting process, Λ_T = ∫_0^T λ dt, M' = N − Λ (⟨M'⟩ = Λ):
(a) εΛ_T = μT + n M'_T − n(S_T − S_0);
(b) ∫_0^T (S_t − ES) dt = (S_0 − S_T)/γ + (N_T − Λ_T)/γ, which is the Poisson-equation identity for û(s) = (s − ES)/γ (L s = −γ(s − ES)).
Proof: dS = −βS dt + dN gives ∫S = (N_T − (S_T − S_0))/β, and Λ_T = μT + α∫S; substitute. ∎ Numerically both hold to 10⁻¹⁴ (relative) in all runs.

## E2. Exponential tails for the integrated intensity — rigorous

Lower tail. If Λ_T ≤ (1−η)μT/ε then by (a) nM'_T ≤ −ημT + n(S_T − S_0). On {|S_T − S_0| ≤ y} put z = (ημT − ny)/n > 0. For θ < 0, e^θ − 1 − θ ≤ θ²/2, so exp(θM'_t − (θ²/2)Λ_t) is a supermartingale and
P(M'_T ≤ −z, Λ_T ≤ V) ≤ exp(−z²/(2V)), V = μT/ε,
i.e. P(Λ_T ≤ (1−η)μT/ε, |ΔS| ≤ y) ≤ exp(−(ημT/n − y)² ε/(2μT)) ≈ exp(−η²μεT/(2n²)).
Upper tail (unit jumps, Bernstein): P(M'_T ≥ z, Λ_T ≤ V') ≤ exp(−z²/(2(V' + z/3))). Rate again of order μεT = (μ/β)γT.
Numerics (η = 0.5): empirical P(Λ_T ≤ 0.5μT/ε) = 0.007, 0, 0 (ε = 0.3) and 0.017, 0, 0 (ε = 0.1) at μεT = 5, 10, 20, against bounds 0.28, 0.08, 0.006 and 0.46, 0.21, 0.046 (valid, 50–1000× loose).

## E3. Exponential tails for occupation functionals (Theorem F with exponential tails) — rigorous modulo integrability

Let h ∈ [0,1] be Λ_h-Lipschitz, m = E_π h. Put û(s) = ∫_0^∞ (P_t h(s) − m) dt. By the monotone-coupling contraction (Problem 3, Lemma C; K = 1) |û(s) − û(s')| ≤ Λ_h|s − s'|/γ, |û(s)| ≤ Λ_h(s + ES)/γ, and Lû = −(h − m). Dynkin gives
∫_0^T (h(S_t) − m) dt = û(S_0) − û(S_T) + 𝓜_T, 𝓜_T = ∫_0^T [û(S_{t−}+1) − û(S_{t−})](dN − λ dt),
with jumps bounded by b_0 = Λ_h/γ and ⟨𝓜⟩_T ≤ b_0²Λ_T. Freedman (Bernstein for compensated point-process integrals with |H| ≤ b_0) and E2 then give, for x = Tm/4 and T so large that the boundary term |û(S_0) − û(S_T)| ≤ b_0(S_0 + S_T + 2ES) ≤ Tm/4 (which holds w.h.p. by the moment bounds of Problem 3 once Tγ ≳ 12κ/m),
P(∫_0^T h ≤ Tm/2) ≤ exp(−x²/(2(b_0²(1+η)μT/ε + b_0x/3))) + P(Λ_T > (1+η)μT/ε) + P(boundary term too large) ≤ 3e^{−cTγ}, c = m²μ / (32 κ²(1+η) β),
where κ := Λ_h·ES (for the trapezoids of `problem05_fisher_concentration.md`, ramp ρ = d/2 = 0.075·ES, so κ = 13.3 and b_0 = κ/μ, independent of ε).
Consequence: the high-probability statement of Theorem F holds with failure probability q for T ≥ T_0(q) = O(log(1/q)/γ) instead of O(1/(γq)).
Numerics: empirical P(∫h ≤ Tm/2) for the bump on I_2 (ε = 0.1, m = 0.047, 300 reps): 0.28, 0.24, 0.14, 0.067 at Tγ = 5, 10, 20, 40 — decay of roughly 0.04–0.056 in log P per unit Tγ. The rigorous rate is c = 3×10⁻⁷ (m = 0.05, κ = 13.3, η = 0.5, μ = β = 1), i.e. 10⁵ smaller: the bound is exponential but quantitatively uninformative (for q = 10⁻³ it requires Tγ ≈ 2×10⁷ against ≈ 170 observed).

## E4. Exponential tail for the Newton martingale — rigorous

In Lemma C1 of `problem05_bootstrap.md`, M_t = ∫_0^t w(dN − λ* du) with w = −2φᵀv/λ_θ, |w| ≤ b_M := 2κ_1, κ_1 = max(1, 1/c⁻ − 1), and ⟨M⟩_t ≤ (4/c⁻)X_t, where X_t is the realized, adapted quantity ∫_0^t (φᵀv)²/λ_θ du (it equals the sum of the vᵀΔFv over completed gaps; no predictable-versus-realized mismatch).
For 0 < θ ≤ 1/b_M, e^{−θw} − 1 + θw ≤ (e−2)θ²w², so Z_t = exp(−θM_t − (e−2)θ²⟨M⟩_t) is a supermartingale; Ville's inequality gives P(∃t: Z_t ≥ e^x) ≤ e^{−x}. If moreover θ ≤ c⁻/(4(e−2)) = c⁻/2.87, then (e−2)θ²⟨M⟩_t ≤ θ X_t, so with probability ≥ 1 − e^{−x}, simultaneously for all t,
−M_t ≤ X_t + x/θ.
Together with Lemma C1 this cancels the −X_t term exactly: with probability ≥ 1 − e^{−x},
sup_{j ≤ n} D_j ≤ D'_0 + x/θ_M + G_n, θ_M = min(1/(2κ_1), c₀⁻/2.87).
(For the box value c₀⁻ = 0.5, κ_1 = 1 and θ_M = 0.174, i.e. 1/θ_M = 5.75; the observed sup|M_j| had median 6–9 and 95 % quantile 12–18 over the whole run, consistent with a deviation of this size.)

## E5. Control of G_n with exponential tails — sketch

G_n ≤ (1+η)S_n + (1+1/η)H_n with S = Σ sᵀF_{i−1}⁻¹s_i, H = Σ hᵀF_{i−1}⁻¹h_i.
- H_n ≤ (ℓ̄/μ_lo)·Σ(aᵀv)²/(λ*λ_θ) = (ℓ̄/μ_lo)(X_n + ½N^A_n) where N^A is a martingale with integrand 2(φᵀv)²/(λ*λ_θ), bounded by 2κ_1²c⁺ and with ⟨N^A⟩ ≤ 4κ_1²c⁺X; E4's argument applies, so H_n ≤ (ℓ̄/μ_lo)(1 + O(1))X_n + O(x) w.p. 1 − e^{−x}.
- S_n: Itô's formula for the stopped vector martingale M^s_u = −∫_0^uφ + Σ jumps gives s_iᵀF⁻¹s_i = a_iᵀF⁻¹a_i/λ*(τ_i)² + (martingale increment). The first term is a sum over events of the predictable function f = φᵀF⁻¹φ/λ*², bounded by ℓ̄/μ, so by Freedman Σ_i f(t_i) ≤ ∫f λ* dt + √(2Vx) + (ℓ̄/μ)x with V ≤ (ℓ̄/μ)∫fλ*dt, and ∫fλ* dt = ∫φᵀF⁻¹φ/λ* dt ≤ c⁺ r(x̄) log det(F_n/F_{i_0}) ≈ d log T. The martingale part has increments controlled by the leverage as well, so it contributes O(ℓ̄√x…).
Expectation of the identity checked: E Σ sᵀF⁻¹s − E Σ aᵀF⁻¹a/λ*² = −0.27 ± 0.16, −0.15 ± 0.14 (ρ_0 = 1) and −0.06 ± 0.05 (ρ_0 = 20), i.e. zero within 1.7σ (an earlier 30-rep run showed a 9 % gap that was heavy-tail noise from the first few high-leverage events).
Status: the algebra is standard but I have not written the bound on the martingale part of S_n with explicit constants; treat the resulting statement G_n ≤ (1+o(1))·d log T + O(x) w.p. 1 − O(e^{−x}) as unproved in detail.

## Theorem C-exp (K = 1; exponential-tail bootstrap)

Under the standing assumptions of `problem05_bootstrap.md` (compact box Θ_c, exact F-norm projection, Assumption M) and with the exponential Theorem F from E3, the iterates satisfy θ_{j+1} ∈ R(δ) for every j with t_j ≥ T_1, with probability ≥ 1 − q, provided
λ_{T_1} δ² ρ_*² > 3·(E D'_0-type bound x_0 + x/θ_M + g_n), with x = ln(c/q), g_n = (1+o(1))·d·log T_n (E5),
i.e. T_1 γ ≳ [x_0 + 5.75·ln(c/q) + d·log T_n] · 4c₀⁺L/(mλ_Q ρ_*² δ²)·(1/γ)·… = O((ln(1/q) + d log T_n)/δ²) in units of 1/γ. The dependence on q is now logarithmic instead of 1/q. The constants (c₀⁺L/(mλ_Q), and the rate c of E3) are as conservative as before.

## What is still not proved

1. **E5 in detail** (the exponential bound for G_n, in particular the martingale part of S_n).
2. **Expected regret with the sharp constant.** Exponential tails make P(bad event at time t) ≲ e^{−ct γ}, but that alone does not give E R_n ≤ (d/2 + o(1))log T_n + O(1): a union bound would need the burn-in region [0, T_c] with T_c ≈ (1/c) log(μT_nψ_max/ε) to be paid at the crude linear rate ψ_max μ/ε, costing O((μ/(βε²)) log T_n) with a large coefficient in front of log T_n. The remedy is an episode/restart argument: split time into alternating periods where the information event holds (tight bound from the restart, costs telescoping log det increments plus a restart cost D) and where it fails (crude rate × duration), and use that failure episodes have exponentially distributed length and expected total duration O(1/(cγ)) independent of T_n. Not carried out.
3. The rigorous constants (rate c of E3) are about 10⁵ times worse than the empirical decay rate; no attempt was made to tighten them (e.g. by bounding the Lipschitz constant of û using the actual state density rather than the crude contraction, or by choosing wider bump intervals).
4. Assumption M, a priori box, K ≥ 2 are unchanged.

Re-derive E2, E3 and E4 yourself: each is a short supermartingale argument that I wrote in one pass; the numerics confirm the exact identities and that the bounds are valid (not violated) but they do not check the constants.

---
## Update (25 Oct 2026)
E5 is now proved (and the expected-regret statement derived) in `problem05_E5_expected_regret.md`; consolidated status in `K1_main_results.md`. The "What is still not proved" items 1 and 2 above are closed there, modulo the inputs listed in its §5. Note also the correction to Lemma 7c / Theorem 7 in its §0.
