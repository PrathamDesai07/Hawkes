# Problem 13, third pass: the four caveats of the second pass

**Status: draft.** Builds on `problem13_second_pass.md` (Theorem TR2). Each section says what is proved (by an argument you should re-derive) and what is only checked numerically or remains conditional. Code: `code/p13c_*.py`, `code/p13d_*.py`.

## 1. Projection

**Claim.** Let Θ_c be convex and Π^F the exact F-norm projection onto it (Problem 4, `project_dense`). If θ* ∈ Θ_c then ‖Π^F(x) − θ*‖_F ≤ ‖x − θ*‖_F (non-expansive in the F-norm, standard Hilbert-space projection). With F_i fixed at step i, the projected step therefore satisfies ‖v_{i+1}‖_{F_i} ≤ ‖unprojected v_{i+1}‖_{F_i}, so the one-step bounds of Lemmas K and N hold unchanged for the potential measured in F_i.
**Two caveats handled.** (a) The potential is U_i = ‖v_i‖²_{F_{i−1}}, and F changes between steps; the change of metric is absorbed by the existing log-det/forgetting argument. (b) For drift, the target θ*_i moves; if θ*_i ∈ Θ_c for all i, the projection is applied to the moved target, and a jump of size D_i enters as an additive term via the dynamic Young inequality, giving V_{i+1} ≤ √((1−η)V_i² + ν_i) + δ_i with ν_i the noise increment and δ_i the drift displacement in the F-norm.
**Diagonal clip.** A coordinate clip is *not* an F-norm projection and is not covered. The earlier K = 2, η = 0.002 failure came from ~1% clipped steps; with the exact projection the bound holds (`p13c_projection_check.py`, `p13c_engine.py`).
**Checks.** Non-expansiveness holds on random draws, Lemma N and the combined bound hold with exact projection. Projected jump-rate slope: −0.47 with large noise, so the jump-rate scaling under projection is *not* sharply confirmed numerically.

## 2. High-probability form of A1′ and A2

Replace the almost-sure leverage bound by the predictable leverage ℓ_i^pred = tr(F_i⁻¹ E[ΔF_i | past]) (F_{i−1}-measurable). Let τ = first i with ℓ_i^pred > ℓ̄ (a predictable stopping time). On {i ≤ τ} the contraction and noise lemmas hold with ℓ̄; stop the recursion at τ, apply the martingale (Freedman) inequality to the stopped process, and add P(τ ≤ n). Result: with probability ≥ 1 − P(τ ≤ n) − δ the bound holds with a log(1/δ) factor.
P(τ ≤ n) is controlled by A4 (window information), which gives ℓ_i^pred ≲ η/(ε²c_Fc⁻(Θ−1)) once the window is full. Numerically (`p13c_hp.py`): with ℓ̄ = 0.5, P(τ ≤ n) = 0 in the window regime and 1 outside it, consistent with the claim. **Not done:** A2's high-probability form still needs an overshoot argument for the stopping time (the process can cross the threshold in one step); I did not close it.

## 3. Step-3 constant

The concentration step of the window argument: truncate g_z = 1/max(λ, z), z = c₁/ε, solve the Poisson equation for the truncated function, and apply Freedman using exact W_ℓ contraction at rate γ. The line-by-line constant C_* is now written out from the truncated Poisson solution, the Freedman variance bound and the small-ball (Gamma majorant, K = 1) estimate. It is valid but astronomically conservative: ≈ 10¹⁰ times the observed constant at k = 2 (real C ≈ 64). Truncation bounds were checked numerically (`p13c_truncation.py`). Treat C_* as a proof of existence of a constant, not as a usable one.

## 4. A4 for K ≥ 2

**Lemma IM (inverse moments, any K).** With the corrected Lyapunov function W^{−a} + ψW^{−a−2}, E[λ^{−a}] ≲ ε^a for a in a range, from a pathwise generator inequality. Checked pathwise (`p13d_lemmaIM.py`) and by MC: K = 2 (β = (1,2), μ = 1.5, ε = 0.2) E[1/λ] = 0.895 ε, E[λ^{−1.2}] = 0.907 ε^{1.2}.
**Reduction.** In the correlation basis, condition on the normalized features, apply Weyl's inequality plus an ε-net to turn window means into a minimum-eigenvalue statement, conditional on:
- **U2**: a uniform kurtosis-type bound for the normalized features;
- **H1**: conditioning of the correlation basis (numerically λ_min ≈ 0.05–0.11; 0.107 at ε = 0.2, K = 1 analogue 0.128).
So **A4 for K ≥ 2 is conditional on U2 and H1, not proved unconditionally.**
**Numerics** (`p13d_K2_window.py`, ε = 0.2, c_win median/10%/1%/min): C = 4: 0.353/0.143/0.068/0.030; C = 16: 0.687/0.430/0.288/0.192; C = 64: 0.840/0.707/0.611/0.568, i.e. the same window law in C = T₀ε as for K = 1 (also checked at ε = 0.4).

## 5. Drifting parameters

Scope: drift in μ only, α fixed, μ_t in a box [μ_lo, μ_hi]. The pair coupling (S^max ≥ S^min driven by the same Poisson measure) contracts the Perron-weighted distance at exactly rate γ for *any* μ_t, because μ_t cancels in the difference: E d(S^max_t, S^min_t) = ωᵀe^{Mt}D₀ regardless of the drift. So the mixing constants carry over, with the box constant κ_I in the information bounds.
**Check** (`p13d_drift_window.py`, K = 2, β = (1,2), ε = 0.4, fast sinusoid drift A = 0.4, P = 60): E d_ω at t = 2, 5, 10, 20 is 0.961, 0.198, 0.0152, 0.0000 against exact 0.957, 0.213, 0.0175, 0.0001 (1500 replicates; sampling error).
Window information under sinusoidal μ-drift (K = 1, A = 0.3, P = 1200), c_win relative to the fixed reference metric: C = 16 median 0.561, 1% 0.188, min 0.106; C = 64 median 0.674, 1% 0.388, min 0.261. Relative to the local stationary information: 0.655/0.249/0.155 and 0.802/0.487/0.384, against 0.66/0.26 (C = 16) and 0.84/0.55 (C = 64) stationary: the drift costs little once the metric is local.
**Not covered:** drift in α (changes γ and the Perron weights; the cancellation argument fails), non-box drift.

## 6. Net status of the four caveats

| caveat | status |
|---|---|
| Projection | Proved for the exact F-norm projection (non-expansive); diagonal clip not covered; projected jump-rate slope numerically noisy |
| High probability | Proved for A1′ via predictable leverage and stopping; A2 overshoot step still open |
| Step-3 constant | Derived line by line; valid but ~10¹⁰ conservative |
| A4, K ≥ 2 | Conditional on U2 and H1; numerics support them at K = 2 |
| A4, drift | μ-drift only (α fixed), box constant; α-drift open |

Remaining open from the earlier passes: the smooth-drift constant is still ~25× loose.
