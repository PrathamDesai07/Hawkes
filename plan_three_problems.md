# Plan: the three items that are still weaker than "Done" (Problems 12, 14, 11)

Method for each item: (1) pin down exactly what is missing; (2) do the mathematics first (write the statement, the proof strategy and the exact generator/algebra computations); (3) only then code: symbolic verification of the algebra, numerical tests of the inequalities, simulations of the consequences; (4) record status honestly. Order: 12 (cleanest, a Lyapunov argument), 14 (recursion + deterministic skeleton), 11 (new structure: observed information, non-consistency).

## Problem 12: middle-tier upper bounds and eigenvector convergence
**Missing:** rigorous upper bounds for the K−1 middle eigenvalues of I = E[aaᵀ/λ]; convergence of the smallest eigenvector to the ridge direction.
**Mathematical idea:** the matrix A_mat = diag(β) − 𝟙αᵀ has real positive eigenvalues γ = a₁ < a₂ < … < a_K with left eigenvectors ℓ_j. For the *fast* modes y_j = ℓ_jᵀS the generator is exactly A y_j = −a_j y_j + μσ_j (σ_j = ℓ_jᵀ𝟙): no coupling to the slow level. A Lyapunov function f_j = y_j²/W (W = Z + z) has generator ≤ −κ_j f_j + C, so E_π[y_j²/W] = O(1), hence E[y_j²/λ] = O(1) uniformly in ε. Courant–Fischer on the span of the y_j's gives the upper bounds; Davis–Kahan with the tangent-bound gap λ₂ ≥ g and the ridge Rayleigh quotient gives the eigenvector statement.
**Code:** sympy check of A f_j; random-state test of the inequality; evaluation of the constants for K = 2, 3 against Monte Carlo; overlap bound vs measured overlaps.

## Problem 14: a tracking theorem for the estimators and the amplitude effect
**Missing:** (a) a theorem on the tracking error of the forgetting-Newton estimator; (b) an analytic explanation of the amplitude-dependent excess of the data-driven preconditioner.
**Mathematical idea (a):** with F_i = (1−η)F_{i−1} + ΔF_i the Newton inequality gives D_i ≤ (1+ρ)[(1−η)D_{i−1} + X_i − 2⟨g_i,v_i⟩ + g_iᵀF_i⁻¹g_i] + (1+1/ρ)‖d_i‖²_{F_i}; take ρ = η/2, use E⟨g,v⟩ ≈ E X (compensator identity, with a drift correction), the leverage control of Theorem R for E g F⁻¹ g ≈ ηd, and F_i ≍ I_ev/η to obtain E‖θ_i − θ*_i‖²_I ≲ 2ηd + 4(‖d‖_I/η)².
**Mathematical idea (b):** in the linear regime the mean error obeys e_{i+1} = (1 − f_i/F̄_i)e_i − d_i; the gain f_i/F̄_i is η times the ratio of the current information to its exponentially lagged average. The information about μ depends strongly on μ (d ln f/d ln μ ≈ −3 or −4), which makes the gain modulation O(1)·A. A deterministic skeleton (no noise) should reproduce the measured excess.
**Code:** deterministic skeleton vs the measured ratios; Monte Carlo check that the theorem's bound dominates the measured MSE and has the right η-scaling.

## Problem 11: the matching estimator rate for k_K < 1
**Missing:** an upper bound matching the information lower bound T^{−(1−k)/2}.
**Mathematical idea:** (i) the error of the MLE is governed by the *observed* information 𝓘_T = ∫dt/λ_t, not its mean: on {𝓘_T ≥ M}, |μ̂ − μ| ≤ C/√M with high probability (monotone score + martingale bound); (ii) the expected information diverges for k < 1 only through rare excursions to the lower boundary: stopping the process when λ first drops below κμ/ε keeps the KL bounded pathwise, and the event {no drop} has probability bounded below, so *no* estimator is consistent at ε = c/T for any k (a stronger and different conclusion than the information dichotomy suggests). 
**Code:** simulate the law of 𝓘_T and of the MLE error for k = 0.6 and 2 at growing T with ε = c/T.
