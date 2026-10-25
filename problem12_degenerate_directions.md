# Problem 12: degenerate directions of the Fisher information for K > 1

Setting: stationary observation (the model chosen for C4), θ = (μ, α_1..α_K), a = (1, S_1..S_K), I(θ*) = E_π[a aᵀ/λ], n = Σα_k/β_k = 1 − ε, k = 2μν₁ with ν₁ = Σα_k/β_k² (Problem 11). Built on Problem 11 (Theorem E, Proposition 3) and Problem 10. Problem 2 is being worked on in another session; where I overlap with it I say so.

**Question (from the problem list):** is the near-singular direction n, μ, or a mixture of the α_k across time scales?

## Answer

None of these alone. The spectrum of I has three tiers, and **which tier is smallest depends on ε and on the β-grid**:

| tier | eigen-direction (in α-coordinates) | size | multiplicity |
|---|---|---|---|
| top (slow mode, the **n-direction**) | ∝ (0, 1/β_k) = ∇_α n (numerics, see below) | ≈ 1.7–1.8/ε (μ = 1.5 runs); rigorous lower bound (μ/ε − 2μ)/\|α\|² | 1 |
| middle (**kernel-shape / time-scale mixtures**) | α-directions w with Σ w_k/β_k = 0 | O(1), independent of ε, set by the β-grid | K − 1 |
| ridge (**μ–n ridge**) | ∝ (1, −εα/μ): μ ↑ by δ, all α_k scaled by (1 − εδ/μ) | ≈ ε/(μ(k−1)) for k > 1, ~ ε^k for k < 1 (Problem 11) | 1 |

**The ridge** keeps the stationary mean rate E λ = μ/ε fixed to first order (μ ↦ μ + δ, ε ↦ ε(1 + δ/μ), so μ/ε is unchanged) and keeps the kernel shape fixed (all α_k scaled together). In (μ, n) coordinates it moves μ by δ and n by only −εδ/μ. So the least informative direction is the mean-rate-preserving trade between μ and n, not n itself: n alone is the *best* determined direction (information ∝ 1/ε).

**As ε → 0 the ridge is always the smallest eigenvalue**, because the middle tier stays O(1). **At moderate ε, a kernel-shape direction can be smaller than the ridge**: for K ≥ 3 with a tight grid the smallest middle eigenvalue is tiny, so the smallest eigen-direction is a mixture of the α_k across time scales. The crossover is at ε* ≈ c_min(K, grid)·μ(k−1), where c_min is the smallest middle eigenvalue (table below).

## What is rigorous

1. **Ridge Rayleigh bound (any K).** For v = (1, −εα/μ), vᵀIv = (1+ε)² m − (ε/μ)(1 + 2ε), m = E_π[1/λ] (Problem 11, Prop 3). With Theorem E (m = O(ε) when k > 1): **λ_min(I) ≤ vᵀIv/|v|² = O(ε)**. With Theorem F (m ≥ cε^{k+η}, k < 1) the same Rayleigh quotient is only an upper bound on λ_min of order m.
2. **Top eigenvalue.** For v = (0, α/|α|): vᵀIv = (E λ − 2μ + μ²m)/|α|² ≥ (μ/ε − 2μ)/|α|², so λ_max(I) ≥ (μ/ε − 2μ)/|α|² ~ 1/ε. (This v is only a test vector; the actual top eigenvector is ∝ (0, 1/β_k), the n-direction, see the numerics.)
3. K = 1: both statements are exact (Problem 10): the eigenvalues are ≈ μ/(εβ²) and ≈ εβ/(μ(2μ−β)) (that is, ε/(μ(k−1)) in the notation above), and the eigenvector of the small one is ∝ (1, −εβ/μ).

## What is heuristic (derived, then checked numerically)

**Middle tier.** Write the state as S_k = ∫e^{−β_k(t−s)}dN_s. For fluctuations of the event stream at n = 1, dN/dt = (martingale increment)/(1 − φ̂(ω)), with φ̂(ω) = Σα_k/(β_k + iω) and martingale intensity E λ. For a direction w with Σ w_k/β_k = 0 (the slow mode cancels) one has wᵀH(ω)/(1 − φ̂(ω)) = −a_w(ω)/D(ω) with a_w = Σ w_k/(β_k(β_k+iω)), D = Σα_k/(β_k(β_k+iω)), which is finite at ω → 0. Because Var(wᵀS) ∝ E λ and the weight 1/λ divides it out, the quadratic form is independent of ε:

wᵀRw ≈ (1/2π)∫ |a_w(ω)|²/|D(ω)|² dω  (R = E[SSᵀ/λ] restricted to the fast directions).

(A first attempt using the pure Poisson shot-noise Gram 1/(β_j+β_k) was wrong by a factor of 2; the clustering factor 1/|1−φ̂|² is what fixes it.) I did not prove this approximation, in particular that the fast fluctuations are asymptotically independent of the level of λ, and I ignored coupling with μ and the slow mode.

## Numerics (stationary Monte Carlo, corrected simulator, μ = 1.5, α_k = (1−ε)β_k/K)

**K = 2, β = (1, 2)** (`code/p12_eigenstructure.py`):

| ε | eigenvalues of I | smallest: overlap with (1, −εα/μ) | λ_min/ε | λ_max·ε | middle vs prediction 0.0667 |
|---|---|---|---|---|---|
| 0.1 | 0.032, 0.065, 18.2 | 0.85 | 0.32 | 1.82 | 0.0647 |
| 0.03 | 0.015, 0.064, 57.2 | 0.99 | 0.51 | 1.72 | 0.0642 |
| 0.01 | 0.0053, 0.066, 171 | 0.9992 | 0.53 | 1.71 | 0.0661 |

The ridge prediction ε/(μ(k−1)) with ν₁ = 0.75, k = 2μν₁ = 2.25, μ = 1.5 is ε/(1.5·1.25) = 0.533ε, which matches the measured 0.53ε at ε = 0.01 (and 0.51ε at 0.03). This uses the conjectured constant m ≈ 2ε/(2μ − 1/ν₁); only m = O(ε) is proved (Theorem E).
The top eigenvector: its overlap with α̂ = α/|α| was 0.7999 (K = 2) and 0.5706 (K = 3) at every ε, which equals cos(α̂, 1/β̂) = 0.800 and 0.571 to three digits. So the top eigenvector is ∝ (0, 1/β_k): the gradient of n = Σα_k/β_k with respect to α, i.e. the pure n-direction. (I did not compute its overlap with 1/β̂ directly.)

**K = 3, β = (1, 2, 4)** (`code/p12_crossover.py`): eigenvalues at ε = 0.03: 0.0055, 0.0225, 0.1652, 60.9; at ε = 0.1: 0.0047, 0.0446, 0.1495, 19.1. Eigenvector overlaps at ε = 0.03: the smallest eigenvector overlaps the predicted smallest fast mode at 0.99 (prediction 0.0064 vs measured 0.0055); the second is the ridge at 0.99 (prediction 0.0287 vs 0.0225); the third overlaps the predicted larger fast mode at 1.00 (prediction 0.187 vs 0.165); the fourth is the slow mode. At ε = 0.1 the pattern is the same (overlaps 1.0, 0.94, 0.96). **So for K = 3 on the (1, 2, 4) grid the smallest eigen-direction is a time-scale mixture, not the ridge, down to ε ≈ 0.007.** The corrected spectral prediction for the upper fast mode is within 2–15 % (K = 2: 0.0667 vs 0.0661; K = 3: 0.176 vs 0.174 at ε = 0.01).

**Grid dependence of the fast block** (n = 1, α_k/β_k = 1/K, geometric grid β_k = r^{k−1}; entries: smallest / largest fast eigenvalue; `code/p12_fastblock.py`):

| r | K = 2 | K = 3 | K = 4 | K = 5 |
|---|---|---|---|---|
| 1.5 | 0.031 / 0.031 | 1.3e-3 / 0.097 | 1.0e-4 / 0.19 | 1.3e-5 / 0.30 |
| 2 | 0.067 / 0.067 | 6.0e-3 / 0.176 | 9.5e-4 / 0.30 | 2.2e-4 / 0.45 |
| 4 | 0.106 / 0.106 | 1.35e-2 / 0.225 | 2.5e-3 / 0.38 | 5.3e-4 / 0.59 |

The smallest middle eigenvalue collapses roughly geometrically with K and with grid compression, so a tight, long dictionary has near-singular shape directions long before ε becomes small. The crossover ε* ≈ c_min·μ(k−1) is therefore tiny for K ≥ 4 or r ≤ 2.

## Consequences for the rest of the project

- Problem 2: the small middle eigenvalues are what Problem 2's condition numbers measure (cond ≈ λ_max/c_min ≈ (1/ε)/c_min once ε < ε*). This matches the Cauchy-Gram growth with K and r found there; the corrected spectral version above is the quantitatively right one at n = 1.
- Problems 5, 6, 9: a preconditioner should treat three scales (1/ε, O(1), ε) separately. The "(d/2)·log T" regret target hides the fact that the ridge (information ε per unit time) is the direction that determines the burn-in, and the shape directions need the grid-dependent constant c_min(K, r).
- Problem 11: the ridge is the direction whose information is bounded iff k > 1; the shape directions have information growing linearly in T (O(1) per unit time), so at ε = c/T they are always learnable for fixed grid.

## Status

- Done: identification of the three tiers and the three eigen-directions, grid dependence, crossover rule.
- Rigorous: the ridge Rayleigh bound (with Theorem E), the top-eigenvalue lower bound, K = 1.
- Not proved: the middle-tier asymptotics, that the eigenvector of the smallest eigenvalue converges to the ridge direction (numerical overlaps only), a lower bound on the middle eigenvalues (the Problem 2 anti-concentration issue), and the behaviour for k ≤ 1 (the ridge eigenvalue is then ~ε^k, so the crossover moves to larger ε).
- Re-derive the Rayleigh quotient for v = (1, −εα/μ) yourself; the spectral derivation is only checked numerically.
