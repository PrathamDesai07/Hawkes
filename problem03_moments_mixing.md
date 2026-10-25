# Problem 3: moment and mixing bounds for the state S as n → 1

Setup: S = (S_1..S_K), dS_k = −β_k S_k dt + dN_t, λ = μ + αᵀS, n = Σα_k/β_k = 1 − ε, κ = β_max/β_min. Generator A.

## Key identity (the one trick everything uses)

Let c_k = α_k/β_k and Z = Σ c_k S_k (so Σc_k = n). Then

**A Z = μ − ελ**, and β_min Z ≤ λ − μ ≤ β_max Z.

Derivation: the continuous drift of Z is −Σα_kS_k = −(λ−μ), and Z jumps by n at rate λ, giving −(λ−μ) + nλ = μ − ελ. So Z has a restoring drift at rate ε·β_min even though no single S_k does.
With one exponential (K=1), Z = (λ−μ)/β and the drift is exactly −εβ(λ−·).

## Lemma 1 (first and second moments, uniform in ε up to explicit powers)

1. E λ = μ/ε (stationary), hence E Z ≤ μ(1−ε)/(εβ_min).
2. A Z² = 2Z(μ − ελ) + n²λ. Using λ ≥ μ + β_min Z in the negative term and λ ≤ μ + β_max Z in the positive one,
   A Z² ≤ −2εβ_min Z² + (2μ(1−ε) + n²β_max) Z + n²μ.
   Taking stationary expectations (A Z² integrates to 0; integrability follows from Lemma 4 below):

   **E Z² ≤ μ / (2β_min ε²) · [ (2μ + n²β_max)/β_min + n²ε ]**, and E λ² ≤ 2μ² + 2β_max² E Z².

Numerical check (exact stationary second moments from the closed linear system for E[S_jS_k], `code/p3_moments_mixing.py`): the bound holds in every case, with E Z² / bound = 1.00 for K=1 as ε → 0 (tight), and 0.47, 0.23, 0.21, 0.10 for β = (1,2), (1,4), (1,2,4), (1,10). The scaling ε²·E λ² → constant is confirmed (for example 1.50 for K=1, μ=β=1, which is the Gamma/CIR value (2βμ+α²)μ/(2β)). The ε^{−2} order is sharp and only the constant carries κ.

## Lemma 1' (all moments p ≥ 1)

Write b = β_min, B = β_max, m_q = E Z^q (stationary). For integer p ≥ 1,
A Z^p = pZ^{p−1}(μ − ελ) + λ Σ_{j=2}^{p} C(p,j) nʲ Z^{p−j}.
Apply λ ≥ μ + bZ to the first term and λ ≤ μ + BZ, n ≤ 1 to the sum:

A Z^p ≤ −pεb Z^p + Σ_{i=1}^{p} a_{p,i} Z^{p−i}, with a_{p,1} = pμ + C(p,2)B and a_{p,i} = μC(p,i) + B·C(p,i+1) for i ≥ 2 (C(p,p+1) = 0).

Stationarity (E A Z^p = 0; integrability is justified by Lemma 4 below, which gives all polynomial moments) yields pεb·m_p ≤ Σ_i a_{p,i} m_{p−i}. Induct on p with m_q ≤ C_q ε^{−q}: each term has m_{p−i} ≤ C_{p−i} ε^{−(p−i)}, and since ε ≤ 1, ε^{−(p−i+1)} ≤ ε^{−p} for every i ≥ 1. Hence

**E Z^p ≤ C_p ε^{−p}**, with C_0 = 1 and C_p = (1/(p·b)) Σ_{i=1}^{p} a_{p,i} C_{p−i}.

The constants depend only on (μ, b, B). Examples with μ = b = B = 1: C_1..C_6 = 1, 2, 5.67, 20.7, 92.1, 484. Since λ − μ ≤ B Z, E λ^p ≤ 2^{p−1}(μ^p + B^p C_p ε^{−p}).
Check (`code/p3b_general_moments.py`, exact stationary moments of S from the closed linear system, p ≤ 4, K ≤ 3, ε ∈ {0.3, 0.1, 0.03}): m_p ε^p / C_p was at most 0.97 everywhere, and approaches 1 only for p = 1 and K = 1 as ε → 0. For K=1, p=4 it is 0.23–0.32. The bound is valid but not tight for p ≥ 2.

## Lemma 4 (exponential moment and tail, uniform in ε)

Let θ ≤ ε/e. Then A e^{θZ} = e^{θZ}[θμ + λ(e^{θn} − 1 − θ)], and e^{θn} − 1 − θ = (e^{θn} − 1 − θn) − εθ ≤ −εθ + θ²n²e^{θn}/2 ≤ −εθ/2 (using θn ≤ 1 so e^{θn} ≤ e, and θ ≤ ε/e). With λ ≥ μ + bZ:

A e^{θZ} ≤ θ e^{θZ} (μ − εbZ/2).

Truncate with ψ_M(z) = e^{θ min(z, M)}: it is bounded, equals e^{θz} for z ≤ M, and A ψ_M ≤ θ e^{θZ}(μ − εbZ/2) on {Z ≤ M} and = 0 on {Z > M} (a jump above M changes nothing, and ψ_M is flat there). Stationarity gives 0 = E Aψ_M(Z). Put z_0 = 4μ/(εb). On {Z > z_0} the bracket is ≤ −μ, and on {Z ≤ z_0} it is ≤ μ, so
0 ≤ μ e^{θz_0} − μ E[e^{θZ}; z_0 < Z ≤ M].
Let M → ∞ (monotone convergence): E[e^{θZ}; Z > z_0] ≤ e^{θz_0}, hence

**E exp(εZ/e) ≤ 2 exp(4μ/(e·β_min))**, independent of ε (take θ = ε/e),
**P(Z ≥ t) ≤ 2e^{4μ/(eβ_min)} · e^{−εt/e}**, and E exp((ε/(eB))(λ − μ)) ≤ the same constant.

So Z and λ are sub-exponential at scale 1/ε, matching the Gamma/CIR scale (the true exponent for K=1 is about 2ε/β for λ, against ε/(eβ) here, so the constant in the exponent is loose by a factor of about 5). In particular all polynomial moments are finite, which justifies E A Z^p = 0 used in Lemma 1'.

Check (`code/p3c_exp_tail.py`, time-weighted simulation, 1.5·10⁵ events, K ∈ {1,2}, ε ∈ {0.2, 0.1}): E exp(εZ/e) ≈ 1.27–1.42 against the bound 8.71 (for K=1, ε=0.1 the Gamma/CIR prediction is about 1.5, consistent). The tail bound is trivial (≥ 1) for t ≲ 6/ε and only informative beyond that. In the runs the simulated tail at t = 4/ε was about 1e-4 to 1e-3, and the event t = 8/ε was never seen. So the bound is correct in scale but loose in constants, and I did not test its far tail.
*Correction (2026-10-25):* the first version of this check reported 2.3–2.7 and tails of 0.03–0.08. Those came from a bug in my thinning loop (rejected proposals did not advance time, so waiting times were too short). Fixed; numbers above are from the corrected script. The proved bounds were not affected.

## Theorem 2 (contraction in a weighted Wasserstein distance, rate ≥ εβ_min)

Metric on states: d(s, s') = Σ_k c_k |s_k − s'_k|.

Proof idea (monotone coupling):
1. Drive all copies with one Poisson random measure: a point at height u fires copy x iff u ≤ λ_x. The linear Hawkes dynamics are monotone, so componentwise s ≤ s' implies S_t ≤ S'_t for all t.
2. For any pair (s, s'), the copies started at s∧s' and s∨s' bracket both S_t and S'_t componentwise, so d(S_t, S'_t) ≤ ΔZ_t := Z(S^∨_t) − Z(S^∧_t), with ΔZ_0 = d(s, s').
3. For the ordered pair, D = S^∨ − S^∧ ≥ 0 decays at rate β_k in each coordinate, and jumps by +1 in every coordinate at rate αᵀD (only the upper copy fires). So d/dt E ΔZ = −Σα_kD_k + n Σα_kD_k = −ε Σα_k E D_k ≤ −εβ_min E ΔZ.

Conclusion: **W_d(P_t(s,·), P_t(s',·)) ≤ e^{−εβ_min t} d(s, s')** for all s, s'. Consequences: unique stationary law, and geometric ergodicity in W_d with rate at least εβ_min and no small-set (minorisation) argument needed. This avoids the usual Foster–Lyapunov difficulty that the level set {Z ≤ R} grows like 1/ε.

Numerical check (`code/p3_coupling_check.py`, K=2, ε=0.2, 3000 replicas): no order violations in any run, and E ΔZ_t agrees with the exact ODE c ᵀ e^{Mt} D_0 to within about 5%. The simulated values sit slightly below the ODE at every time point; I did not quantify the Monte Carlo error, so I don't know whether that is noise or a small bias.

## Proposition 3 (the rate is exactly of order ε, and the true gap is bounded above)

Means obey E[S_t | S_0] − m = e^{Mt}(S_0 − m) with M = −diag(β) + 1αᵀ. The decay rate γ of the slowest mode solves

**Σ_k α_k / (β_k − γ) = 1**, with 0 < γ < β_1 (real root; it exists because Σα_k/β_k = n < 1), and γ = ε / Σ_k(α_k/β_k²) + O(ε²).

The corresponding left eigenvector gives a linear function ℓ(S) with E[ℓ(S_t)|S_0] − Eℓ = e^{−γt}(ℓ(S_0) − Eℓ), so correlations decay exactly like e^{−γt}. Hence the L² spectral gap is at most γ and the relaxation time is at least 1/γ.
Together with Theorem 2:

**εβ_min ≤ (contraction rate) ≤ γ ≈ ε / Σ(α_k/β_k²) ≤ ε β_max / n.**

For K=1 both ends equal εβ (exact). Numerically γ/ε is 1.00 (K=1), 1.3 (β=1,2), 1.6 (β=1,4), 1.8 (β=1,10), 0.6 (β=0.3,1,3), and ε/Σ(α/β²) matches γ to about 1% at ε = 0.01.

## What this means for the rest of the project

- **Moment bounds:** E λ^p ≲ ε^{−p} for every p (Lemma 1'), plus a sub-exponential tail at scale 1/ε (Lemma 4), with constants depending only on μ, p, β_min, β_max. These are what the martingale concentration in Problem 5 needs for the gradient bounds (Freedman needs conditional variance bounds, which the moments give).
- **Mixing time:** 1/γ ≍ 1/(εβ), independent of μ. This is different from the problem list's guess of burn-in ~ 1/(με) in time. That guess may refer to information accumulation (T μ (1−n) ≳ 1) rather than relaxation of the state, so Problem 6 should treat the two burn-in notions separately. In events rather than time, relaxation takes about μ/(ε²β) events, because the event rate is μ/ε.
- **Constants carry κ = β_max/β_min.** The bounds are not uniform over the β-grid, in the same direction as the K-dependence found in Problem 2.

## Status and open points

- Proved here: Lemma 1 (p = 1, 2), Lemma 1' (all p, by induction), Lemma 4 (exponential moment and tail, which also settles the integrability step), Theorem 2, Proposition 3 (the exact rate equation; the first-order expansion of γ is computed, not formally derived). Numerics agree with every bound. This meets the "Done when" of Problem 3 (explicit ε-dependence of the moment bounds and of the contraction rate), subject to your re-derivation.
- Not done: whether the true L² spectral gap equals γ (only an upper bound is shown); constants uniform over the grid (they carry κ and, for p-th moments, grow like the C_p recursion); a bound on sup over a time horizon [0, T] of λ_t (Lemma 4 is a one-time stationary tail, so a union bound over time or events needs an extra argument, probably an exponential supermartingale or a Palm-measure version); behaviour when started far from stationarity with ε ≈ 1/T; sharper constants in the tail exponent.
- Re-derive the key identity A Z = μ − ελ and the coupling argument yourself, as the problem list asks. The numerics support them but are not a proof.
