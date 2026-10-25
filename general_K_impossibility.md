# Impossibility for the stationary record, any K (bypassing TV(π_θ, π_θ′))

Written 25 Oct 2026. **Result in one line:** I did **not** obtain TV(π_θ, π_θ′) = O(δ) for K ≥ 2 (it is probably true; see "What was tried"), but the statistical consequence it was needed for — no consistent estimator of μ from a stationary record with unobserved past when k_K > 1 — **can be proved for every K by a different route that never compares the two stationary laws.** This is my proof, written in one pass; re-derive it before relying on it. Code: `code/verify_survival.py`, `code/verify_ghost_hit.py`.

## Setting (as in `problem11_dichotomy.md` §4b)
α_k > 0, β_1 < … < β_K, n = Σα_k/β_k = 1 − ε, b = β_1, B = β_K; γ ∈ (0,b) the root of Σα_k/(β_k − γ) = 1; ω_k = α_k/(β_k − γ) (Σω_k = 1), Z = ωᵀS, ν̃ = Σω_k/β_k, γ = ε/ν̃, k_K = μ/c₀ with c₀ = n/(2ν̃). Assume k_K > 1. Ridge pair θ = (μ, α), θ′ = (μ + δ, α(1 − εδ/μ)) (all α_k scaled), u = θ′ − θ. Q_θ, Q_θ′ are the laws of the event record on [0,T] of the stationary processes (the state S₀ unobserved), ε = c/T.

## What was tried first: the direct route, and why it does not extend
For K = 1 the stationary law is accessible through its explicit characteristic function (`theoremA_step2_and_stationary_tv.md`). For K ≥ 2 the characteristic function of S ∈ R^K_+ solves a first-order PDE solved by a K-dimensional flow ẇ = −V(w), V(u) = Bu − α(1 − e^{−|u|}), φ(u) = exp(−μ∫₀^∞(1 − e^{−|w(s)|})ds); it has no closed form, there is no scalar Gamma comparator for the joint law (the fast coordinates fluctuate on a different scale from the slow mode Z), and a Fourier–Sobolev bound in R^K would require uniform complex-flow estimates that I did not manage. Comparing the stationary laws through their Poisson-cluster (shot-noise) representation fails too: the two Lévy intensities have different immigration rates over the infinite past, so the KL between the underlying Poisson random measures is infinite even though the laws of S₀ are close. **TV(π_θ, π_θ′) for K ≥ 2 remains unproved.**

## The route that works: run the θ′-dynamics from π_θ for a time L = K₀/γ, then couple to the stationary θ′-process

Let R_L be the law of the record on [0,T] of the process with S_{−L} ~ π_θ and θ′-dynamics on [−L, T].

**Step 1 (KL, same initial law).** P_θ restricted to [−L,T] and the process R_L restricted to [−L,T] share the initial law π_θ, so their KL is exactly (L+T)·r′ with r′ = E_{π_θ}[λ log(λ/λ′) − λ + λ′] ≤ ρ·uᵀI(θ)u, ρ = max(μ/μ′, 1/(1 − εδ/μ)). For the ridge direction uᵀIu = δ²[(1+ε)²m − (ε/μ)(1+2ε)] ≤ δ²(1+ε)²m with m = E_π[1/λ] ≤ C_Eε by Theorem E (this uses k_K > 1; verified in `verification_thmA_E_F.md`). With L = K₀/γ and γ = ε/ν̃, εL = K₀ν̃ and εT = c:
 (L+T)·r′ ≤ ρ(1+ε)²C_E(K₀ν̃ + c)δ².
By Pinsker and data processing, **TV(Q_θ, R_L) ≤ δ·[ρ(1+ε)²C_E(K₀ν̃ + c)/2]^{1/2}.** No comparison of π_θ with π_θ′ occurs.

**Step 2 (R_L versus the stationary θ′-record, by superposition).** For a linear Hawkes process the intensity from time −L on is μ′ + g_s(t) + Σ_iα′(t − t_i), where s = S_{−L} and g_s(t) = Σ_kα′_ks_ke^{−β_k(t+L)} is deterministic given s. Hence (Hawkes–Oakes cluster representation with a deterministic time-varying baseline) the process equals, in law, the union of (E) the empty-start process driven by immigrants of rate μ′ after −L and (X_s) independent "ghost" immigrants of intensity g_s(t) together with their clusters. R_L uses s ~ π_θ, the stationary θ′-process uses s′ ~ π_θ′; E can be shared. So the two records coincide on [0,T] unless X_s or X_{s′} has an event in [0,T]:
 **TV(R_L, Q_θ′) ≤ P(X_s hits [0,T]) + P(X_{s′} hits [0,T]).**

**Step 3 (ghost-hit probability, needs the cluster survival bound).** Let q(r) = P(a cluster rooted at 0 has an event after time r). A ghost immigrant born at τ = t + L ∈ [0, L) reaches [0,T] only if its cluster has an event after relative time L − τ. A union bound over the Poisson ghost immigrants gives
 P(X_s hits | s) ≤ ∫₀^{L+T} h_s(τ)p(τ)dτ, h_s(τ) = Σα′_ks_ke^{−β_kτ}, p(τ) = min(1, q(L − τ)) for τ < L and 1 on [L, L+T].

**Lemma S (survival probability, any K).** Let Y_k(s) = (α_k/β_k)e^{−β_ks} + α_k∫₀^se^{−β_k(s−v)}q(v)dv. Then q(s) = 1 − exp(−|Y(s)|), |Y| = ΣY_k, and Y solves
 Y′ = −BY + α(1 − e^{−|Y|}), Y(0) = (α_k/β_k).
Moreover ζ = Σ_kY_k/(β_k − γ) satisfies ζ′ = −γζ − r(|Y|), r(x) = x − 1 + e^{−x} ≥ 0, and for s ≥ 1/γ
 **q(s) ≤ C′γe^{−γs}, C′ = (B − γ)/(c(1 − e^{−1})), c = c_r(b − γ)², c_r = r(R₀)/R₀², R₀ = (B − γ)/(b − γ).**
*Proof.* Offspring of the root form a Poisson process of intensity α(r) = Σα_ke^{−β_kr}; a sub-cluster rooted at r has no event after s with probability 1 − q(s − r) for r ≤ s, and any offspring at r > s is itself an event after s. Hence 1 − q(s) = exp(−ᾱ(s) − ∫₀^sα(r)q(s−r)dr), with ᾱ(s) = Σ(α_k/β_k)e^{−β_ks}; the exponent is |Y(s)|, and differentiating the definition of Y_k gives the ODE. For ζ: using Σα_k/(β_k − γ) = 1 and β_k/(β_k − γ) = 1 + γ/(β_k − γ), ζ′ = −Σβ_kY_k/(β_k − γ) + q = −|Y| − γζ + (1 − e^{−|Y|}) = −γζ − r(|Y|). Since Y_k ≥ 0, (b − γ)ζ ≤ |Y| ≤ (B − γ)ζ; ζ is non-increasing, so |Y| ≤ (B − γ)ζ(0) < (B − γ)/(b − γ) = R₀; r(x)/x² is decreasing (it is h(x)/x² from `theoremA_step2_and_stationary_tv.md`), so r(|Y|) ≥ c_r|Y|² ≥ c ζ². Then v = 1/ζ satisfies v′ ≥ γv + c, so ζ(s) ≤ ζ₀e^{−γs}/(1 + (cζ₀/γ)(1 − e^{−γs})) ≤ (γ/(c(1 − e^{−1})))e^{−γs} for s ≥ 1/γ, and q ≤ |Y| ≤ (B − γ)ζ. ∎
(The weaker first-moment bound q(s) ≤ E[#events after s] ≤ (Σα_k/γ)e^{−γs}, from m(t) ≤ (Σα_k)e^{−γt}, loses a factor of order 1/γ² = O(ε⁻²) and would force L ≈ (2/γ)ln(1/ε), turning the cost in Step 1 into a ln(1/ε) factor; the quadratic damping term −r(|Y|) is what removes it.)

Applying Lemma S for τ ≤ L − 1/γ: ∫₀^{L−1/γ}h_s(τ)C′γe^{−γ(L−τ)}dτ ≤ C′γe^{−γL}Σα_ks_k/(β_k − γ) = C′γe^{−γL}Z(s) (up to the factor α′/α = 1 + O(εδ)). The remaining part, τ ≥ L − 1/γ, is at most e^{−b(L−1/γ)}G(s) with G(s) = Σα_ks_k/β_k ≤ (λ − μ)/b. By Lemma G (E_πZ = μ/γ, E_πλ = μ/ε):
 **E P(X_s hits) ≤ C′μe^{−γL}(1 + O(εδ)) + e^{−b(L−1/γ)}μ/(bε),**
and the second term is ≤ e^{−bK₀/(2γ)} ≪ δ (b/γ ≍ b ν̃/ε → ∞). The same bound holds for X_{s′} with s′ ~ π_θ′ (E_{π_θ}Z′ ≤ (1 + O(εδ))E_{π_θ}Z).

## Theorem V (any K, stationary record, unobserved past)
Let k_K > 1 and ε = c/T → 0. There exist ε₀, δ₁, C (depending on the dictionary, μ, c) such that for ε ≤ ε₀ and 0 < δ ≤ δ₁,
 **TV(Q_θ, Q_θ′) ≤ Cδ(ln(1/δ) + 1)^{1/2}**
(take K₀ = ln(1/δ) in Steps 1–3: Step 1 gives ≤ Cδ(K₀ + c)^{1/2}; Steps 2–3 give ≤ 2C′μe^{−K₀(1−O(δ))} + o(δ) = O(δ)). Hence for δ small enough TV ≤ 1/2 and, by Le Cam, every estimator satisfies max_{θ,θ′}P(|μ̂ − μ_i| ≥ δ/2) ≥ 1/4: **μ is not consistently estimable from a stationary record with unobserved past when k_K > 1, with or without α known, for every K.** For K = 1 this is weaker than the O(δ) bound of `theoremA_step2_and_stationary_tv.md` (a ln(1/δ) factor), but it needs neither a closed-form transform nor the comparison TV(π_θ, π_θ′).

## Numerical checks (K = 1, 2, 3)
- **Lemma S derivation, flow vs cluster simulation (60 000 clusters):** K = 1, ε = 0.1: s = 0.5, 2, 5, 10, 20, 40: MC 0.5225, 0.3666, 0.2048, 0.0964, 0.0287, 0.0035 vs flow 0.524, 0.3693, 0.2057, 0.0958, 0.0285, 0.0035. K = 2 (β = (1,2)): MC 0.4938, 0.316, 0.1586, 0.0653, 0.0148, 0.001 vs flow 0.4964, 0.3185, 0.1599, 0.0648, 0.0147, 0.001.
- **ζ′ = −γζ − r(|Y|)** with the analytic derivative: max abs difference 2.3 × 10⁻¹³ (K = 2, ε = 0.05).
- **Survival bound:** the Bernoulli inequality holds on the whole grid for K = 1, 2, 3 and ε = 0.2, 0.05. The observed ratio q(s)/(γe^{−γs}) for s ≥ 1/γ is at most 1.96 / 2.96 (K = 1), 1.50 / 2.28 (K = 2), 1.28 / 1.97 (K = 3), against the proof constants C′ = 5.4 / 4.5, 19.2 / 12.6, 78.1 / 40.5. The first-moment bound is 12–223× larger than the true q at s = 2/γ, growing as ε decreases.
- **End-to-end ghost-hit probability** (K = 2, β = (1,2), ε = 0.05, μ = 1.5, T = 1/ε, 3000 stationary states s ~ π_θ, direct simulation of the ghost immigrants and their clusters): P(X_s hits [0,T]) = 0.646, 0.279, 0.104, 0.0437 for K₀ = 1, 2, 3, 4 (L = K₀/γ), i.e. decay ≈ e^{−1} per unit of K₀ as e^{−γL} predicts, with an O(1) constant (≈ 2.4). The exact expected number of hitting ghost clusters computed from the flow is 1.297, 0.358, 0.1207, 0.0432 (an upper bound on the probability, as it must be). E_πZ = 22.42 vs μ/γ = 22.63.

## Limits
- Needs k_K > 1 and Theorem E (hence ε small in terms of the dictionary); constants C′, C_E are explicit but not optimised, and ε₀ is not explicit.
- The superposition/cluster representation with a deterministic time-varying baseline is standard but I did not write a proof.
- The ln(1/δ) factor is probably an artefact; whether TV(Q_θ, Q_θ′) = O(δ) for K ≥ 2 (and whether TV(π_θ, π_θ′) = O(δ)) is open.
- Nothing here covers k_K ≤ 1 (no impossibility is expected there: information diverges).
- This is a statement about the **stationary** experiment. In the empty-start experiment (the actual streaming setting) the μ-information diverges like ln T for every k and the relevant lower bound is logarithmic, not inconsistency (Problem 10 §4).
