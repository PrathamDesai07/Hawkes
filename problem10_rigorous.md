# Problem 10: proofs, cold-start information, and what remains open

Companion to [problem10_fisher_k1.md](problem10_fisher_k1.md). It addresses four caveats from that note. Notation: K = 1, λ = μ + αS, dS = −βS dt + dN, α = (1−ε)β, a = (1, S), k₀ = 2μ/β, c(u) = εβu + h(αu), h(x) = e^{−x} − 1 + x.

---

## 1. Derivation of the reduction (written out so you can re-derive it)

**Fisher information of the observed path.** Log-likelihood on [0,T]: ℓ_T(θ) = ∫ log λ_t(θ) dN_t − ∫ λ_t(θ) dt, with λ_t(θ) = θᵀa_t and a_t = (1, S_t) predictable and *independent of θ* (β is fixed). Score: ∇ℓ_T = ∫ (a_t/λ_t)(dN_t − λ_t dt), a martingale. Its predictable quadratic variation is ∫ (a aᵀ/λ²)·λ dt. Hence

**I_T(θ) = E_θ ∫₀^T a_t a_tᵀ / λ_t dt**, (integrable because 1/λ ≤ 1/μ and E|a|² < ∞ by the moment bounds of Problem 3).

Stationary observation: I_T = T·I(θ), with I(θ) = E_π[a aᵀ/λ].

**Reduction for K = 1.** Put m = E_π[1/λ]. The generator on f(λ) = λ is A λ = −β(λ−μ) + αλ = βμ − βελ; stationarity gives E λ = μ/ε. Since S = (λ−μ)/α:

- I_μμ = E[1/λ] = m.
- I_μα = E[S/λ] = E[(λ−μ)/(αλ)] = (1 − μm)/α.
- I_αα = E[S²/λ] = E[(λ−μ)²/(α²λ)] = (E λ − 2μ + μ²m)/α².
- det I = (1/α²)[m(Eλ − 2μ + μ²m) − (1−μm)²] = (m·Eλ − 1)/α².

## 2. Derivation of the Laplace-transform equation

Let φ(u) = E_π e^{−uλ}, u ≥ 0. For f_u(λ) = e^{−uλ} (bounded, C¹) the generator of the piecewise-deterministic process is
A f_u = e^{−uλ}[uβ(λ−μ) + λ(e^{−uα} − 1)].
It is bounded for u > 0 (λe^{−uλ} is bounded), so stationarity gives E_π A f_u = 0. With E[λe^{−uλ}] = −φ'(u) (dominated convergence, E λ < ∞):

0 = −uβμ φ − (uβ + e^{−uα} − 1) φ', i.e. **φ'(u) = −uβμ φ(u)/c(u)**, c(u) = uβ + e^{−uα} − 1 = εβu + h(αu).

For u > 0, c(u) ≥ εβu > 0, and s/c(s) → 1/(εβ) as s → 0, so Γ(u) = ∫₀^u s/c(s) ds is finite. The linear ODE gives φ(u) = φ(u₀) exp(−βμ[Γ(u) − Γ(u₀)]); letting u₀ ↓ 0 with φ(u₀) → 1 (φ is a Laplace transform):

**φ(u) = exp(−βμ Γ(u))**, and by Tonelli (1/λ = ∫₀^∞ e^{−uλ} du): **m = ∫₀^∞ φ(u) du.**

Independent check (`code/p10e_derivation_checks.py`, simulated K=1 path, ε = 0.2, μ = β = 1, 4·10⁵ events, not using any formula above): E λ = 5.007 (theory 5.000); φ(u) at u = 0.05, 0.2, 0.5, 1, 2 agrees with the ODE to within 0.5% (e.g. 0.4177 vs 0.4184); the simulated matrix I = ((0.2730, 0.9088), (0.9088, 5.124)) against the reduction with exact m: ((0.2739, 0.9076), (0.9076, 5.115)).

---

## 3. Caveat 1 resolved: the asymptotics are proved

**Theorem A (K = 1, stationary, μ and β fixed, ε → 0).** Let m_ε = E_π[1/λ].
1. If k₀ > 1: m_ε/ε → 2/(2μ−β).
2. If k₀ = 1: m_ε/(ε ln(1/ε)) → 2/β.
3. If k₀ < 1: m_ε/ε^{k₀} → C = (2/β)^{k₀} ∫₀^∞ u^{−k₀} F₀(u) du, where F₀(u) = exp(−βμ ∫₀^u D₀(s) ds), D₀(s) = (β²s²/2 − h(βs))/(h(βs)·β²s/2).

Consequently (from §1) I_μμ = m_ε, I_αα = μ/(εα²)(1 + o(1)), I_μα → 1/β, and det·β² → β/(2μ−β) when k₀ > 1.

**Proof.** Let a = α²/(2β), k_ε = 2βμ/α² (so k_ε → k₀ and k_ε − k₀ = O(ε)).

*Step 1 (Gamma factorisation).* G(s) = s/c(s). Write G(s) = 1/(εβ + α²s/2) + D_ε(s) with
D_ε(s) = (α²s²/2 − h(αs)) / ((εβs + h(αs))(εβ + α²s/2)).
Since h(x) ≤ x²/2, D_ε ≥ 0. And βμ ∫₀^u ds/(εβ + α²s/2) = k_ε ln(1 + au/ε). Therefore
**φ_ε(u) = P_ε(u) F_ε(u)**, P_ε(u) = (1 + au/ε)^{−k_ε}, F_ε(u) = exp(−βμ ∫₀^u D_ε) ∈ (0,1], decreasing.
In particular φ_ε ≤ P_ε, the Laplace transform of a Gamma law. (Checked numerically in `code/p10c_lemma_checks.py`, together with the bounds below.)

*Step 2 (bounds on D_ε).* Write x = αs. Use x²/2·(1 − x/3) ≤ h(x) ≤ x²/2 and x²/2 − h(x) ≤ x³/6.
- For x ≤ 1.5: denominator ≥ (x²/4)(αx/2), numerator ≤ x³/6, so D_ε ≤ 4/(3α).
- For x ≥ 1.5, D_ε is bounded by a constant (numerator ≤ x²/2, denominator ≥ h(1.5)αx/2).
- For x ≥ 1: x²/2 − h(x) ≥ 0.132 x² (the ratio is increasing in x and equals 0.132 at x = 1) and, for ε ≤ 1/4, the denominator is ≤ αx². So D_ε(s) ≥ 0.13/β for s ≥ 1/α, and hence **F_ε(u) ≤ exp(−0.13 μ (u − 1/α))** uniformly in small ε.
- Pointwise, D_ε → D₀, so by dominated convergence F_ε → F₀ locally uniformly. Also, for u ≤ 1.5/α the bound D_ε ≤ 4/(3α) gives F_ε(u) ≥ exp(−4βμu/(3α)); this is the form used in Step 3 for small u, which is all that is needed there.

*Step 3, case k₀ > 1.* m_ε/ε = ∫₀^∞ (1 + av)^{−k_ε} F_ε(εv) dv. As ε → 0, F_ε(εv) ≥ exp(−4βμεv/(3α)) → 1, and the integrand is ≤ (1 + βv/4)^{−(1+k₀)/2} for small ε, which is integrable. Dominated convergence gives m_ε/ε → ∫₀^∞ (1 + βv/2)^{−k₀} dv = (2/β)/(k₀ − 1) = 2/(2μ − β).

*Step 3, case k₀ < 1.* m_ε = ε^{k_ε} ∫₀^∞ (ε + au)^{−k_ε} F_ε(u) du. The integrand is ≤ (au)^{−k_ε} F_ε(u): integrable near 0 uniformly (k_ε ≤ k' < 1) and, by Step 2, exponentially decaying for large u. It converges pointwise to (βu/2)^{−k₀}F₀(u). Dominated convergence gives m_ε/ε^{k_ε} → C, and ε^{k_ε − k₀} → 1.

*Step 3, case k₀ = 1 (β = 2μ).* Here k_ε = (1−ε)^{−2} > 1, and k_ε − 1 ~ 2ε. Fix δ ∈ (0,1). On [δ, ∞): the integrand is ≤ (au/ε)^{−k_ε}F_ε(u), so the contribution is O(ε). On [0, δ]: exp(−4βμδ/(3α)) ≤ F_ε ≤ 1 and
∫₀^δ (1 + au/ε)^{−k_ε} du = (ε/(a(k_ε−1)))[1 − (1 + aδ/ε)^{−(k_ε−1)}] = (ε/a) ln(1/ε)(1 + o(1)),
because (1 + aδ/ε)^{−(k_ε−1)} = 1 − 2ε ln(1/ε)(1 + o(1)). Hence m_ε = (ε/a) ln(1/ε)(1 + O(η(δ)) + o(1)) with η(δ) → 0 as δ → 0, and 1/a → 2/β. ∎

Numerical confirmation: `code/p10_fisher_k1_exact.py` (k₀ = 2: m/ε → 1.998 at ε = 1e-4; k₀ = 0.6: predicted C = 6.052, m/ε^{k₀} = 5.93 and rising) and `code/p10c_lemma_checks.py` (k₀ = 1: the local slope d(m/ε)/d ln(1/ε) is 2.0007, 2.0001, 2.0000, 2.0000 for ε = 1e-5 … 1e-8, against the predicted 2/β = 2).

---

## 4. Caveat 2 resolved: total information over [0, T] is not T·I(θ*) from a cold start

Two facts, one rigorous and one structural, change the conclusion of the earlier note.

**Stationary observation.** If λ_0 ~ π then I_T = T·I(θ*) exactly (§1). So with ε = c/T the μ-information is bounded for k₀ > 1 (→ 2c/(2μ−β)), log-divergent at k₀ = 1, and polynomially divergent for k₀ < 1. This is what the previous note's heuristic computed.

**Cold start (S_0 = 0).** The process is not stationary on [0, T] when ε ≈ 1/T.

*Exact transient formula (K = 1).* Let M_T = ∫₀^T E[1/λ_t] dt and L_T = ∫₀^T E λ_t dt. The same reduction as §1 holds pathwise in expectation:
I_μμ = M_T, I_μα = (T − μM_T)/α, I_αα = (L_T − 2μT + μ²M_T)/α², det = (M_T L_T − T²)/α²,
with L_T = (μ/ε)T − (μ/ε − μ)(1 − e^{−βεT})/(βε) (from E λ_t = μ/ε − (μ/ε − μ)e^{−βεt}).
The Laplace transform φ_t(u) = E e^{−uλ_t} solves ∂_tφ + c(u)∂_uφ = −βμ u φ, a first-order linear PDE solved by characteristics du/ds = c(u). With u_T(u₀) the flow of c for time T from u₀ (so Θ(u_T) − Θ(u₀) = T, Θ' = 1/c), φ_0 = e^{−λ_0u}, and φ = exp(−βμΓ) the stationary transform:

**M_T = ∫₀^∞ du₀ · e^{−λ₀u₀} / (c(u₀) φ(u₀)) · ∫_{u₀}^{u_T(u₀)} φ(u) du.**

(Derivation: φ_t(u) = φ_0(u₀)·φ(u)/φ(u₀) along the characteristic; change variables u → u₀ with du = (c(u)/c(u₀)) du₀, then dt = du/c(u).)
Validation (`code/p10b_transient.py`): formula vs Monte Carlo (4000 paths) for μ = β = 1: ε = 0.2, T = 5: 2.777 vs 2.802; T = 20: 7.309 vs 7.351; ε = 0.5, T = 10: 6.274 vs 6.292.

*Rigorous lower bound, any K.* From Problem 3, A Z = μ − ελ with Z = Σ(α_k/β_k)S_k and λ − μ ≤ β_max Z. From a cold start Z_0 = 0 and d E Z_t/dt = μ − ε E λ_t ≤ μ, so E Z_t ≤ μt, so **E λ_t ≤ μ(1 + β_max t)**. By Jensen, E[1/λ_t] ≥ 1/E λ_t, hence

**I_T(μ,μ) ≥ ln(1 + β_max T)/(μ β_max)** for every K and every ε.

So for μ with α known, cold-start information always diverges at least like ln T, including when k₀ > 1.

*Efficient information (α unknown), K = 1, rigorous.* eff_T := det/I_αα = (M L − T²)/(L − 2μT + μ²M) ≥ (M L − T²)/(L + μ²M) ≥ M(1 − x)(1 − y), x = T²/(ML), y = μ²M/L. With T = c/ε: L_T = (μ/ε²)(c − (1 − e^{−βc})/β)(1 + O(ε)), so T²/L_T → c²/(μ(c − (1−e^{−βc})/β)) is bounded, M_T → ∞ (by the bound above), so x → 0; and y ≤ μT/L = O(ε) → 0. Hence **eff_T ≥ (1 − o(1)) ln(1 + βT)/(μβ) → ∞.**

*Scaling (numerics, K = 1, exact formula, β = 1, T = 2/ε).* d M_T/d ln T between successive rows:

| k₀ | ε = 3e-2 → 3e-4: M_T from … to … | observed behaviour | conjecture |
|---|---|---|---|
| 2 (μ=1) | 8.45 → 18.08; slope 2.17, 2.11, 2.06, 2.03 | tends to 2.0 | M_T ~ [2/(β(2μ−β))] ln T |
| 1 (μ=0.5) | 25.5 → 93.6; slope 11.2, 13.6, 15.9, 18.2 | slope grows roughly linearly in ln T; slope/(2 ln T at the midpoint) = 1.18, 1.15, 1.14, 1.11, drifting down | M_T ~ a·ln² T with a ≈ 1 (the coefficient is not pinned down) |
| 0.6 (μ=0.3) | 61.9 → 549; M_T/T^{0.4}: 11.5, 13.2, 14.6, 15.5, 16.2 | slowly rising | M_T ~ C' T^{1−k₀} |

The k₀ > 1 coefficient has a heuristic derivation: for 1 ≪ t ≪ 1/ε the CIR-limit gives E[1/λ_t] ≈ 2/((2βμ − α²)t), and ∫ dt/t gives the logarithm. It is a conjecture; I did not prove an upper bound on M_T.
Efficient information tracks M_T (e.g. k₀ = 2: 14.57 vs M_T = 18.08 at T = 6667; their difference stays O(1)).

**Consequence.** The "bounded vs divergent" dichotomy of the previous note holds only for stationary observation. For a cold start the baseline information diverges for every k₀, and the dichotomy becomes *logarithmic vs polynomial* divergence, with the same threshold k₀ = 1. Which observation model contribution C4 means has to be decided before Problem 11 is stated as a theorem.

---

## 5. Caveat 3: K ≥ 2

What is done for K ≥ 2:
- Fisher information formula (§1) holds for any K.
- Rigorous cold-start lower bound on the α-known μ-information: ln(1 + β_max T)/(μβ_max) (§4).
- Monte Carlo of the cold-start information matrix (`code/p10d_coldstart_mc.py`; validated against the exact K=1 formula to about 2%): for K = 2, β = (1, 2), α_k = (1−ε)β_k/2, μ = 1, ε = 0.01 the efficient μ-information is 1.06, 1.85, 2.54, 3.27, 4.20, 6.05, 8.12 at T = 2, 5, 10, 20, 40, 80, 160, against K = 1 (exact) 1.07, 1.75, 2.44, 3.30, 4.33, 5.56, 7.03. Same order and same slow growth, so the K = 1 picture is not an artefact. The MC has 50 paths, so T ≥ 80 is noisy; this is consistent with, not proof of, logarithmic growth.

What is not done for K ≥ 2:
- The stationary asymptotics: the Laplace transform now has K arguments and no closed-form ODE solution, so the Gamma-majorant argument of Theorem A has no direct analogue. Problem 2's conditioning results (λ_min(I) ~ c_K·ε) are the numerical substitute.
- A rigorous lower bound on the *efficient* information for K ≥ 2 (needs anti-concentration of S_t over dyadic time windows, the same obstacle as in Problem 2).

## 6. Errata (found while doing this)

- **Monte Carlo thinning bug.** In some of my earlier scripts the thinning loop discarded rejected proposals instead of advancing time (waiting times came out too short). It affected my first versions of `p10b_transient.py` (caught immediately: it disagreed with the exact formula) and `p3c_exp_tail.py` (the Problem 3 exponential-tail check: corrected numbers are now E exp(εZ/e) ≈ 1.27–1.42, and the Problem 3 note was corrected). The same pattern was in the first version of `p2_identifiability.py` (Problem 2 Monte Carlo). That script has since been fixed by another session; I did not re-verify the Problem 2 tables from the fixed script, so check that note before relying on its numbers. The K = 1 exact (ODE) results, on which the Problem 1, 3 and 10 conclusions rest, never used a simulation.
- Everything in this note that is simulated (`p10b_transient`, `p10d_coldstart_mc`, `p10e_derivation_checks`) uses the corrected loop and was validated against an exact quantity before use.

## 7. Status

- Caveat 1 (asymptotics not proved): **resolved**, Theorem A.
- Caveat 2 (stationarity): **resolved**, and it overturned the earlier heuristic for k₀ > 1 (§4).
- Caveat 3 (K ≥ 2): **partly resolved** (definition, rigorous μ-information bound, MC); stationary asymptotics and rigorous efficient-information bound open.
- Caveat 4 (you re-derive): §1 and §2 are written out line by line so you can check them; both were also verified by independent simulation. The step I would check hardest is Step 2 of Theorem A (the explicit constants 4/(3α) and 0.13/β).
- Open: proofs of the cold-start coefficients (conjectured in §4).
