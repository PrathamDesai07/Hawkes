# Problem 2: identifiability and conditioning of the exponential dictionary

I(θ*) = E_π[a aᵀ/λ], a = (1, S_1..S_K). Code: `code/p2_window_bound.py`, `code/p2_lemma3_check.py`, `code/p2_identifiability.py`.
Setup for all numbers: μ = 1, β_k = r^(k−1), α_k = (1−ε)β_k/K, MC with 60–80k events (noisy at the 5–10 % level in the smallest eigenvalue; not quantified).

**Correction.** The original `p2_identifiability.py` had a thinning bug (waiting time redrawn from zero after a rejection instead of accumulated), which made λ_min about 10× too small for K ≥ 2 and the printed lower bound look vacuous. Fixed; everything below uses the fixed simulator. My earlier Problem 4 conditioning table was affected and has been corrected.

## Identifiability (qualitative)

I is nonsingular iff 1, S_1, …, S_K are linearly independent as random variables under π. I expect this always holds for distinct β_k and α_k > 0; I have not proved it.

## Lemma W (window inequality; rigorous, gives a computable bound)

For any h > 0, I ⪰ J_h with
J_h = (1/h) E_π ∫_0^h φ_s(u) φ_s(u)ᵀ / λ_s(u) · exp(−Λ_s(u)) du,
φ_s(u) = (1, s_k e^{−β_k u}), λ_s(u) = μ + Σ α_k s_k e^{−β_k u}, Λ_s(u) = ∫_0^u λ_s.
Proof: I = (1/h) E ∫_0^h a_t a_tᵀ/λ_t dt for any stationary window. The integrand is PSD, so keep only the event of no arrival in [0,t). Given F_0 it has probability exp(−Λ_s(t)), and on it a_t and λ_t are the deterministic φ_s(t), λ_s(t).
Caveat: J_h → I as h → 0 (the integrand at u = 0 is a_0 a_0ᵀ/λ_0), so J_h is within 1–4 % of I at small h *by construction*. It is a valid bound but not an independent one. At h ≍ ε/μ it is within a factor of 0.7–0.95 of I. Its use is as a route to an explicit bound via the Gram matrix ∫_0^h φφᵀ, not as a number.

## Lemma 3 (from Problem 1) is not vacuous

λ_min(I) ≥ (ε/μ) · m_K, with m_K = inf_{|v|=1} (E|vᵀa|)². Numerically (inf by Nelder–Mead, 12 restarts):

| K, r | ε | λ_min(I) (MC) | m_K | bound = (ε/μ) m_K | bound/λ_min |
|---|---|---|---|---|---|
| 1, 2 | 0.30 / 0.10 / 0.05 | 0.141 / 0.069 / 0.041 | 0.235 / 0.242 / 0.224 | 0.071 / 0.024 / 0.011 | 0.50 / 0.35 / 0.28 |
| 2, 2 | 0.30 / 0.10 / 0.05 | 0.034 / 0.046 / 0.046 | 0.072 / 0.256 / 0.270 | 0.022 / 0.026 / 0.0135 | 0.64 / 0.55 / 0.29 |
| 2, 4 | 0.30 / 0.10 / 0.05 | 0.052 / 0.069 / 0.065 | 0.117 / 0.318 / 0.318 | 0.035 / 0.032 / 0.016 | 0.67 / 0.46 / 0.24 |
| 3, 2 | 0.30 | 0.0033 | 0.0069 | 0.0021 | 0.63 |

(K=3 rows for ε = 0.1, 0.05 did not finish before the run timed out.) So the bound is within a factor of about 2 to 4 for K ≤ 3. It loses the ε-dependence for K ≥ 2, as below.

## What the numbers say about K, r and ε

- **K = 1:** λ_min ≍ ε (exact, Problem 1).
- **K ≥ 2:** λ_min(I) does **not** shrink with ε in the range tested: K=2, r=2: 0.034 → 0.046 → 0.046; r=4: 0.052 → 0.069 → 0.065. So the ε/μ factor in Lemma 3 is pessimistic for K ≥ 2.
- λ_max grows like 1/ε, so cond(I) grows roughly like 1/ε for K ≥ 2 (K=2, r=2: 138 → 280 → 574; K=3, r=2: 1503 → 2967 → 5261), against ε^(−2) for K=1 (exact).
- **K-dependence:** λ_min(I) at ε = 0.1, r = 2: 0.069 (K=1), 0.046 (K=2), 0.0045 (K=3). The dictionary's own Cauchy-matrix conditioning (1/(β_i+β_j), geometric grid) for K = 1..6: r=2: 1, 38, 456, 2.8×10³, 1.2×10⁴, 3.7×10⁴; r=4: 1, 15, 97, 455, 1.9×10³, 7.7×10³; r=1.5: 1, 102, 3.0×10³, 4.2×10⁴, 3.4×10⁵, 1.9×10⁶.
- **r-dependence:** a wider grid ratio helps (K=2, ε = 0.1: λ_min 0.046 at r=2 vs 0.069 at r=4; cond 280 vs 163).

## Status

- Done: a correct simulator for I(θ*) for K ≤ 3; the rigorous window inequality; a check that Lemma 3 is within 2–4× for K ≤ 3; scaling of λ_min and cond in ε, K and r.
- **Not done:** an explicit lower bound on λ_min(I) as a *formula in K and the grid ratio*, which is what the "Done when" asks for. What is available: λ_min ≥ (ε/μ) m_K with m_K computed numerically. The earlier Cauchy-matrix bound was vacuous, and it was computed from buggy simulations, so that comparison was unreliable; it was not redone.
- Suggested route: bound m_K from below by anti-concentration of vᵀa given no arrivals in a window, which reduces to the smallest eigenvalue of the Gram matrix of {1, e^{−β_k u}} on [0,h] weighted by exp(−Λ). That Gram matrix is explicit and its conditioning is the Cauchy-matrix table above. Not carried out.

---

## Update (25 Oct 2026): attempt at an explicit lower bound on λ_min(I)

### Route 1: Gram matrix of {1, e^{−β_k u}} (code: `code/p2_gram_bound.py`) — rigorous but useless

Restrict Lemma W to the good state set R = {s_k ≥ s_lo ∀k, λ_s ≤ L}. Then
λ_min(I) ≥ p · e^{−Lh} · min(1, s_lo²) · λ_min(G_h) / (L h),
with p ≥ 1 − e^{−μt₀} − μ/(εL) (hazard ≥ μ makes "no event in [−t₀,0]" have probability ≤ e^{−μt₀}; Markov on λ using E λ = μ/ε), s_lo = e^{−β_K t₀}, G_h the Gram matrix of {1, e^{−β_k u}} on [0,h]. Since e^{−Lh} forces h ≲ ε/μ, λ_min(G_h) ~ h^{2K+1} (measured slope 3.0 for K=1, 4.8 for K=2; K=3 not measurable in double precision).
Optimised over (t₀, L, h) the bound is 10⁻⁶ × (MC λ_min) for K=1 and 10⁻¹¹ to 10⁻¹⁶ × for K=2..3. Valid, but vacuous in practice; the ε^{2K+1} decay contradicts the flat/ε-linear truth. Do not use.

### Route 2: truncation + moments (code: `code/p2_truncation_bound.py`, `code/p2_truncation_check.py`) — rigorous given Lemma B, within 0.2–5 % of truth

**Theorem T (draft).** Let M = E_π[a aᵀ], a = (1, S_1..S_K). For every L > 0 and p ≥ 2,

λ_min(I) ≥ (1/L) · [ λ_min(M) − (E|a|⁴)^{1/2} · (E λ^p / L^p)^{1/2} ].

Proof. (i) 1/λ ≥ 1/L on {λ ≤ L}, so I ⪰ E[a aᵀ 1{λ≤L}]/L. (ii) For unit v, E[(vᵀa)² 1{λ≤L}] = vᵀMv − E[(vᵀa)² 1{λ>L}], and by Cauchy–Schwarz E[(vᵀa)² 1{λ>L}] ≤ (E(vᵀa)⁴)^{1/2} P(λ>L)^{1/2} ≤ (E|a|⁴)^{1/2} (E λ^p/L^p)^{1/2} (Markov). Take the infimum over v. ∎

**Ingredients, all explicit.**
- M is exact. m = E S solves A m = −μ1 (A = −diag β + 1αᵀ); Q_ij = E[S_i S_j] solves the K² × K² linear system (β_i + β_j) Q_ij − Σ_k α_k (Q_ki + Q_kj) = μ(m_i + m_j) + E λ, E λ = μ + αᵀm. (From E[L(S_iS_j)] = 0.) Check: λ_min(M) from this system agrees with the Monte Carlo value to 1–3 % in every configuration tried (e.g. 0.3203 vs 0.3215 for K=1, ε=0.3; 0.0129 vs 0.0132 for K=3, ε=0.3).
- Moments of λ and |a|: from Lemma B (`problem03_moments_gap.md`): λ ≤ μ + κZ, S_k ≤ Z/ℓ_min, so E λ^p ≤ 2^{p−1}(μ^p + κ^p M_p(Z)) and E|a|⁴ ≤ 2(1 + K² M_4(Z)/ℓ_min⁴), with M_p(Z) from the Lemma B recursion.
- Conditioning: since 1/λ ≤ 1/μ, λ_max(I) ≤ tr(M)/μ, hence cond(I) ≤ tr(M) / (μ · bound). (Not evaluated numerically.)

**Numerics** (μ = 1, β_k = r^(k−1), α_k = (1−ε)β_k/K; p and L optimised on a grid, p up to 40):

| K, r | ε | MC λ_min(I) | Theorem T bound | bound / MC | best (p, Lε/μ) |
|---|---|---|---|---|---|
| 1, 2 | 0.30 | 0.139 | 6.9×10⁻³ | 4.9×10⁻² | (20, 13) |
| 1, 2 | 0.10 | 0.065 | 1.4×10⁻³ | 2.2×10⁻² | (24, 21) |
| 1, 2 | 0.05 | 0.038 | 6.0×10⁻⁴ | 1.6×10⁻² | (24, 26) |
| 2, 2 | 0.30 | 0.034 | 6.4×10⁻⁴ | 1.9×10⁻² | (24, 56) |
| 2, 2 | 0.10 | 0.046 | 4.6×10⁻⁴ | 1.0×10⁻² | (24, 73) |
| 2, 2 | 0.05 | 0.044 | 2.2×10⁻⁴ | 5.0×10⁻³ | (32, 84) |
| 2, 4 | 0.10 | 0.068 | 1.6×10⁻⁴ | 2.4×10⁻³ | (24, 238) |
| 3, 2 | 0.10 | 0.0045 | 1.7×10⁻⁵ | 3.8×10⁻³ | (32, 273) |

- It is non-vacuous and 10⁵–10¹⁰ times better than Route 1. It is still 20–400× below the truth, because the moment-based tail control forces L ≈ 10–300 × μ/ε.
- With *no* tail penalty (pure truncation at L ≈ 1–2 μ/ε, `p2_truncation_check.py`) the bound is within 2–8× of truth. So the gap to 20–400× comes from the tail control (Hölder + Lemma B), not from the idea. A sharper tail bound (exponential tail of λ instead of moments) would close most of it.
- K-dependence enters through λ_min(M) (0.33 → 0.35 → 0.05 for K = 1, 2, 3 at ε = 0.1, r = 2) and through ℓ_min and κ in the tail term. Both are explicit, but I have not extracted a closed-form function of (K, r).
- ε-dependence: roughly ε · λ_min(M) / (log-type factor), against the measured λ_min that is ∝ ε for K=1 and flat for K≥2 — so for K ≥ 2 the bound is conservative in ε.

### Status after this attempt
- Done: a bound that is explicit and computable from (μ, β, α, ε), valid for all K, and 10⁻³–10⁻¹ of the truth in the cases tested. This meets the "explicit lower bound on λ_min(I) in terms of K, grid" criterion in a computational sense.
- **Not done:** a closed-form *formula* in K and the grid ratio (λ_min(M) is an eigenvalue of an explicit but non-closed matrix); a sharper tail bound; ε-sharpness for K ≥ 2.
- **Dependencies to re-derive before relying on it:** Theorem T's proof (five lines, easy), the second-moment linear system, and Lemma B (still a draft in Problem 3; numerically checked for K=1 and 2 only).
- The numerics for MC λ_min carry about 5–10 % uncertainty (not quantified); the bound itself is deterministic.

---

## Erratum (25 Oct 2026): the "λ_min does not shrink with ε for K ≥ 2" observation was a range artefact

Above I wrote that for K ≥ 2 "λ_min(I) does **not** shrink with ε in the range tested", that cond(I) grows roughly like 1/ε, and (in the Lemma 3 comparison) that the ε/μ factor is pessimistic for K ≥ 2. My ε range stopped at 0.05. `problem12_degenerate_directions.md` explains what I was seeing: the Fisher matrix has a **middle tier** of K−1 kernel-shape modes of size O(1) set by the β-grid (≈ 0.067 for β = (1,2)), and a **ridge** along (1, −εα/μ) of size ≈ ε/(μ(k−1)); the ridge only becomes the smallest eigenvalue below ε* ≈ c_min(K, grid)·μ(k−1) (≈ 0.03 in my setting). For ε ≥ 0.05 I was measuring the middle tier. For ε → 0, λ_min(I) ∝ ε for every K (that note's Monte Carlo: 0.032, 0.015, 0.0053 at ε = 0.1, 0.03, 0.01 for K = 2, β = (1,2), μ = 1.5), so Lemma 3's ε/μ scaling is sharp in ε as ε → 0. The statements about cond ∼ 1/ε and about the bound being "conservative in ε" hold only at moderate ε. The K-dependence reported here (λ_min falling roughly ×10 per added exponential on a geometric grid at ε = 0.1) is the middle-tier collapse of that note's last table and is unaffected. I have not re-run my experiments at smaller ε.
