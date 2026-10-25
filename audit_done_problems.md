# Audit of the problems marked Done (1, 3, 10, 11, 12, 14), 2026-10-25

Scope: re-check of proof steps, independent re-runs on fresh random parameters, source checks, cross-session consistency. **This is an audit, not a certificate of correctness:** the proofs have not been formally verified, the numerical checks are finite samples, and several results are explicitly heuristic or numerical (flagged in each note).

## Errors found and fixed during this audit

| where | error | effect | fix |
|---|---|---|---|
| `problem14_drift_model.md`, tracker row 14 | Hall–Willett Theorem 1 was read as a **squared** deviation Σ‖·‖² (the text extraction had dropped a subscript). The LaTeX source shows Σ‖λ_{t+1} − Φ_t(λ_t)‖₂, the ℓ₂ **norm** (a path length relative to known dynamics). | The literature comparison was wrong: their variation measure is the same type (path length) as the bounded-variation class chosen in Problem 14, not a different one. The decision itself is unaffected. | Paragraph and tracker note rewritten from the LaTeX source. |
| `problem10_rigorous.md`, Theorem A, step 2 | (a) For x ≥ 1.5 the bound on D_ε used "denominator ≥ h(1.5)αx/2", which gives a bound growing linearly in x, not a constant. (b) "denominator ≤ αx² for ε ≤ 1/4": the crude product bound is 1.11αx² at ε = 1/4 (true for ε ≤ 1/5). | Both claims are true (numerically den/(αx²) ≤ 0.67 at ε = 1/4, and the correct argument via h(x)/x increasing gives D_ε ≤ 2.1/α), so Theorem A is unaffected; the proof text was sloppy. | Wording corrected, ε ≤ 1/5 stated for the crude proof, numerics noted. |

Errors found and fixed earlier in this work, for completeness: a thinning-loop bug in some Monte Carlo scripts (first versions of `p2_identifiability.py`, `p3c_exp_tail.py`, `p10b_transient.py`), a non-causal sampling artifact in `p14_mech.py` and `p14_k2.py`, a wrong noise coefficient in a random-walk prediction, and the retracted "c_lag ≈ 1.8" claim in Problem 14.

## Independent re-checks that passed

- **Problem 1.** Counterexample and Lemma 1 (re-derived: for b = c·a the exp-concavity constant is 1/(cλ − 1)²); Lemma 2 follows from λ_θ/λ* ≤ max ratios (one line); Lemma 3 is Cauchy–Schwarz; the K = 1 table is reproduced from the exact formulas (`p1_check.py`).
- **Problem 3.** On 30 random parameter sets with K = 1..4: E Z² bound holds in all, and γ = ε/ν̃ holds to 1e-8 in all. Lemma 4's pathwise inequality A e^{θZ} ≤ θ e^{θZ}(μ − εbZ/2) holds on 150 000 random states (minimum slack 1.2e-4). Moment recursion p ≤ 4 reproduced against exact moments.
- **Problem 10.** The efficient-information formula agrees with the inverse of the exact K = 1 information matrix to 9e-16 on 8 random settings; the Jensen floor m ≥ ε/μ and the Theorem C bound hold on 25 and 23 random settings (max m/bound = 0.997); E λ_t ≤ μ(1 + β_max t) holds in 240 of 240 random cold-start checks for K = 1..4; the k = 1 slope is 2.0000 down to ε = 1e-8.
- **Problem 11.** Lemma H (corrected Lyapunov function) holds pathwise in 800 000 random states over 40 parameter sets, K = 1..4, ε between 0.005 and 0.5, after switching to a numerically stable evaluation. **Note:** a first, naive evaluation reported 570 "violations"; these were floating-point cancellation at very large states (the generator is a difference of O(10⁶) terms), not real violations, which is why the stable form is the one to trust.
- **Problem 12.** The claim that the top eigenvector is ∝ (0, 1/β_k) had only been inferred from an overlap with α̂; a direct check now gives |cos| = 0.9996 for K = 2 and for K = 3.
- **Problem 14.** Scripts reproduce; the time-varying-gain filter, transfer-function and decomposition results are as reported; the referenced code files all exist.
- All scripts referenced by the Done notes exist and run.

## Cross-session consistency

`problem03_moments_gap.md` (from another session) proves the contraction rate is exactly γ using the Perron weights; my `problem03_moments_mixing.md` proves rate ≥ εβ_min with weights α_k/β_k. They agree (γ = ε/ν̃ ≥ εβ_min since ν̃ ≤ 1/β_min); mine is the weaker statement. The two notes should be merged. The Perron-weight Lyapunov variable also reappears in Problem 11 (Lemma G).

## What this audit does not establish

- The proofs (Problem 1 lemmas, Problem 3 Theorem 2 and Lemmas 1'/4, Problem 10 Theorem A and the transient formula, Problem 11 Theorems B, E, F and Lemma H, Problem 12's rigorous parts) were checked by re-derivation and numerics only, not by an independent mathematician or a proof assistant. The notes ask you to re-derive the key steps yourself; that remains necessary.
- Heuristic or numerical-only parts, as flagged in the notes: Problem 10 cold-start coefficients; Problem 11 matching upper bound, K ≥ 2 divergence of the α-unknown information, k_K = 1; Problem 12 middle-tier asymptotics and eigenvector convergence; Problem 14 mechanism of the amplitude effect, K ≥ 4 and ε < 0.4 for K ≥ 2.
- Solo (ICASSP 2018) has not been read (closed access).
- "Done" means a draft with a stated scope, not a finished theorem in every case; Problem 11 in particular has no matching upper bound.
