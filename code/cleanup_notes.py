"""One-off cleanup of stale statements in problem14_drift_model.md and a sentence in box_cost_and_lower_bound.md."""
import io

def edit(path, pairs):
    s = io.open(path, encoding='utf8').read()
    for a, b in pairs:
        assert a in s, (path, a[:70])
        s = s.replace(a, b)
    io.open(path, 'w', encoding='utf8', newline='\n').write(s)

edit('problem14_drift_model.md', [
 # header: conventions + units
 ("Target form for Problem 13 (average tracking MSE of a constant-gain estimator with gain η):\n**η σ²(ε) + Δ²/η² + J D²/(η T)**, with σ²(ε) = r/(2·eff) in event units",
  "Target form for Problem 13 (average tracking MSE of a constant-gain estimator with gain η):\n**η σ²(ε) + Δ²/η² + J D²/(η T)** (the jump term is in event units and up to the factor ½ of the tanh formula below: in time units it is ρD²/(2ηr), ρ = J/T per unit time), with σ²(ε) = r/(2·eff) in event units"),
 # path-length constants
 ("avg bias² ≤ min(Δ_max V_T/(η²T), V_T²/(2ηT)).",
  "avg bias² ≤ min(Δ_max V_T/(η²T), V_T²/(ηT)) (the second bound is what the Young/ℓ¹–ℓ^∞ argument gives; an earlier version of this note had V_T²/(2ηT), whose factor ½ is not derived)."),
 ("and ≤ V_T²/(2ηT).", "and ≤ V_T²/(ηT) (Σe² ≤ ‖e‖_∞‖e‖₁ ≤ V_T·V_T/η)."),
 # noise convention
 ("- Noise: stationary variance ≈ η σ²/2.",
  "- Noise: stationary variance ≈ η σ_ξ²/2 for noise variance σ_ξ² per step. **Convention:** the coefficient called σ²(ε) in the target form is a(ε) = MSE/η = σ_ξ²/2 = r/(2·eff), i.e. σ_ξ² = r/eff."),
 # F2: retract stale claims
 ("- The excess for α unknown is **about 1.6–2.0, the same at K = 1 and K = 2, and independent of ε, η and of the amplitude** (ε = 0.4, η = 0.004: ratios 1.76 at A = 0.5 and 1.82 at A = 0.2 with 64 replicates, split-half 0.0075 / 0.0079). So it is not nonlinearity in A.",
  "- ~~The excess for α unknown is about 1.6–2.0, the same at K = 1 and K = 2, and independent of ε, η and of the amplitude.~~ **RETRACTED (see R0, R1):** the mean-square ratio is too noisy a statistic (the same algorithm at the same settings gave 1.2–1.9 across seed sets), the K = 2 numbers came from the peek-ahead script, and the decomposition in R1 shows an amplitude-dependent excess (adaptive/fixed = 1.32, 1.48, 1.75 at A = 0.2, 0.35, 0.5) that is absent with a fixed preconditioner. The hypotheses (i) and (ii) rejected below remain rejected."),
 ("**I did not identify the mechanism**; empirically the lag term should carry a constant c_lag ≈ 1.8 ± 0.2 when α is estimated (c_lag ≈ 1.2 when α is known): Δ²/η² → c_lag·Δ²/η² in the C5 bound, until the cause is understood.",
  "~~I did not identify the mechanism; empirically the lag term should carry a constant c_lag ≈ 1.8 ± 0.2 ...~~ **RETRACTED (R0, R1):** there is no constant c_lag. The lag term is the response of the time-varying-gain filter (R2), plus an amplitude-dependent excess of the data-driven preconditioner (R1; not derived analytically)."),
 # F3: superseded bias ratios
 ("Bias ratios to the filter formula: 1.64 and 1.87 (table above). So the K = 1 conclusions (noise coefficient r/(2 eff), bias excess ≈ 1.6–2) carry over to K = 2 at the single (β, ε, μ) tested.",
  "Bias ratios to the filter formula: 1.64 and 1.87 (table above) — **SUPERSEDED (R0, R3):** these came from the peek-ahead script `p14_k2.py`. With causal sampling the K = 2 response follows the filter to 6 % (β = (1,2), η = 0.002). The noise-coefficient part of this paragraph stands: r/(2 eff) carries over to K = 2."),
 # caveat-table row
 ("| K ≥ 2, α unknown not run | **done at one setting** (K = 2, β = (1,2), μ = 1.5, ε = 0.4): noise coefficient and the bias excess agree with K = 1; K ≥ 3, other grids, other ε not done |",
  "| K ≥ 2, α unknown not run | **done** (see R3): noise coefficient r/(2 eff) agrees with K = 1 (3–13 %); the earlier bias-excess statement is superseded; K = 3 on (1,2,4) and (1,3,9) and K = 2 on (1,2), (1,4) at ε = 0.4 done, K ≥ 4 and ε < 0.4 not run |"),
 ("| still open | mechanism of c_lag; K ≥ 3 and other grids; drift shapes other than a sinusoid, jump class at other ε, k, K; Solo (paper not accessible) |",
  "| still open | analytic derivation of the amplitude effect (R1); K ≥ 4 and ε < 0.4 for K ≥ 2; Solo (paper not accessible) (the earlier items c_lag, K = 3, other drift shapes and the jump class are addressed in R1–R5) |"),
])

edit('box_cost_and_lower_bound.md', [
 ("the log T slope is unaffected (E R_n − ln T_n is constant in n for each s, from the table in the earlier circularity note).",
  "I measured E R_n at one horizon only (n = 30 000), so the table does not by itself show that the log T slope is unaffected; the theory says the slope (the leading constant) is box-free."),
])
print("ok")
