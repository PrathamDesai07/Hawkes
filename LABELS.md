# Label map: resolving the clashing theorem names across the notes

Several notes reuse the same label for different results. Rather than renaming inside the notes (which would break their cross-references), use these unique names when citing. The same names are the ones used in `latex/hawkes_master.tex`.

| note | label used there | unique name | what it is |
|---|---|---|---|
| `problem01_curvature_lemma.md` | Lemma 1 / 2 / 3 | Curv-1 (rank-one Hessian), Curv-2 (expected curvature), Curv-3 (λ_min bound) | |
| `problem03_moments_mixing.md` | Lemma 1, 1′, 4; Theorem 2; Prop. 3 | Mom-1, Mom-p, ExpMom, Contr-W (rate εβ_min), Rate-γ | |
| `problem03_moments_gap.md` | Lemma A, B, C | Lyap-Zℓ, Mom-B (weaker constants than Mom-p), Contr-ℓ (exact rate γ) | |
| `problem10_rigorous.md` | Theorem A | **Fisher-asymptotics (K = 1)** | m_ε/ε → 2/(2μ−β) etc. |
| `problem11_dichotomy.md` | Theorem B | **Imp-S0 (K = 1)** | impossibility given the initial state |
| `problem11_dichotomy.md` | Theorem C | Lyap-C (superseded by Theorem E) | region 2μ > n²B only |
| `problem11_dichotomy.md` | Lemma G, H; Theorems E, F; Cor. D, E′ | unchanged | |
| `problem02_identifiability.md` | Theorem T (draft) | **Trunc-λmin** | truncation bound on λ_min(I) |
| `theoremA_step2_and_stationary_tv.md` | Theorem T, Lemma F | **TV-Gamma (K = 1)**, **TV-H1** | stationary-law comparison |
| `thmT_explicit_constants.md` | Theorem T′, Φ | **TV-Gamma-explicit** | explicit constants |
| `general_K_impossibility.md` | Theorem V, Lemma S | unchanged | any-K impossibility |
| `problem05_circularity.md` | Theorem B | **Box-K1** | ratio set and projection |
| `problem07_leverage.md` | Theorem 7 | **Cond-Reg** (corrected: L2 replaced by truncation) | |
| `problem05_fisher_concentration.md` | Lemmas F1, F2; Theorem F | **Occ-Var**, **Two-Int**, **Conc-J** | |
| `problem05_bootstrap.md` | Theorem C | **Boot-poly** (superseded) | |
| `problem05_exponential_tails.md` | E1–E4; E5 sketch; Theorem C-exp | **Boot-exp** | |
| `problem05_E5_expected_regret.md` | E5a–E5d, Prop. E5, Theorem R | **Tail-G**, **Reg-sharp** | |
| `K2_regret_chain.md` | Theorem FK, Theorem R^(K) | **Conc-J-K** (replaces Assumption M), **Reg-sharp-K** | |
| `box_cost_and_lower_bound.md` | Prop. A1, Prop. B | **Box-cost**, **Minimax-LB** | |

Rule of thumb: whenever a note says "Theorem B", "Theorem C" or "Theorem T" it needs the file name to be unambiguous.
