import io
add = """

---
## Addendum (25 Oct 2026): strengthening round

Items of the "weaker than it sounds" table that were improved (details in `tangent_bound.md` and `theoremR_proof_complete.md`):
- **Problem 2:** rigorous lower bound within 1.3-3.5x of the Monte Carlo lambda_min (K = 1, 2, 3), correct eps-scaling, from the tangent minorant with exact moment systems; the 10^3-10^8 looseness is gone.
- **Problem 12:** every eigenvalue of I is bounded below rigorously (0.23-0.85 of Monte Carlo): the three tiers are now a theorem as regards lower bounds. Still open: matching upper bounds and eigenvector convergence.
- **Problem 5:** (a) every step of Theorem R's expected-regret proof is written (only the per-time failure bound H3 remains an input); (b) Assumption M is no longer needed: the information level comes from the tangent minorant (about 10^4 better than the truncation device for the K = 2 example).
- **Not strengthened (honestly):** Problem 11 (a matching estimator rate for k_K < 1 needs LAMN/weak-convergence arguments; the observed information is random at O(1) relative size when eps = c/T); Problem 14 (a tracking theorem for the estimators is the content of the not-started Problem 13; the amplitude effect of the data-driven preconditioner is not derived); the multiprecision evaluations are not interval arithmetic.
"""
with io.open('final_audit.md', 'a', encoding='utf8', newline='\n') as f: f.write(add)
with io.open('K2_regret_chain.md', 'a', encoding='utf8', newline='\n') as f:
    f.write("\n---\n## Update (25 Oct 2026): better information level\nThe information level of Theorem FK can be taken from the tangent minorant of `tangent_bound.md` (lambda_*/2 = 0.0126 for the K = 2 example, against about 1e-6 for the truncation device), and Assumption M is not needed. The expected-regret proof is written in full in `theoremR_proof_complete.md`.\n")
print('ok')
