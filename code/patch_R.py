import io

def sub(path, a, b):
    s = io.open(path, encoding='utf8').read()
    assert a in s, (path, a[:70])
    io.open(path, 'w', encoding='utf8', newline='\n').write(s.replace(a, b, 1))

sub('theoremR_proof_complete.md',
    "E log det(F_n/F_0) ≤ d ln T_n + log det(I(θ\\*)/c⁻) + C_4.",
    "E log det(F_n/F_0) ≤ d ln T + log det(I(θ\\*)/c⁻) + C_4, where the horizon is a deterministic time T and F_n is the matrix after the last gap completed before T (so E J = E∫₀^T aaᵀ/λ\\* dt, no random stopping time enters; n = N_T is then random, and the statement is E R_T ≤ …, which is the same statement as for n events up to the O(1) boundary gap).")
print('ok')
