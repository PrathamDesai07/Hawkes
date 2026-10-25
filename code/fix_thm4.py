import io
p = 'problem14_tracking_theorem.md'
s = io.open(p, encoding='utf8').read()
a = s.index(" mean_i E‖θ_{i+1} − θ*_{i+1}‖²_{I_ev} ≤ (2/c_F²)")
b = s.index("and in particular, for the idealised window")
new = (" mean_i E‖θ_{i+1} − θ*_{i+1}‖²_{I_ev} ≤ (2/c_F²)·(1+η/2)(1+η′)/(1−B)·ηd + (2/c_F)(1+η/2)·mean_i EΞ_i + (2C_F/c_F)(1+η/2)(1+2/η)·mean_i‖d_i‖²_{I_ev}/η + (transient),\n"
       "where the second-to-last coefficient (1+2/η)/η ≈ 2/η²,\n")
s = s[:a] + new + s[b:]
io.open(p, 'w', encoding='utf8', newline='\n').write(s)
print('ok')
