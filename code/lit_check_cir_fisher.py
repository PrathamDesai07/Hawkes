"""Cross-check against Ben Alaya, Kebaier, Tran (arXiv 1708.07070) Theorem 2.1 (subcritical CIR, LAN): for dX = (a - bX)dt + sqrt(2 sigma X) dB (Feller: a >= sigma),
   I(a,b) = (1/(2 sigma)) [[ b/(a - sigma), -1 ], [ -1, a/b ]]   (information per unit time).
 Map from K=1 Hawkes (diffusion approximation  d lam = (beta mu - beta eps lam) dt + beta sqrt(lam) dB):  a = beta*mu, b = beta*eps, 2 sigma = beta^2.
 My exact Hawkes Fisher matrix in (mu, alpha) -> (a, b) = (beta mu, beta - alpha):  I_aa = I_mm/beta^2, I_ab = -I_ma/beta, I_bb = I_aa(alpha).
 Also the critical-case coefficient:  Cor 2.3(i): I(a0) = 1/(2 sigma (a0 - sigma)) per unit log-time  ->  Fisher info for mu = beta^2 I(a) = 2/(beta (2 mu - beta)) = the Problem-10 cold-start coefficient A."""
import math
from p10_fisher_k1 import exact_m
beta = 1.0
for mu in (1.0, 0.8):
    sig = beta**2/2; a = beta*mu
    print(f"mu={mu} (k0={2*mu/beta}), beta=1:  CIR (a,b,sigma) = ({a}, beta*eps, {sig}); Feller a>=sigma: {a >= sig}")
    for eps in (0.1, 0.01, 0.001):
        al = (1 - eps)*beta; m = exact_m(mu, al, beta); El = mu/eps
        Imm = m; Ima = (1 - mu*m)/al; Iaa = (El - 2*mu + mu*mu*m)/al**2
        H = (Imm/beta**2, -Ima/beta, Iaa); b = beta*eps
        C = (b/(2*sig*(a - sig)), -1/(2*sig), a/(2*sig*b))
        rel = [abs(h - c)/abs(c) for h, c in zip(H, C)]
        print(f"   eps={eps:6g}: Hawkes (I_aa, I_ab, I_bb) = ({H[0]:.5g}, {H[1]:.5g}, {H[2]:.5g});  CIR Thm 2.1 = ({C[0]:.5g}, {C[1]:.5g}, {C[2]:.5g});  rel.diff = ({rel[0]:.2e}, {rel[1]:.2e}, {rel[2]:.2e})")
    print(f"   critical-case coefficient: beta^2/(2 sigma (a - sigma)) = {beta**2/(2*sig*(a - sig)):.4f}   vs  2/(beta(2mu-beta)) = {2/(beta*(2*mu - beta)):.4f}")
