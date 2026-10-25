#!/bin/bash
# Problem 6 completion: stationary information I(theta*) vs eps for k0 <= 1 and K = 2, 3.
cd "$(dirname "$0")"; mkdir -p p6_out3
for mu in 0.3 0.5 0.8 1 2; do for eps in 0.3 0.1 0.03 0.01 0.003; do
  t=$(awk -v e=$eps -v m=$mu 'BEGIN{printf "%d", 4e6*e/m}')
  ./p6_moments.exe K=1 mu=$mu beta=1 eps=$eps tmax=$t reps=8 > p6_out3/I_k1_mu${mu}_eps${eps}.txt
done; done
for eps in 0.3 0.1 0.03 0.01 0.003; do
  t=$(awk -v e=$eps 'BEGIN{printf "%d", 3e6*e/1.5}')
  ./p6_moments.exe K=2 mu=1.5 beta=1,2 eps=$eps tmax=$t reps=8 > p6_out3/I_k2_mu1.5_eps${eps}.txt
  t=$(awk -v e=$eps 'BEGIN{printf "%d", 3e6*e/2}')
  ./p6_moments.exe K=3 mu=2 beta=1,2,4 eps=$eps tmax=$t reps=8 > p6_out3/I_k3_mu2_eps${eps}.txt
done
echo done > p6_out3/DONE3
