#!/bin/bash
# Problem 6, second sweep: potential defect and exit time T_exit(delta) (K=1 over eps, mu incl. k0<=1; K=2,3 long horizon).
cd "$(dirname "$0")"; mkdir -p p6_out2
for eps in 0.3 0.1 0.03 0.01; do for prior in scaled fixed; do
  ./p6_scaling.exe K=1 mu=1 beta=1 eps=$eps reps=32 tmax=6000 tmin=10 start=truth prior=$prior rho0=1 seed=2000 > p6_out2/k1_mu1_eps${eps}_${prior}.txt
done; done
for mu in 0.3 0.5 0.6 1.5 2 4; do for eps in 0.1 0.03; do
  ./p6_scaling.exe K=1 mu=$mu beta=1 eps=$eps reps=32 tmax=6000 tmin=10 start=truth prior=scaled rho0=1 seed=3000 > p6_out2/k1_mu${mu}_eps${eps}_scaled.txt
done; done
for eps in 0.3 0.1 0.03; do
  ./p6_scaling.exe K=2 mu=1.5 beta=1,2 eps=$eps reps=32 tmax=30000 tmin=10 start=truth prior=scaled rho0=1 seed=4000 > p6_out2/k2_eps${eps}_scaled.txt
  ./p6_scaling.exe K=3 mu=2 beta=1,2,4 eps=$eps reps=24 tmax=30000 tmin=10 start=truth prior=scaled rho0=1 seed=5000 > p6_out2/k3_eps${eps}_scaled.txt
done
echo done > p6_out2/DONE
