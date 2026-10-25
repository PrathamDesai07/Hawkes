#!/bin/bash
# Problem 6 sweep driver: K=1 over eps, start in {far,truth}, prior in {scaled,fixed}; K=2; mu-dependence.
cd "$(dirname "$0")"; mkdir -p p6_out
for eps in 0.3 0.1 0.03 0.01; do
 for start in far truth; do for prior in scaled fixed; do
  ./p6_scaling.exe K=1 mu=1 beta=1 eps=$eps reps=32 tmax=3000 start=$start prior=$prior rho0=1 > p6_out/k1_eps${eps}_${start}_${prior}.txt
 done; done
done
for mu in 0.6 2; do for eps in 0.1 0.03; do
  ./p6_scaling.exe K=1 mu=$mu beta=1 eps=$eps reps=32 tmax=3000 start=truth prior=fixed rho0=1 > p6_out/mu${mu}_eps${eps}_truth_fixed.txt
done; done
for eps in 0.3 0.1 0.03; do for prior in scaled fixed; do
  ./p6_scaling.exe K=2 mu=1.5 beta=1,2 eps=$eps reps=24 tmax=3000 start=near prior=$prior rho0=1 > p6_out/k2_eps${eps}_near_${prior}.txt
done; done
echo done > p6_out/DONE
