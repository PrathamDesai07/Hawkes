#!/bin/bash
# Problem 21: accuracy on identical simulated paths (all methods), plus timing.
cd "$(dirname "$0")"; mkdir -p p21_out
for eps in 0.3 0.1 0.03; do
  ./p21_baselines.exe mode=acc K=1 mu=1 beta=1 eps=$eps n=100000 reps=64 seed=21000 > p21_out/acc_k1_eps${eps}.csv
done
./p21_baselines.exe mode=acc K=2 mu=1.5 beta=1,2 eps=0.1 n=100000 reps=32 seed=22000 > p21_out/acc_k2_eps0.1.csv
./p21_baselines.exe mode=time K=1 eps=0.1 n=200000 > p21_out/time_k1.txt
./p21_baselines.exe mode=time K=2 mu=1.5 beta=1,2 eps=0.1 n=200000 > p21_out/time_k2.txt
echo done > p21_out/DONE
