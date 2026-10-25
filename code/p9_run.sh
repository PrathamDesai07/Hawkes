#!/bin/bash
# Problem 9 sweep: Bayes regret of Newton over the stability band, L_emp, I_G benchmark; two-point checks.
cd "$(dirname "$0")"; mkdir -p p9_out
./p9_bayes.exe mode=bayes eps=0.3  mu=1 reps=300 tmax=3000  tmin=3  seed=100 > p9_out/b_mu1_eps0.3.txt
./p9_bayes.exe mode=bayes eps=0.1  mu=1 reps=300 tmax=3000  tmin=3  seed=200 > p9_out/b_mu1_eps0.1.txt
./p9_bayes.exe mode=bayes eps=0.03 mu=1 reps=200 tmax=10000 tmin=3  seed=300 > p9_out/b_mu1_eps0.03.txt
./p9_bayes.exe mode=bayes eps=0.01 mu=1 reps=120 tmax=10000 tmin=3  seed=400 > p9_out/b_mu1_eps0.01.txt
./p9_bayes.exe mode=bayes eps=0.003 mu=1 reps=60 tmax=10000 tmin=3  seed=500 > p9_out/b_mu1_eps0.003.txt
for mu in 0.3 0.6 2; do for eps in 0.1 0.03; do
 ./p9_bayes.exe mode=bayes eps=$eps mu=$mu reps=200 tmax=3000 tmin=3 seed=600 > p9_out/b_mu${mu}_eps${eps}.txt
done; done
for eps in 0.3 0.1 0.03; do for d in 0.05 0.15; do
 ./p9_bayes.exe mode=twopoint eps=$eps mu=1 reps=64 tmax=2000 delta=$d seed=700 > p9_out/tp_eps${eps}_d${d}.txt
done; done
echo done > p9_out/DONE
