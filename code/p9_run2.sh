#!/bin/bash
# Problem 9: valid two-point configurations (alpha-perturbation delta <= 0.4 eps keeps theta' inside the stable set).
cd "$(dirname "$0")"; mkdir -p p9_out
for eps in 0.3 0.1 0.03; do for f in 0.4 0.1; do
 d=$(awk -v e=$eps -v f=$f 'BEGIN{printf "%.4f", e*f}')
 ./p9_bayes.exe mode=twopoint eps=$eps mu=1 reps=96 tmax=3000 delta=$d seed=800 > p9_out/tpv_eps${eps}_f${f}.txt
done; done
echo done > p9_out/DONE2
