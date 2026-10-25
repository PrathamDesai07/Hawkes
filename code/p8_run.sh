#!/bin/bash
# Problem 8 sweep: regret notions vs parameter-error notions (K=1 over eps and start; K=2).
cd "$(dirname "$0")"; mkdir -p p8_out
for eps in 0.3 0.1 0.03; do for st in near truth; do
  ./p8_link.exe K=1 mu=1 beta=1 eps=$eps reps=64 n=20000 start=$st seed=12000 > p8_out/k1_eps${eps}_${st}.txt
done; done
for eps in 0.3 0.1; do
  ./p8_link.exe K=2 mu=1.5 beta=1,2 eps=$eps reps=48 n=20000 start=near seed=13000 > p8_out/k2_eps${eps}_near.txt
done
echo done > p8_out/DONE
