#!/bin/bash
# Problem 19 sweep. Output: p19_out/*.csv (raw per-replicate rows, see p19_sim.cs).
cd "$(dirname "$0")"; mkdir -p p19_out
./p19_sim.exe K=1 mu=1 beta=1 eps=0.3  reps=200 tmax=30000 start=near seed=1000 > p19_out/k1_eps0.3_near.csv
./p19_sim.exe K=1 mu=1 beta=1 eps=0.1  reps=200 tmax=30000 start=near seed=2000 > p19_out/k1_eps0.1_near.csv
./p19_sim.exe K=1 mu=1 beta=1 eps=0.03 reps=160 tmax=30000 start=near seed=3000 > p19_out/k1_eps0.03_near.csv
./p19_sim.exe K=1 mu=1 beta=1 eps=0.01 reps=96  tmax=20000 start=near seed=4000 > p19_out/k1_eps0.01_near.csv
./p19_sim.exe K=1 mu=1 beta=1 eps=0.1  reps=160 tmax=30000 start=far  seed=5000 > p19_out/k1_eps0.1_far.csv
./p19_sim.exe K=1 mu=1 beta=1 eps=0.03 reps=96  tmax=30000 start=far  seed=6000 > p19_out/k1_eps0.03_far.csv
./p19_sim.exe K=1 mu=1 beta=1 eps=0.01 reps=64  tmax=20000 start=far  seed=7000 > p19_out/k1_eps0.01_far.csv
./p19_sim.exe K=2 mu=1.5 beta=1,2 eps=0.3  reps=96 tmax=30000 start=near seed=8000 > p19_out/k2_eps0.3_near.csv
./p19_sim.exe K=2 mu=1.5 beta=1,2 eps=0.1  reps=96 tmax=30000 start=near seed=9000 > p19_out/k2_eps0.1_near.csv
./p19_sim.exe K=2 mu=1.5 beta=1,2 eps=0.03 reps=64 tmax=20000 start=near seed=9500 > p19_out/k2_eps0.03_near.csv
echo done > p19_out/DONE
