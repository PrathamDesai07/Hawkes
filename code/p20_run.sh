#!/bin/bash
# Problem 20: fits on the Deribit BTC option-trade stream (events.bin from p20_prep.py).
cd "$(dirname "$0")"; mkdir -p p20_out
./p20_fit.exe mode=select tr0=0 tr1=30 va0=30 va1=60 "dicts=0.03;0.3;3;30;0.1,1;1,10;0.3,3;3,30;0.1,1,10;0.03,0.3,3;0.3,3,30;0.1,1,10,100" > p20_out/select.txt
./p20_fit.exe mode=daily d0=0 d1=365 "dicts=1;0.3,3;0.1,1,10;0.1,1,10,100" > p20_out/daily_13_20.txt
./p20_fit.exe mode=daily d0=0 d1=365 h0=0 h1=24 "dicts=0.1,1,10" > p20_out/daily_0_24.txt
./p20_fit.exe mode=online d0=0 d1=365 etamult=5 "dicts=0.1,1,10" > p20_out/online_eta5.txt
./p20_fit.exe mode=online d0=0 d1=365 etamult=1 "dicts=0.1,1,10" > p20_out/online_eta1.txt
./p20_fit.exe mode=online d0=0 d1=365 etamult=0 "dicts=0.1,1,10" > p20_out/online_growing.txt
./p20_fit.exe mode=rolling d0=339 d1=349 wlen=1 step=1 "dicts=0.1,1,10" > p20_out/rolling_oct2025.txt
./p20_fit.exe mode=rolling d0=318 d1=328 wlen=1 step=1 "dicts=0.1,1,10" > p20_out/rolling_calm.txt
./p20_fit.exe mode=rolling d0=119 d1=126 wlen=1 step=1 "dicts=0.1,1,10" > p20_out/rolling_mar2025.txt
echo done > p20_out/DONE
# pooled batch MLE over the first 30 / 120 / 365 days (for the growing-window comparison)
for r in "0 30" "0 120" "0 365"; do set -- $r; ./p20_fit.exe mode=select tr0=$1 tr1=$2 va0=$2 va1=$(( $2 + 1 > 365 ? 365 : $2 + 1 )) "dicts=0.1,1,10" | tail -1 | sed "s/^/pooled days [$1,$2): /"; done > p20_out/pooled.txt
