#!/bin/zsh
# w6 identification step: independent high-precision reruns of W_6(k), k=6..15 (d=5), one process at a time, nice 10
# usage: run_eval_hi.sh DPS
cd "<project>/research/08_theory/WTS_weight_k_structure/w6"
PY="<project>/.venv/bin/python"
DPS=$1
for k in 6 7 8 9 10 11 12 13 14 15; do
  nice -n 10 "$PY" s02_eval.py $k 5 --dps $DPS >> eval_log_$DPS.txt 2>&1
done
echo "eval$DPS chain done" >> eval_log_$DPS.txt
