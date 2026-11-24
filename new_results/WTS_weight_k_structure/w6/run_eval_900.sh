#!/bin/zsh
# dps 900 run of W_6(k), k=6..11 (5 points for the quartic + 1 check), one process at a time, nice 10
cd "<project>/research/08_theory/WTS_weight_k_structure/w6"
PY="<project>/.venv/bin/python"
for k in 6 7 8 9 10 11; do
  nice -n 10 "$PY" s02_eval.py $k 5 --dps 900 >> eval_log_900.txt 2>&1
done
echo "eval900 chain done" >> eval_log_900.txt
