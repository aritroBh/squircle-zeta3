#!/bin/zsh
# weight-4 and weight-5 data at dps 400 for the lower-weight template tests: W_4 (d=3) k=4..8, W_5 (d=4) k=5..10; one process at a time, nice 10
cd "<project>/research/08_theory/WTS_weight_k_structure/w6"
PY="<project>/.venv/bin/python"
for k in 4 5 6 7 8; do
  nice -n 10 "$PY" s02_eval.py $k 3 --dps 400 >> eval_log_w4_400.txt 2>&1
done
for k in 5 6 7 8 9 10; do
  nice -n 10 "$PY" s02_eval.py $k 4 --dps 400 >> eval_log_w4_400.txt 2>&1
done
echo "w4/w5 400 chain done" >> eval_log_w4_400.txt
