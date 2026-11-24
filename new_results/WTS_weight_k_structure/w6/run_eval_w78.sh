#!/bin/zsh
# bonus: W_7 (d=6) for k=7..15 and W_8 (d=7) for k=8..15, dps 150, one process at a time, nice 10
cd "<project>/research/08_theory/WTS_weight_k_structure/w6"
PY="<project>/.venv/bin/python"
for k in 7 8 9 10 11 12 13 14 15; do
  nice -n 10 "$PY" s02_eval.py $k 6 --dps 150 >> eval_log_w78.txt 2>&1
done
for k in 8 9 10 11 12 13 14 15; do
  nice -n 10 "$PY" s02_eval.py $k 7 --dps 150 >> eval_log_w78.txt 2>&1
done
echo "w78 chain done" >> eval_log_w78.txt
