#!/bin/zsh
# lower-weight guide data: W_4 (d=3) and W_5 (d=4) for k=4..15, dps 150, one process at a time, nice 10
cd "<project>/research/08_theory/WTS_weight_k_structure/w6"
PY="<project>/.venv/bin/python"
for k in 4 5 6 7 8 9 10 11 12 13 14 15; do
  nice -n 10 "$PY" s02_eval.py $k 3 4 --dps 150 >> eval_log_lower150.txt 2>&1
done
echo "lower chain done" >> eval_log_lower150.txt
