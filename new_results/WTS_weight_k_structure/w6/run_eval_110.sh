#!/bin/zsh
# independent re-run at higher working precision (mp.mp.dps = 110), one process at a time
cd "<project>/research/08_theory/WTS_weight_k_structure/w6"
PY="<project>/.venv/bin/python"
nice -n 10 "$PY" s02_eval.py 5 5 4 3 2 1 --dps 110 >> eval_log_110.txt 2>&1
for k in 6 7 8 9 10 11 12 13 14 15; do
  nice -n 10 "$PY" s02_eval.py $k 5 --dps 110 >> eval_log_110.txt 2>&1
done
echo "eval110 chain done" >> eval_log_110.txt
