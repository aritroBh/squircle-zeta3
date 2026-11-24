#!/bin/zsh
# sequential, niced high-precision runs of W_4(k); resumable (w4_eval.py skips existing json)
cd "<project>/research/08_theory/WTS_weight_k_structure/w4/"
PY="<project>/.venv/bin/python"
for dps in 300 360; do
  for k in 4 5 6 7 8 3; do
    nice -n 10 $PY w4_eval.py 4 $k $dps > results/log_w4_k${k}_dps${dps}.txt 2>&1
  done
done
echo done > results/id_run_hp.done
