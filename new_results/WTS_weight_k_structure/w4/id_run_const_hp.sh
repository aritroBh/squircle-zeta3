#!/bin/zsh
cd "<project>/research/08_theory/WTS_weight_k_structure/w4/"
PY="<project>/.venv/bin/python"
for dps in 600 560 800; do
  nice -n 10 $PY id_const.py $dps > results/log_idconst_dps${dps}.txt 2>&1
done
echo done > results/id_run_const_hp.done
