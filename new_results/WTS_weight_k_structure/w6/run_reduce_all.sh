#!/bin/zsh
# sequential, one process at a time, nice 10
cd "<project>/research/08_theory/WTS_weight_k_structure/w6"
PY="<project>/.venv/bin/python"
for k in "$@"; do
  nice -n 10 "$PY" s01_reduce.py $k >> reduce_log.txt 2>&1
done
echo "reduce chain done: $@" >> reduce_log.txt
