#!/bin/zsh
# one process at a time, each under 20 minutes, resumable
cd "<project>/research/08_theory/WTS_weight_k_structure/theory"
PY="<project>/.venv/bin/python"
for k in 2 3 4 5 6 7 8 9 10 11 12 13 14 15; do
  if [ ! -f kern_$k.pkl ]; then
    nice -n 10 timeout 1200 $PY -B run_kernel_k.py $k >> kernels.log 2>&1 || echo "k=$k FAILED/TIMEOUT" >> kernels.log
  fi
done
echo ALLDONE >> kernels.log
