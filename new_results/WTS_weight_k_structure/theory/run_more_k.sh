#!/bin/zsh
cd "<project>/research/08_theory/WTS_weight_k_structure/theory"
PY="<project>/.venv/bin/python"
for k in 16 17 18; do
  if [ ! -f kern_$k.pkl ]; then
    nice -n 10 timeout 1200 $PY -B run_kernel_k.py $k >> kernels_more.log 2>&1 || echo "k=$k FAILED/TIMEOUT" >> kernels_more.log
  fi
done
echo ALLDONE >> kernels_more.log
