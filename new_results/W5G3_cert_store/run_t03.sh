#!/bin/bash
cd "$(dirname "$0")"
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 VECLIB_MAXIMUM_THREADS=4 MKL_NUM_THREADS=4
PY="<project>/.venv/bin/python"
date +%s > t03_start.epoch
$PY t03_solve.py pp > t03_pp.out 2>&1 &
A=$!
$PY t03_solve.py pm > t03_pm.out 2>&1 &
B=$!
wait $A $B
$PY t03_solve.py mp > t03_mp.out 2>&1 &
A=$!
$PY t03_solve.py mm > t03_mm.out 2>&1 &
B=$!
wait $A $B
date +%s > t03_end.epoch
echo done > t03_done.flag
