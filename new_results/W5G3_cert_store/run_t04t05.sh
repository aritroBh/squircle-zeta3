#!/bin/bash
cd "$(dirname "$0")"
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 VECLIB_MAXIMUM_THREADS=4 PYTHONDONTWRITEBYTECODE=1
PY="<project>/.venv/bin/python"
date +%s > t04_run2_start.epoch
RESUME=1 $PY t04_assemble.py > t04.out 2>&1 || { echo FAIL04 > run_fail.flag; exit 1; }
date +%s > t04_run2_end.epoch
rm -f t05.log t05.out
$PY t05_w4.py > t05.out 2>&1 || { echo FAIL05 > run_fail.flag; exit 1; }
date +%s > t05_end.epoch
echo done > run_done.flag
