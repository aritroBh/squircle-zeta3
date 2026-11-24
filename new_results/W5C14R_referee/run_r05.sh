#!/bin/sh
PY="<project>/.venv/bin/python"
export OMP_NUM_THREADS=4 VECLIB_MAXIMUM_THREADS=4 OPENBLAS_NUM_THREADS=4
for a in "4194301 +" "4194301 -" "4194287 +" "4194287 -"; do set -- $a; "$PY" r05_span.py $1 $2 > r05_$1_$2.out 2>&1 || { echo FAIL $a; exit 1; }; done
echo ALLDONE
