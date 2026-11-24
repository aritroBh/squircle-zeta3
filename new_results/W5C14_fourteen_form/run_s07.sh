#!/bin/sh
PY="<project>/.venv/bin/python"
for s in pp pm mp mm; do OMP_NUM_THREADS=2 "$PY" s07_dixon.py $s > s07_$s.log 2>&1 || exit 1; done
echo ALLDONE
