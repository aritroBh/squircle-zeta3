#!/bin/bash
PY='<project>/.venv/bin/python'
cd "$(dirname "$0")"
for job in "prove_tops.py" "prove_E3b34.py" "prove_S12.py" "prove_S34.py" "prove_adj34.py 3 0" "prove_adj34.py 3 1" "prove_adj34.py 4 2" "prove_adj34.py 4 1" "prove_adj34.py 4 0"; do
  echo "##### $job  start $(date +%T)"
  /usr/bin/time -p nice -n 10 "$PY" -B $job 2>&1 | tail -25
  echo "##### $job  end $(date +%T)"
done
