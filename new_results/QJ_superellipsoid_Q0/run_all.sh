#!/bin/sh
# Runs the whole QJ job in order (about 1 minute; 09 is the slow direct check of Q0).
PY="<project>/.venv/bin/python"
cd "$(dirname "$0")"
$PY 00_reference.py > 00_reference.log 2>&1
$PY 01_J.py   > 01_J.log 2>&1
$PY 02_P2c.py > 02_P2c.log 2>&1
$PY 03_P2a.py > 03_P2a.log 2>&1
$PY 04_P1.py  > 04_P1.log 2>&1
$PY 09_Q0_direct.py 75 k > 09_Q0_k3.log 2>&1
$PY 09_Q0_direct.py 65 s > 09_Q0_simplex.log 2>&1
$PY 05_reduce.py > 05_reduce.log 2>&1
if grep -l "FAIL\|Traceback" *.log; then echo "PROBLEM in the logs above"; else echo "all logs clean"; fi
