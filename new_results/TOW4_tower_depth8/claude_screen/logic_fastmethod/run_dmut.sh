#!/bin/bash
# usage: run_dmut.sh TAG PATCH DELTA S
cd "<project>/research/08_theory/TOW4_tower_depth8/claude_screen/logic_fastmethod"
PY="<project>/.venv/bin/python"
export PYTHONPATH="$PWD/lib"
out=$(nice -n 10 timeout 900 "$PY" -B dmut/$1_$2.py $3 $4 2>&1)
if echo "$out" | grep -q ": True"; then v="TRUE"; elif echo "$out" | grep -q ": False"; then v="FALSE"; else v="ERR($(echo "$out" | tail -1 | cut -c1-60))"; fi
echo -e "$2\t($3,$4)\t$1\t$v"
