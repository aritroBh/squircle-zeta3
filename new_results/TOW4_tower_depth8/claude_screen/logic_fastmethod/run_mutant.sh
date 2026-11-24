#!/bin/bash
# usage: run_mutant.sh MUTNAME DELTA S MUTKEY DROPC
cd "<project>/research/08_theory/TOW4_tower_depth8/claude_screen/logic_fastmethod"
PY="<project>/.venv/bin/python"
export PYTHONPATH="$PWD/lib"
export MUTKEY="$4"; export DROPC="$5"
out=$(nice -n 10 timeout 900 "$PY" -B mut/$1.py $2 $3 2>&1)
if echo "$out" | grep -q "exact numerator == 0 : True"; then
  echo -e "$1\t($2,$3)\tUNDETECTED (printed True)" 
elif echo "$out" | grep -q "exact numerator == 0 : False"; then
  echo -e "$1\t($2,$3)\tdetected: numerator != 0"
else
  echo -e "$1\t($2,$3)\tdetected: $(echo "$out" | tail -1 | cut -c1-120)"
fi
