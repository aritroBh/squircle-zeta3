#!/bin/bash
cd "<project>/research/08_theory/TOW4_tower_depth8/claude_screen/logic_fastmethod"
PY="<project>/.venv/bin/python"
export PYTHONPATH="$PWD/lib"
out=$(nice -n 10 timeout 900 "$PY" -B emut/$1.py $2 2>&1)
if echo "$out" | grep -q "numerator == 0 : True"; then v="UNDETECTED(True)"; elif echo "$out" | grep -q "numerator == 0 : False"; then v="detected(num!=0)"; else v="detected(ERR $(echo "$out" | tail -1 | cut -c1-70))"; fi
echo -e "$1\tdelta=$2\t$v"
