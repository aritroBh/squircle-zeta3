#!/bin/zsh
# Planner chain: scan the remaining D=21 values, collect every bad index (error > 1e-1150), rerun the LLL on fixed values.
cd "<project>/research/08_theory/ZQ5_zagier_weight5/planner"
PY="<project>/.venv/bin/python"
nice -n 10 $PY planner_check_d21_vals3.py $(cat d21_remaining.txt) > planner_check_d21_vals3.out 2>&1
BAD=$(cat planner_check_d21_vals2.out planner_check_d21_vals3.out | $PY -c "
import sys,re
b=[]
for l in sys.stdin:
    m=re.match(r'\s*(\d+) x=\S+ \|mine-stored\|=(\S+)',l)
    if m and 'e' in m.group(2) and int(m.group(2).split('e')[1])>-1150: b.append(m.group(1))
print(','.join(b))")
echo "bad indices: $BAD" > planner_d21_rerun_fixed.out
nice -n 10 $PY planner_d21_rerun_fixed.py "$BAD" 830 1000 >> planner_d21_rerun_fixed.out 2>&1
echo CHAIN_DONE >> planner_d21_rerun_fixed.out
