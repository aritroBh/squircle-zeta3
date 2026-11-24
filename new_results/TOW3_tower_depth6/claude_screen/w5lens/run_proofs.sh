#!/bin/zsh
PY='<project>/.venv/bin/python'
cd "<project>/research/08_theory/TOW3_tower_depth6/claude_screen/w5lens"
for a in 5 6; do
  /usr/bin/time -p nice -n 10 $PY -B prove_E3b56.py $a > out_E3_$a.txt 2>&1
done
for ds in "5 0" "5 1" "5 2" "6 0" "6 1" "6 2" "6 3"; do
  set -- ${=ds}
  /usr/bin/time -p nice -n 10 $PY -B prove_V56.py $1 $2 > out_V_$1_$2.txt 2>&1
done
for ds in "5 2" "5 1" "5 0" "6 3" "6 2" "6 1" "6 0"; do
  set -- ${=ds}
  /usr/bin/time -p nice -n 10 $PY -B prove_adj56.py $1 $2 > out_adj_$1_$2.txt 2>&1
done
echo done > out_ALLDONE.txt
