#!/bin/zsh
PY="<project>/.venv/bin/python"
cd "<project>/research/08_theory/TOW4_tower_depth8/claude_screen/w6lens"
for a in "7 3" "8 4" "8 3" "7 2" "8 2" "7 1" "8 1" "7 0" "8 0"; do
  set -- ${=a}
  echo "adj delta=$1 s=$2"
  /usr/bin/time -p nice -n 10 $PY -B prove_adj78_fast.py $1 $2
done
echo ALLDONE
