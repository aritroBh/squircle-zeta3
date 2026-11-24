#!/bin/zsh
PY="<project>/.venv/bin/python"
cd "<project>/research/08_theory/TOW4_tower_depth8/claude_screen/w6lens"
for dl in 7 8; do nice -n 10 $PY -B prove_E3b78.py $dl; done
for ds in "7 3" "8 4" "8 3" "7 2" "8 2"; do echo "adj $ds"; nice -n 10 $PY -B prove_adj78_fast.py $ds; done
