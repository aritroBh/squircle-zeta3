#!/bin/bash
cd "<project>/research/08_theory/WTS_weight_k_structure/w5"
PY="<project>/.venv/bin/python"
for f in F8 F8s F20; do nice -n 10 $PY id5_search.py $f 1200 c3 c2 c1 c0 delta > id5_$f.log 2>&1; done
nice -n 10 $PY id5_search.py F20s 1400 c3 c2 c1 c0 delta > id5_F20s.log 2>&1
nice -n 10 $PY id6_blockrel.py 1400 > id6_blockrel.log 2>&1
nice -n 10 $PY id9_diag.py 3000 nosqrt2 -- P0 Q0 A22_0 A31_0 A40_0 > id9_diag.log 2>&1
echo finished > run_logs.done
