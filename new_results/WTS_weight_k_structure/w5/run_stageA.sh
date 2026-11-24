#!/bin/bash
cd "<project>/research/08_theory/WTS_weight_k_structure/w5"
PY="<project>/.venv/bin/python"
# Stage A: no Li2(b) atoms, Li4/Li5 only at 1/2 (+zeta5); with sqrt2 multiples; D=3000
nice -n 10 $PY id7_big.py results/cm_D3300_m11.json 3000 c2 noli2 li4half li5half > id7_c2_D3000_stageA.log 2>&1
nice -n 10 $PY id7_big.py results/cm_D3300_m11.json 3000 c1 noli2 li4half li5half > id7_c1_D3000_stageA.log 2>&1
