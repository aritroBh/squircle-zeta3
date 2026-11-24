#!/bin/bash
cd "<project>/research/08_theory/WTS_weight_k_structure/w5"
PY="<project>/.venv/bin/python"
while kill -0 57814 2>/dev/null; do sleep 10; done
# Stage A+: no Li2(b) atoms, Li4 only at 1/2, Li5 at all 8 points (+zeta5); with sqrt2 multiples
nice -n 10 $PY id7_big.py results/cm_D3300_m11.json 3250 c2 noli2 li4half > id7_c2_D3250_stageAplus.log 2>&1
