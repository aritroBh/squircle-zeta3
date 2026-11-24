#!/bin/bash
cd "<project>/research/08_theory/CHN5_p5_w12_exact"
PY=<project>/.venv/bin/python
for start in 170 340 510 680 850; do
  cnt=170
  if [ $start -eq 850 ]; then cnt=74; fi
  nice -n 10 "<project>/.venv/bin/python" run_eval12.py $start $cnt >> run_all.log 2>&1
done
echo "ALL CHUNKS DONE" >> run_all.log
