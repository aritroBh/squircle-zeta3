#!/bin/zsh
# sequential evaluation of W_{d+1}(k); one process at a time, nice 10
cd "<project>/research/08_theory/WTS_weight_k_structure/w6"
PY="<project>/.venv/bin/python"
for k in 5 6 7 8 9 10; do
  nice -n 10 "$PY" s02_eval.py $k 5 4 3 2 1 >> eval_log.txt 2>&1
done
for k in 3; do nice -n 10 "$PY" s02_eval.py $k 3 2 1 >> eval_log.txt 2>&1; done
for k in 4; do nice -n 10 "$PY" s02_eval.py $k 4 3 2 1 >> eval_log.txt 2>&1; done
# extension points (extra k beyond the requested 5..10) once their finals exist
for k in 11 12 13 14 15; do
  while [ ! -f finals_$k.json ]; do sleep 10; done
  nice -n 10 "$PY" s02_eval.py $k 5 2 1 >> eval_log.txt 2>&1
done
echo "eval chain done" >> eval_log.txt
