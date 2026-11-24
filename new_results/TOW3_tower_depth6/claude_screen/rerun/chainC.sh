#!/bin/zsh
PY='<project>/.venv/bin/python'
cd "$(dirname "$0")"
run() {
  local tag="$1"; shift
  local t0=$(date +%s)
  echo "=== START $tag $(date +%H:%M:%S)" >> chainC.log
  timeout 1200 nice -n 10 "$PY" -B "$@" > "out_${tag}.txt" 2>&1
  local rc=$?
  local t1=$(date +%s)
  echo "=== END $tag rc=$rc wall=$((t1-t0))s :: $(tail -n 3 out_${tag}.txt | tr '\n' '|')" >> chainC.log
}
for ds in "1 0" "2 0" "2 1" "3 0" "3 1" "4 0" "4 1" "4 2"; do
  set -- ${=ds}
  run "V_$1_$2" prove_V56.py $1 $2
done
for dl in 1 2 3 4; do run "E3_$dl" prove_E3b56.py $dl; done
for ds in "3 0" "3 1" "4 2" "4 1" "4 0"; do
  set -- ${=ds}
  run "ADJ_$1_$2" prove_adj56.py $1 $2
done
