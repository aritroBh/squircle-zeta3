#!/bin/zsh
PY='<project>/.venv/bin/python'
cd "$(dirname "$0")"
run() {
  local tag="$1"; shift
  local t0=$(date +%s)
  echo "=== START $tag $(date +%H:%M:%S)" >> chain_sym3.log
  timeout 1100 nice -n 10 "$PY" -B "$@" > "sym3_${tag}.txt" 2>&1
  local rc=$?
  local t1=$(date +%s)
  echo "=== END $tag rc=$rc wall=$((t1-t0))s :: $(tail -n 2 sym3_${tag}.txt | tr '\n' '|')" >> chain_sym3.log
}
run "E3_3" prove_E3b56.py 3
run "E3_4" prove_E3b56.py 4
for ds in "3 1" "4 2" "3 0" "4 1" "4 0"; do
  set -- ${=ds}
  run "ADJ_$1_$2" prove_adj56.py $1 $2
done
for ds in "3 0" "3 1" "4 0" "4 1" "4 2"; do
  set -- ${=ds}
  run "V_$1_$2" prove_V56.py $1 $2
done
