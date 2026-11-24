#!/bin/zsh
PY='<project>/.venv/bin/python'
cd "$(dirname "$0")"
run() {
  local tag="$1"; shift
  local t0=$(date +%s)
  echo "=== START $tag $(date +%H:%M:%S)" >> chainB.log
  timeout 1200 nice -n 10 "$PY" -B "$@" > "out_${tag}.txt" 2>&1
  local rc=$?
  local t1=$(date +%s)
  echo "=== END $tag rc=$rc wall=$((t1-t0))s :: $(tail -n 3 out_${tag}.txt | tr '\n' '|')" >> chainB.log
}
for ds in "6 3" "6 2" "6 1" "6 0"; do
  set -- ${=ds}
  run "ADJ_$1_$2" prove_adj56.py $1 $2
done
