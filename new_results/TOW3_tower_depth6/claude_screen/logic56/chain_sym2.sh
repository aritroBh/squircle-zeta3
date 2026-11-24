#!/bin/zsh
PY='<project>/.venv/bin/python'
cd "$(dirname "$0")"
run() {
  local tag="$1"; shift
  local t0=$(date +%s)
  echo "=== START $tag $(date +%H:%M:%S)" >> chain_sym2.log
  timeout 1150 nice -n 10 "$PY" -B "$@" > "sym_${tag}.txt" 2>&1
  local rc=$?
  local t1=$(date +%s)
  echo "=== END $tag rc=$rc wall=$((t1-t0))s :: $(tail -n 2 sym_${tag}.txt | tr '\n' '|')" >> chain_sym2.log
}
for ds in "6 1" "5 0" "6 0"; do
  set -- ${=ds}
  run "ADJ_$1_$2" prove_adj56.py $1 $2
done
