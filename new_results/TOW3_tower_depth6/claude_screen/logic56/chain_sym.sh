#!/bin/zsh
PY='<project>/.venv/bin/python'
cd "$(dirname "$0")"
run() {
  local tag="$1"; shift
  local t0=$(date +%s)
  echo "=== START $tag $(date +%H:%M:%S)" >> chain_sym.log
  timeout 1100 nice -n 10 "$PY" -B "$@" > "sym_${tag}.txt" 2>&1
  local rc=$?
  local t1=$(date +%s)
  echo "=== END $tag rc=$rc wall=$((t1-t0))s :: $(tail -n 2 sym_${tag}.txt | tr '\n' '|')" >> chain_sym.log
}
for ds in "5 0" "5 1" "5 2" "6 0" "6 1" "6 2" "6 3"; do
  set -- ${=ds}
  run "V_$1_$2" prove_V56.py $1 $2
done
run "E3_5" prove_E3b56.py 5
run "E3_6" prove_E3b56.py 6
for ds in "5 2" "6 3" "6 2" "5 1"; do
  set -- ${=ds}
  run "ADJ_$1_$2" prove_adj56.py $1 $2
done
