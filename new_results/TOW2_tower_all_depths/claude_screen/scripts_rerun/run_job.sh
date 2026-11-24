#!/bin/bash
# usage: run_job.sh NAME script args...
PY='<project>/.venv/bin/python'
name="$1"; shift
cd "$(dirname "$0")"
s=$(date +%s)
{ echo "START $name $(date +%T)"; nice -n 10 "$PY" -B "$@"; echo "EXIT $?"; } > "out_$name.log" 2>&1
e=$(date +%s)
echo "WALL $((e-s)) s" >> "out_$name.log"
