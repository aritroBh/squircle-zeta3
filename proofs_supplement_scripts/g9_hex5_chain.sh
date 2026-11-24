#!/bin/sh
# Re-run of the notes' pipeline for Theorem 26 (steps S1 to S6), in a scratch directory G9WORK.
# usage: G9WORK=/some/empty/dir sh g9_hex5_chain.sh   (run from the folder scripts/)
set -e
PY=<project>/.venv/bin/python
for s in g9_hex5_s1_engine g9_hex5_s2_atoms g9_hex5_s3_jloc g9_hex5_s4_jtab g9_hex5_s5_blocks g9_hex5_s6_assembly; do
  t0=$(date +%s)
  "$PY" $s.py > "$G9WORK/$s.log" 2>&1
  t1=$(date +%s)
  echo "$s done in $((t1-t0)) s"
done
