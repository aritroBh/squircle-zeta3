#!/bin/sh
# W4RX: full rerun (about one minute).  Independent scripts: evalG, check_a, check_b, check_cd,
# classify_rows, frame_table (import only w4rx_core / w4rx_rel).  The two scripts marked RERUN
# import the original code and are used only for comparison.
set -e
cd "$(dirname "$0")"
PY=../../../.venv/bin/python
export PYTHONDONTWRITEBYTECODE=1
mkdir -p data logs
$PY evalG.py 900 600 data/G_values_N900.json   > logs/evalG_N900.log
$PY evalG.py 1100 660 data/G_values_N1100.json > logs/evalG_N1100.log
$PY extract_original_rows.py  > logs/extract_original_rows.log      # RERUN of original code (comparison data)
$PY check_a.py      > logs/check_a.log
$PY check_b.py      > logs/check_b.log
$PY check_cd.py     > logs/check_cd.log
$PY classify_rows.py > logs/classify_rows.log
$PY frame_table.py  > logs/frame_table.log
$PY reg_identity_check.py > logs/reg_identity_check.log
$PY rerun_original_verifier.py > logs/rerun_original_verifier.log  # RERUN of original verifier
grep -h "VERDICT\|exact agreement\|exact equality\|identical ordered" logs/*.log
