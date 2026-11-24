#!/bin/sh
# re-run every step and keep its printed output (used for result.md)
PY="./.venv/bin/python"
cd "research/10_frontier/E1_edge_law"
"$PY" s01_edge_sym.py > s01_log.txt 2>&1
"$PY" s03_beta_table.py > s03_log.txt 2>&1
"$PY" s04_exact_local.py > s04_log.txt 2>&1
"$PY" s08_families.py > s08_log.txt 2>&1
"$PY" s14_reductions.py > s14_log.txt 2>&1
"$PY" s12_basis_example.py > s12_log.txt 2>&1
"$PY" s15_family_check.py > s15_log.txt 2>&1
"$PY" s10_geometric_forms.py > s10_log.txt 2>&1
"$PY" s09_theta_check.py > s09_log.txt 2>&1
"$PY" s07_pslq_H.py 80 120 > s07_log.txt 2>&1
"$PY" s05_fd.py > s05_log.txt 2>&1
"$PY" s02_validate.py 30 > s02_log.txt 2>&1
"$PY" s06_generic_hex.py 30 > s06_log_30.txt 2>&1
"$PY" s11_degeneration.py > s11_log.txt 2>&1
"$PY" s11b_degeneration_small.py > s11b_log.txt 2>&1
echo ALLDONE > batch_done.txt
