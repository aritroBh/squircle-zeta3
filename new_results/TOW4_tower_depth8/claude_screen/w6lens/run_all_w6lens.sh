#!/bin/zsh
# Re-run of every weight-six-lens check (all exact arithmetic; each < 20 min on a loaded machine).
# Run from this folder.  Items refer to the lens: (1) Lemma 3 hypothesis, (2) no +-i / const, (3) 5-term identity,
# (4) WT6S tables, (5) k=6 diagonal feed.
PY="<project>/.venv/bin/python"
cd "$(dirname "$0")"
run() { echo "== $*"; nice -n 10 $PY -B "$@"; }
run chk_kernel.py        # (2)  ipos/ineg/const of degree-5 kernels k=5..21
run chk_quartic.py       # (3)  K6=5K7-10K8+10K9-5K10+K11, windows k=7..21, k=5 off, degree exactly 4
run chk_wt6s.py          # (4)  e5m_binomial_coeffs.json 60/60, deg5_residue_polynomials.json 192/192
run chk_tower6.py        # (1)  degree-6 towers depth 0..8 (TOWN1 k<=52 and WTS tower_k.pkl k<=14), (3,3) k=6 off
run chk_feed.py          # Lemma 2 rebuilt: degree-5 kernels from degree-6 towers == stored kernels, k=6..21
run chk_feed2.py         # (5)  feed of Delta_33 into the degree-5 kernels at k=6 (full discrepancy => zero function)
run chk_mech.py          # (5)  F(DL1) = -DL2, DL2(i)=0, DL2'(i) != 0
run chk_e3p3.py          # (5)  TOW4 E3 formula at p=3, delta=0..8  vs  true Laurent coefficients orders 8..0
run chk_k52.py           # quartic law k=6..52 from own tower-built kernels
run chk_w6num.py         # W6(k) = sum res*N (150 digits) vs sum C(k-6,m-1)E5m, k=6..21; k=5 offset
run chk_forced.py        # own forced recursion: closing F(delta,top+1)=0, delta<=8, 396 (d,m) points
run chk_forced2.py       # own forced family vs degree-6 off-diagonal towers (1269 cells) + negative control
run chk_forced3.py       # (3,3): a_delta(6,3) at k=6 and ehat-Newton at k=7..52 (423 cells)
run chk_cache.py         # forced_F78.pkl == own recursion
run chk_negctl.py        # E3 certificate negative controls
run prove_E3b78.py 7     # TOW4 E3 numerators (copied script)
run prove_E3b78.py 8
for a in "7 3" "8 4" "8 3" "7 2" "8 2" "7 1" "8 1" "7 0" "8 0"; do set -- ${=a}; run prove_adj78_fast.py $1 $2; done
