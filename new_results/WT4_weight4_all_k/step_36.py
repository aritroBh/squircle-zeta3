import json, os
work='<project>/research/08_theory/WT4_weight4_all_k'
res = json.load(open(os.path.join(work,'hp_results.json')))
chk = json.load(open(os.path.join(work,'check_tsv_vs_quad.json')))
d1v = res['D1@900']['value']; d2v = res['D2@900']['value']
md = f"""# WT4: the weight-4 part W_4(k) of the squircle coefficient M_k for all k (Kimi K3)

**Status: in progress. Verdicts so far: kernels/TSVs EXACT (controls passed); D1, D2 values at 860 digits (two-precision agreement ~1e-810); closed forms and certificates pending.**

With W_4(k) = W4'' + (k-4) D1 + binom(k-4,2) D2 (kernel polynomial law, KSPAN), this job computes
D1 and D2 exactly and proves closed forms.

## Task 1 (DONE): exact kernels, deltas, G-expressions

`wt4_A_kernels.py` (fastred reduce_k_fast on binteg_3..14): degree-3 keys (2,1),(3,0) for k=3..14.
- CONTROL: k=4 kernels == M4a omega4.json (2,1),(3,0) **EXACT MATCH**.
- Deltas (coefficients p+q*sqrt2; D1/D2 have q=0 everywhere, poles in {{0,a,-a,-b}}):
  - D1: K_21 = 1/v + (3/4)/(v-a) - 2/(v+a) - (3/4)/(v+b);  K_30 = (-1/12)/v - (1/2)/(v-a) + (5/6)/(v+a) - (1/6)/(v+b)
  - D2: K_21 = (1/2)/v + (1/2)/(v-a) - 1/(v+a) - (1/2)/(v+b)  == Rw2 (the weight-2 kernel!) ;  K_30 = (-1/3)/(v-a) + (1/3)/(v+a)
- QUADRATIC LAW K^(k) = K4 + (k-4) D1 + binom(k-4,2) D2 for k=7..14, both keys: **EXACT True**.

`wt4_B_gexpr.py` (w4_symbol.log_product): d1.tsv ({420} terms... see file), d2.tsv; CONTROL rebuilt W4pp G-expression
== M4a gpl_w4pp.tsv EXACTLY (432 terms). Alphabet {{0, +-a, +-b}} throughout.

## Task 1 check (DONE): TSV vs quadrature, dps 130 (`wt4_C_check.py`, flint ctx.prec=700)

| target | |G - quad| | flint rad |
|---|---|---|
| D1 | {chk['D1']['diff']} | {chk['D1']['flint_rad']} |
| D2 | {chk['D2']['diff']} | {chk['D2']['flint_rad']} |
| W4pp control | {chk['W4pp(control)']['diff']} | {chk['W4pp(control)']['flint_rad']} |

Control quad(W4pp) - M4a quad_90 W4pp = {chk['control_quad_vs_m4a']} (90-digit reference limit).

## Task 2 (DONE): D1, D2 to 860 digits (`wt4_D_hp.py`, tanh-sinh, two working dps 900/850)

| constant | value (first 100 digits) | |@900-@850| |
|---|---|---|
| D1 | {d1v[:100]} | 7.742e-812 |
| D2 | {d2v[:100]} | 1.498e-810 |

err_est 1.0e-908 both. Files d1_860.txt, d1_810.txt, d2_860.txt, d2_810.txt.

## Tasks 3-5: pending (LLL closed forms; M4b2-style certificates; k=3 statement)
"""
open(os.path.join(work,'result.md'),'w').write(md)
print('result.md written')
open(os.path.join(work,'progress.log'),'a').write('result.md initial version written\n')
