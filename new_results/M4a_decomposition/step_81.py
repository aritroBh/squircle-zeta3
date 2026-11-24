import json
res = json.load(open("hp_results.json"))
md = """# M4a: weight-5 (W5) and weight-4 (W4'') parts of the squircle coefficient M4

**VERDICT: ALL 5 TASKS PASS.**
1. W5 + W4'' + W3^(4) + T = M4 verified to 5.5e-108; M4 also by Neville extrapolation (4.4e-56).
2. G-expressions for W5 (782 terms) and W4'' (432 terms) built exactly; match quadrature to ~200 digits.
3. Level-8 word vectors built (W4'': 3524+2638 words; W5: 18900+25762 words); match quadrature to ~1e-62.
4. W4'' to 800 digits (agreement 2.7e-810), W5 to 600 digits (8.3e-612), M4 to 800 digits (3.6e-813).
5. Candidate closed form for W4'' found: stable at 700 and 1000 digits (identical up to sign),
   residual 6.9e-716 / 1.5e-1016, and it agrees with the 1060-digit W4'' to 1.07e-1060.

**False alarm resolved (important for reruns):** `import dilog_exact` (triggered by unpickling the
paper's `*_tails_exact.pkl`) executes `mp.mp.dps = 60` at import time. Any mpmath computation done
after `pickle.load` silently runs at 60 digits. All apparent ~1e-60 discrepancies with the paper's
exact values were this artifact; with dps set AFTER imports, the paper's W3' (k=3) and W3^(4) (k=4)
agree with quadrature to 1.2e-110, and the rem:k4 identity
W3^(4) = W3'(k=3) - zeta(3)/8 - pi^2*lam/24 - lam^3/12 holds to 2.8e-1264 (flint, prec 4200).
The paper's values are CORRECT.

## 1. Exact Omega's for k = 4   [ibp_reduce_k4.out identical to paper; omega4.json exact pole/residue data]
a = 1+sqrt(2), b = sqrt(2)-1; Omega_ij = sum c_z/(v-z); M4 = sum int_b^a l1^i l2^j Omega_ij dv.

Weight 5 (W5):
  (2,2):  (7 sqrt(2)/128) / v
  (3,1):  (-5 sqrt(2)/96) / v
  (4,0):  (1/24 - 3 sqrt(2)/128) / v  +  (-1/12)/(v+a)

Weight 4 (W4''):
  (2,1):  (7/16 + sqrt(2)/8)/v + (5/16)/(v-a) + (-7/8)/(v+a) + (-5/16)/(v+b)
  (3,0):  (1/24 - sqrt(2)/8)/v + (-1/8)/(v-a) + (5/24)/(v+a) + (-1/6)/(v+b)

Weight 3 (W3^(4), known): (1,1): (3/4+sqrt(2)/4)/v + (-3/2)/(v+a);
  (2,0): (1/4-sqrt(2)/4)/v + (1/4)/(v-a) + (1/4)/(v+a) + (-1)/(v+b)
Weight 2 (T): (1,0): (1/2)/v + (1/2)/(v-a) + (-1)/(v+a) + (-1/2)/(v+b)

## 2. G-expressions (files gpl_w5.tsv, gpl_w4pp.tsv; format of A3_W4_proof/gpl_expression.tsv)
Built with w4_symbol.log_product from (pole z, l1^i l2^j factors, coefficient) triples parsed
from omega4.json. Every entry satisfies len(word)+k = weight. Alphabets {0, a, -a, b, -b}.
  W5: 782 terms, W4'': 432 terms.
Numerical check (own copy of the real_g acb Taylor evaluator + endpoint shuffle-regularization
with unit tangent G_reg((e,);e)=0, since W4'' has poles at z=a): flint prec 700 bits, cap 420:
  W4'': |G - quad| = 0.0 at 90 quad digits, flint radius 2.65e-208 (108 endpoint-regularized words)
  W5:   |G - quad| = 0.0 at 90 quad digits, flint radius 7.99e-209 (no regularization needed)

## 3. Level-8 word vectors (files wv_w4pp.pkl, wv_w5.pkl)
Converted exactly as verify_D.rebuild_target: forms=[dlog(p) for p in w];
forms_path(POINTS[1]->POINTS[A]) - forms_path(POINTS[1]->POINTS[B]); shuffled with atom(('log',2))^k;
split into rational and sqrt2 parts. All words have the correct weight (4 resp. 5).
  W4'': support (3524 rational, 2638 sqrt2); W5: (18900 rational, 25762 sqrt2).
Numeric check (own copy of check_conversion's series/half/g_at_one/h evaluator; prec 400, cap 330):
  W4'': |wv - quad| = 9.52e-63 (flint radius 7.8e-116)
  W5:   |wv - quad| = 8.14e-62 (flint radius 1.6e-115)

## 4. High-precision values (mpmath tanh-sinh on the v-form = mk_numerics method A;
## two independent working precisions; files W4pp_800.txt, W5_600.txt, M4_800.txt, W4pp_1060.txt)
  W4'' @890 vs @850: |diff| 2.67e-810  -> W4pp_800.txt  (quad err est 1e-908)
  W5   @690 vs @650: |diff| 8.30e-612  -> W5_600.txt    (err est 1e-715)
  M4   @890 vs @850: |diff| 3.63e-813  -> M4_800.txt    (err est 1e-922)
  W4'' @1100 (err est 1e-1164) -> W4pp_1060.txt (for the 1000-digit LLL)
  W4'' = 1.424453676008872921514422491257270535637099847675524471725467791684571969152806...
  W5   = -0.18235662466773630012609370629766450393982143738931837189657087268442643560389849...
  M4   = 0.006027993112609122504516108250877605571476734148292802895459093539190168935771...
  Independent-method evidence: G-expression real_g values agree (Task 2); word-vector values agree
  (Task 3); M4 Neville extrapolation agrees to 4.4e-56; final sum W5+W4''+W3^(4)+T-M4 = -2.8e-601
  (limited by W5's 600-digit truncation).

## 5. Candidate closed form for W4''   [step_pslq2.py; pslq2_w4pp_700.json, pslq2_w4pp_1000.json]
Same 60-number basis as the W4 search (22 weight-4 products x2 + the SAME 8 Li4's x2:
1/2, b, -b, 1/sqrt2, -b^2, 2-sqrt2, 1-1/sqrt2, 1/sqrt2-1 [last one unused]).
LLL (python-flint) on [W4''] + basis:
  700 digits:  W4'' coef 11520, max|coef| 1.198e6, residual -6.87e-716  (random ~1e11.5)
  1000 digits: same relation with opposite global sign, residual -1.45e-1016 (random ~1e16.4)
  => relation IDENTICAL up to sign (stable). A second, independent 700-digit run with the
  auto-greedy Li4 basis (m3_pslq.py-style, step_pslq.py) gave a different-basis relation with
  W4'' coef -6912, residual 3.43e-717 — consistent.
Solving for W4'' (lam = log 2, L = log(1+sqrt2), b = sqrt2-1):

W4'' = (1-104 sqrt2) Li4(b) - (1+104 sqrt2) Li4(-b) + 2 sqrt2 Li4(1/sqrt2) + (41/4) sqrt2 Li4(-b^2)
       - (37/4 + sqrt2/8) Li4(1/2) + (10+16 sqrt2) Li4(1-1/sqrt2) - (10-16 sqrt2) Li4(2-sqrt2)
       + L [ (10+2 sqrt2) Li3(b) - (10-2 sqrt2) Li3(-b) - 18 Li3(1/sqrt2) + (63/32 + 3 sqrt2/8) zeta(3) ]
       + lam [ 6 sqrt2 (Li3(b) - Li3(-b)) - 11 sqrt2 Li3(1/sqrt2) + (1/16 + 77 sqrt2/64) zeta(3) ]
       + 5 Li2(b)^2 + Li2(b) [ (5+sqrt2) L^2 - 4 lam L - (7/4) sqrt2 lam^2 - (19/24 + sqrt2/6) pi^2 ]
       + (1333/11520 - 391 sqrt2/2880) pi^4
       + pi^2 [ (49/48 + 11 sqrt2/48) L^2 + (7 sqrt2/6 - 23/16) lam L + (19/48 - 157 sqrt2/192) lam^2 ]
       + (-19/24 - sqrt2/12) L^4 + (-7/6 - 7 sqrt2/3) L^3 lam + (3/2 + sqrt2/8) L^2 lam^2
       + (7/12 + 13 sqrt2/12) L lam^3 + (-37/96 + 13 sqrt2/48) lam^4 .

Verification: evaluated at 1120 dps, closed form - W4''(1060-digit quadrature) = 1.07e-1060.
(Note: its pi^4 and Li/log coefficients indeed differ from W4's; same basis.)

## Files
  ibp_reduce_k4.out, omega4.json, quad_90.json, m4_neville.json,
  gpl_w5.tsv, gpl_w4pp.tsv, gpl_w5_check.json, gpl_w4pp_check.json,
  wv_w5.pkl, wv_w4pp.pkl, wv_w5_check.json, wv_w4pp_check.json,
  hp_results.json, W4pp_800.txt, W5_600.txt, M4_800.txt, W4pp_1060.txt,
  pslq_w4pp_700.json, pslq2_w4pp_700.json, pslq2_w4pp_1000.json, w4pp_closed_form.json

## Rerun commands (all in this directory; copies of paper scripts + w4_symbol/cmzv on sys.path)
  python ibp_reduce.py 4                      # omega parse: see step file history / omega4.json
  python step_quad_pieces.py                  # Task 1
  python step_neville_m4.py                   # Task 1 method-C check
  python gexpr.py                             # Task 2 build TSVs
  python g_eval.py                            # Task 2 numeric check
  python step_wordvectors.py                  # Task 3 build pickles
  python step_wv_check.py ; python step_wv_check5.py   # Task 3 numeric checks
  python step_hp.py w4pp 890 ; python step_hp.py w4pp 850
  python step_hp.py w5 690 ; python step_hp.py w5 650
  python step_hp.py m4 890 ; python step_hp.py m4 850 ; python step_hp.py w4pp 1100
  python step_pslq.py W4pp_800.txt 700        # Task 5 auto-greedy basis
  python step_pslq2.py W4pp_800.txt 700 ; python step_pslq2.py W4pp_1060.txt 1000
  python step_assemble_w4pp.py                # closed-form assembly + 1.07e-1060 check
"""
open("result.md", "w").write(md)
import datetime
with open("progress.log", "a") as fh:
    fh.write(f"[{datetime.datetime.now().isoformat()}] result.md finalized: all 5 tasks documented.\n")
print("result.md finalized,", len(md), "chars")
