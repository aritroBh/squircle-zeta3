import json, datetime
q = json.load(open("quad_90.json"))
nev = json.load(open("m4_neville.json"))
md = f"""# M4a: weight-5 (W5) and weight-4 (W4'') parts of the squircle coefficient M4

**VERDICT (so far): Task 1 PASS.** The decomposition M4 = W5 + W4'' + W3^(4) + T is verified
numerically at 108 digits (residual {q['resid_exact_sum']}); M4 independently confirmed by Neville
extrapolation (method C) at 55 digits. Tasks 2-5 in progress.

**Note on a false alarm resolved:** `import dilog_exact` (triggered by unpickling the paper's
`*_tails_exact.pkl`) executes `mp.mp.dps = 60` at import time. Any mpmath computation done after
`pickle.load` silently runs at 60 digits. All apparent ~1e-60 "discrepancies" with the paper's exact
values were this artifact; with dps set AFTER all imports/loads, the paper's reduced W3' (k=3) and
W3^(4) (k=4) agree with quadrature to 1.2e-110 and with the rem:k4 identity
W3^(4) = W3'(k=3) - zeta(3)/8 - pi^2 lam/24 - lam^3/12 to 2.8e-1264 (flint, prec 4200).

## Exact Omega's for k = 4 (from ibp_reduce.py 4, re-run; output identical to paper's
## outputs/ibp_reduce_k4.out; parsed to pole/residue data omega4.json with exact reconstruction check)

a = 1+sqrt(2), b = sqrt(2)-1.  Omega_{{ij}} = sum_z c_z/(v-z):

Weight 5 (i+j=4):
  (2,2):  7 sqrt(2)/128 * 1/v
  (3,1):  -5 sqrt(2)/96 * 1/v
  (4,0):  (1/24 - 3 sqrt(2)/128) * 1/v  +  (-1/12)/(v+a)

Weight 4 (i+j=3):
  (2,1):  (7/16 + sqrt(2)/8)/v + (5/16)/(v-a) + (-7/8)/(v+a) + (-5/16)/(v+b)
  (3,0):  (1/24 - sqrt(2)/8)/v + (-1/8)/(v-a) + (5/24)/(v+a) + (-1/6)/(v+b)

Weight 3 (i+j=2) [known, m4_tails_exact]:
  (1,1):  (3/4 + sqrt(2)/4)/v + (-3/2)/(v+a)
  (2,0):  (1/4 - sqrt(2)/4)/v + (1/4)/(v-a) + (1/4)/(v+a) + (-1)/(v+b)

Weight 2 (i+j=1) = T integrand:
  (1,0):  (1/2)/v + (1/2)/(v-a) + (-1)/(v+a) + (-1/2)/(v+b)

## Task 1 numerics (mpmath tanh-sinh, 90 dps working; file quad_90.json)

  W5   = {q['W5']}
  W4'' = {q['W4pp']}
  W3^4 = {q['W34_quad']}   (= paper's exact value to 1.16e-110)
  T    = {q['T']}   (= (pi^2+6 lam^2)/24 exactly)
  M4   = {q['M4']}   (full k=4 integrand; quad err est 1e-147)

  W5 + W4'' + W3^4(exact) + T - M4 = {q['resid_exact_sum']}   [< 50 digits: PASS at ~108 digits]
  W5 + W4'' + W3^4(quad)  + T(quad) - M4 = {q['resid_quad_sum']}

  Method C (Neville extrapolation of the perimeter, 12 points p=400..819200):
  M4_neville = {nev['M4_neville']}
  M4_neville - M4_quad = 4.37e-56   [independent check: PASS at 55 digits]

## Rerun commands so far
  python ibp_reduce.py 4 > ibp_reduce_k4.out     # then step that wrote omega4.json
  python step_quad_pieces.py
  python step_neville_m4.py
"""
open("result.md", "w").write(md)
with open("progress.log", "a") as fh:
    fh.write(f"[{datetime.datetime.now().isoformat()}] result.md written (Tasks 0-1 documented).\n")
print("result.md written")
