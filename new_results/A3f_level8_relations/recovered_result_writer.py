result = r'''# A3f: level-8 relations toward closing Conjecture 4 (W4) — status report

Work dir: `research/08_theory/A3f_level8_relations`. Venv: `<home>)/.venv/bin/python`.
Engine/data: copied from `A3c_claude_elementary/Q` into `Q/` (import only; originals untouched).

## Status table

| step | goal | status | result |
|---|---|---|---|
| 0 audit | reproduce 3/10/37, check all rows | PROVED | exact-Q quotients 3/10/37 (fmpq rref + primes 1000000009, 1000000021, 999999937); weight-4 quotient 142 (rank 5042) confirmed mod second prime 1000000021; D residue 141/142 nonzero confirmed; ALL relation rows (6/88/1321/17625) numerically zero at ~60 digits (worst 2.9e-57) |
| 1 weight 2 | 10 -> 9 | PROVED | missing relation = Rogers five-term at (x,y)=(sqrt2/2, sqrt2-1) (+4 equivalent instances); quotient exactly 9 = Deligne-Goncharov bound; LLL hint (325 digits, residual 1.5e-329) = 3 x five-term row |
| 2 weight 3 | 37 -> 27 | PARTIAL (34) | lifts of the new weight-2 relation add exactly 3 -> quotient exactly 34 over Q; 9 three-term (T3) + 9 inversion (I3) one-variable Li3 instances add NOTHING (already in span); 7 relations still missing; Kummer K3 two-variable equation: symbol solution space computed (4-dimensional), genuine part identified, not yet verified/specialized |
| 3 weight 4 | 142 -> enough | NOT STARTED (audit only) | still 142 (upper bound, mod 2 primes); D residue 141/142 nonzero in current quotient |
| 4 D reduction | D in relation span | OPEN | not yet tested at weight 4 |

Conjecture 4 is NOT yet proved. Honest bottom line: weight 2 closed at the DG bound; weight 3 at 34 (target 27); weight 4 untouched at 142 (target 81); c0,c1 still open (certified |c0|<1e-322, |c1|<1e-321 from W4 motivic/ball-arithmetic work, reused).

## The new weight-2 relation (exact certificate, weight 2 CLOSED)

Rogers five-term relation for L(z) = Li2(z) + (1/2) log z log(1-z), valid for 0<x,y<1
(Rogers 1907; classical, provable by elementary differentiation; cf. Zagier, "The dilogarithm function"):
    L(x) + L(y) = L(xy) + L(u) + L(v),  u = x(1-y)/(1-xy), v = y(1-x)/(1-xy).
Specialization (all five arguments in (0,1) and in the level-8 cross-ratio set):
    x = sqrt2/2, y = sqrt2 - 1,  xy = 1 - sqrt2/2,  u = 2 - sqrt2,  v = 3 - 2*sqrt2.
Equivalent instances (each alone raises rank 54->55): (sqrt2/2, 2sqrt2-2), (sqrt2-1, sqrt2/2),
(sqrt2-1, sqrt2-1), (2sqrt2-2, sqrt2/2).
As a level-8 MPV row (exact, integer coefficients; via cmzv.atom, atoms previously checked to 1e-51):
  +H(1,2) -H(1,3) -H(1,5) +H(1,6) -H(2,0) +H(2,1) -H(2,2) +H(2,3) -H(2,4) +2H(2,6) -2H(2,7)
  +H(3,1) +H(3,2) -2H(3,3) +H(3,7) -H(4,2) +2H(4,5) +H(4,6) -2H(4,7) +H(5,0) -2H(5,1)
  +H(5,4) +3H(5,5) -4H(5,6) +H(5,7) +H(6,2) +2H(6,4) -3H(6,5) -2H(6,6) +3H(6,7) -2H(7,2)
  +3H(7,3) -2H(7,4) -H(7,5) +3H(7,6) -2H(7,7) = 0.
Certificate check: adding this row to the 88 standard rows raises the exact-Q rank from 54 to 55,
so the quotient drops from 10 to 9 = DG bound. The LLL hint on the 10 free-basis values,
  [0,-6,-1,7,3,-9,-3,21,-5,5]  (words (6,6),(6,7),(7,0)..(7,7); residual 1.47557e-329 at 330 digits),
equals 3 x this row in the quotient (exact normal-form comparison).
Files: weight2_new_relation.{json,pickle}, weight2_fiveterm.pickle (20 five-term + 5 reflection rows),
weight2_values.json (all 64 word values, 330 digits).

## Weight-3 details (34 = 27 + 7)

- lower basis = audited weight-1 basis + weight-2 basis + the new relation; make_relations(3, lower):
  1328 rows, exact-Q rank 542 (mod 1000000009/1000000021 agree) -> quotient 34.
  The new lifts contribute exactly 3 (= weight-1 quotient dimension).
- 9 T3 (Landen) and 9 I3 (inversion) one-variable Li3 instances on cross-ratios: all numerically zero
  (worst 1.7e-57) but rank contribution 0 — they are consequences of standard relations + lifts.
- Kummer two-variable trilogarithm (K3), ansatz = the 20 Kummer monomial arguments m_i(x,y) of K4:
  the condition "sum_i c_i (1-m_i) x m_i x m_i is a totally symmetric tensor" has a 4-dimensional
  solution space (free parameters c16..c19 = coefficients of Li3(x), Li3(y), Li3(-x/xi), Li3(-y/eta)).
  The x<->y-symmetric subfamily has VANISHING symbol and evaluates to ~0 (1e-60) — these are
  combinations of the (already standard) one-variable identities. The genuine two-variable K3
  equations are the x<->y-ANTISYMMETRIC solutions (c16=-c17, c18=-c19), whose symmetric tensors
  integrate to antisymmetric log-cubics G(x,y); derivation computed but verification/specialization
  not finished (step_11_k3_symbol.py, step_12_k3_derive.py). Then: specialize K3 at S-unit (x,y),
  convert to MPV rows, expect quotient 34 -> 27.

## Evidence table (all numbers printed by my runs)

| check | digits/precision | residual / result |
|---|---|---|
| w1 rows (6) numeric | acb prec=200 (~60d) | worst 1.0868029e-59 |
| w2 rows (88) numeric | acb prec=200 | worst 7.9810019e-59 |
| w3 rows (1321) numeric | acb prec=200 | worst 4.7792757e-58 |
| w4 rows (17625) numeric | acb prec=200 | worst 2.8852122e-57 |
| exact-Q ranks w1/2/3 | exact (fmpq rref) | 4/54/539 -> quotients 3/10/37; 3 primes agree |
| w4 rank | mod 1000000021 (Codex used 1000000009) | rank 5042, quotient 142; D residue 141/142 |
| w2 values (64 words) | prec=1100 bits, cap=2048 | ball radius <= 1.6e-329 |
| known w2 rows at 330 digits | mpmath 335 dps | worst 5.0e-330 |
| LLL missing relation | 330 digits | residual 1.47557e-329; next vectors ~63-digit coeffs (kernel 1-dim) |
| 25 five-term/reflection rows | acb prec=200 | worst 1.1611553e-58; exact rank 54->55, quotient 9 |
| 17 T3/I3 rows | acb prec=200 | worst 1.7205439e-57; rank contribution 0 |
| probes | 330 digits | Li2(sqrt2-1) NOT in Q-span{pi^2,log^2 2,log^2 a,log2 log a} (>1e-256); Cl2(pi/4) NOT in span{G, pi log2, pi log a} (>1e-240) |

## Pipeline bug found and fixed (matters for any high-precision reuse of Q/check_conversion.py)
`check_conversion.py` bakes `N = ctx.cap` and `roots` at import time (cap 320, prec 200), and the
kind-1 series converge only like |1/(2(1-zeta8))| ~ 0.654^n, so cap=320 yields only ~1e-28 accuracy
there. Fix before raising precision: set ctx.prec, ctx.cap, CC.N, and recompute CC.roots AFTER import
(step_05_weight2_values.py). With cap=2048: residuals ~1e-330 as expected.

## Rerun commands
    PY='<home>)/.venv/bin/python'
    cd research/08_theory/A3f_level8_relations
    "$PY" -B step_02_audit.py            # exact-Q quotients 3/10/37 + mod-3-prime cross-check
    "$PY" -B step_03_numeric_w123.py     # 60-digit check of all w1-3 rows
    "$PY" -B step_04_audit_w4.py         # w4: 60-digit check of 17625 rows + 2nd-prime rank 5042/142
    "$PY" -B step_05_weight2_values.py   # 64 word values at ~325 digits
    "$PY" -B step_06_weight2_lll.py      # LLL -> missing relation (hint)
    "$PY" -B step_07_weight2_fiveterm.py # Rogers five-term instances -> quotient 9
    "$PY" -B step_08_weight2_certificate.py  # exact certificate: rank 55, LLL hint = 3x row
    "$PY" -B step_09_weight3_lifts.py    # weight-3 quotient 34 with new lifts (exact)
    "$PY" -B step_10_weight3_onevar.py   # T3/I3 Li3 instances add 0
    "$PY" -B step_11_k3_symbol.py        # K3 symbol solution space (4-dim)
    "$PY" -B step_12_k3_derive.py        # K3 candidate combinations + numerics

## What remains (ordered)
1. Verify the antisymmetric K3 solutions numerically + degeneration, specialize at S-unit (x,y)
   points, convert to MPV rows -> expect quotient 27 at weight 3 (or run LLL at weight 3 as a guide).
2. Weight 4: rebuild relations with lower = {w1 basis, w2 basis + new relation, w3 basis + new
   relations}; add Kummer K4 instances (B_kummer_li4 library, already project-verified) specialized
   so all 20 arguments are level-8 cross-ratios, plus I4/D4; mod-p rank first (2 primes), exact after.
3. After each batch, test D (both parts) in the relation span; if residue 0 exactly, Conjecture 4 is
   proved and c0 = c1 = 0 follows since value(D) = (c0 + sqrt2 c1) pi^4 = 0 with pi^4 != 0.
'''
open('result.md','w').write(result)
open('progress.log','a').write(
 'step_09: weight-3 quotient 34 (exact Q, 2 primes) with new weight-2 relation lifted (+3).\n'
 'step_10: T3/I3 one-variable Li3 instances add 0 rank (already spanned); worst numeric 1.7e-57.\n'
 'step_11/12: Kummer K3 symbol solution space = 4-dim; symmetric part symbol-zero (known); antisymmetric part = genuine K3, derived not verified.\n'
 'result.md written: weight2 CLOSED (quotient 9), weight3 at 34/27, weight4 at 142/81, D reduction OPEN.\n')
print('result.md and progress.log written;', len(result), 'chars')
