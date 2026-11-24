step09 = r'''import re
txt115 = open('p5_run115.txt').read()
txt150 = open('p5_run150.txt').read()
p5_115 = re.search(r'P\^\(5\)\(1\) = ([0-9.]+)', txt115).group(1)
p5_150 = re.search(r'P\^\(5\)\(1\) = ([0-9.]+)', txt150).group(1)
q_150 = re.search(r'P5/sqrt3 = ([0-9.]+)', txt150).group(1)
common = 0
for a,b in zip(p5_115.replace('.',''), p5_150.replace('.','')):
    if a==b: common+=1
    else: break
print("common leading digits of P5 between dps=115 and dps=150 runs:", common)

md = f"""# HEX5: independent numeric check of the regular hexagon's order-five coefficients

**Verdict: CONFIRMED.**

The coefficient of zeta(5) in P^(5)(1) for the regular hexagon is **2835*sqrt(3)/32** and that of
zeta(2)zeta(3) is **1395*sqrt(3)/32**, exactly as claimed (paper Prop. 8 + full locality).
Both coefficients were recovered *numerically* (no symbolic input from the paper's edge law):
a stable 18-dimensional integer relation for P^(5)(1)/sqrt(3) on the prescribed weight<=5 basis
was found by PSLQ at two independent precisions (125 and 145 digits, identical integer vector),
with residuals at the expected level ~ maxcoeff * 10^(-digits).

## Objects and method

- Body: K_p = {{ sum_i |<n_i,x>|^p <= 1 }}, n_i = (cos(i pi/3), sin(i pi/3)), i = 1,2,3 (weights 1).
- P(p) = perimeter of K_p; Taylor expansion about p = 1: P(1+d) = sum_k P^(k)(1)/k! d^k.
- Engine: "Route A" of ../CX3k_review_draft8/polygon_engine.py — exact formal Taylor series in d of the
  polar arc-length integrand sqrt(r^2+r'^2), r = S^(-1/(1+d)), S(theta) = sum_i |<n_i,x(theta)>|^(1+d),
  integrated coefficientwise by tanh-sinh quadrature over the kink wedges. Series truncation N = 7 (exact to d^5).
- Speed-up (validated, used only for the high-precision runs): the integrand is pi/3-periodic and symmetric
  about each wedge midpoint (verified numerically to 1e-40, step_03.py), so ONE wedge (pi/6, pi/2) with
  paired evaluation, multiplied by 6, reproduces the full sum. This optimized engine (step_04.py) reproduces
  P^(0..3)(1) of the original engine exactly and P^(5) to all digits the original reaches (step_08.py).

## Task 1: sanity (P'''(1)) and P^(5)(1) to >= 80 digits

Paper's closed form: P'''(1) = (sqrt3/8)(16 log^3 2 + 3 pi^2 log 2 - 96 log^2 2 + 96 log 2 - 9 zeta(3) - 6 pi^2).

| quantity | value | check |
|---|---|---|
| P'''(1) engine (step_01, dps=60) | -5.145455017135615837171681131123975071115553566419741699 | matches paper |
| P'''(1) paper closed form | -5.145455017135615837171681131123975071115553566419741721 | abs diff 2.24e-53, rel 4.4e-54 |
| P(1) | 3.464101615137754587054892683011744733886... | = 2*sqrt(3) (hexagon perimeter) |

P^(5)(1), two independent quadrature runs (optimized engine; tanh-sinh level-doubling until consecutive
levels agree below eps = 10^-(dps-5)):

| run | final level diff (k<=5) | P^(5)(1) |
|---|---|---|
| dps=115 (step_04) | 5.3e-115 (level 6->7; prev level diff 6.7e-88) | {p5_115} |
| dps=150 (step_04) | 7.2e-151 (level 6->7; prev level diff 6.7e-88) | {p5_150} |

The two runs share {common} leading digits of P^(5)(1). Convergence is doubly exponential in the level
(diffs 1.3e-18, 2.1e-41, 6.7e-88, 5.3e-115 / 7.2e-151), i.e. the last two quadrature levels agree far
below the target precision. Also cross-checked against the *unmodified* original engine at dps=50
(step_08): 656.415361968... agreeing to all ~10 digits that run converges to.

## Task 2: PSLQ on the weight <= 5 basis

Basis (17 elements, in this order): {{1, log2, pi^2, log2^2, z3, pi^2 log2, log2^3, pi^4, z3 log2,
pi^2 log2^2, log2^4, z5, z2 z3, pi^4 log2, z3 log2^2, pi^2 log2^3, log2^5}} with z2 z3 = zeta(2)zeta(3)
= pi^2 zeta(3)/6. PSLQ vector: [P^(5)(1)/sqrt(3), basis...] (dimension 18), mpmath.pslq,
tol = 10^-(dps-25), maxcoeff = 1e9.

| PSLQ precision | result | residual | max coeff | residual/maxcoeff |
|---|---|---|---|---|
| dps=100 | spurious small relation (rejected) | 2.7e-73 | 36116 | 7.6e-78 >> 1e-100*maxcoeff*1e-4; fails acceptance |
| dps=125 | **relation R (below)** | 7.34e-122 | 61440 | 1.2e-127 ~ maxcoeff*10^-125*2 |
| dps=145 | **same relation R** | 4.81e-142 | 61440 | 7.8e-148 |

Relation R (identical integer vector at dps=125 and dps=145; residual rechecked independently at dps=170:
9.59e-146; reconstruction from R reproduces the stored 150-digit value of P^(5)(1)/sqrt(3) to 149.7
relative digits, step_07.py):

  128 * P^(5)(1)/sqrt(3) =
       30720 log2 - 4920 pi^2 - 61440 log2^2 + 6240 pi^2 log2 + 30720 log2^3 + 448 pi^4
       + 5760 z3 log2 - 1920 pi^2 log2^2 - 5120 log2^4 + 11340 z5 + 5580 z2 z3
       - 123 pi^4 log2 - 1440 z3 log2^2 + 160 pi^2 log2^3 + 256 log2^5

(coefficients of 1 and z3 are exactly 0; Li4(1/2), Li5(1/2) are NOT needed - a relation exists without them).

Equivalently:

  P^(5)(1) = (sqrt(3)/128) * [ 30720 log2 - 4920 pi^2 - 61440 log2^2 + 6240 pi^2 log2 + 30720 log2^3
      + 448 pi^4 + 5760 zeta(3) log2 - 1920 pi^2 log2^2 - 5120 log2^4 + 11340 zeta(5)
      + 5580 zeta(2) zeta(3) - 123 pi^4 log2 - 1440 zeta(3) log2^2 + 160 pi^2 log2^3 + 256 log2^5 ]

Subtraction variant (step_07): fixing the claimed zeta terms and PSLQ-ing the remainder
R' = P^(5)(1)/sqrt(3) - (2835/32) z5 - (1395/32) z2z3 on the remaining 15 basis elements gives exactly
the restricted relation R at dps=125 (residual 7.34e-122) - the same relation, confirming consistency.

## Task 3: the zeta(5) and zeta(2)zeta(3) coefficients

| coefficient | this computation | claimed | match |
|---|---|---|---|
| zeta(5) in P^(5)(1) | sqrt(3) * 11340/128 = **2835 sqrt(3)/32** (= 88.59375 sqrt3) | 2835 sqrt(3)/32 | EXACT |
| zeta(2)zeta(3) in P^(5)(1) | sqrt(3) * 5580/128 = **1395 sqrt(3)/32** (= 43.59375 sqrt3) | 1395 sqrt(3)/32 | EXACT |

## Evidence / convergence summary

| check | digits / residual |
|---|---|
| P'''(1) vs paper closed form | abs diff 2.24e-53 (dps=60 run; engine internal level diff 3.35e-56) |
| P^(5)(1): last two TS levels, dps=115 | 5.32e-115 |
| P^(5)(1): last two TS levels, dps=150 | 7.19e-151 |
| P^(5)(1): dps=115 vs dps=150 runs | {common} common leading digits |
| PSLQ relation stability | identical 18-dim integer vector at dps=125 and dps=145 |
| PSLQ residual (independent recheck at dps=170) | 9.59e-146 (maxcoeff 61440) |
| reconstruction of P^(5)(1)/sqrt3 from the relation | 149.7 relative digits |

## Files

- step_01.py : sanity, P'''(1) via original engine vs paper closed form (also P^(0..2)(1))
- step_02.py : driver used for first high-precision attempts with the original engine (superseded)
- step_03.py : verifies pi/3 periodicity and wedge-midpoint symmetry of the integrand series
- step_04.py : optimized engine (one wedge x6, paired evaluation, per-level logging); P^(0..5)(1) at dps=50/115/150
- step_05.py : head-to-head original vs optimized per-level behaviour (diagnostic)
- step_06.py : PSLQ of P^(5)(1)/sqrt3 on the 17-element weight<=5 basis at dps=100/125/145
- step_07.py : independent dps=170 residual + reconstruction check; subtraction-variant PSLQ
- step_08.py : unmodified original engine at dps=50, kmax=5 (cross-check of P^(5))
- step_09.py : writes this report
- p5_run115.txt, p5_run150.txt : raw high-precision outputs (P^(0..5)(1), P5/sqrt3)
"""
open('result.md','w').write(md)
print("result.md written,", len(md), "chars")
'''
open('step_09.py','w').write(step09)
import subprocess
r = subprocess.run(['<project>/.venv/bin/python','step_09.py'],capture_output=True,text=True)
print(r.stdout); print(r.stderr[-300:])
