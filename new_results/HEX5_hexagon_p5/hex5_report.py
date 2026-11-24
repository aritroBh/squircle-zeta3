import re
txt115 = open('p5_run115.txt').read()
txt150 = open('p5_run150.txt').read()
p5_115 = re.search(r'P\^\(5\)\(1\) = ([0-9.]+)', txt115).group(1)
p5_150 = re.search(r'P\^\(5\)\(1\) = ([0-9.]+)', txt150).group(1)
common = 0
for a,b in zip(p5_115.replace('.',''), p5_150.replace('.','')):
    if a==b: common+=1
    else: break
print("common leading digits of P5 (dps=115 vs dps=150):", common)

md = f"""# HEX5: independent numeric check of the regular hexagon's order-five coefficients

**Verdict: CONFIRMED.**

The coefficient of zeta(5) in P^(5)(1) for the regular hexagon is **2835*sqrt(3)/32** and that of
zeta(2)zeta(3) is **1395*sqrt(3)/32**, exactly as claimed (paper Prop. 8 + full locality).
Both were recovered purely numerically: a stable 18-dimensional integer relation for P^(5)(1)/sqrt(3)
on the prescribed weight<=5 basis was found by PSLQ, with the identical integer vector at two
independent precisions (125 and 145 digits) and residuals at the expected level ~ maxcoeff*10^(-digits).

## Objects and method

- Body: K_p = {{ sum_i |<n_i,x>|^p <= 1 }}, n_i = (cos(i pi/3), sin(i pi/3)), i=1,2,3, weights 1.
- P(p) = perimeter of K_p; Taylor expansion about p=1: P(1+d) = sum_k P^(k)(1)/k! d^k.
- Engine: "Route A" of ../CX3k_review_draft8/polygon_engine.py - exact formal Taylor series in d of the
  polar arc-length integrand sqrt(r^2+r'^2), r = S^(-1/(1+d)), S(th) = sum_i |<n_i,x(th)>|^(1+d),
  integrated coefficientwise by tanh-sinh quadrature over the kink wedges. Series order N=7 (exact in d to d^5).
- Validated speed-up used for the high-precision runs (step_03.py): the integrand series is pi/3-periodic
  and symmetric about the wedge midpoint (both verified to 1e-40), so ONE wedge (pi/6, pi/2) with paired
  evaluation, times 6, equals the full sum. The optimized engine (step_04.py) reproduces P^(0..3)(1) of the
  original engine exactly and P^(5)(1) to every digit the original reaches (step_08.py cross-check).

## Task 1: sanity Pd3(1), and P^(5)(1) to >= 80 digits

Paper closed form: Pd3(1) = (sqrt3/8)(16 log^3 2 + 3 pi^2 log2 - 96 log^2 2 + 96 log2 - 9 zeta(3) - 6 pi^2).

| quantity | value | check |
|---|---|---|
| Pd3(1), original engine, dps=60 (step_01) | -5.145455017135615837171681131123975071115553566419741699 | matches paper |
| Pd3(1), paper closed form | -5.145455017135615837171681131123975071115553566419741721 | abs diff 2.24e-53 (rel 4.4e-54) |
| P(1) | 3.464101615137754587054892683011744733886... | = 2 sqrt(3), the hexagon perimeter |

P^(5)(1), two independent runs of the optimized engine (tanh-sinh with level doubling until consecutive
levels agree below eps = 10^-(dps-5)):

| run | last-two-levels diff (k<=5) | P^(5)(1) |
|---|---|---|
| dps=115 (step_04) | 6.66e-88 (level 5->6), then 5.32e-115 (level 6->7, stop) | {p5_115} |
| dps=150 (step_04) | 6.66e-88 (level 5->6), then 7.19e-151 (level 6->7, stop) | {p5_150} |

The two runs share {common} leading digits of P^(5)(1). Convergence between successive quadrature levels is
doubly exponential (diffs 1.3e-18, 2.1e-41, 6.7e-88, then 5.3e-115 / 7.2e-151). Cross-check with the
*unmodified* original engine at dps=50, kmax=5 (step_08): P^(5)(1) = 656.415361968..., agreeing with the
high-precision value 656.41536196762... to all ~10 digits that low-precision run converges to.

## Task 2: PSLQ on the weight <= 5 basis

Basis (17 elements, this order): 1, log2, pi^2, log2^2, z3, pi^2 log2, log2^3, pi^4, z3 log2, pi^2 log2^2,
log2^4, z5, z2 z3, pi^4 log2, z3 log2^2, pi^2 log2^3, log2^5, where z2 z3 = zeta(2)zeta(3) = pi^2 zeta(3)/6.
PSLQ vector [P^(5)(1)/sqrt(3), basis...] (dim 18), mpmath.pslq, tol=10^-(dps-25), maxcoeff=1e9 (step_06).

| PSLQ precision | outcome | residual | max coeff |
|---|---|---|---|
| dps=100 | spurious relation, REJECTED (residual/maxcoeff = 7.6e-78, far above maxcoeff*10^-100 = 3.6e-96) | 2.74e-73 | 36116 |
| dps=125 | **relation R** | 7.34e-122 | 61440 |
| dps=145 | **relation R again, identical vector** | 4.81e-142 | 61440 |

Acceptance criterion residual ~ maxcoeff*10^(-digits): at dps=125, maxcoeff*10^-125 = 6.1e-121 vs residual
7.3e-122 (ok); at dps=145, 6.1e-141 vs 4.8e-142 (ok). The dps=100 "relation" fails the same test - 100
digits are simply insufficient for an 18-dim lattice with coefficients ~6e4; the stable 125/145-digit
relation is the accepted one.

Relation R (step_06; independently re-verified at dps=170 in step_07: residual 9.59e-146, and the value
reconstructed from R matches the stored 150-digit P^(5)(1)/sqrt(3) to 149.7 relative digits):

  128 * P^(5)(1)/sqrt(3) =
      30720 log2 - 4920 pi^2 - 61440 log2^2 + 6240 pi^2 log2 + 30720 log2^3 + 448 pi^4
      + 5760 z3 log2 - 1920 pi^2 log2^2 - 5120 log2^4 + 11340 z5 + 5580 z2 z3
      - 123 pi^4 log2 - 1440 z3 log2^2 + 160 pi^2 log2^3 + 256 log2^5

(coefficients of 1 and of z3 are exactly 0; Li4(1/2) and Li5(1/2) are NOT needed - a stable relation
exists without them, so they were never added).

Equivalently:
  P^(5)(1) = (sqrt(3)/128) * [ 30720 log2 - 4920 pi^2 - 61440 log2^2 + 6240 pi^2 log2 + 30720 log2^3
      + 448 pi^4 + 5760 zeta(3) log2 - 1920 pi^2 log2^2 - 5120 log2^4 + 11340 zeta(5)
      + 5580 zeta(2) zeta(3) - 123 pi^4 log2 - 1440 zeta(3) log2^2 + 160 pi^2 log2^3 + 256 log2^5 ].

Subtraction variant (step_07): with the claimed zeta terms fixed,
R' = P^(5)(1)/sqrt(3) - (2835/32) z5 - (1395/32) z2z3 PSLQ-ed on the remaining 15 basis elements
(dimension 16) returns exactly the restriction of R at dps=125 (residual 7.34e-122).

## Task 3: the zeta(5) and zeta(2)zeta(3) coefficients

| coefficient | this computation | claimed | match |
|---|---|---|---|
| zeta(5) | 11340/128 * sqrt(3) = **2835 sqrt(3)/32** (= 88.59375 sqrt(3)) | 2835 sqrt(3)/32 | EXACT |
| zeta(2)zeta(3) | 5580/128 * sqrt(3) = **1395 sqrt(3)/32** (= 43.59375 sqrt(3)) | 1395 sqrt(3)/32 | EXACT |

## Evidence table

| check | value |
|---|---|
| Pd3(1) vs paper closed form (step_01) | abs diff 2.24e-53; engine last-level diff 3.35e-56 |
| P^(5)(1) last two TS levels, dps=115 run | 5.32e-115 |
| P^(5)(1) last two TS levels, dps=150 run | 7.19e-151 |
| P^(5)(1) dps=115 vs dps=150 | {common} common leading digits |
| PSLQ relation stability | identical 18-dim integer vector at dps=125 and dps=145 |
| PSLQ residual, independent recheck at dps=170 (step_07) | 9.59e-146 (maxcoeff 61440) |
| reconstruction of P^(5)(1)/sqrt(3) from R | matches stored value to 149.7 relative digits |
| original engine dps=50 cross-check of P^(5) (step_08) | agrees to 10 digits (its own accuracy) |

## Files that produced this

- step_01.py  sanity: Pd3(1) via the original engine vs the paper's closed form (dps=60)
- step_02.py  first high-precision driver for the original engine (superseded by step_04)
- step_03.py  numerical verification of the pi/3 periodicity and wedge-midpoint symmetry
- step_04.py  optimized engine (one wedge x6, paired evaluation, per-level log); P^(0..5)(1) at dps=50/115/150
- step_05.py  head-to-head per-level diagnostic, original vs optimized (also exposed the 15-digit
              sqrt(3)/2 initialization pitfall that was fixed in step_04)
- step_06.py  PSLQ of P^(5)(1)/sqrt(3) on the 17-element basis at dps=100/125/145
- step_07.py  independent dps=170 residual + reconstruction check; subtraction-variant PSLQ
- step_08.py  unmodified original engine at dps=50, kmax=5 (cross-check)
- hex5_report.py  writes this result.md
- p5_run115.txt, p5_run150.txt  raw outputs (P^(0..5)(1) and P5/sqrt3 at both precisions)
(step_09.py..step_28.py in the directory are the tool's own auto-saved transcript blocks.)
"""
open('result.md','w').write(md)
print("result.md written,", len(md), "chars")
