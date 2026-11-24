import os
os.chdir("<project>/research/08_theory/OCT5D3_octagon_p5_closed")
md = open("result.md").read()
md += """
## TASK 3 RESULT (closure ACCEPTED) + lower weights CLOSED

### Lower weights (all verified by direct evaluation at 1300 digits, residuals ~1e-1298)
- Wt1 = (4800-3360 sqrt2) lam + (2880-1920 sqrt2) L            [also PROVED exactly]
- Wt2 = (360-280 sqrt2) pi2 + (-20880+14760 sqrt2) lam^2 + (-19200+13440 sqrt2) lam L
        + (-7680+5280 sqrt2) L^2 + (-3840+2880 sqrt2) Li2(b)    [LLL identical 800/1500, h=2.09e4]
- Wt3 = (-4320+2880 sqrt2) Li3(b) + (15120-10800 sqrt2) Li3(-b) + (11520-8640 sqrt2) Li3(1/sqrt2)
        + (-17985/2+6915 sqrt2) z3 + (6480-4560 sqrt2) L^3 + (17280-12240 sqrt2) lam L^2
        + (34200-24120 sqrt2) lam^2 L + (34410-24330 sqrt2) lam^3 + (5760-4320 sqrt2) Li2(b) lam
        + (5760-4320 sqrt2) Li2(b) L + (1370-980 sqrt2) pi2 lam + (-900+660 sqrt2) pi2 L
        [LLL identical 900/1500, h=6.88e4; 12-atom basis, Lchi3 eliminated via eq:chi3]
- Wt4 = 29 terms (paper R-basis: 7 Li4 at {b,-b,-b^2,1/sqrt2,1/2,1-1/sqrt2,2-sqrt2}, Li3*log, Li2(b)*.,
        pi2*., z3*log, pure-log^4) — coefficients in closed_forms_lower.pkl
        [LLL identical 900/1500, residuals 2.45e-895 / 1.52e-1495, h=1.37e6]

### Weight-5 remainder R = Wt5 - classical(x0) - C: CLOSED (NUMERICAL)
- Basis: 51 atoms = 21 planner invisible (zeta2 x 12 weight-3, zeta3 x 5, zeta4 x 2, zeta5, Lchi5)
  + K32 (the Li32 symbol-kernel constant, PROVED dim-1 kernel) + 7 mixed Li32xclassical kernel
  constants (PROVED kernel of 61 atoms: 29 = 22 classical + K32 + 6-7 mixed) + 22 classical
  symbol-kernel combos.
- LLL at 1100 digits: relation rel1100 (height 2.13e9, residual 1.26e-1090 at working precision).
- LLL at 1500 digits: relation rel1500 (height 2.2e9, residual 1.47e-1490 at working precision).
- rel1500 VERIFIED by direct evaluation at 900 (resid 1.17e-894), 1100, 1500 (1.4e-1494).
- rel1500 - rel1100 lies in the VERIFIED 26-dim atom value-kernel (rank check exact)
  => the two relations coincide modulo exact atom-relations. ACCEPTED.
- Coefficients are representative/basis-dependent (as expected); in the RREF representative:
  z5:        82937921/2944 - 56662239/2944 sqrt2
  Lchi5:     -400/3 + 39752/3 sqrt2
  z2*z3:     -16159155/2944 - 8342511/2944 sqrt2
  K32:       -3272927/1104 + 7011289/2944 sqrt2
  (full relation in rel1500_frac.pkl)
"""
open("result.md","w").write(md)
print("result.md updated, length", len(md))
