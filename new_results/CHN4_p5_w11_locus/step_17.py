summary = """# CHN4 working summary (auto-updated after each stage)

## Goal
Does V10+ ∩ {f11=0} ∩ {W12 irred coords=0} have positive points? (f11 = zeta(3,3,5) coord of c11.)
Expected answer: NO => P5 has irreducibles at every positive exponent vector (through W12).

## Stage 1 (DONE): setup + validation
- Copied engine files from ../CHN3_p5_locus (point_engine, mzv_exact, kimilib, engine5, common, C4.txt, Q7p.txt, subres1_A/B, R_a1a2a4, forms10.pkl, test_pt_w1112.pkl).
- Validated point_engine vs test_pt_w1112.pkl at (1,2,3,4,5):
  - W11 full coords (reg+m0) match reference EXACTLY (all 9).
  - W12 reg coords match reference for irreducible indices 7,8,11 and the CLASS idx9+(14/9)idx10.
- KEY: m0 (Gamma-product) part contributes EXACTLY 0 to all irreducible content:
  W11 idx7 (z3 z(3,5)), idx8 (z(3,3,5)); W12 idx7 (z2^2 z(3,5)), idx8 (z2 z(3,7)), idx11 (z(2,1,1,8)), and the CLASS idx9+(14/9)idx10 (individually idx9,10 get m0 but it cancels in the combination).
  => use fast regular-only engine (coords_reg) for all irreducible content.
- Reversal symmetry of f11 (W11 idx8) and idx7 confirmed at a generic Fraction point.
- Timing: W11 ~2.2 s/point exact Fraction; W12 ~8.3 s/point.

## Planner anchor values (planner_v10_w12.out), V10 points a3=1:
- (1,1,1,1.293668784620243,0.07051999939180706): f11 = z(3,3,5) = -4.2895814 (NONZERO)
- (1,1,1,5.196691946779181,6.530403653425654): f11 = -10747354.0 (NONZERO)
- Both f11 NEGATIVE at the two known positive V10 points.

## Next: Stage 2 - exact f11 form
- test divisibility by a2a3a4(Sum a); interpolate quotient Q11 (deg 7, 330 monomials, 170 w/ reversal).
"""
open('summary.md','w').write(summary)
print("summary.md written")
