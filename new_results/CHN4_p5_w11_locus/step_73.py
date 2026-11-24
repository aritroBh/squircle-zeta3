summary = """# CHN4 working summary (auto-updated)

## Goal
Does V10+ ∩ {f11=0} ∩ {W12 irred coords=0} have positive points? f11 = zeta(3,3,5) coord of c11.
**ANSWER (numerical): NO.** On the W11-irreducible-free locus Gamma, W12 irreducibles are always present.

## Stage 1 (DONE): setup + validation
- point_engine validated vs test_pt_w1112.pkl (W11 full coords match exactly; W12 irreducibles match).
- m0 (Gamma-product) contributes EXACTLY 0 to all irreducible content (W11 idx7,8; W12 idx7,8,11,CLASS). Use regular-only engine.
- Timing: W11 ~2.2s/pt exact; W12 ~8s/pt exact, ~35s/pt mpf@120d.

## Stage 2 (DONE): exact f11 form   [EXACT]
- f11 divisible by a2*a3*a4*(Sum a) (vanishes on each hyperplane).
- Interpolated quotient Q11 (deg 7, reversal-invariant, 170 orbits) from 170 rank-170 pts;
  VERIFIED f11 = a2a3a4(Sum a)Q11 EXACTLY on 48 held-out pts (0 mismatches).
- **Q11 = -(1/10) P7**, P7 integer, primitive, IRREDUCIBLE over Q, reversal-symmetric, deg 7, 209 terms.
  P7(1,1,1,1,1) = -5104. Matches planner f11 at both V10 pts. Files: Q11.txt, P7.txt, Q11_poly.pkl.
- Rank-167 curiosity = reversal-duplicate points (reversal-inv polys can't tell a from rev a). Not a bug.
- RUNNING: 693-point rigorous proof (D=f11-L*Q11 = 0 on rank-693 reversal-inv deg-11 det set).

## Stage 3 (DONE): sign of Q11 on V10+   [NUMERICAL]
- Robust V10 solver: exact a4 isolation via G.intervals() + mp bisection + common-root a5 (C4 & Q7').
  CAUTION: a5=-B/A gives SPURIOUS lifts at degenerate a4 (e.g. a4=1 double root); must use common-root a5.
- **Q11 CHANGES SIGN on V10+**: confirmed at 60 digits:
  Q11=-0.0495 at (1/100,1/100); Q11=-140429 at (1,1); Q11=+0.222 at (1/100,3).
- Wide survey [1e-3,1e3]^2: 4057 cells show BOTH signs across sheets. Gamma is extensive (not boundary-confined).
- => Gamma = V10+ ∩ {f11=0} is a REAL curve with positive points. P5 irreducible-free thru W11 EXISTS on Gamma.

## Stage 4 (DONE, core): W12 on Gamma   [NUMERICAL]
- Gamma traced (Q11=0 on V10+): branch with a2~a4, a5~1.0-1.4, a1 in [0.0005,0.2].
- mpf engine @120d at 6 Gamma pts: at EVERY pt f11 & z3z35 vanish (~1e-49..1e-72) AND
  ALL W12 irreducibles NONZERO:
    a1=0.0005: z2^2z35=-1.80 z2z37=-9.67 CLASS=35.5 z2118=15.0
    a1=0.01:   -46.72 -250.29 CLASS=919.26 z2118=389.34
    a1=0.1:    -3596 -19265 CLASS=70757 z2118=29968
    (etc; max|W12 irred| from 35 to 70757, never 0)
- => no Gamma point is irreducible-free through W12.

## Logic of the answer (NO):
- Off V10 (C4 or Q7' !=0): irreducible at W8 (z(3,5)) or W10 (z(3,7)).   [EXACT, CHN3]
- On V10 \\ Gamma (f11 !=0): irreducible at W11 (z(3,3,5)).  (z3z(3,5) vanishes on V10: control)
- On Gamma (V10 and f11=0): W12 irreducibles present.   [NUMERICAL]
- Hence every positive a has an irreducible MZV through W12.

## Remaining / TODO
- Finish 693-pt divisibility proof (rigor upgrade for f11=L*Q11).
- Control: z3z(3,5) coord of c11 is a multiple of C4 (test on {C4=0}).
- More Gamma coverage (larger a1, other components) for robustness.
- Boundary limits of Q11/Gamma (a_i -> 0/inf).
"""
open('summary.md','w').write(summary)
print("summary.md updated")
