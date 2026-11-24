summary = r"""# CHN4 final summary

## Question
Is there any positive exponent vector a=(a1..a5) (open orthant) at which the 5-chain P5 volume has
NO irreducible MZV through weight 12? I.e. is V10+ ∩ {f11=0} ∩ {W12 irreducible coords=0} nonempty
in the open orthant?  (f11 = zeta(3,3,5) coordinate of c11.)

## ANSWER: NO.  [part EXACT, part NUMERICAL]
Every positive point carries an irreducible MZV by weight 12:
- off V10 (C4 or Q7' != 0): irreducible at W8 (z(3,5)) or W10 (z(3,7)).         [EXACT]
- on V10 \ Gamma (f11 != 0): irreducible at W11 (z(3,3,5)).                     [EXACT]
- on Gamma (C4=Q7'=f11=0, a_i>0): W12 irreducibles present.                     [NUMERICAL]
The ONLY irreducible-free-through-W12 points lie on the BOUNDARY (some a_i=0), never in the open orthant.

## The exact W11 form f11   [EXACT, PROVED]
**f11 = zeta(3,3,5) coord of c11 = -(1/10) a2 a3 a4 (a1+..+a5) P7(a)**,
P7 = irreducible / Q, primitive, reversal-symmetric, homogeneous deg 7, 209 terms (P7.txt; Q11=f11/(a2a3a4 Sum a)=-P7/10, Q11.txt).
- Divisibility by a2a3a4(Sum a) PROVED: D=f11-a2a3a4(Sum a)Q11 = 0 on a rank-693 reversal-invariant
  deg-11 determining set (run_proof_div.py; det11 rank=693 verified mod 2 primes). 0 nonzeros.
- P7 irreducible (sympy factor_list). Matches planner f11 at both V10 points.

## Control   [EXACT]
W11 zeta3*zeta(3,5) coord = a2a3a4(Sum a)*C4*H3 (H3 deg 3, 21 terms, all-positive coeffs).
=> z3z(3,5) vanishes on V10, so on V10 the ONLY W11 irreducible is f11. (run_ctrl_z3z35.py; H7 interp
verified on 20 fresh pts, H7/C4 exact division remainder 0.)

## V10 = irreducible-free locus through W10   [EXACT, from forms10.pkl]
- W8 z(3,5) = a2a3a4(Sum a) C4 / 10
- W10 z2 z(3,5) = -a2a3a4(Sum a)(a2+a4)(a1+a3+a5) C4 / 10   (vanishes at C4=0)
- W10 z(3,7) = a2a3a4(Sum a) Q7' / 224                      (vanishes at Q7'=0)
So through W10, positive irreducible-free points are EXACTLY V10={C4=0,Q7'=0}.

## Gamma = V10+ ∩ {f11=0} is a NONEMPTY real curve   [NUMERICAL]
Q11 (=-P7/10) CHANGES SIGN on V10+ (robust V10 solver: exact a4 isolation via G.intervals + mp
bisection + COMMON-ROOT a5; a5=-B/A gives spurious lifts at degenerate a4, do NOT use it blindly).
- Q11=-0.0495 at (a1,a2)=(1/100,1/100); Q11=+0.222 at (1/100,3); Q11=-140429 at (1,1). [60 digits]
- Sign change on ONE continuous sheet (a1=1/100: Q11 crosses 0 between a2=2.725 and 2.750).
- Wide survey [1e-3,1e3]^2: 4057 (a1,a2) cells show BOTH signs across sheets. Gamma is extensive.

## W12 kills Gamma in the open orthant   [NUMERICAL]
mpf engine @120d at 8 interior Gamma points (a1=5e-5..2.0, two regions). Controls (W8 z(3,5), W10
z(3,7), W11 z3z(3,5), f11) vanish ~1e-49..1e-81 (confirms on Gamma). ALL four W12 irreducibles
(z2^2 z(3,5), z2 z(3,7), CLASS z(3,9)+(14/9)z(5,7), z(2,1,1,8)) are NONZERO at every point:
  normalized max|W12 irred| / max|c12| = 6.1e-8 .. 5.6e-3, always > 0.
The ratio -> 0 (linearly) only as a1->0 toward the boundary stratum; at every a_i>0 point it is > 0.

## Boundary structure (the "hidden" locus)   [EXACT, Fraction engine]
- On {a1=0, a2=a4, a5=a3}: C4=Q7'=Q11=0 IDENTICALLY, and ALL W12 irreducibles vanish EXACTLY
  (z2^2z35=z2z37=CLASS=z2118=0 at t=1,2,3, a3=1,2). Same for reversal {a5=0,a2=a4,a1=a3}.
- Faces {a2=0},{a3=0},{a4=0}: ALL irreducibles (W8..W12) vanish at all tested points (a2a3a4
  divides every irreducible coordinate). {a1=0}/{a5=0} OFF-stratum: irreducibles present.
=> Closed-orthant irreducible-free-through-W12 locus = faces a2,a3,a4=0 + strata
  {a1=0,a2=a4,a5=a3} and {a5=0,a2=a4,a1=a3}. (The a1=0/a5=0 strata are P4-chain degenerations:
  P5(a1=0) = square x P4(a2,a3,a4,a5), and {a2=a4,a5=a3} is P4's alternating locus.)

## Verdict labels
- EXACT (computer algebra): f11=-(1/10)a2a3a4(Sum a)P7 [PROVED]; V10=W10 locus; W11 control; boundary vanishing.
- NUMERICAL (mpmath): Q11 sign change / Gamma nonempty; W12 nonzero on interior Gamma.
- CONJECTURAL: MZV independence (only to read "coefficient !=0" as "irreducible genuinely present").

## Key files
Q11.txt, P7.txt (f11 quotient); Q11_poly.pkl; proof_div.log (693-pt proof); ctrl_z3z35.log;
H7_z3z35.pkl, H7_div_C4_quotient.pkl (control); gamma_final.log/pkl (W12 on Gamma);
gamma_final_pts.pkl; det11_points.pkl (rank-693 det set); interp_points.pkl; signmap.pkl.
Scripts: run_f11_evals.py, run_proof_div.py, run_ctrl_z3z35.py, run_gamma_final.py, survey_lib.py.
"""
open('summary.md','w').write(summary)
print("summary.md written,", len(summary), "chars")
