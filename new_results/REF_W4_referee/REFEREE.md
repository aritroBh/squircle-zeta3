# Referee report: Theorem thm:W4 ("W_4 up to pi^4"), paper/ieee/main.tex l.968-980

Verdict: PROVED AFTER REPAIRS. The mathematics holds. The repairs are to the written sketch: it leaves out
one ingredient and is loose in one place. No change to the statement is needed. One optional clarification is suggested.

## 1. Statement check (own code: statement_check.py, statement_check.log)
- W_4 was computed from its definition by mpmath tanh-sinh at 90 digits. R (T3) was transcribed by hand from ieee/main.tex.
- W_4 = 0.79578505983224081132035981962972676254652741387102124725072771573639...
- W_4 - R = -1.0e-88, U - R0 = 5.1e-89, V - R1 = -1.1e-88: they agree to the full 88 digits of the run.
- My transcription of T3 matches conjecture_terms() (29 terms, conjecture_terms_dump.txt) exactly, term by term.
  So the exact scripts, which parse paper/squircle_note/main.tex, check the same R that ieee/main.tex prints.

## 2. Proof check
| Step | Kind | Re-run / check |
|---|---|---|
| W_4 = sum of 148 GPLs G_1(w;a)-G_1(w;b), basepoint 1, letters {0,+-a,+-b} | exact algebra (shuffle) | re-run PASS. I checked the decomposition of l, m by hand ((1-b)(1+a)=2) |
| Full symbol of W_4 = symbol of R (16 words over {2,a}) | exact | re-run PASS |
| Delta_{3,1}(W4-R)=0 (pi^2 kept in left factor), each prefactor part | exact (flint, Q) + classical Li2/Li3 identities | re-run PASS. Telescoped and literal Goncharov cuts agree |
| Delta_{1,3}(W4-R)=0 (right factor mod pi^2) | exact + same identities | re-run PASS |
| Delta_{2,2}=0 from coassociativity + injectivity of Delta_{1,1} on weight 2 | exact (det 2) | re-run PASS. Also automatic: Ext^1(Q(0),Q(2))=0 for totally real F |
| Reduced coaction zero => U-R0, V-R1 in Q*pi^4 | theorem (by hand) | see below |
| Bounds |c0|<1e-322, |c1|<1e-321 | certified Arb | re-run at 1100 bits/EPS 2^-1200: output identical to planner_run_1100.txt |

The motivic step is correct. Objects: motivic iterated integrals on P^1 minus {0,+-a,+-b,infinity} over F=Q(sqrt2),
with an ordinary basepoint at 1 and unit-tangent basepoints at a and b, plus Li_n^m at real points of F inside (-1,1).
All of these lie in MT(F) (Deligne-Goncharov 2005). In fact they lie in MT(O_F[1/2]), because every difference is
(sqrt2)^j a^m. There are no roots of unity, no level 8 and no Galois conjugation. The real embedding
sqrt2 -> +1.414 is used throughout, and all paths are real segments that cross no letter.
Argument: P^m(MT(F)) is O(U^dR) tensor Q[L^m] (Brown 2013, the cited [Brown]). If the (3,1), (2,2) and (1,3) components vanish,
then D^u is primitive in O(U^dR)_4. Primitives of weight 4 correspond to Ext^1(Q(0),Q(4)) = K_7(F) tensor Q, which is 0 because r_2=0 (Borel).
So D^m is U-invariant, which gives D^m in Q (L^m)^4 and D = q*16*pi^4 with q rational. This holds unconditionally, for each prefactor
part separately, because every relation used has rational coefficients. So c0 and c1 are individually rational.
No period conjecture is used.

The value-argument concern does not apply here. I checked all 16 relations used in the reductions numerically to 120 digits (check_relations.py):
inversion, the three-term relation, duplication and Spence-Kummer at (b,1/sqrt2), (2-sqrt2,1/sqrt2), (1/sqrt2,2b) and (-1,b^2).
All of their Li3 arguments are real and < 1 (range -5.83 to 0.854), so no branch is ever chosen.
The antiderivatives in polylog_exact.py (P1, P2, J1, J2) are real on each region. They assert that no path crosses s=0 or s=1,
and I checked their derivatives by hand, including the s>1 branch written via r=1/s.

## 3. Certified bound
w4_certified.py bounds c0=(U-R0)/pi^4 and c1=(V-R1)/pi^4 separately, so it bounds the individual c_i and not only the combination.
- Assumptions: Arb/python-flint is correct (acb.integral with the analytic flag forwarded to every log; acb polylog for |z|<1).
- The endpoint tails are bounded analytically by (1/eps) int_0^eps (log(1/s)+C)^k ds = sum C(k,j)(w+C)^j (k-j)!. I checked this formula.
- Cosmetic: the line "tail rad 0.000e+00" is a float underflow in printing (about 1e-353). The arb value itself is used, so the bound is unaffected.
- The bound does not and cannot give c_i=0. The paper correctly leaves c_i=0 to Theorem thm:W4R.

## 4. Repairs (text only)
R1 (main.tex, proof sketch of thm:W4, 2nd paragraph). The sketch omits that the coproduct factors are reduced with classical
Li2/Li3 functional equations and exact real antiderivatives. These identities must hold for motivic periods, not only for numbers.
They do: they specialise functional equations on M_{0,4} and M_{0,5}, whose constants lie in MT(Z), at real points with all
arguments <1. Proposed sentence to add after "...vanish, separately for the rational and the sqrt2 parts.":
"The reduction of the coproduct factors uses exact antiderivatives and the inversion, three-term, duplication and four
Spence--Kummer relations for Li_3, all at real arguments below 1; these are specialisations of functional equations
and hold for motivic periods."
R2 (same sketch, first paragraph). "Write W_4-R as a sum of iterated integrals on P^1 minus {0,+-a,+-b,infinity}" is loose.
R also contains Li_4 at 1/2, 1/sqrt2 and 2-sqrt2, which are naturally integrals on P^1 minus {0,1,infinity}. Suggested wording:
"Write W_4 as a sum of iterated integrals on P^1 minus {0,+-a,+-b,infinity} (base point 1), and W_4-R as a combination of
motivic periods of mixed Tate motives over Q(sqrt2)".
R3 (optional, statement). Make the per-part meaning explicit, since the individual bounds rest on it: "Write R=R_0+sqrt2 R_1 by
splitting the coefficients, and W_4 = U + sqrt2 V likewise. Then U-R_0=c_0 pi^4 and V-R_1=c_1 pi^4 with c_0,c_1 rational."
(Rationality makes c_0 and c_1 unique, so the statement as printed is not wrong.)
No sentence in the statement claims more than is proved.

## 5. Housekeeping (not affecting the theorem)
- result.md (A3_W4_proof) still says "Conjecture 4 remains open". That refers to c_i=0 and is now superseded by thm:W4R.
- The exact checks themselves were not re-derived with independent code; they were re-run and passed. The independent checks were
  the numerical value of the statement, the exact identity of T3 with the checked R, the 16 relations, the antiderivatives by hand,
  and the motivic argument.
