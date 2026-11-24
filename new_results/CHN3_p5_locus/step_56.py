summary = '''# CHN3: the 5-chain P5 — where do irreducible MZVs vanish to all orders?

Work dir: `<project>/research/08_theory/CHN3_p5_locus/`.
Copied from CHN2 (read-only sources, copied not edited): engine5.py, kimilib.py, mzv_exact.py, mzv_num.py, mzv_vals.py,
symengine.py, engine.py, forms10.pkl, forms8_W8.pkl, fft_check.pkl, numcheck_pt1/2.pkl; from CHN1: w12_planes.py/.json.

## Step 1 — exact geometry of V10 = {C4=0} ∩ {Q7'=0} (DONE, EXACT)

C4 (quartic, 55 terms) and Q7' (sextic, 182 terms) extracted from forms10.pkl:
P8 = a2a3a4(SA)C4/10, P10(zeta(3,7)) = a2a3a4(SA)Q7'/224, P10(z2 z(3,5)) = -a2a3a4(SA)(a2+a4)(a1+a3+a5)C4/10.
Both irreducible over Q; both reversal-invariant (subs simultaneous=True!); C4(1..1) = -232.
Degrees: C4 deg 4 in a1,a5 and deg 3 in a2,a3,a4; Q7' deg 6 in a1,a5, deg 5 in a2,a3,a4.
Leading coeffs in a5 are CONSTANTS (4 and 64) => resultant in a5 loses nothing.

Dehomogenise a3=1 (positive points have a3>0; nothing positive lost).
**Res_a5(C4,Q7') = 16 * a1 * (a1+a3) * G(a1,a2,a3,a4)**, G irreducible over Q of degree 22 (2138 terms dehomogenized).
Files: R_a1a2a4.txt (dehomogenized), R_factors.txt.

Projection of V10 to (a1,a2,a3,a4)-space = {a1=0} ∪ {a1+a3=0} ∪ {G=0}:
- {a1=0}: boundary face (carries {a1=0,a5=0} stratum and the P4 curve {a1=0,C3=Q5=0}) — NOT in open positive orthant.
- {a1+a3=0}: no positive points at all.
- {G=0}: the ONLY component meeting the open positive orthant.

On {G=0} (and A != 0), V10 is the graph **a5 = -B/A**, A,B the deg-1 subresultant coeffs in (a1,a2,a4)
(subres1_A.txt 692 terms, subres1_B.txt 748 terms; subres0.txt == R exactly). So V10's positive part is birational
to the irreducible degree-22 surface {G=0}: a SINGLE Q-irreducible component, nonempty in the positive orthant.

EXACT certified positive point: (1,1,1,a4*,a5*), a4* = unique real root of G(1,1,a4) in (1,2)
(isolated to 2^-66 via dup_refine_real_root), a5* = -B(1,1,a4*)/A(1,1,a4*) certified in
(0.07051999939180706,...) > 0 by exact interval arithmetic (Fraction); A(1,1,a4*) ~ -2.03e15 != 0.
Second positive point: a4* in (5,6) (~5.1967), a5* ~ 6.5304. (a4*=1 double root gives a5 in {0,-1}, not positive;
a4~3.4362 gives a5~ -2.0832, not positive.)

Labels: EXACT (sympy resultant/factor/subresultants; engine forms not certified) + NUMERICAL (nroots values).

## Step 2 — engine to W=11,12 (IN PROGRESS)
run_engine12.py / run_eng_diag.py: build_engine(12) is slow (pure-Python Fraction stuffle arithmetic, load ~140).
Checkpointed: eng12.pkl, vw12.pkl (pending).
'''
open('summary.md','w').write(summary)
print("summary written")
import time
time.sleep(60)
print(open('engdiag.log').read())
import os
print([f for f in ['eng12.pkl','vw12.pkl'] if os.path.exists(f)])
