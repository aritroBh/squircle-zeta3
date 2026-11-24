import os, time
base='<project>/research/08_theory'
W=base+'/CHN5_p5_w12_exact'
os.chdir(W)
summary = """# CHN5: an exact weight-12 cut of the P5 curve Gamma

## Question
Does Gamma+ (open orthant) contain a point with A = CLASS = 0?
A = coef z2^2 z(3,5) of c12; CLASS = coef z(3,9) + (14/9) coef z(5,7); B = coef z2 z(3,7); Z = coef z(2,1,1,8).

## VERDICT (current): NO — A < 0 on ALL of Gamma+, so no interior point is irreducible-free through W12.
- EXACT: A, B, CLASS, Z are each a2a3a4(Sum a) * (reversal-inv degree-8 quotient) [PROVED, 924-pt rank-924].
- EXACT: on V10, B = (75/14) A and Z = -(25/3) A [PROVED, ideal membership in (C4,Q7')].
  => on V10 the W12 irreducible content is spanned by A and CLASS alone.
- NUMERICAL: A = L*QA < 0 and CLASS = L*QC > 0 at every interior Gamma point (traced a1/a3 = 1e-11..1e9+).

## Step 1: exact forms [EXACT, PROVED]
- det12: 924-pt determining set for reversal-inv deg-12, rank 924 mod two primes (nmod_mat). 255-pt deg-8 interp subset.
- m0 contributes EXACTLY 0 to A,B,CLASS,Z at all 924 pts.
- D = coord - a2a3a4(Sum a)*Q8 == 0 at ALL 924 pts for all four coords => PROVED identities.
  Quotient term counts (reversal-orbit coeffs): QA 186, QB 182, QC 230, QZ 115. Reversal-inv, degree 8.

## Step 2: proportionalities [EXACT, PROVED]
- QB - (75/14) QA = U1*C4 + V1*Q7'  (U1 deg4, 50 terms; V1 deg2, 6 terms)  [ideal membership, verified]
- QZ + (25/3) QA = U2*C4            (U2 deg4, 50 terms; V2 = 0)            [ideal membership, verified]
- (LP: no sign-definite representative of QA/QC mod (C4,Q7') on the WHOLE orthant; sign must be argued on V10.)

## Step 3: A on Gamma+ [NUMERICAL]
- Gamma = {C4=Q7'=P7=0}, a3=1. Traced via V10-sheet Q11 sign changes + Newton (mpf), a1 = 1e-11..1e9+.
- ALL interior points: QA < 0, QC > 0 (L > 0 => A < 0, CLASS > 0). No exceptions (119+ pts).
- Stratum artifacts at a1=1 ({a5=0,a2=a4,a1=a3} curve, QA=QC=0 there) excluded; genuine near-stratum
  interior branch has QA<0.
- Small-a1 end: approaches stratum pt (0, 8/3, 1, 8/3, 1); QA/S^8 ~ -8.16201e-6*a1 (linear, negative slope),
  QC/S^8 ~ +1.606e-4*a1. Clean, no crossing.
- Large-a1 end: Gamma extends to a1>=1e9 (a3->0 face); QA/S^8 -> ~ -1.89e-6, -2.33e-5 (negative plateaus).
- min |QA|/S^8 over interior ~ 1e-7..1.9e-6 (never 0 in interior; ->0 only at the a1=0 boundary stratum).
- Cross-checks: 11 independent dense solves agree (QA<0, QC>0).
- Caveat: branch discovery via scanning is not a certified completeness proof; resultant elimination in progress.

## Files
det12_points.pkl, Q8_quotients.pkl, eval12_*.pkl (924 pts), interp12.log, ideal.log,
gamma_trace.pkl/log (numpy, superseded), gamma_trace2.pkl/log (mpf), progress.log, run_*.py
"""
open('summary.md','w').write(summary)
print("summary.md updated")
import time
time.sleep(60)
print(open('resultants.log').read() if os.path.exists('resultants.log') else "(no log yet)")
