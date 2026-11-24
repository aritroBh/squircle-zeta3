import os, time
base='<project>/research/08_theory'
W=base+'/CHN5_p5_w12_exact'
os.chdir(W)
summary = """# CHN5: an exact weight-12 cut of the P5 curve Gamma (work in progress)

## Question
Does Gamma+ (open orthant) contain a point with A = CLASS = 0?
A = coef z2^2 z(3,5) of c12; CLASS = coef z(3,9) + (14/9) coef z(5,7).

## Status so far
- Setup: engine copied from CHN4, validated. Basis W12 confirmed:
  idx7=A (z2^2 z(3,5)), idx8=B (z2 z(3,7)), idx9=z(3,9), idx10=z(5,7), idx11=Z (z(2,1,1,8)).
  CLASS = idx9 + (14/9)*idx10. m0 contributes EXACTLY 0 to A,B,Z,CLASS (idx9,10 individually nonzero).
- Reversal invariance of all 12 coords confirmed at 2 generic points.
- det12: 924-pt determining set for reversal-inv deg-12 polys; rank 924 mod 2 primes (nmod_mat).
  255-pt deg-8 interp subset, rank 255. [det12_points.pkl]
- Hyperplane sanity: A,B,CLASS,Z vanish EXACTLY at 4 generic pts each on a2=0,a3=0,a4=0,Sa=0.
- Numerical (from CHN4 gamma_final.pkl): B=(75/14)A and Z=-(25/3)A hold to ~50 digits at all 8
  Gamma points; CLASS/A NONCONSTANT (-19.68..-16.14) -> A and CLASS are 2 independent functions on V10.
- Gamma tracer built: solves C4=Q7'=P7=0 per a1 (dense-guess Newton); finds MULTIPLE branches
  (2 at a1=0.01, 4 at a1=1.0044, 5 at a1=2.0 incl. near-boundary ones).

## Running
- 924-pt exact W12 eval campaign (g=vw+vwr, m0; ~7.2s/pt). ETA ~19:40.
- Then: interpolate quotients QA,QB,QC,QZ (deg 8, 255 orbits); verify D=0 on 924 pts (PROOF);
  factor quotients; ideal membership for proportionalities; Gamma trace with QA,QC signs.
"""
open('summary.md','w').write(summary)
print("summary.md written")
import pickle
print("chunks so far:", [f for f in os.listdir('.') if f.startswith('eval12_') and f.endswith('.pkl')])
