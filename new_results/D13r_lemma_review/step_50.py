work='research/08_theory/D13r_lemma_review'
result = r"""# D13r — Independent adversarial review of the lemmas behind the certified decay theorem

Reviewer mandate: try to BREAK each lemma. All numbers printed by code in this directory
(step_NN.py files, mpmath @60 dps / sympy exact / python-flint only to reproduce code lines).

## Headline findings

1. **REAL BUG (final cell)**: `step_10.py`, `step_10_J96.py`, `step_10_J192.py` compute
   `ecw = (rt2/pi).sqrt()*e^{-c_w/2}` = **0.4069446** = sqrt(sqrt(2)/pi)*e^{-c_w/2},
   but the tail identity e^{-c2 d^2} = y*sqrt(2/pi)*e^{-c_w/2} needs **0.4839414** = sqrt(2/pi)*e^{-c_w/2}.
   Off by 2^{-1/4} = 0.840896 (the comment `# sqrt(2/pi)` shows the intent; `rt2` is sqrt(2), not 2).
   TM/TR/TLm final-cell pieces are 2^{1/4} too small. **TM is VIOLATED**: certified TM contribution
   (2/pi)(pi/2)^{3/2}*TMreg = **0.1289691** < truth (2/pi)*TM_true*N^{3/2}*e^{theta} = **0.1441235**
   at N = 1e10 (erfc-exact, 60 dps); also 0.1289691 < 0.1324781 at N = 1e12. TR/TLm buggy values still
   dominate their truths (0.29596 > 0.11345, 1.61775 > 0.09084). Corrected final-cell C = **49.94909**
   (was 49.56260). Corrected TM bound 0.15342 > truth 0.14412 (margin 6.4%, tightest at junction).
2. **The main theorem survives with the SAME constant**: overall C = max(194.3809 [100,400],
   617117.5536 [400,1e10], 49.95 corrected final) = **617117.6**, dominated by the sweep, which
   evaluates e^{-c2 d^2} directly (no ecw). Even the buggy 49.56 >> final-cell truth ~2.
3. Two harmless derivation errors (numbers still valid): (a) step_37/38's A1m phase box
   [0, z^2/(1+z^2)] does NOT contain the required A1 = (1-arctan(zz)/zz)/y (needs A1' <= y z^2/(1+z^2),
   false for y<1/3) — but the full (y,sigma,A1) z-box contains z_true at every one of 41x4+21x7
   sampled N in [100, 1.5e10] (sigma-sweep covers the phase with ~8x margin); (b) the identity
   e^{RePhi0+theta} = e^{-(N/2)S(y0)} is not exact (difference (pi sqrt2/(4t))(rho-1)^2 > 0) — but
   the h3-box PROVABLY contains the truth for all t>0.
4. Smaller items: D13h's breg is not the exact regular form of (phi+pi/4)/((pi/rt2)y) (exact form
   derived here, residual 5.7e-44), but min(breg_h - Btrue) = 0.08333 > 0 and |Btrue| <= 1/6 < breg_abs
   ~ 0.25, so d6 valid; final-cell TLfN misses a (pi/2)^{3/2} factor (value ~1e-19229, irrelevant);
   TL0's constant derivation is questionable (my derivation gives 2x the coded value) but the coded
   value bounds the truth with ratio 0.389/0.382/0.381 at N=100/400/1000 (margin 2.6x).

## Verdict table

| lemma | verdict | evidence (worst margin, digits) |
|---|---|---|
| L1 |gamma'|,|gamma''| bounds | HOLDS (proof + numeric) | series+integral-test proof; 1659 pts @60dps: worst |g'|/B=0.80094308518 (lim 0.8009431025), |g''|/B=0.62499997038 (lim 0.625); psi'/psi'' ratios ->1^- |
| L2 lambdaR regular form | HOLDS (proof + numeric) | hand algebra; 900 pts @60dps, worst residual 9.5e-53 |
| L3 Gh decreasing; |g| decreasing; beta bound | HOLDS (proved) | sympy: Gh'=-(rho/rt2)(2e^3+4rt2e^2t+4et^2+rt2t^3)/g^2 exactly, manifestly <0; Re[b/((n+w)(n+w+3/2))]=(n^2+4nt+1.5n+3t+2t^2)/(rt2 D)>0 (resid 1.6e-58); beta worst ratio 0.7071 |
| L4 Phi3b/Phi4b + monotonicity | HOLDS (proved) | Phi''/'''/'''' sympy residuals 0; 2500 pts worst ratios 0.9999997/0.9999996 (tight at inf, never exceed); monotonicity via n'/n <= m/x < d'/d (sympy numerators all-positive) |
| L5 ReX2=th(M1-L2); logcf decreasing | HOLDS (proof + numeric) | 500 pts worst residual 8.0e-46 (rel ~1e-54); logcf' numerator h=u-(1+u)log(1+u), h'=-log(1+u)<0; M1-L2 max -0.25 / -0.209 |
| L6 Binet | HOLDS (proof + numeric) | |mu(z)|<=1/(12|z|cos(arg z)) proved ((r^2-s^2t^2)^2 sin^2 th >=0); ray worst ratio 0.7071=1/rt2<1; real-z E two-sided: no violation dense grid; final-cell ReE bound: max(ReE-bound)=-7.4e-26 dense incl t->0 |
| L7 final-cell factorisations + tail id | HOLDS (with caveats 3a/3b) | all residuals <1e-45 @60dps; tail id 3.8e-60; z=MV/L~ residual 2.9e-56; leading-term identity exact |
| L8 assembly | HOLDS-NUMERICALLY, MINUS the ecw bug | coverage [0,1]U[1,x0/2]U[x0/2,x0-d]U[x0+-d]U[x0+d,oo) complete (d<=x0/3); W1 oddness + Taylor decomposition hand-verified; TM IBP vs erfc: crude dominates (0.368/0.316, 0.297/0.282, 0.258/0.250); signed window/sc = 2.19/2.00/2.08 (bounded) vs certified W-pieces; TL0 ratio 0.389/0.382/0.381; contiguity of 29+43 cells + final cell verified; sweep cells 0/20/41/42 rerun exact match (nsub=4) |

## Overall verdict on the theorem

The STATEMENT |M_k^main + (sqrt(pi)/N)e^{-sqrt(pi N)}cos sqrt(pi N)| <= C N^{-3/2} e^{-sqrt(pi N)}
for all N >= 100 with C = 617117.6 **stands** (dominated by the unaffected sweep; truth max is ~2).
The final-cell certificate as written is **not a valid proof** of its stated 49.5626 (TM piece
provably smaller than the truth, 0.12897 < 0.14412 at N=1e10); the one-character fix
(sqrt(2/pi) instead of sqrt(sqrt(2)/pi)) repairs it with corrected C = 49.9491.

## Files

step_01.py ... step_40.py (this directory): exploration (01-12); L1 numerics (13-15); L2 + Gh
sympy (16); |g| decreasing + beta (17); L4 (18); L5 (19); L6 (20); L7 factorisations (21);
mS inexactness + box membership proof (22-23); D13h step_08/09 review (24-25); A1 phase-box bug
discovery + geometry scan (26-27); D13h breg check (28); exact Btrue regular form + z identity +
leading term (29-31); final-cell breg validity + sympy decomposition (32); sweep z-box containment
(33); corrected identity confirmations (34); breg margin (35); final-cell C rerun (36); TM violation
dissection: true vs certified, pi/2 bookkeeping (37-38); ecw bug isolated: (rt2/pi).sqrt()=0.6709382
not 0.7978846 (39); TM exact at 1e8..1e14 + TR/TLm truths + bug in all three step_10 files +
corrected C (40); y0 enclosure + saddle monotonicity sympy + sweep spot-checks (41-42).
progress.log: per-lemma log.
"""
open(work+'/result.md','w').write(result)
with open(work+'/progress.log','a') as f:
    f.write("FINAL: main theorem C=617117.6 stands (sweep-dominated); final-cell ecw bug found (TM violated 0.12897<0.14412 @N=1e10), corrected C=49.9491; A1m phase-box and mS/breg derivation slips documented as harmless. result.md written.\n")
print("result.md written,", len(result), "chars")
