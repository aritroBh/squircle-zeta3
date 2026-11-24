result = r'''# TOW5 — the tower law at v = i: **PROVED TO DEPTH 10** (all d, all k); the weight-SEVEN polynomial law PROVED for every k >= 7; Task E: reduction verified at delta = 9, 10 (no overclaim)

## VERDICT

**PROVED TO DEPTH 10.** The tower law "the depth-delta tower coefficient u_{pq}^{(k)}(delta) at v = +i
is a polynomial in k of degree floor(delta/2)" is proved for depths **delta <= 10**, all keys (p,q) of
all degrees d, all k (off-diagonal: k >= d; diagonal: k >= 2p+1 with k = 2p exceptional), extending
TOW4 (delta <= 8). The proof is **FIT-FREE**: the e-family is the FORCED ascending recursion
F(delta,0) := rho(delta)/e00 = a_delta, F(delta,s+1) := Lambda^{(delta)}(s) - F(delta,s), built
entirely from the PROVED initial data a_delta (Bell assembly from the proved R_1..R_5 — **R_5 is proved
here**, closing the last gap) and the proved tau constants (extended to n = 10). The Newton identities
N(delta,s) for s < floor(delta/2) hold by construction; the single remaining identity per depth is the
closing condition F(delta, floor(delta/2)+1) = 0 (TOW3 M1 pre-screen equivalence), PROVED as an exact
rational identity in (d,m): **F(9,5) = 0 and F(10,6) = 0** — together with the adjacent-key identities
(fold + e-hat, exact identities in p) and blindness E3 at depths 9 and 10.

**Consequence (Task D): the weight-seven polynomial law for every k >= 7.** Tower law at degree 7 to
depth 10 = 2*7-4 (proved) =>(WTS Lemma 3)=> kernel law K(6): degree-6 kernel residues are polynomials
in k of degree <= 5. **The range is k >= 7, not merely k >= 8: degree 7 is odd, so there is NO
diagonal key of degree 7** (all four keys (7,0),(6,1),(5,2),(4,3) are off-diagonal; and j = i-1 never
holds for i+j = 6, so no degree-6 work item has a diagonal fold feed at degree 7 either). The extra
level k = d+1 = 7 that Lemma 3 leaves unproved for even n is therefore INSIDE the proved range — this
answers the Task-D question (there is no degree-7 diagonal key; nothing is exceptional at k = 7).
=>(WTS Lemma 2, fixed real-linear image)=> W_7(k) is a quintic in k for k >= 7. The +-i simple-pole
coefficients are residues, hence quintics for k >= 7, and vanish at k = 8..13 (6 points > degree 5),
hence identically for k >= 7. The level k = 7 is settled EXACTLY: recorded degree-6 kernels at k = 7
equal the quintic fitted from k = 8..13 at every pole, every QSI component, all four keys
(**56/56 components, 0 bad**). Hence
  **W_7(k) = sum_{m=1..6} C(k-7, m-1) E_{6,m} = E_{6,1} + (k-7) E_{6,2} + C(k-7,2) E_{6,3}
    + C(k-7,3) E_{6,4} + C(k-7,4) E_{6,5} + C(k-7,5) E_{6,6} for every k >= 7**
(constants E_{6,1}..E_{6,6} determined by the weight-7 kernel data). **k = 6 is OFF the quintic**:
kern_6 degree-6 residues differ from the quintic at **18/56 components** (the initial kernel, as at
every lower weight). Weights polynomial in k for ALL k >= w: w = 2,3,4,5 (previous), 6 (TOW4),
**7 (here)**; w >= 8 open (needs delta >= 12).

**Task E: the all-depth closing condition is REDUCED to the master x = -1 identity (TOW4's proved
equivalence) and that identity is now verified symbolically at delta = 9 (7.4s) and delta = 10 (52.7s);
the all-depth proof is NOT achieved (same obstruction as TOW4; no overclaim).**

## Task A — R_5, rho(9), rho(10) closed forms, PROVED

**R_5 = S^{(5)}(t*)/S(t*) = N_5/(3 Pi_5), PROVED for all (d,m)** (`prove_S5.py`): TOW2's prove_S34
method extended to l = 5 (126-term finite sum, proved method). The new Stirling tail
  P_5(x) = S(x, x-5) = x(x-1)(x-2)(x-3)(x-4)^2 (x-5)^2 (3x^2-23x+38)/11520
is PROVED by the same recurrence (P_r(x+1) = (x+1-r) P_{r-1}(x) + P_r(x), exact polynomial identities,
r <= 5) + anchors (P_r(a) = 0, a = 0..r), hence P_r(n) = S(n,n-r) for all n >= 0 (induction via the
Stirling recurrence). The identity total(5) == N_5/(3 Pi_5) is proved symbolically (sp.cancel of the
difference is exactly 0; the gcd-reduced denominator is exactly 3 Pi_5 — note sp.together alone
under-cancels; N_5 has degree 10 = 2l in d, 35 terms, leading -3072 i d^10). CONFIRMED by direct exact
differentiation of S(t) at **24/24 integer keys d = 6..10 and 67/67 keys at d in {12,15,18,22,26,30}**.

a_9 = Bell_9(L_1..L_9), a_10 = Bell_10(L_1..L_10) (a_9 via proved R_1..R_4; a_10 via sl_5 from the
proved R_5 by exact cumulant algebra). The L cache is verified EQUAL to the proved tow2lib-style L(r)
**symbolically** (sqrt2-aware, difference 0) for **r = 1..10** (`verify_L.py`). tau extended to n = 10:
tau1(9) = -99/256 - 297i/256, tau2(9) = 99/256 - 297i/256, tau1(10) = 1 - 99i/256, tau2(10) = 1 +
99i/256 (land in Q(i), sqrt2 cancels; exact vs defs, difference 0, plus 40-digit numeric ~ 5e-41;
conj(tau1_n) = (-1)^n tau2_n holds at n = 8, 9, 10).

  a_9  = N_9  / (5806080 (2d-9)(2d-7)(2d-5)(2d-3))           [N_9: deg 13 in d, 40 terms, lc -256 i d^13]
  a_10 = N_10 / (116121600 (2d-11)(2d-9)(2d-7)(2d-5)(2d-3))  [N_10: deg 15 in d, 51 terms, lc -1024 d^15]
(Pi-parts Pi_4, Pi_5, matching a_delta ~ Pi_{floor(delta/2)}; full numerators in a9_closed.pkl /
a10_closed.pkl / proof_depth10.md.)

Exact check vs step-B data: rho(delta)/rho(0) = a_delta, rho(0) = e00 (1+[p>q] inside e00; initial
level k0 = d for ALL keys; sqrt2-part of ratio asserted 0): **425 cells each for delta = 9, 10 (all
keys 2d-4 >= delta, d <= 40), 0 bad.** The rho(10) check is a 425-fold independent confirmation of the
proved R_5.

## Task B — exact tower data at depths 9, 10

TOWN1's trunc_engine copied and run at D = 10 (exact Q(sqrt2,i); validated against TOWN1 at
delta <= 8: 162/162 rho cells, 36/36 tower cells): (i) rho initial data, degrees 7..40
(`data_rho_D10_d7_40.pkl`, 11.7s); (ii) full towers, levels k = 7..24, degrees 7..12
(`data_towers_D10_K24.pkl`, 6.2s; 1001 cells: 91 keys with depth 9, 91 with depth 10).

## Task C — Newton identities at depths 9, 10, PROVED

**Non-adjacent (d,m symbolic).** FORCED family from rho alone (fit-free); N(delta,s) for
s < floor(delta/2) hold by construction; the closing conditions are the only further identities:
**F(9,5) = 0 and F(10,6) = 0, sp.together + numerator exactly 0** (`forced910.py`; full family
delta <= 10 built and cached in 54s). Sanity: F == tow3lib.e (delta <= 6 all s); F(top) == zigzag
(delta <= 10). New closed forms: F(9,4) = -16i(d-6)/Pi_4 (zigzag top, A_5 = 16), F(10,5) = -61/Pi_5
(A_6 = 61); F(9,1..3), F(10,1..4) in proof_depth10.md / forced_F910.pkl.

**Adjacent (p symbolic).** The 11 identities (9,0..4) and (10,0..5) hold with the diagonal
ehat-conversion and the fold, exact numerator 0 (`prove_adj910_fast.py`; times 0.1s--92s). The
adjacent TOP identities (9,4), (10,5) are NEW.

**Blindness E3.** The feed of Delta_pp(delta') (delta' <= delta) into the child (p,p-1) at depth
delta (DL2-direct + fold) vanishes identically at delta = 9, 10, exact numerator 0
(`prove_E3b910.py`, 13.9s/22.9s).

**Full-law confirmation vs step-B data:** **990 cells, 0 bad** (all keys degree 7..12,
delta <= min(10,2d-4), levels k = 7..24) **+ 11 diagonal-initial-level cells k = 2p, 0 bad**
(`check_full910.py`).

## Task D — the weight-seven law for every k >= 7, PROVED

(see the VERDICT). Exact certificates (`verify_W7.py`, WTS kern_k.pkl ground-truth exact kernels,
degree-6 keys (6,0),(5,1),(4,2),(3,3), real poles {0, a = 1+sqrt2, -a, -b = 1-sqrt2}): W7a k = 7 ==
quintic(8..13) **56/56**; W7b 6th differences vanish k = 7..13 (56/56); W7c +-i simple-pole
coefficients vanish k = 7..13; W7d k = 6 off the quintic **18/56** (initial kernel). The paper's
range may be changed to "W_w(k) polynomial for every k >= w" for w = 2..7.

## Task E — toward all depths (the master x = -1 specialisation at delta = 9, 10)

TOW4's PROVED reduction (closing at delta ⟺ the finite alternating identity a_delta = feed(U(-1)),
U(dp) = sum_s (-1)^s F(dp,s) = U-hat(dp, x = -1) ⟺ the master equation at x = -1 ⟺ the value
identity V(delta)) is unchanged; the master identity feed(U(-1)) - a_delta == 0 is now verified
symbolically at **delta = 9 (7.4s) and delta = 10 (52.7s)**. The structural evidence continues: at
delta = 9, 10, U(delta) has exactly the same denominator as a_delta and coincides in the top-3
leading numerator terms (first differing at d^10 resp. d^12, i.e. 3 below the top, as at delta = 7, 8).
No closed form for U(delta) was found and the identity remains self-referential at general delta (the
dp = delta term re-introduces U(delta) at the parents), so the all-depth proof is NOT achieved — same
obstruction as TOW4, now confirmed two depths higher. No overclaim: what is proved is the equivalence
and its validity at delta = 9, 10.

## Evidence table (all checks EXACT; residual = numerator/remainder EXACTLY 0)

| lemma | statement | certificate | residual | time |
|---|---|---|---|---|
| TAU910 | tau1/tau2(9),(10) in Q(i), exact vs defs; conj symmetry | step (tau910.pkl) | diff 0 (+40-digit numeric ~5e-41) | 0.2s |
| R5 | R_5 = N_5/(3 Pi_5), all (d,m); P_5 proved | prove_S5.py | cancel'd diff 0; direct-diff 24/24 + 67/67 | 10s |
| L | L cache == tow2lib-style L(r), r = 1..10, symbolic | verify_L.py | sqrt2-aware diff 0 | 161s |
| RHO | a_9,a_10 Bell; rho(9)/e00==a_9, rho(10)/e00==a_10 | step (rho check) | 425 cells each, 0 bad | 7s/54s |
| CLOSE | closing F(9,5)=0, F(10,6)=0 in (d,m) | forced910.py | numerator 0 | 54s total |
| CLOSE-sanity | F==e (delta<=6 all s); F(top)==zigzag (delta<=10) | step | diff 0 | 6s |
| ADJ | 11 adjacent identities (9,0..4),(10,0..5) in p | prove_adj910_fast.py | numerator 0 | 0.1--92s each |
| E3 | blindness delta = 9, 10 in p | prove_E3b910.py | numerator 0 | 13.9s/22.9s |
| LAW | tower law delta<=10 vs step-B data | check_full910.py | 990 cells 0 bad + 11 diag-init 0 bad | 16s |
| W7 | k=7 vs quintic(8..13); quinticity 7..13; +-i vanish; k=6 off | verify_W7.py | 56/56 comp; k=6: 18/56 off | 0.2s |
| MASTER | closing ⟺ master x=-1 coefficient, delta = 9, 10 | step (Task E) | feed(U)-a_delta == 0 | 7.4s/52.7s |

Digit counts: every certificate is exact rational arithmetic over Q(i) (sympy) or exact Fraction
Q(sqrt2,i) arithmetic (data checks); residuals are exactly 0, no floating point anywhere except the
independent 40-digit mpmath confirmation of tau(9), tau(10), additionally proved symbolically
(difference 0). Data cells: rho 2 x 425; full law 990 + 11; weight-7 4 keys x 56 components x 4 checks.
Domains: depth delta needs 2d-4 >= delta (delta = 9, 10 => d >= 7); adjacent (p,p-1) needs
4p-6 >= delta (delta = 9, 10 => p >= 4); k >= d off-diagonal, k >= 2p+1 diagonal stable, k = 2p
exceptional; all feed/Pi denominators are nonzero odd integers in these domains, so numerator-0
rational identities prove the identities throughout.

## Files (research/08_theory/TOW5_tower_depth10/)

- `tow5lib.py` (R_5 via N5_proved; sl_5; tau to n = 10; L cache to r = 10; arho delta <= 10; forced
  loader), `prove_S5.py` (P_5 + R_5 = N_5/(3 Pi_5), PROVED + direct check), `verify_L.py` (L cache
  r = 1..10 symbolic), `forced910.py` (forced recursion + closing proofs, saves `forced_F910.pkl`),
  `prove_adj910_fast.py` (11 adjacent identities), `prove_E3b910.py` (blindness 9, 10),
  `check_full910.py` (990+11-cell law), `verify_W7.py` (Task D), `verify_tow5.py` (per-lemma
  PASS/FAIL), `N5_proved.pkl`, `tau910.pkl`, `L910_candidate.pkl`, `a9_closed.pkl`, `a10_closed.pkl`,
  `forced_F910.pkl`, `data_rho_D10_d7_40.pkl`, `data_towers_D10_K24.pkl`, `result.md`,
  `proof_depth10.md`, `progress.log`.
- Byte-inherited (copied, provenance preserved): `tow2lib.py`, `tow3lib.py`, `tow4lib.py`,
  `L_simplified.pkl`, `L_parsed.pkl`, `check_data.py`, `arho_simplified.pkl`, `prove_S34.py` (TOW2/4);
  `trunc_engine.py`, `driver.py` (TOWN1); `fastred.py` (WTS); `forced_F78.pkl` (TOW4). Data: WTS
  `theory/kern_*.pkl` (read-only).

**Honest boundaries:** the tower law at depths >= 11 is unproved (rho(11) needs R_1..R_5 = proved,
reachable; rho(12) needs R_6, not proved); the all-depth closing condition remains REDUCED to the
master x = -1 identity but NOT proved for all delta (verified at delta = 9, 10); the weight-eight law
(needs delta <= 12) is open.
'''
work = '<project>/research/08_theory/TOW5_tower_depth10'
open(work+'/result.md','w').write(result)
print("result.md written:", len(result), "chars")
with open(work+'/progress.log','a') as f:
    import time
    f.write(f"{time.strftime('%H:%M:%S')} step: proof_depth10.md + result.md written (A-D complete, E recorded)\n")
