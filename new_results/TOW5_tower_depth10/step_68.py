proof = r'''# TOW5 — proof document: the tower law at depths 9 and 10 (PROVED); the weight-seven law (PROVED for every k >= 7)

Inherited from TOW2 (delta <= 4), TOW3 (delta <= 6), TOW4 (delta <= 8), all PROVED and re-verified
by referees (TOW2R, CX8, TOW3R, TOW4R, CX10): the setup; Lemma T; Lemmas S1/S4/S3/S5 (R_1..R_4 =
N_l/Pi_l for all (d,m); `prove_S12.py`, `prove_S34.py`); the Bell theorem rho(delta)/e00 = a_delta =
[y^delta] exp(sum_{r>=1} L_r y^r) with L_r built (Faa di Bruno, blocks of size >= 2 since t'(i) = 0)
from sl_l, l <= floor(r/2), plus explicit algebraic pieces; the zigzag top coefficients; the TOW1
reduction E0 + {N(delta,s)} + E3 ==> tower law; the tau constants in Q(i); the adjacent machinery
(fold + ehat; e00_pp anti-real so Dc = -conj(Delta); e00_{p,p-1} real — CX8 #2 carried forward); the
FORCED-recursion equivalence (TOW3 M1): F(delta,0) := a_delta, F(delta,s+1) := Lambda^{(delta)}(s) -
F(delta,s) from rho alone reproduces every closed form, and the tower law at depth delta is equivalent
to the single closing condition F(delta, floor(delta/2)+1) = 0 together with E3. The TOW4R/CX10 errata
are applied from the start (F(7,3) = +5i(d-5)/Pi_3; adjacent domain delta = 7..10 is p >= 4; at even
delta the closing identity is blind to the seed a_delta, so rho(delta) at even delta rests on the proved
initial-data derivation; k = w-1 is an exact finite kernel check; verify scripts check return codes).
Referee normalisations (CX8/TOW3R): everything is normalised by rho(0) = e00 (off-diagonal) resp.
rho_pp(0) (diagonal); the off-diagonal e vs the stable diagonal e-hat are kept distinct; uniqueness is
by finite descending induction in degree at fixed k; zigzag = sec + tan (A_n positive Euler zigzag
numbers, sech x + tanh x with absolute values in the code).

Domains (stated explicitly): depth delta needs 2d-4 >= delta (delta = 9, 10 => d >= 7); adjacent
(p,p-1) needs 4p-6 >= delta (delta = 9, 10 => p >= 4); k >= d off-diagonal; k >= 2p+1 diagonal stable,
k = 2p exceptional. In these domains every feed/Pi denominator is a nonzero odd integer, so
numerator-0 rational identities prove the identities throughout.

## 1. R_5, then rho(9), rho(10) — PROVED closed forms (Task A)

**R_5 = S^{(5)}(t*)/S(t*) = N_5/(3 Pi_5), PROVED for all (d,m)** (`prove_S5.py`). Method = TOW2's
`prove_S34.py` extended to l = 5: the normalized derivative is a finite sum (126 terms,
n = 0..5, r = 0..n, a = n..5, c = 0..5-a) of explicit rational functions of (d,m) over Q(i)
(TOW2 proof_all_depths.md Sec. 1, proved method). The single new ingredient is the Stirling tail
  P_5(x) = S(x, x-5) = x(x-1)(x-2)(x-3)(x-4)^2 (x-5)^2 (3x^2 - 23x + 38)/11520,
PROVED by the same recurrence + anchors as TOW2: (i) P_r(x+1) = (x+1-r) P_{r-1}(x) + P_r(x) as exact
polynomial identities for r <= 5 (symbolic, difference 0); (ii) P_r(a) = 0 for a = 0..r; hence
P_r(n) = S(n, n-r) for all integers n >= 0 (both sides satisfy a_r(n+1) = (n+1-r) a_{r-1}(n) + a_r(n)
by the Stirling recurrence S(n+1,k) = k S(n,k) + S(n,k-1) at k = n+1-r, and agree at n = 0; induction).
The identity total(5) == N_5/(3 Pi_5) is PROVED symbolically (sp.cancel of the difference is exactly 0;
the gcd-reduced denominator of the total is exactly 3 Pi_5 — sp.together alone under-cancels, the gcd
reduction gives 3 Pi_5; the factor 3 is the P_5 denominator 11520 = 2^8.3^2.5 residue). N_5 has degree
10 = 2l in d, 35 terms, leading coefficient -3072 i d^10 (matching 1/(2^5 5!) of the Stirling tail).
CONFIRMED independently: direct exact differentiation of S(t) = sum_j C_j sum_s Stir2(p,s) Stir2(q,j-s)
t^{2(j-s)} (1-t)^{2s} (t^2+(1-t)^2)^{d-j} at t* = (1+i)/2 vs R_5(d,m=p): **24/24 integer keys
d = 6..10 and 67/67 keys at d in {12,15,18,22,26,30}, residual 0**.

  N_5 = -3072 i d^10 - 7680 d^9 + 64000 i d^9 + 15360 d^8 m + 153600 d^8 - 542720 i d^8
        - 307200 d^7 m - 30720 i d^7 m - 1246720 d^7 + 2408448 i d^7 + 30720 i d^6 m^2
        + 2478080 d^6 m + 537600 i d^6 m + 5304064 d^6 - 5980416 i d^6 + 46080 d^5 m^2
        - 537600 i d^5 m^2 - 10403328 d^5 m - 3640320 i d^5 m - 12756736 d^5 + 8249088 i d^5
        - 30720 d^4 m^3 - 614400 d^4 m^2 + 3624960 i d^4 m^2 + 24581120 d^4 m + 12057600 i d^4 m
        + 17382400 d^4 - 6060480 i d^4 + 409600 d^3 m^3 + 30720 i d^3 m^3 + 2795520 d^3 m^2
        - 11942400 i d^3 m^2 - 33091840 d^3 m - 20305920 i d^3 m - 12473440 d^3 + 2752352 i d^3
        - 15360 i d^2 m^4 - 1858560 d^2 m^3 - 230400 i d^2 m^3 - 5018880 d^2 m^2 + 20136960 i d^2 m^2
        + 23950400 d^2 m + 16060800 i d^2 m + 3663136 d^2 - 1734112 i d^2 - 7680 d m^4 + 115200 i d m^4
        + 3345920 d m^3 + 337920 i d m^3 + 2989440 d m^2 - 16060800 i d m^2 - 7326272 d m
        - 4350720 i d m - 18624 d + 846912 i d + 3072 m^5 - 168960 i m^4 - 1992960 m^3
        + 4350720 i m^2 + 37248 m

**sl_5 = (log S)^{(5)}(t*)** by exact cumulant algebra from the PROVED R_1..R_5:
sl_5 = 24 R_1^5 - 60 R_1^3 R_2 + 20 R_1^2 R_3 + 30 R_1 R_2^2 - 5 R_1 R_4 - 10 R_2 R_3 + R_5; its
reduced denominator is 3 (2d-11)(2d-9)(2d-7)(2d-5)^2 (2d-3)^5 (only (2d-odd) factors; the (d-j)
factors of the raw total cancel in the cumulant combination).

**a_9 = Bell_9(L_1..L_9), a_10 = Bell_10(L_1..L_10)** (a_9 uses sl_l, l <= 4 = proved R_1..R_4;
a_10 uses sl_5 = proved R_5). The L cache (L_simplified 1..8 + L910_candidate 9,10) is verified
EQUAL to the proved tow2lib-style L(r) construction (tow2lib.Lr with sl_5) **symbolically** in
Q(i,sqrt2)(d,m) (sqrt2-aware reduction: numerator sqrt2 -> w, both w-parity parts exactly 0) for
**r = 1..10** (`verify_L.py`; r = 9: 82s, r = 10: 61s) — closing the r = 9, 10 gap, so rho(9), rho(10)
rest entirely on proved ingredients. tau constants extended to n = 10 (land in Q(i), sqrt2 cancels;
verified EXACTLY vs the tow2lib definitions, difference 0, plus 40-digit numeric ~ 5e-41; the symmetry
conj(tau1_n) = (-1)^n tau2_n holds at n = 8, 9, 10):
  tau1(9) = -99/256 - 297 i/256,  tau2(9) = 99/256 - 297 i/256;
  tau1(10) = 1 - 99 i/256,        tau2(10) = 1 + 99 i/256.

  a_9  = N_9  / (5806080 (2d-9)(2d-7)(2d-5)(2d-3))          [N_9: deg 13 in d, 40 terms, lc -256 i d^13]
  a_10 = N_10 / (116121600 (2d-11)(2d-9)(2d-7)(2d-5)(2d-3)) [N_10: deg 15 in d, 51 terms, lc -1024 d^15]
(Pi-parts Pi_4, Pi_5 = prod_{j=1}^l (2d-2j-1), matching a_delta ~ Pi_{floor(delta/2)}; the higher
(2d-odd) powers in sl_5's denominator cancel in the Bell assembly. Full numerators in a9_closed.pkl /
a10_closed.pkl; leading terms N_9 = -256 i d^13 + 12288 i d^12 - 4608 d^11 - 273920 i d^11 + ...,
N_10 = -1024 d^15 + 61440 d^14 + (-1726720 + 23040 i) d^13 + ....)

*Exact check vs step-B data (trunc_engine, D = 10, exact Q(sqrt2,i)).* rho(delta)/rho(0) = a_delta,
rho(0) = e00 (the 1+[p>q] factor inside e00; initial level k0 = d for ALL keys; the sqrt2-part of the
ratio asserted 0 in every cell): **425 cells each for delta = 9, 10, all keys with 2d-4 >= delta,
d <= 40; 0 bad.** (`data_rho_D10_d7_40.pkl`; the 425-cell rho(10) check is also a 425-fold independent
confirmation of the PROVED R_5 through sl_5 and L(10).)

## 2. The e-family at depths 9, 10 — the FORCED family (fit-free), and the closing conditions (Task C, non-adjacent)

The e-family is the FORCED ascending family, no fit anywhere:
  F(delta,0) := a_delta (proved above);  F(delta,s+1) := Lambda^{(delta)}(s) - F(delta,s),
  Lambda^{(delta)}(s) = i (2d-1) sum_{dp=0}^{delta} [ tau1(delta-dp) F(dp,s)(d+1,m+1)
                       - tau2(delta-dp) F(dp,s)(d+1,m) ] / (2d-dp-1).
The Newton identities N(delta,s), e(delta,s)+e(delta,s+1) = Lambda^{(delta)}(s), hold BY CONSTRUCTION
for s < floor(delta/2). The remaining identity at depth delta is the closing condition
F(delta, floor(delta/2)+1) = 0 (equivalently the top equation F(delta,top) = Lambda^{(delta)}(top)).

*Theorem (closing conditions, all (d,m)).* **F(9,5) = 0 and F(10,6) = 0 as exact rational identities
in (d,m) over Q(i)**: `forced910.py`, sp.together + numerator exactly 0 (built on the PROVED delta <= 8
family; F(9,1) 2.4s, F(10,1) 4.4s; closing checks F(9,5), F(10,6) numerator 0; full family delta <= 10
cached in forced_F910.pkl). Hence the non-adjacent Newton identities at depths 9 (s = 0..4) and 10
(s = 0..5) all hold.

*Sanity (the M1 claim reproduced).* F(delta,s) == tow3lib.e(delta,s) for delta <= 6, all s (exact
difference 0); F(delta, floor(delta/2)) == zigzag top for delta <= 10 (exact difference 0).

The new closed forms (denominators exactly Pi_{floor(delta/2)} times small integers):
  F(9,1) = [deg 10 in d, 40 terms] / (40320 Pi_4)
  F(9,2) = ( -16 i d^7 + 528 i d^6 + (-80 - 7752 i) d^5 + 160 d^4 m + (1680 + 65160 i) d^4 - 3360 d^3 m
             + (-13660 - 335369 i) d^3 + (27320 - 240 i) d^2 m + (50100 + 1042137 i) d^2 + 240 i d m^2
             + (-100200 + 1440 i) d m + (-68040 - 1774488 i) d - 1440 i m^2 + 136080 m + 1246500 i ) / (240 Pi_4)
  F(9,3) = ( -20 i d^4 + 420 i d^3 + (-30 - 3457 i) d^2 + 60 d m + (180 + 12768 i) d - 360 m - 16956 i ) / (12 Pi_4)
  F(9,4) = -16 i (d-6) / Pi_4    (zigzag top, (-1)^5 A_5 (d-6)/Pi_4, A_5 = 16)
  F(10,1) = [deg 12 in d, 54 terms] / (645120 Pi_5)
  F(10,2) = [deg 9 in d, 42 terms] / (2880 Pi_5)
  F(10,3) = ( -80 d^6 + 2640 d^5 + (-37496 + 240 i) d^4 - 480 i d^3 m + (290096 - 4320 i) d^3
              + 8640 i d^2 m + (-1271513 + 26508 i) d^2 + (-240 - 53016 i) d m + (2937427 - 53532 i) d
              + 240 m^2 + 107064 i m - 2733468 ) / (96 Pi_5)
  F(10,4) = ( -16 d^3 + 288 d^2 + (-1777 + 8 i) d - 16 i m + 3554 ) / Pi_5
  F(10,5) = -61 / Pi_5           (zigzag top, (-1)^5 A_6/Pi_5, A_6 = 61)
(Full numerators of F(9,1), F(10,1), F(10,2) in forced_F910.pkl; tops proved == zigzag by the sanity
check above.)

## 2b. Adjacent-key identities at depths 9, 10 — PROVED (Task C, adjacent)

*Theorem (all p).* For the adjacent key (p, p-1) (d = 2p-1) and (delta,s) in { (9,0),(9,1),(9,2),
(9,3),(9,4), (10,0),(10,1),(10,2),(10,3),(10,4),(10,5) }, the identity e(delta,s) + e(delta,s+1) =
Lambda^{(delta)}(s) holds with the diagonal ehat-conversion and the fold, as an exact rational identity
in p over Q(i): `prove_adj910_fast.py`, exact numerator 0 via the lcm-denominator polynomial check
(cache cancelled ehat/Ediag fractions, then clear the lcm denominator and expand the exact polynomial
numerator; method validated at delta <= 8 by TOW4/TOW4R). Timings: (9,0) 50s, (9,1) 24s, (9,2) 10s,
(9,3) 2.6s, (9,4) 0.3s, (10,0) 86s, (10,1) 42s, (10,2) 20s, (10,3) 6.2s, (10,4) 1.0s, (10,5) 0.1s.
The adjacent TOP identities (9,4) [odd, s=4] and (10,5) [even, s=5] are NEW (TOW4 proved adjacent tops
at (7,3),(8,4)). Sanity: reproduces the delta <= 8 adjacent identities (same script family).

## 2c. Blindness E3 at delta = 9, 10 — PROVED (Task C, E3)

*Theorem (all p).* The feed of Delta_pp(delta') = rho_pp(delta') - sum_s (-1)^s ehat_pp(delta',s)
(delta' <= delta) into the child (p, p-1) at depth delta (DL2-direct + fold) vanishes identically:
`prove_E3b910.py`, exact numerator 0 (delta = 9: 13.9s; delta = 10: 22.9s; method validated at
delta <= 8). e00_pp is anti-real, so Dc = -conj(Delta); e00_{p,p-1} is real.

## 3. Conclusion: the tower law at depths <= 10 — PROVED

By TOW1's reduction theorem: E0 (PROVED, WT2b) + the identities N(delta',s) for delta' <= 10
(Section 2, all key types; the diagonal law follows from the off-diagonal parent via the
ehat-conversion) + blindness E3(delta') for delta' <= 10 (Section 2c) ==> the tower law holds at
depths delta <= 10 for all keys and all k (off-diagonal k >= d; diagonal k >= 2p+1 with k = 2p
exceptional). The descending finite-induction uniqueness argument (TOW2R/CX8) gives U = u: for each
fixed integer k, start at degree k with the true initial data and descend; at off-diagonal initial
keys the forced family agrees with rho (Section 1), at stable keys the Newton identities make it
satisfy T1, at diagonal initial keys the actual rho is retained and at k = 2p E3 gives the stable
extrapolation; the triangular recursion has a unique solution. Exact confirmation vs step-B data:
**990 cells, 0 bad** (all keys of degree 7..12, delta <= min(10, 2d-4), levels k = 7..24;
off-diagonal via the forced family, diagonal via ehat; plus **11 diagonal-initial-level cells k = 2p
vs a_delta(2p,p), 0 bad**). `check_full910.py`, data `data_towers_D10_K24.pkl` (engine validated
against TOWN1 at delta <= 8: 162/162 rho cells, 36/36 tower cells).

## 4. Task D — the weight-seven law for every k >= 7, PROVED

Specializing to degree 7: the tower law to depth 10 = 2*7-4 implies, by WTS Lemma 3 (PROVED; tower
law for degree n = 7 at depths [0, 2*6-2] = [0,10]), the kernel law K(6): every degree-6 kernel
residue is a polynomial in k of degree <= 5. **The range is k >= 7, not merely k >= 8: degree 7 is
odd, so there is NO diagonal key of degree 7** (the four keys (7,0),(6,1),(5,2),(4,3) are all
off-diagonal; the kernel feeds (WTS Lemma 2) of every degree-6 work item involve only these
off-diagonal degree-7 towers — W_{ij} = -(i+1) Gm_{i+1,j} DL1 - [i>j](j+1) Gm_{i,j+1} DL2 -
[j=i-1] i F(Gm_ii DL1), and j = i-1 never holds for i+j = 6, so no fold term and no diagonal (4,4)
feed at degree 7 either). The proved tower law for off-diagonal keys holds from k = d = 7, so the
single extra level k = d+1 = 7 that Lemma 3 leaves "covered by the data, not proved" for even n is
here INSIDE the proved range. (This answers the Task-D question: there is no degree-7 diagonal key;
nothing is exceptional at k = 7.) W_7(k) is a k-independent fixed real-linear image of these residues
(WTS Lemma 2, PROVED; the fold makes the map real-linear in the +i tower coordinates and their
conjugates, CX8 #10 carried forward), hence a quintic in k for k >= 7. The +-i simple-pole
coefficients are residues, hence quintics for k >= 7; they vanish at k = 8..13 (6 points > degree 5),
hence identically zero for k >= 7.

Exact certificates (`verify_W7.py`, WTS `kern_k.pkl` ground-truth exact kernels, degree-6 keys
(6,0),(5,1),(4,2),(3,3), real poles {0, a = 1+sqrt2, -a, -b = 1-sqrt2}, QSI components):
- W7a: the recorded degree-6 kernels at **k = 7** equal the quintic fitted from k = 8..13 at every
  pole, every QSI component, all four keys: **56/56 components, 0 bad** (confirms the sharpened
  Lemma-3 range k >= 7 exactly).
- W7b: 6th finite differences of every residue-component vanish on k = 7..13 (56/56; quinticity).
- W7c: the +-i simple-pole coefficients (ipos/ineg) vanish at k = 7..13 (all four keys).
- W7d: **k = 6 is OFF the quintic**: kern_6 degree-6 residues differ from the quintic(k = 8..13) at
  **18/56 components** (the initial kernel, mirroring k = w-1 exceptional at every lower weight).
Hence
    W_7(k) = sum_{m=1..6} C(k-7, m-1) E_{6,m}
           = E_{6,1} + (k-7) E_{6,2} + C(k-7,2) E_{6,3} + C(k-7,3) E_{6,4} + C(k-7,4) E_{6,5} + C(k-7,5) E_{6,6}
    for every k >= 7
(a quintic in k; constants E_{6,1}..E_{6,6} determined by the weight-7 kernel data, k = 7..12). Weights
polynomial in k for ALL k >= w: w = 2,3,4,5 (previous jobs), 6 (TOW4), **7 (here)**; w >= 8 open
(needs delta >= 12).

## 5. Task E — toward all depths (the master x = -1 specialisation at delta = 9, 10)

See result.md. (Recorded as attempted; no overclaim.)
'''
work = '<project>/research/08_theory/TOW5_tower_depth10'
open(work+'/proof_depth10.md','w').write(proof)
print("proof_depth10.md written:", len(proof), "chars")
