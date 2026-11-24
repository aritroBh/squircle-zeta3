# TOW4 — proof document: the tower law at depths 7 and 8 (PROVED); the weight-six law (PROVED); all-depth attempt (Task D)

Inherited from TOW2 (delta <= 4) and TOW3 (delta <= 6), all PROVED there and re-verified by referees
(TOW2R, CX8, TOW3R): setup; Lemma T; Lemmas S1/S4/S3/S5 (R_1..R_4 = S^{(l)}(t*)/S(t*) = N_l/Pi_l for
all (d,m); `prove_S12.py`, `prove_S34.py`); the Bell theorem rho(delta)/e00 = a_delta = [y^delta]
exp(sum_{r>=1} L_r y^r) with L_r built (Faa di Bruno, only blocks of size >= 2 since t'(i) = 0) from
sl_l, l <= floor(r/2), plus explicit algebraic pieces; the zigzag top coefficients and the all-s proof
of the top equations for non-adjacent keys (`prove_tops.py`); the reduction (TOW1): E0 + {N(delta,s)} +
E3 ==> tower law; the tau constants tau1_n, tau2_n in Q(i), n <= 7; the adjacent machinery (fold +
ehat; **e00_pp is anti-real** so the actual conjugate in e00-units is Dc = -conj(Delta); e00_{p,p-1}
is real — CX8 #2 carried forward); and the **FORCED-recursion equivalence** (TOW3 claude_screen M1
pre-screen, `s03_forced.py`): the ascending recursion F(delta,0) := a_delta, F(delta,s+1) :=
Lambda^{(delta)}(s) - F(delta,s) from rho alone reproduces every closed form, and the tower law at
depth delta is equivalent to the single **closing condition** F(delta, floor(delta/2)+1) = 0 together
with E3. Domains (stated explicitly): depth delta needs 2d-4 >= delta (delta = 7, 8 => d >= 6);
adjacent (p,p-1) needs 4p-6 >= delta (p >= 3); k >= d off-diagonal; k >= 2p+1 diagonal stable,
k = 2p exceptional. In these domains every feed/Pi denominator is a nonzero odd integer.

## 1. rho(7), rho(8) — PROVED closed forms (Task A)

Ingredients: a_7 = Bell_7(L_1..L_7), a_8 = Bell_8(L_1..L_8) use sl_l, l <= 4 (R_1..R_4 PROVED). The
simplified cache L(1..8) (`L_simplified.pkl`, TOW2 via TOW3) is here verified EQUAL to the proved
tow2lib.L(r) construction **symbolically** (exact difference == 0) for **r = 1..8** (`verify_tow4.py`
lemma L; r = 8: 9.9s) — this closes the r = 7, 8 gap left at TOW3 (TOW3R's step_04 did r <= 6), so
rho(7), rho(8) rest entirely on proved ingredients. tau constants are extended to n = 8:
tau1(8) = -169/128 + 41 i/128, tau2(8) = -169/128 - 41 i/128 (sqrt2 cancels, they land in Q(i));
verified EXACTLY against the tow2lib definitions (difference 0) and to 40 digits numerically;
the symmetry conj(tau1_n) = (-1)^n tau2_n holds at n = 8.

  a_7 = N_7 / (40320 (2d-7)(2d-5)(2d-3)),  N_7 =
    64 i d^10 - 2048 i d^9 + (672 + 29808 i) d^8 - 1344 d^7 m + (-17472 - 259952 i) d^7
    + 34944 d^6 m + (196560 + 1503964 i) d^6 + (-393120 + 6720 i) d^5 m + (-1241520 - 6034952 i) d^5
    - 6720 i d^4 m^2 + (2483040 - 114240 i) d^4 m + (4758978 + 17057277 i) d^4
    + 114240 i d^3 m^2 + (-9514596 + 734160 i) d^3 m + (-11077248 - 33707943 i) d^3
    + (-10080 - 734160 i) d^2 m^2 + (22137696 - 2113440 i) d^2 m + (14484330 + 44731447 i) d^2
    + 6720 d m^3 + (50400 + 2113440 i) d m^2 + (-28968660 + 2293200 i) d m + (-8162700 - 35621565 i) d
    - 33600 m^3 - 2293200 i m^2 + 16325400 m + 12303900 i
  a_8 = N_8 / (645120 (2d-9)(2d-7)(2d-5)(2d-3)),  N_8 = [deg 12 in d, 57 terms; see `a8_closed.pkl`
    and `verify_tow4.py` output; leading terms 256 d^12 - 10752 d^11 + (208384 - 3584 i) d^10 + ...]

(Denominator Pi-parts Pi_3, Pi_4 = prod_{j=1}^l (2d-2j-1), matching a_delta ~ Pi_{floor(delta/2)}.)

*Exact check vs TOWN1 towers_D8_K52.pkl.* rho(delta)/rho(0) = a_delta, with rho(0) = e00 (the
1+[p>q] factor is inside e00; the initial level is k0 = d for ALL keys — for diagonal keys this is the
exceptional level k = 2p whose rho_pp is the true initial data; the sqrt2-part of the ratio is asserted
0 in every cell): **717 cells each for delta = 7, 8, all keys with 2d-4 >= delta, d <= 52; 0 bad.**
(`verify_tow4.py` lemma rho.) Sanity: reproduces TOW3's a_5, a_6 checks.

## 2. The e-family at depths 7, 8 — the FORCED family (fit-free), and the closing conditions (Task B, non-adjacent)

The e-family is the FORCED ascending family, no fit anywhere:
  F(delta,0) := a_delta (proved above);  F(delta,s+1) := Lambda^{(delta)}(s) - F(delta,s),
  Lambda^{(delta)}(s) = i (2d-1) sum_{dp=0}^{delta} [ tau1(delta-dp) F(dp,s)(d+1,m+1) - tau2(delta-dp) F(dp,s)(d+1,m) ] / (2d-dp-1).
The Newton identities N(delta,s), e(delta,s)+e(delta,s+1) = Lambda^{(delta)}(s), hold BY CONSTRUCTION for
s < floor(delta/2). The remaining identity at depth delta is the closing condition F(delta,
floor(delta/2)+1) = 0 (equivalently the top equation F(delta,top) = Lambda^{(delta)}(top)).

*Theorem (closing conditions, all (d,m)).* F(7,4) = 0 and F(8,5) = 0 as exact rational identities in
(d,m) over Q(i): `forced78.py`, sp.together + numerator exactly 0 (delta = 7: 0.1s; delta = 8: 0.0s;
the full family delta <= 8 built and cached in 47.9s). Hence the non-adjacent Newton identities at
depths 7 (s = 0..3) and 8 (s = 0..4) all hold.

*Sanity (the M1 claim reproduced).* The forced family coincides with the previously fitted/proved
forms: F(delta,s) == tow3lib.e(delta,s) for delta <= 6, all s (exact difference 0); and F(delta,
floor(delta/2)) == zigzag top for delta <= 8 (exact difference 0).

The new closed forms (denominators exactly Pi_{floor(delta/2)} times small integers):
  F(7,1) = ( 16 i d^7 - 416 i d^6 + (80 + 4840 i) d^5 - 160 d^4 m + (-1360 - 32520 i) d^4 + 2720 d^3 m
             + (8980 + 135069 i) d^3 + (-17960 + 240 i) d^2 m + (-26960 - 342004 i) d^2 - 240 i d m^2
             + (53920 - 1200 i) d m + (30300 + 479005 i) d + 1200 i m^2 - 60600 m - 278550 i ) / (480 Pi_3)
  F(7,2) = ( 4 i d^4 - 68 i d^3 + (6 + 455 i) d^2 - 12 d m + (-30 - 1375 i) d + 60 m + 1500 i ) / (6 Pi_3)
  F(7,3) = -2 i (d-5) / Pi_3   (zigzag top, A_4 = 5)
  F(8,1) = ( 64 d^9 - 2304 d^8 + (38064 - 480 i) d^7 + 960 i d^6 m + (-378384 + 12960 i) d^6 - 25920 i d^5 m
             + (2485836 - 149520 i) d^5 + 2880 d^4 m + 299040 i d^4 m + (-11129496 + 940800 i) d^4
             - 2880 d^3 m^2 + (-43200 - 1881600 i) d^3 m + (33687881 - 3386790 i) d^3 + 43200 d^2 m^2
             + (219600 + 6773100 i) d^2 m + (-65764401 + 6557550 i) d^2 + (-219600 + 1440 i) d m^2
             + (-372960 - 13115100 i) d m + (74061780 - 5276460 i) d - 960 i m^3 + 372960 m^2 + 10552920 i m
             - 35953740 ) / (5760 Pi_4)
  F(8,2) = ( 16 d^6 - 432 d^5 + (5032 - 48 i) d^4 + 96 i d^3 m + (-32056 + 720 i) d^3 - 1440 i d^2 m
             + (116233 - 3684 i) d^2 + (48 + 7368 i) d m + (-223115 + 6228 i) d - 48 m^2 - 12456 i m
             + 173022 ) / (48 Pi_4)
  F(8,3) = ( 20 d^3 - 300 d^2 + (1549 - 10 i) d + 20 i m - 2592 ) / (4 Pi_4)
  F(8,4) = 16 / Pi_4   (zigzag top, (-1)^4 A_5/Pi_4 with A_5 = 16)
[The tops are proved == zigzag by the sanity check above; no separate numerator claim is made here.]

## 2b. Adjacent-key identities at depths 7, 8 — PROVED (Task B, adjacent)

*Theorem (all p).* For the adjacent key (p, p-1) (d = 2p-1) and (delta,s) in { (7,0),(7,1),(7,2),(7,3),
(8,0),(8,1),(8,2),(8,3),(8,4) }, the identity e(delta,s) + e(delta,s+1) = Lambda^{(delta)}(s) holds with
the diagonal ehat-conversion and the fold, as an exact rational identity in p over Q(i):
`prove_adj78_fast.py`, exact numerator 0 via the lcm-denominator polynomial check. Timings:
(7,0) 56s, (7,1) 20s, (7,2) 6s, (7,3) 0.8s, (8,0) 97s, (8,1) 40s, (8,2) 15s, (8,3) 3s, (8,4) 0.3s.
(The method — cache cancelled ehat/Ediag fractions, then clear the lcm denominator and expand the exact
polynomial numerator — was validated against TOW3's slow sp.simplify method on (7,2),(8,3),(8,4),(7,3),
identical True.) Sanity: reproduces the delta <= 6 adjacent identities. The adjacent TOP identities
(7,3) [odd, s=3] and (8,4) [even, s=4] are NEW (TOW3 proved adjacent tops at (5,2),(6,3)).

## 2c. Blindness E3 at delta = 7, 8 — PROVED (Task B, E3)

*Theorem (all p).* The feed of Delta_pp(delta') = rho_pp(delta') - sum_s (-1)^s ehat_pp(delta',s)
(delta' <= delta) into the child (p, p-1) at depth delta (DL2-direct + fold) vanishes identically:
`prove_E3b78.py`, exact numerator 0 (delta = 7: 23.6s; delta = 8: 34.9s; method validated on
delta = 5, 6: True). e00_pp is anti-real, so Dc = -conj(Delta); e00_{p,p-1} is real.

## 3. Conclusion: the tower law at depths <= 8 — PROVED

By TOW1's reduction theorem: E0 (PROVED, WT2b) + the identities N(delta',s) for delta' <= 8 (Section 2,
all key types; the diagonal law follows from the off-diagonal parent via the ehat-conversion) +
blindness E3(delta') for delta' <= 8 (Section 2c) ==> the tower law holds at depths delta <= 8 for all
keys and all k (off-diagonal k >= d; diagonal k >= 2p+1 with k = 2p exceptional). The descending
finite-induction uniqueness argument (TOW2R/CX8) gives U = u: for each fixed integer k, start at degree
k with the true initial data and descend; at off-diagonal initial keys the forced family agrees with rho
(Section 1), at stable keys the Newton identities make it satisfy T1, at diagonal initial keys the actual
rho is retained and at k = 2p E3 gives the stable extrapolation; the triangular recursion has a unique
solution. Exact confirmation vs TOWN1 towers: **49953 cells, 0 bad** (all keys d <= 52,
delta <= min(8, 2d-4), k <= min(d+8, 52); off-diagonal via the forced family, diagonal via ehat; plus
222 diagonal-initial-level cells k = 2p vs a_delta(2p,p), 0 bad). `check_full78.py`.

## 4. Task C — the weight-six law for every k >= 6, PROVED

Specializing to degree 6: the tower law to depth 8 = 2*6-4 implies, by WTS Lemma 3 (PROVED; tower law
for degree n = 6 at depths [0, 2*5-2] = [0,8] for k >= 7), the kernel law K(5): every degree-5 kernel
residue is a polynomial in k of degree <= 4 for k >= 7. (At k >= 7 all degree-6 keys, including the
diagonal (3,3) stable from k = 2p+1 = 7, satisfy the tower law; k = 6 is (3,3)'s exceptional level and
is outside Lemma 3's range.) W_6(k) is a k-independent fixed real-linear image of these residues
(WTS Lemma 2, PROVED; the fold makes the map real-linear in the +i tower coordinates and their
conjugates, CX8 #10 carried forward), hence a quartic in k for k >= 7. The +-i simple-pole coefficients
are residues, hence quartics for k >= 7; they vanish at k = 7..11 (5 points > degree 4), hence
identically zero for k >= 7. The level k = 6 is settled EXACTLY: the recorded degree-5 kernels at k = 6
(`kern_6.pkl`, ground truth from the actual reduction) equal the quartic fitted from k = 7..11 at every
pole {0, a, -a, -b}, every QSI component, all three keys (5,0),(4,1),(3,2): **48/48 components, 0 bad**.
(Independent confirmation: 5th differences of every residue-component vanish on k = 7..13, 48/48; the
+-i poles vanish at k = 6..13.) Hence
    W_6(k) = sum_{m=1..5} C(k-6,m-1) E_{5,m}   for every k >= 6
(a quartic in k; constants E_{5,1}..E_{5,5} as determined by WT6S). k = 5 is exceptional, as always —
confirmed: kern_5 degree-5 residues differ from the quartic(k = 7..11) at 15/48 components (mirroring
WT5's P(4) != W_5(4)). This upgrades WT6S's exact-data law to a proved law for all k >= 6 (paper
prop:w6k).

## 5. Task D — toward all depths: the closing condition F(delta, floor(delta/2)+1) = 0

(Recorded as attempted; see result.md for the verdict. No overclaim: what is proved is stated with its
certificate; what is reduced/conjectured is labelled.)


## Referee errata (TOW4R + CX10, both PROOF COMPLETE; recorded by the M1 session 2026-10-28; the text above is unchanged)
F(7,3) = +5i(d-5)/Pi_3 (printed -2i(d-5)/Pi_3 is wrong; code correct); adjacent domain at delta = 7, 8 is p >= 4; at even delta the closing identity is blind to the seed a_delta, so rho(8) rests on the proved initial-data derivation; k = 6 is an exact finite kernel check (E3 does not cover the k = 6 kernel feed); Task D: the finite telescoping identity is correct, the broader equivalence has a repairable gap, x-indexing off by one; verify_tow4.py should check return codes.
