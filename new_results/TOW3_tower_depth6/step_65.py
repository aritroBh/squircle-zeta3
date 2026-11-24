work = '<project>/research/08_theory/TOW3_tower_depth6'
proof = r"""# TOW3 — proof document: the tower law at depths 5 and 6 (PROVED), weight-5 law (PROVED), Task D

Inherited from TOW2 (all PROVED there): setup, Lemma T, Lemmas S1/S4/S3/S5
(R_1..R_4 = S^{(l)}(t*)/S(t*) = N_l/Pi_l for all (d,m); `prove_S12.py`, `prove_S34.py`),
the Bell theorem rho(delta)/e00 = a_delta = [y^delta] exp(sum_{r>=1} L_r y^r) with L_r built
(Faa di Bruno, only blocks of size >= 2 since t'(i) = 0) from sl_l, l <= floor(r/2), plus
explicit algebraic pieces; the zigzag top coefficients and the all-s proof of the top equations
for non-adjacent keys (`prove_tops.py`); the reduction (TOW1): E0 + {N(delta,s)} + E3 ==>
tower law; the tau-constants tau1_n, tau2_n in Q(i), n <= 7 (`check_tau`, re-run here: True).

## 1. rho(5), rho(6) — PROVED closed forms (Task A)

Ingredients: a_5 = Bell_5(L_1..L_5) uses sl_l, l <= 2 (R_1, R_2 PROVED);
a_6 = Bell_6(L_1..L_6) uses sl_l, l <= 3 (R_1, R_2, R_3 PROVED). The simplified cache
L(1..8) (`L_simplified.pkl`, TOW2) was verified EQUAL to the proved L(r) construction of
tow2lib (exact random-point check at r = 5,6; the L(r) formula is itself the proved assembly).

  a_5 = ( -16 i d^7 + 304 i d^6 + (-80-2520 i) d^5 + 160 d^4 m + (1040+11840 i) d^4
          - 2080 d^3 m + (-5140-34069 i) d^3 + (10280-240 i) d^2 m + (11440+60081 i) d^2
          + 240 i d m^2 + (-22880+960 i) d m + (-9600-59920 i) d - 960 i m^2 + 19200 m + 25200 i )
        / (480 (2d-5)(2d-3))
  a_6 = ( -64 d^9 + 1728 d^8 + (-20976+480 i) d^7 - 960 i d^6 m + (150432-10080 i) d^6
          + 20160 i d^5 m + (-702876+89040 i) d^5 + (-2880-178080 i) d^4 m + (2221452-423840 i) d^4
          + 2880 d^3 m^2 + (34560+847680 i) d^3 m + (-4760069+1146870 i) d^3
          - 34560 d^2 m^2 + (-138960-2293260 i) d^2 m + (6677868-1670730 i) d^2
          + (138960-1440 i) d m^2 + (186480+3341460 i) d m + (-5524995+1019460 i) d
          + 960 i m^3 - 186480 m^2 - 2038920 i m + 1984500 )
        / (5760 (2d-7)(2d-5)(2d-3))

Exact check vs TOWN1 towers_D8_K52.pkl: rho(delta)/e00 = a_delta (the factor 1+[p>q] is
already inside e00 = rho(0)): **720 cells each for delta = 5, 6, all keys with 2d-4 >= delta,
d <= 52; 0 bad.** Sanity: TOW2-proved a_4: 723/723. (QSI division; sqrt2-part of ratio = 0
asserted in every cell.)

## 2. The Newton identities at depths 5, 6 — PROVED (Task B, non-adjacent)

The full e-family at depths <= 6 (e(delta,0) = a_delta; zigzag tops; sub-tops):

  e(2,1) = -1/Pi_1
  e(3,1) = i (d-3)/Pi_1
  e(4,1) = (d^3 - 9 d^2 + (113/4) d - i d/2 + i m - 29)/Pi_2                       [TOW2, PROVED]
  e(4,2) = 2/Pi_2  (zigzag top)                                                   [TOW2, PROVED]
  e(5,1) = ( -4 i d^4 + 52 i d^3 + (-6-269 i) d^2 + 12 d m + (24+638 i) d - 48 m - 552 i )/(12 Pi_2)
  e(5,2) = -2 i (d-4)/Pi_2  (zigzag top, A_3 = 2)                                  [top, PROVED all s]
  e(6,1) = ( -16 d^6 + 336 d^5 + (-3064+48 i) d^4 - 96 i d^3 m + (15424-576 i) d^3
             + 1152 i d^2 m + (-44629+2364 i) d^2 + (-48-4728 i) d m + (68999-3228 i) d
             + 48 m^2 + 6456 i m - 43398 )/(96 Pi_3)
  e(6,2) = ( -4 d^3 + 48 d^2 + (-199+2 i) d - 4 i m + 268 )/(2 Pi_3)
  e(6,3) = -5/Pi_3  (zigzag top, A_4 = 5)                                          [top, PROVED all s]

Provenance of the three new sub-top forms: e(5,1), e(6,2) were FITTED in TOW2 (`e_fits.pkl`;
CANDIDATES there); e(6,1) is FITTED here from the TOWN1 towers (60 keys, exact linear algebra,
deg-6 numerator over Pi_3). Each is VERIFIED against the exact towers at **414 fresh keys
(d <= 40), 0 bad**, and — decisively — each is PROVED by the symbolic verification of EVERY
Newton identity it enters (below; this is the required upgrade from fit to proved).

*Theorem (non-adjacent identities, all (d,m)).* For every key with p >= q+2 and
(delta,s) in { (5,0),(5,1),(5,2),(6,0),(6,1),(6,2),(6,3) }:
    e(delta,s) + e(delta,s+1) = Lambda^{(delta)}(s).
Proof: exact rational identities in (d,m) over Q(i); `prove_V56.py`, numerator 0:
  (5,0) [V5] 0.8s, (5,1) 0.2s, (5,2) 0.0s, (6,0) [V6] 7.3s, (6,1) 0.6s, (6,2) 0.1s, (6,3) 0.0s.
The two value identities V(5) = (5,0), V(6) = (6,0) are exactly TOW2's stated gap.

## 3. Adjacent keys and blindness — PROVED (Task B)

*Theorem (adjacent, all p).* For the keys (p,p-1) and the seven (delta,s) above, the same
identities hold with the diagonal ehat-conversion and the fold: exact rational identities in p
over Q(i); `prove_adj56.py`, remainder 0: (5,0) 228s, (5,1) 30s, (5,2) 1.6s, (6,0) 668s,
(6,1) 128s, (6,2) 17s, (6,3) 0.8s. Sanity: reproduces TOW2's (3,0),(3,1),(4,0),(4,1),(4,2).
The adjacent TOP identities (5,2) [odd, s=2] and (6,3) [even, s=3] are NEW (TOW2 proved
adjacent tops only at (3,1),(4,2)).

*Theorem (blindness E3 at delta = 5, 6, all p).* The feed of
Delta_pp(delta') = rho_pp(delta') - sum_s (-1)^s ehat_pp(delta',s) (delta' <= delta) into the
child (p,p-1) at depth delta (DL2-direct + fold) vanishes identically: `prove_E3b56.py`,
numerator 0: delta = 5 [9.4s], delta = 6 [20.9s]. Sanity: reproduces delta = 3, 4.

## 4. Conclusion and Task C — the weight-five law for every k, PROVED

By TOW1's reduction theorem: E0 (PROVED, WT2b) + the identities N(delta,s) for delta <= 6
(Sections 2-3, all key types; the diagonal law follows from the off-diagonal parent via the
ehat-conversion) + blindness E3(delta) for delta <= 6 (Section 3) ==>
**the tower law holds at depths delta <= 6 for all keys and all k** (off-diagonal k >= d;
diagonal k >= 2p+1 with k = 2p exceptional).
Exact confirmation vs TOWN1 towers: **38939 cells, 0 bad** (all keys d <= 52, k <= min(d+8,52),
delta <= 6; off-diagonal via e-family, diagonal via ehat; initial levels by rho).

*Weight five (Task C).* Specializing to degree 5: the tower law to depth 6 = 2*5-4 implies,
by WTS Lemma 3 (PROVED there), the kernel law K(4): every degree-4 kernel residue is a
polynomial in k of degree <= 3 for k >= 6. W_5(k) is a k-independent linear image of these
residues (WTS Lemma 2, PROVED), hence a cubic in k for k >= 6. The level k = 5 = d+1 is not
covered by Lemma 3; it is settled EXACTLY here: the recorded degree-4 kernels at k = 5
(`kern_5.pkl`) equal the cubic fitted from k = 6..9, for all three keys (4,0),(3,1),(2,2),
every pole {0, a, -a, -b}, every QSI component: **40/40, 0 bad**. (Independent confirmation:
4th differences of every residue-component vanish on k = 6..13, 40/40.) Hence
    W_5(k) = W5'' + (k-5) E1 + C(k-5,2) E2 + C(k-5,3) E3    for every k >= 5
(k = 4 exceptional, as always; WT5 measured P(4) - W_5(4) != 0). The constants E1, E2, E3,
W5'' are WT5's (E3 classical proved-form candidate there; E1, E2 depth-two exact + classical
numerical). This upgrades WT5's exact-data statement (5 <= k <= 14) to a proved law: the
polynomiality itself is now proved for every k >= 5.

## 5. Task D — sub-top tower coefficients and the M_d map (in progress)
"""
with open(f'{work}/proof_depth6.md','w') as f: f.write(proof)
print('proof_depth6.md updated')
