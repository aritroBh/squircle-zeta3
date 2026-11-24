work = '<project>/research/08_theory/TOW4_tower_depth8'
result = r"""# TOW4 — the tower law at v = i: **PROVED TO DEPTH 8** (all d, all k); the weight-six polynomial law PROVED for every k >= 6; all-depth attempt (Task D)

## VERDICT

**PROVED TO DEPTH 8.** The tower law "the depth-delta tower coefficient u_{pq}^{(k)}(delta) at v = +i is a
polynomial in k of degree floor(delta/2)" is proved for depths **delta <= 8**, all keys (p,q) of all degrees d,
all k (off-diagonal: k >= d; diagonal: k >= 2p+1, k = 2p exceptional), extending TOW3 (delta <= 6). The proof
is FIT-FREE: the e-family is the FORCED ascending recursion F(delta,0) := rho(delta)/e00 = a_delta,
F(delta,s+1) := Lambda^{(delta)}(s) - F(delta,s), built entirely from the PROVED initial data a_delta
(Bell assembly from the proved R_1..R_4) and the proved tau constants. The Newton identities N(delta,s) for
s < floor(delta/2) hold by construction; the single remaining identity at each depth is the closing condition
F(delta, floor(delta/2)+1) = 0 (TOW3 M1 pre-screen equivalence), here PROVED as an exact rational identity in
(d,m) for delta = 7 (F(7,4) = 0) and delta = 8 (F(8,5) = 0), together with the adjacent-key identities
(fold + ehat-conversion, exact identities in p) and blindness E3 at depths 7 and 8.

**Consequence (Task C): the weight-six polynomial law for every k >= 6.** Tower law at degree 6 to depth
8 = 2*6-4 (proved here) ==(WTS Lemma 3, PROVED)==> kernel law K(5): degree-5 kernel residues are polynomials
in k of degree <= 4 for k >= 7 (the degree-6 diagonal key (3,3) is stable for k >= 7; k = 6 is its exceptional
level and is outside Lemma 3's range) ==(WTS Lemma 2, k-independent fixed real-linear image)==> W_6(k) is a
quartic in k for k >= 7. The +-i simple-pole coefficients are residues, hence quartics for k >= 7; they vanish
at k = 7..11 (5 points > degree 4), hence identically for k >= 7. The level k = 6 is settled EXACTLY: the
recorded degree-5 kernels at k = 6 equal the quartic fitted from k = 7..11 at every pole, every QSI component,
all three keys (5,0),(4,1),(3,2) (**48/48 components, 0 bad**; quarticity on k = 7..13 independently confirmed,
48/48; the +-i poles vanish at k = 6..13). Hence
    W_6(k) = sum_{m=1..5} C(k-6,m-1) E_{5,m}  =  E_{5,1} + (k-6) E_{5,2} + C(k-6,2) E_{5,3} + C(k-6,3) E_{5,4} + C(k-6,4) E_{5,5}
    for every k >= 6
(k = 5 exceptional, as always — confirmed: kern_5 degree-5 residues differ from the quartic at 15/48
components, mirroring WT5's P(4) != W_5(4); constants E_{5,m} as determined by WT6S). This upgrades WT6S's
exact-data law to a proved law for all k >= 6 (paper prop:w6k). Weights polynomial in k for ALL k >= w:
w = 2 (TOW1), 3 (WT3+TOW1), 4 (TOW2+k=4 exact, TOW2R/CX8), 5 (TOW3+k=5 exact), **6 (here, k=6 exact)**;
w >= 7 open (needs delta >= 10).

## What is proved — every step in writing (all indices) or exact symbolic code (numerator/remainder 0)

| item | statement | status |
|---|---|---|
| tau constants n = 8 | tau1(8) = -169/128 + 41i/128, tau2(8) = conj; conj(tau1_n) = (-1)^n tau2_n | PROVED (exact vs tow2lib def, diff 0; 40-digit numeric) |
| L(1..8) cache = proved construction | L_simplified == tow2lib.L(r) symbolically, r = 1..8 (R_1..R_4 proved TOW2) | PROVED, diff 0 (r=8: 9.9s); closes TOW3's r<=6 gap |
| rho(7), rho(8) closed forms | a_7 = N_7/(40320 Pi_3), a_8 = N_8/(645120 Pi_4), Bell from proved L | PROVED + exact check **717 cells each, 0 bad** |
| closing conditions (non-adjacent) | F(7,4) = 0, F(8,5) = 0, exact rational identities in (d,m) | PROVED, numerator 0 (forced78.py, total 47.9s) |
| forced family = true e-family | F == tow3lib.e (delta<=6, all s); F(delta,top) == zigzag top (delta<=8) | PROVED (sanity; the M1 claim reproduced) |
| adjacent identities (7,0..3),(8,0..4) | fold + ehat-conversion, exact identities in p | PROVED, exact numerator 0 (prove_adj78_fast.py; 9/9) |
| blindness E3 at delta = 7, 8 | Delta_pp feed into child (p,p-1) vanishes identically | PROVED all p, numerator 0 (23.6s, 34.9s) |
| tower law delta <= 8 | above + TOW1 reduction theorem | PROVED; exact check **49953 cells, 0 bad** (+ 222 diagonal-initial) |
| weight-6 law (Task C) | K(5) for k >= 7 (proved) + k = 6 exact check 48/48 | PROVED for all k >= 6 |

## Key structures

- a_7, a_8: proof_depth8.md Section 1 (explicit numerators N_7, N_8 over Q(i)); denominators Pi_3, Pi_4.
- The FORCED family F(delta,s) for delta <= 8 is cached (forced_F78.pkl, cancel'd single fractions);
  denominators are exactly Pi_{floor(delta/2)} times small integers; tops are the zigzag family.
- Adjacent/e3 machinery byte-inherited from TOW3 (prove_adj56/prove_E3b56 structure), with e00_pp
  anti-real (Dc = -conj(Delta)) and e00_{p,p-1} real (CX8 #2 carried forward).

## Domains (stated explicitly, per TOW3R wording fix)

Depth delta needs 2d-4 >= delta (delta = 7, 8 => d >= 6); adjacent (p,p-1) needs 4p-6 >= delta (p >= 3);
k >= d off-diagonal; k >= 2p+1 diagonal stable, k = 2p exceptional. All feed/Pi denominators are nonzero
odd integers in these domains; the exact numerator-0 checks are rational identities over Q(i), valid
wherever denominators are nonzero.

## Files (research/08_theory/TOW4_tower_depth8/)

- `tow4lib.py` (tau to n = 8, arho delta <= 8, forced family loader), `forced78.py` (forced recursion +
  closing proofs), `prove_adj78_fast.py` (adjacent), `prove_E3b78.py` (blindness), `check_full78.py`
  (49953-cell confirmation), `check_data.py` (referee QSI arithmetic, copied), `forced_F78.pkl`,
  `a7_closed.pkl`, `a8_closed.pkl`, `L_simplified.pkl` (copied), `result.md`, `proof_depth8.md`,
  `progress.log`, `verify_tow4.py`.

## Task D — all-depth attempt

See Section D of proof_depth8.md: the closing condition F(delta, floor(delta/2)+1) = 0 for ALL delta,
via the forced recursion / TOW2 master equation. (Status recorded there; no overclaim.)
"""
open(work + '/result.md','w').write(result)
print("wrote result.md", len(result))
