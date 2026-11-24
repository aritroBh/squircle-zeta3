# TOW5 — the tower law at v = i: **PROVED TO DEPTH 10** (all d, all k); the weight-SEVEN polynomial law PROVED for every k ≥ 7

## VERDICT

**PROVED TO DEPTH 10.** The tower law "the depth-δ tower coefficient u_{pq}^{(k)}(δ) at v = +i is a polynomial in k of degree ⌊δ/2⌋" is proved for depths **δ ≤ 10**, all keys (p,q) of all degrees d, all k (off-diagonal: k ≥ d; diagonal: k ≥ 2p+1 with k = 2p exceptional), extending TOW4 (δ ≤ 8). The proof is **FIT-FREE**: the e-family is the FORCED ascending recursion F(δ,0) := ρ(δ)/e00 = a_δ, F(δ,s+1) := Λ^{(δ)}(s) − F(δ,s), built entirely from the PROVED initial data a_δ (Bell assembly from the proved R_1..R_5 — **R_5 is proved here**, closing the last gap) and the proved τ constants (extended to n = 10). The Newton identities N(δ,s) for s < ⌊δ/2⌋ hold by construction; the single remaining identity per depth is the closing condition F(δ, ⌊δ/2⌋+1) = 0 (TOW3 M1 equivalence), PROVED as an exact rational identity in (d,m): **F(9,5) = 0 and F(10,6) = 0** — together with the adjacent-key identities (fold + ê, exact identities in p) and blindness E3 at depths 9, 10. Referee errata (TOW3R/TOW4R/CX8/CX10) applied throughout (normalisation by ρ(0)=e00; off-diagonal e vs stable diagonal ê distinct; finite descending-induction uniqueness; adjacent domain δ=9,10 is p ≥ 4; even-δ closing is blind to the seed, so ρ(10) rests on the proved initial-data derivation; verify checks return codes; zigzag = sec + tan).

**Consequence (Task D): the weight-seven polynomial law for every k ≥ 7.** Tower law at degree 7 to depth 10 = 2·7−4 (proved) ⇒(WTS Lemma 3)⇒ kernel law K(6): degree-6 kernel residues are polynomials in k of degree ≤ 5. **The range is k ≥ 7, not merely k ≥ 8: degree 7 is odd, so there is NO diagonal key of degree 7** (all four keys (7,0),(6,1),(5,2),(4,3) are off-diagonal; and j = i−1 never holds for i+j = 6, so no degree-6 work item has a diagonal fold feed at degree 7 either). The extra level k = d+1 = 7 that Lemma 3 leaves unproved for even n is therefore INSIDE the proved range — this answers the Task-D question (nothing is exceptional at k = 7). ⇒(WTS Lemma 2, fixed real-linear image)⇒ W_7(k) is a quintic in k for k ≥ 7. The ±i simple-pole coefficients are residues, hence quintics for k ≥ 7, vanish at k = 8..13 (6 points > degree 5), hence identically for k ≥ 7. The level k = 7 is settled EXACTLY: recorded degree-6 kernels at k = 7 equal the quintic fitted from k = 8..13 at every pole, every QSI component, all four keys (**56/56, 0 bad**). Hence
  **W_7(k) = Σ_{m=1..6} C(k−7, m−1) E_{6,m} = E_{6,1} + (k−7)E_{6,2} + C(k−7,2)E_{6,3} + C(k−7,3)E_{6,4} + C(k−7,4)E_{6,5} + C(k−7,5)E_{6,6} for every k ≥ 7**
(constants E_{6,1}..E_{6,6} determined by the weight-7 kernel data). **k = 6 is OFF the quintic**: kern_6 degree-6 residues differ at **18/56 components** (the initial kernel, as at every lower weight). Weights polynomial in k for ALL k ≥ w: w = 2,3,4,5 (previous), 6 (TOW4), **7 (here)**; w ≥ 8 open (needs δ ≥ 12).

**Task E: the all-depth closing condition remains REDUCED to the master x = −1 identity (TOW4's proved equivalence), now verified symbolically at δ = 9 (7.4s) and δ = 10 (52.7s); all-depth proof NOT achieved (same obstruction; no overclaim).**

## Task A — R_5, ρ(9), ρ(10), PROVED

**R_5 = S⁽⁵⁾(t*)/S(t*) = N_5/(3 Π_5), PROVED for all (d,m)** (`prove_S5.py`): TOW2's prove_S34 method extended to l = 5 (126-term finite sum, proved method). The new Stirling tail P_5(x) = S(x,x−5) = x(x−1)(x−2)(x−3)(x−4)²(x−5)²(3x²−23x+38)/11520 is PROVED by the same recurrence (P_r(x+1) = (x+1−r)P_{r−1}(x) + P_r(x), exact polynomial identities, r ≤ 5) + anchors (P_r(a)=0, a=0..r), hence P_r(n) = S(n,n−r) ∀n ≥ 0 (induction via the Stirling recurrence). The identity total(5) == N_5/(3 Π_5) is proved symbolically (sp.cancel of the difference is exactly 0; the gcd-reduced denominator is exactly 3 Π_5 — note sp.together alone under-cancels; N_5 has degree 10 = 2l in d, 35 terms, leading −3072 i d^10). CONFIRMED by direct exact differentiation of S(t): **24/24 integer keys d = 6..10 and 67/67 keys at d ∈ {12,15,18,22,26,30}**.

a_9 = Bell_9(L_1..L_9), a_10 = Bell_10(L_1..L_10) (a_9 via proved R_1..R_4; a_10 via sl_5 from the proved R_5 by exact cumulant algebra). L cache verified EQUAL to the proved tow2lib-style L(r) **symbolically** (√2-aware, difference 0) for **r = 1..10** (`verify_L.py`). τ extended to n = 10 (land in Q(i), √2 cancels; exact vs defs, difference 0, +40-digit numeric ~5e-41; conj(τ1_n) = (−1)ⁿτ2_n holds): τ1(9) = −99/256 − 297i/256, τ2(9) = 99/256 − 297i/256; τ1(10) = 1 − 99i/256, τ2(10) = 1 + 99i/256.

  a_9  = N_9  / (5806080 (2d−9)(2d−7)(2d−5)(2d−3))           [N_9: deg 13, 40 terms, lc −256 i d^13]
  a_10 = N_10 / (116121600 (2d−11)(2d−9)(2d−7)(2d−5)(2d−3))  [N_10: deg 15, 51 terms, lc −1024 d^15]

Exact check vs step-B data: ρ(δ)/ρ(0) = a_δ (ρ(0) = e00; initial level k0 = d; √2-part asserted 0): **425 cells each for δ = 9, 10 (all keys 2d−4 ≥ δ, d ≤ 40), 0 bad** (the ρ(10) check is a 425-fold independent confirmation of R_5).

## Task B — exact tower data at depths 9, 10

TOWN1's trunc_engine copied and run at D = 10 (exact Q(√2,i); validated vs TOWN1 at δ ≤ 8: 162/162 ρ cells, 36/36 tower cells): (i) ρ initial data, degrees 7..40 (`data_rho_D10_d7_40.pkl`); (ii) full towers, levels k = 7..24, degrees 7..12 (`data_towers_D10_K24.pkl`, 1001 cells); (iii) extended towers, levels k = 7..30, degrees 7..18 (`data_towers_D10_K30_deg18.pkl`, 16291 cells).

## Task C — Newton identities at depths 9, 10, PROVED

**Non-adjacent (d,m symbolic).** FORCED family from ρ alone (fit-free); the closing conditions: **F(9,5) = 0 and F(10,6) = 0, sp.together + numerator exactly 0** (`forced910.py`; family δ ≤ 10 built in 54s). Sanity: F == tow3lib.e (δ ≤ 6 all s); F(top) == zigzag (δ ≤ 10): F(9,4) = −16i(d−6)/Π_4 (A_5 = 16), F(10,5) = −61/Π_5 (A_6 = 61); other F(9,·), F(10,·) in proof_depth10.md / forced_F910.pkl.

**Adjacent (p symbolic).** The 11 identities (9,0..4) and (10,0..5) hold with the diagonal ê-conversion and the fold, exact numerator 0 (`prove_adj910_fast.py`; 0.1s–92s each). The adjacent TOP identities (9,4), (10,5) are NEW.

**Blindness E3.** The feed of Δ_pp(δ′) (δ′ ≤ δ) into the child (p,p−1) at depth δ (DL2-direct + fold) vanishes identically at δ = 9, 10, exact numerator 0 (`prove_E3b910.py`, 13.9s/22.9s).

**Full-law confirmation:** **990 cells + 11 diagonal-initial, 0 bad** (degree 7..12, k = 7..24; `check_full910.py`); **extended 16225 cells + 66 diagonal-initial, 0 bad** (degree 7..18, k = 7..30; `check_full910_ext.py`).

## Task D — the weight-seven law for every k ≥ 7, PROVED

(see the VERDICT). Exact certificates (`verify_W7.py`, WTS kern_k.pkl ground-truth exact kernels, degree-6 keys (6,0),(5,1),(4,2),(3,3), real poles {0, a = 1+√2, −a, −b = 1−√2}): W7a k = 7 == quintic(8..13) **56/56**; W7b 6th differences vanish k = 7..13 (56/56); W7c ±i simple-pole coefficients vanish k = 7..13; W7d k = 6 off the quintic **18/56** (initial kernel). The paper's range may be changed to "W_w(k) polynomial for every k ≥ w" for w = 2..7.

## Task E — toward all depths (no overclaim)

TOW4's PROVED reduction (closing at δ ⟺ the finite alternating identity a_δ = feed(U(−1)), U(dp) = Σ_s (−1)^s F(dp,s) = Û(dp,−1) ⟺ the master equation (1+x)Û = ρ + xΛ(Û) at x = −1 ⟺ the value identity V(δ)) is unchanged; the master identity feed(U(−1)) − a_δ == 0 is now verified symbolically at **δ = 9 (7.4s) and δ = 10 (52.7s)** (verify_tow5.py MASTER). Structural evidence continues: at δ = 9, 10, U(δ) has exactly the same denominator as a_δ and coincides in the top-3 leading numerator terms (first differing at d^10 resp. d^12, 3 below the top, as at δ = 7, 8).

*New structural observation (verified symbolically, δ = 2..10, `taskE_defect.py`):* the defect U(δ) − a_δ has a universal prefactor — δ even: (−1)^{δ/2+1}·2δ(δ−1)·Q_δ; δ odd: (−1)^{(δ−1)/2}·2δ(δ−1)·i(d−(δ+3)/2)·Q_δ, with Q_δ a polynomial (exact division, remainder 0, all δ = 2..10; deg_d Q_δ = 3(⌊δ/2⌋−1)). So U(δ) = a_δ + prefactor·Q_δ has a definite closed form per δ; identifying Q_δ in general would give an independent grip on U(δ) = Û(δ,−1). Not achieved: the identity remains self-referential at general δ (the dp = δ term re-introduces U(δ) at the parents), so the all-depth proof is NOT achieved — same obstruction as TOW4, now confirmed two depths higher.

## Evidence table (all checks EXACT; residual = numerator/remainder EXACTLY 0)

| lemma | statement | certificate | residual | time |
|---|---|---|---|---|
| TAU910 | τ1/τ2(9),(10) in Q(i), exact vs defs; conj symmetry | verify_tow5.py (+tau910.pkl) | diff 0 (+40-digit ~5e-41) | 0.2s |
| R5 | R_5 = N_5/(3 Π_5), all (d,m); P_5 proved | prove_S5.py | cancel'd diff 0; direct-diff 24/24 + 67/67 | 10s |
| L | L cache == tow2lib-style L(r), r = 1..10, symbolic | verify_L.py | √2-aware diff 0 | 161s |
| RHO | a_9,a_10 Bell; ρ(9)/e00==a_9, ρ(10)/e00==a_10 | verify_tow5.py | 425 cells each, 0 bad | 7s/54s |
| CLOSE | closing F(9,5)=0, F(10,6)=0 in (d,m) | forced910.py | numerator 0 | 54s total |
| CLOSE-sanity | F==e (δ≤6 all s); F(top)==zigzag (δ≤10) | verify_tow5.py | diff 0 | 6s |
| ADJ | 11 adjacent identities (9,0..4),(10,0..5) in p | prove_adj910_fast.py | numerator 0 | 0.1–92s each |
| E3 | blindness δ = 9, 10 in p | prove_E3b910.py | numerator 0 | 13.9s/22.9s |
| LAW | tower law δ≤10 vs step-B data | check_full910.py / check_full910_ext.py | 990 + 11, 0 bad; extended 16225 + 66, 0 bad | 16s/18s |
| W7 | k=7 vs quintic(8..13); quinticity 7..13; ±i vanish; k=6 off | verify_W7.py | 56/56 comp; k=6: 18/56 off | 0.2s |
| MASTER | closing ⟺ master x=−1 coefficient, δ = 9, 10 | verify_tow5.py | feed(U)−a_δ == 0 | 7.4s/52.7s |

Digit counts: every certificate is exact rational arithmetic over Q(i) (sympy) or exact Fraction Q(√2,i) arithmetic (data checks); residuals exactly 0, no floating point except the independent 40-digit mpmath confirmation of τ(9), τ(10), additionally proved symbolically (difference 0). Data cells: ρ 2×425; full law 990+11 (core) and 16225+66 (extended); weight-7 4 keys × 56 components × 4 checks. Domains: depth δ needs 2d−4 ≥ δ (δ = 9, 10 ⇒ d ≥ 7); adjacent (p,p−1) needs 4p−6 ≥ δ (δ = 9, 10 ⇒ p ≥ 4); k ≥ d off-diagonal, k ≥ 2p+1 diagonal stable, k = 2p exceptional; all feed/Π denominators are nonzero odd integers in these domains, so numerator-0 rational identities prove the identities throughout. **`verify_tow5.py`: 11/11 lemmas PASS in 623s < 40 min** (log `verify_tow5.log`).

## Files (new_results/TOW5_tower_depth10/)

- `tow5lib.py` (R_5 via N5_proved; sl_5; τ to n = 10; L cache to r = 10; arho δ ≤ 10; forced loader), `prove_S5.py` (P_5 + R_5 = N_5/(3 Π_5), PROVED + direct check), `verify_L.py`, `forced910.py` (closing proofs → `forced_F910.pkl`), `prove_adj910_fast.py`, `prove_E3b910.py`, `check_full910.py`, `check_full910_ext.py`, `verify_W7.py`, `verify_tow5.py`, `taskE_defect.py`, `N5_proved.pkl`, `tau910.pkl`, `L910_candidate.pkl`, `a9_closed.pkl`, `a10_closed.pkl`, `forced_F910.pkl`, `data_rho_D10_d7_40.pkl`, `data_towers_D10_K24.pkl`, `data_towers_D10_K30_deg18.pkl`, `result.md`, `proof_depth10.md`, `progress.log`.
- Byte-inherited (copied, provenance preserved): `tow2lib.py`, `tow3lib.py`, `tow4lib.py`, `L_simplified.pkl`, `L_parsed.pkl`, `check_data.py`, `arho_simplified.pkl`, `prove_S34.py` (TOW2/4); `trunc_engine.py`, `driver.py` (TOWN1); `fastred.py` (WTS); `forced_F78.pkl` (TOW4). Data read-only: WTS `theory/kern_*.pkl`.

**Honest boundaries:** tower law at depths ≥ 11 unproved (ρ(11) needs only R_1..R_5 = proved, reachable; ρ(12) needs R_6, not proved); the all-depth closing condition remains REDUCED to the master x = −1 identity but NOT proved for all δ (verified at δ = 9, 10); the weight-eight law (needs δ ≤ 12) is open.