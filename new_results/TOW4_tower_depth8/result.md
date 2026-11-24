# TOW4 — the tower law at v = i: **PROVED TO DEPTH 8** (all d, all k); the weight-six polynomial law PROVED for every k ≥ 6; Task D: a PROVED reduction (no overclaim)

## VERDICT

**PROVED TO DEPTH 8.** The tower law "the depth-δ tower coefficient u_{pq}^{(k)}(δ) at v = +i is a polynomial in k of degree ⌊δ/2⌋" is proved for depths **δ ≤ 8**, all keys (p,q) of all degrees d, all k (off-diagonal: k ≥ d; diagonal: k ≥ 2p+1 with k = 2p exceptional), extending TOW3 (δ ≤ 6). The proof is **FIT-FREE**: the e-family is the FORCED ascending recursion F(δ,0) := ρ(δ)/e00 = a_δ, F(δ,s+1) := Λ^{(δ)}(s) − F(δ,s), built entirely from the PROVED initial data a_δ (Bell assembly from the proved R_1..R_4) and the proved τ constants (extended here to n = 8). The Newton identities N(δ,s) for s < ⌊δ/2⌋ hold by construction; the single remaining identity per depth is the closing condition F(δ, ⌊δ/2⌋+1) = 0 (TOW3 M1 pre-screen equivalence), PROVED as an exact rational identity in (d,m): **F(7,4) = 0 and F(8,5) = 0** — together with the adjacent-key identities (fold + ê-conversion, exact identities in p) and blindness E3 at depths 7 and 8.

**Consequence (Task C): the weight-six polynomial law for every k ≥ 6.** Tower law at degree 6 to depth 8 = 2·6−4 (proved) ⇒(WTS Lemma 3)⇒ kernel law K(5): degree-5 kernel residues are polynomials in k of degree ≤ 4 for k ≥ 7 (degree-6 diagonal key (3,3) stable for k ≥ 7; k = 6 is its exceptional level, outside Lemma 3) ⇒(WTS Lemma 2, fixed real-linear image)⇒ W_6(k) is a quartic in k for k ≥ 7. The ±i simple-pole coefficients are residues, hence quartics for k ≥ 7, and vanish at k = 7..11 (5 points > degree 4), hence identically for k ≥ 7. The level k = 6 is settled EXACTLY: recorded degree-5 kernels at k = 6 equal the quartic fitted from k = 7..11 at every pole, every QSI component, all three keys (5,0),(4,1),(3,2) (**48/48 components, 0 bad**). Hence
  **W_6(k) = Σ_{m=1..5} C(k−6,m−1) E_{5,m} = E_{5,1} + (k−6)E_{5,2} + C(k−6,2)E_{5,3} + C(k−6,3)E_{5,4} + C(k−6,4)E_{5,5} for every k ≥ 6**
(k = 5 exceptional, confirmed: kern_5 residues differ from the quartic at 15/48 components, mirroring WT5's P(4) ≠ W_5(4); constants E_{5,m} as determined by WT6S). Paper prop:w6k. Weights polynomial in k for ALL k ≥ w: w = 2,3,4,5 (previous jobs), **6 (here)**; w ≥ 7 open (needs δ ≥ 10).

**Task D: a PROVED reduction, all-depth proof NOT achieved (obstruction identified).** See below; no overclaim.

## Task A — ρ(7), ρ(8) closed forms, PROVED

a_7 = Bell_7(L_1..L_7), a_8 = Bell_8(L_1..L_8) (use sl_l, l ≤ 4 = proved R_1..R_4). The simplified cache L(1..8) is verified EQUAL to the proved tow2lib.L(r) construction **symbolically** (exact difference 0) for **r = 1..8** — closing the r = 7,8 gap left at TOW3 (TOW3R did r ≤ 6). τ extended to n = 8: τ1(8) = −169/128 + 41i/128, τ2(8) = conjugate (√2 cancels; verified exactly vs the tow2lib definitions, difference 0, plus 40-digit numeric; conj(τ1_n) = (−1)ⁿτ2_n holds).

  a_7 = N_7/(40320 (2d−7)(2d−5)(2d−3))   [N_7: deg 10 in d, 41 terms, leading 64i d¹⁰ − 2048i d⁹ + (672+29808i)d⁸ + …]
  a_8 = N_8/(645120 (2d−9)(2d−7)(2d−5)(2d−3))   [N_8: deg 12 in d, 57 terms, leading 256 d¹² − 10752 d¹¹ + …]
  (Π-parts Π_3, Π_4, matching a_δ ~ Π_{⌊δ/2⌋}; full numerators in proof_depth8.md §1 and a7_closed.pkl/a8_closed.pkl.)

Exact check vs TOWN1 towers_D8_K52.pkl: ρ(δ)/ρ(0) = a_δ, ρ(0) = e00 (the 1+[p>q] factor is inside e00; initial level k0 = d for ALL keys — for diagonal keys the exceptional level k = 2p, whose ρ_pp is the true initial data; √2-part of ratio asserted 0): **717 cells each for δ = 7, 8 (all keys 2d−4 ≥ δ, d ≤ 52), 0 bad.**

## Task B — Newton identities at depths 7, 8, PROVED

**Non-adjacent (d,m symbolic).** The FORCED family F is defined from ρ alone (fit-free); N(δ,s) for s < ⌊δ/2⌋ hold by construction; the closing conditions are the only further identities: **F(7,4) = 0 and F(8,5) = 0, sp.together + numerator exactly 0** (forced78.py; full family δ ≤ 8 built in 47.9s). Sanity (M1 claim reproduced): F == tow3lib.e (δ ≤ 6, all s, difference 0); F(δ,⌊δ/2⌋) == zigzag top (δ ≤ 8, difference 0). New closed forms (denominators exactly Π_{⌊δ/2⌋} × small integers): F(7,1) = N_{71}/(480 Π_3), F(7,2) = N_{72}/(6 Π_3), F(7,3) = −2i(d−5)/Π_3 (zigzag), F(8,1) = N_{81}/(5760 Π_4), F(8,2) = N_{82}/(48 Π_4), F(8,3) = N_{83}/(4 Π_4), F(8,4) = 16/Π_4 (zigzag, A_5 = 16). Full numerators in proof_depth8.md §2.

**Adjacent (p symbolic, fold + ê-conversion).** For the adjacent key (p, p−1), all (δ,s) ∈ {(7,0..3),(8,0..4)}: e(δ,s)+e(δ,s+1) = Λ^{(δ)}(s) as exact rational identities in p over Q(i) — **9/9, exact numerator 0** via an lcm-denominator polynomial check with cached cancelled ê/Ediag fractions (prove_adj78_fast.py; method validated against TOW3's slow sp.simplify method on 4 cases, identical True). Adjacent TOP identities (7,3) [odd s=3] and (8,4) [even s=4] are NEW. Timings: (7,0) 56s, (7,1) 20s, (7,2) 6s, (7,3) 0.8s, (8,0) 97s, (8,1) 40s, (8,2) 15s, (8,3) 3s, (8,4) 0.3s. Sanity: reproduces all δ ≤ 6 adjacent identities. (This crushes the depth-6 bottleneck: (6,0) took 668s at TOW3; here (8,0) = 97s.)

**Blindness E3 (all p).** The feed of Δ_pp(δ′) = ρ_pp(δ′) − Σ_s(−1)^s ê_pp(δ′,s) (δ′ ≤ δ) into the child (p,p−1) at depth δ (DL2-direct + fold) vanishes identically: **numerator 0, δ = 7 (23.6s), δ = 8 (34.9s)** (prove_E3b78.py; validated on δ = 5,6). e00_pp anti-real ⇒ Dc = −conj(Δ); e00_{p,p−1} real (CX8 #2 carried forward).

**Conclusion.** By TOW1's reduction theorem (E0 + all N(δ′,s), δ′ ≤ 8, both key types + E3), the tower law holds at depths δ ≤ 8 for all keys and all k; the descending finite-induction uniqueness argument (TOW2R/CX8) gives U = u. Exact confirmation vs TOWN1: **49953 cells, 0 bad** (all keys d ≤ 52, δ ≤ min(8,2d−4), k ≤ min(d+8,52); off-diagonal via F, diagonal via ê) + **222 diagonal-initial cells k = 2p vs a_δ(2p,p), 0 bad** (check_full78.py).

## Task D — toward all depths: a PROVED reduction (closing ⟺ master x = −1 ⟺ value identities)

**Theorem (reduction, PROVED in writing, all indices, and verified symbolically at δ = 7, 8).** For each δ, given the zigzag top equations (PROVED all s, TOW2) and closing at all δ′ < δ, the following are equivalent, and the tower law at ALL depths is equivalent to any of them for all δ:
1. the closing condition F(δ, ⌊δ/2⌋+1) = 0;
2. the **finite alternating identity** a_δ = i(2d−1) Σ_{dp=0}^{δ} [ τ1(δ−dp)·U(dp)(d+1,m+1) − τ2(δ−dp)·U(dp)(d+1,m) ] / (2d−dp−1), where U(dp) := Σ_{s=0}^{⌊dp/2⌋} (−1)^s F(dp,s) = Û(dp, x=−1) (a FINITE terminating sum at each dp ≤ δ, given closing at lower depths);
3. the **master equation at x = −1**, depth-δ coefficient: ρ = Λ(Û(−1)) (TOW2's (1+x)Û = ρ + xΛ(Û) specialized at x = −1, with the x-indexing as corrected by CX8 #7: x⁰ = value identities, x^{s≥1} = slopes, x^{top+1} = tops);
4. the value identity V(δ) (TOW2's descending formulation).
Proof: (1)⟺(2) is the telescoping of the recursion, F(δ,top+1) = Σ_{s=0}^{top}(−1)^{top−s}Λ(δ,s) + (−1)^{top+1}a_δ = (−1)^{top}[Σ_s(−1)^sΛ(δ,s) − a_δ], and Σ_s(−1)^sΛ(δ,s) = feed of Σ_s(−1)^sF(dp,s) = feed of U. (2)=(3) is the definition of the x = −1 specialization. (1)⟺(4) is the ascending/descending duality. Certificate: **feed(U(−1)) − a_δ ≡ 0 symbolically at δ = 7 (4.8s), δ = 8 (23.1s)** (verify_tow4.py lemma MASTER).

**Structural evidence.** U(δ) = Û(δ,−1) has exactly the same denominator as a_δ (480Π₂, 5760Π₃, 40320Π₃, 645120Π₄ at δ = 5,6,7,8) and coincides with a_δ in all leading numerator terms (e.g. at δ = 7 both begin 64i d¹⁰ − 2048i d⁹ + (672+29808i)d⁸ − 1344 d⁷m + …, diverging only at d⁷). This makes the reduction a sharp, testable statement per δ.

**Obstruction (honest).** The all-depth proof reduces to the single scalar identity (2) for all δ — equivalently to proving that the formal x = −1 specialization of the master equation TERMINATES (CX8's "a generating identity alone does not prove the required Newton termination"). No closed form for U(δ) was found, and the identity resists a direct induction on δ: the dp = δ term re-introduces U(δ) at the parents, so (2) at depth δ is self-referential unless an independent grip on U(δ) = Û(δ,−1) exists. I did NOT prove the closing condition for all δ; what is proved is the equivalence above and its validity at δ = 7, 8.

## Evidence table (all checks EXACT; residual = numerator/remainder EXACTLY 0)

| lemma | statement | certificate | residual | time |
|---|---|---|---|---|
| TAU8 | τ1(8),τ2(8) ∈ Q(i), exact vs defs; conj symmetry | verify_tow4.py | diff 0 (+40-digit numeric) | <1s |
| L | L cache == tow2lib.L(r), r = 1..8, symbolic | verify_tow4.py | diff 0 | r=8: 9.9s |
| RHO | a_7,a_8 Bell; ρ(7)/e00==a_7, ρ(8)/e00==a_8 | verify_tow4.py | 717 cells each, 0 bad | 3s/8s |
| CLOSE | closing F(7,4)=0, F(8,5)=0 in (d,m) | forced78.py | numerator 0 | 47.9s total |
| CLOSE-sanity | F==e (δ≤6 all s); F(top)==zigzag (δ≤8) | verify_tow4.py | diff 0 | — |
| ADJ | 9 adjacent identities (7,0..3),(8,0..4) in p | prove_adj78_fast.py | numerator 0 | 0.3–97s |
| E3 | blindness δ = 7, 8 in p | prove_E3b78.py | numerator 0 | 23.6s/34.9s |
| LAW | tower law δ≤8 vs TOWN1 | check_full78.py | 49953 cells 0 bad + 222 diag-init 0 bad | 27s |
| W6 | k=6 vs quartic(7..11); quarticity 7..13; ±i vanish 6..13 | verify_tow4.py | 48/48 comp 0 bad | 10s |
| MASTER | closing ⟺ master x=−1 coefficient, δ = 7,8 | verify_tow4.py | feed(U)−a_δ ≡ 0 | 4.8s/23.1s |

Digit counts: every certificate is exact rational arithmetic over Q(i) (sympy) or exact Fraction Q(√2,i) arithmetic (data checks); residuals are exactly 0, no floating point anywhere except the independent 40-digit mpmath confirmation of τ(8), which is additionally proved symbolically (difference 0). Data cells: ρ 2×717; full law 49953 + 222; weight-6 3×48. Domains: depth δ needs 2d−4 ≥ δ (δ = 7,8 ⇒ d ≥ 6); adjacent needs 4p−6 ≥ δ (p ≥ 3); k ≥ d off-diagonal, k ≥ 2p+1 diagonal stable, k = 2p exceptional; all feed/Π denominators are nonzero odd integers in these domains, so numerator-0 rational identities prove the identities throughout.

## Files (new_results/TOW4_tower_depth8/)

- `tow4lib.py` (τ to n = 8; arho δ ≤ 8; forced-family loader), `forced78.py` (forced recursion + closing proofs, saves `forced_F78.pkl`), `prove_adj78_fast.py` (adjacent, fast exact method), `prove_E3b78.py` (blindness), `check_full78.py` (49953-cell law), `check_data.py` (referee QSI arithmetic, copied), `verify_tow4.py` (10/10 lemmas PASS, 473s < 30 min), `a7_closed.pkl`, `a8_closed.pkl`, `L_simplified.pkl` (copied from TOW3), `result.md`, `proof_depth8.md`, `progress.log`, `verify_tow4.log`.
- Byte-inherited from TOW3 (copied, provenance preserved): `tow2lib.py`, `tow3lib.py`, `L_parsed.pkl`, `prove_V56/adj56/E3b56/tops.py`. Data: TOWN1 `towers_D8_K52.pkl`, WTS `theory/kern_*.pkl` (read-only).
- Step scripts that produced the results (this session's step_NN.py in the same folder): step_27 (a_7,a_8 Bell assembly), step_31/32 (ρ 717-cell checks), step_34 (L cache r=1..8 symbolic), step_39/41 (τ(8)), step_43–49 (forced recursion forced78.py + closing), step_52 (F==e, F(top)==zigzag sanity), step_56–63 (adjacent fast method + 9/9 batch), step_64 (E3 δ=7,8), step_65–67 (check_full78 49953 cells), step_69–73 (Task C weight-6 + k=5 exceptional), step_78 (master x=−1 at 7,8), step_80–82 (verify_tow4.py 10/10).

**Honest boundaries:** tower law at depths ≥ 9 unproved (ρ(9) reachable with R_1..R_4, but ρ(10) needs R_5, not proved); the all-depth closing condition is REDUCED to the master x=−1 identity (2) but NOT proved for all δ; the weight-seven law (needs δ ≤ 10) is open.