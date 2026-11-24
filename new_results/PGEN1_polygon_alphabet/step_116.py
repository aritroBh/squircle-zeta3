result = r"""# PGEN1: the trilogarithm alphabet of a general polygon, and the cobracket of the far kink lines

Labels: **PROVED** = exact rational/symbol algebra (machine-checked identities); **NUMERICAL** = high-precision checks with digits stated; **OPEN** = not settled.

Data inputs (read-only): OCT9d `short_form.pkl`, `TW3.pkl` (generic octagon, symbol-proved); OCT9b `sol_B123.pkl`, `sol_B0.pkl`, `Sunit_lists.pkl` (OCT9c basis); OCT6 `result.md` 34-term table (regular octagon); DEC1 `result.md` 15-term table (decagon); DOD5 `dod5_tallies.pkl` (12-gon, symbol-proved).

## 0. Setup

On an edge e (t from corner to corner) the terms of the weighted sum are A0*t, A1*(1-t) (corner kink lines) and B_j(1+kappa_j t) (far kink lines); far line j crosses the line of e at t = p_j = -1/kappa_j. The weight-3 part of the edge is a combination of iterated integrals G(a1,a2,a3;1) with letters in {0,1,p_1,...,p_s}: polylogarithms on M_{0,s+3} with marked points {0,1,inf,p_1,...,p_s}.
delta(N3) is computed as the (2,1)-cobracket: for a presentation N3 = sum c_i Li_3(z_i) modulo products,
  delta(N3) = sum_i c_i (z_i ^ (1-z_i)) (x) z_i  in  ^2(F^x) (x) F^x (over Q).
Over Q and over the real quadratic fields Q(sqrt2), Q(sqrt3), Q(sqrt5) the pre-Bloch group p_2(F) (x) Q injects into ^2 F^x (Borel: B_2(F) (x) Q = 0 since r_2 = 0), so the wedge element decides delta = 0 faithfully, and any two symbol-exact presentations give the same wedge element (a Li_3 functional equation has delta in B_2 (x) F^x = 0). All wedge elements below are exact.

## 1. Task 1 -- the alphabet rule

**Rule (planner's conjecture):** each edge needs Li_3 only at cross-ratios of its points {0,1,inf,p_j}, plus products of lower weight and zeta(3).

### 1(a). thm:hexagon (s = 1). PROVED.
With p = -1/kappa = 1/(1-r), the cross-ratios of {0,1,inf,p} form a single S_3-orbit {1-r, r, 1/r, 1/(1-r), (r-1)/r, r/(r-1)}. The theorem's two arguments x = 1-r and x/(x-1) are both in this orbit (checked for r in {2, 1/2, 3, 1/3, 5/2, 7/4}, exact). At r = 2 and r = 1/2 the two points become -1 and 1/2, whose wedges vanish for ANY coefficients (delta = 0): this is exactly why prop:nonlocal's non-local zeta(3) is meaningful (canonical) at r = 2, 1/2.

### 1(b). The generic rational octagon (s = 2 per edge). PROVED (symbol level).
kappa values (oct9c_common.K12), marked points, and cross-ratio sets CR_e:

| edge | kappa_1, kappa_2 | p_1, p_2 | |CR_e| | letters of CR_e | rank span{s(z): z in CR_e} mod D | TW3_e in span? |
|---|---|---|---|---|---|---|---|
| 0 | 3/2, -7/12 | -2/3, 12/7 | 30 | {2,3,5,7} | 10 (of 20) | YES |
| 1 | 1/3, 1/9 | -3, -9 | 30 | {2,3,5} | 8 (of 8) | YES |
| 2 | 5/3, -5/9 | -3/5, 9/5 | 30 | {2,3,5} | 8 (of 8) | YES |
| 3 | 1/4, 1/24 | -4, -24 | 18 | {2,3,5} | 5 (of 8) | YES |

- "TW3_e in span" is an exact rational rank computation (PROVED): each edge's weight-3 symbol (OCT9d's TW3, itself symbol-proved) lies in span{S(Li_3(z)) : z in CR_e} + D, D = shuffle-decomposables. We extracted explicit cross-ratio-only presentations of all four edges (exact symbol identities): 25, 18, 18, 14 nonzero Li_3 terms respectively, saved in CR_presentations.pkl.
- The 8 shared OCT9d arguments and the 20 edge-0 arguments are NOT themselves in CR_e (they were chosen to minimize invisible constants); some follow by duplication: -9/16 = -(3/4)^2 with 3/4 in CR_1, 4/9 = (2/3)^2 with 2/3 in CR_1; 7-4sqrt3 = (2-sqrt3)^2 for the 12-gon.
- **The dimensions 8 and 20 explained.** For a letter space V of dimension n, the weight-3 symbols modulo products form Q_3 = V^{(x)3}/D with dim Q_3 = (n^3-n)/3 (Witt's formula: the degree-3 part of the free Lie algebra on n letters). Edges 1-3 have letters in {2,3,5} (n = 3): dim 8. Edge 0 has letters in {2,3,5,7} (n = 4): dim 20. The prime 7 enters through p_2 = 12/7 (the ratio 5/12). The CR_e alone span 10, 8, 8, 5 of these dimensions; the 99 {2,3,5}-units span the full 8, and the {2,3,5,7}-units span the full 20 (OCT9b/OCT9d). So "8 and 20" are the Witt dimensions of the edge's letter alphabet; the edge's own cross-ratios are enough for its target (proved), though they need not span the whole alphabet space.

### 1(c). Regular octagon, decagon, 12-gon. PROVED (exact algebra).
- **Regular octagon (s = 2):** far ratios r = sqrt2 and sqrt2/2; p_1 = 2+sqrt2, p_2 = -(1+sqrt2). CR_e = 15 points: orbit(sqrt2) = {sqrt2, 1-sqrt2, sqrt2/2, 1-sqrt2/2, 2+sqrt2, -1-sqrt2}, orbit(-sqrt2/2) = {-sqrt2/2, 1+sqrt2/2, -sqrt2, sqrt2-1, 2-sqrt2, 1+sqrt2}, orbit(1/2) = {1/2, 2}. All six Li_3 points of the closed form (1-sqrt2, 1-sqrt2/2, sqrt2-1, 2-sqrt2, sqrt2/2, -sqrt2/2) lie in CR_e; the rational orbit {1/2, 2} is where log 2 comes from. Letters {2, 1+sqrt2} (n = 2): dim Q_3 = 2.
- **Decagon (s = 2 effective; the third non-local term is constant):** p_1 = -phi, p_2 = phi^2. CR_e = the single 6-point orbit {+-phi, +-1/phi, phi^2, 1/phi^2}, which contains both Li_3 points used (+-1/phi). Letters collapse to span{phi}: every z and 1-z is +-phi^k.
- **12-gon (s = 4):** p_1 = (3+sqrt3)/2, p_2 = 4+2sqrt3, p_3 = -(3+2sqrt3), p_4 = -(1+sqrt3)/2. |CR_e| = 63, all over the 3-letter space {1+sqrt3, sqrt3, 2+sqrt3} with (1+sqrt3)^2 = 2(2+sqrt3) (n = 3, dim Q_3 = 8). All 8 DOD5 Li_3 points are in CR_e or duplication: 7 directly in CR_e (with full S_3-orbits), and 7-4sqrt3 = (2-sqrt3)^2 is the square of a cross-ratio. The span of {S(Li_3(z)) : z in CR_e} mod products is the FULL 8-dimensional Q_3, and all 8 DOD5 Li_3 symbols lie in it (exact ranks, PROVED).

### 1(d). Symbol-level proof for s <= 4 far kink lines. PROVED (exact computation).
Symbol recursion (verified numerically at weights 2 and 3): with the chain (a_0 = 0, a_1,...,a_n, a_{n+1} = 1),
  S(I(0;a_1,...,a_n;1)) = sum_i S(I(0;...a_i removed...;1)) (x) (a_{i+1}-a_i)/(a_i-a_{i-1}).
For s in {1,2,3,4} take formal letters {p_i, q_i = 1-p_i, d_{ij} = p_i-p_j} (n = 2,5,9,14 letters) and all convergent distinct-letter words G(a1,a2,a3;1) (a_1 != 0, a_3 != 1):

| s | n | dim Q_3 = (n^3-n)/3 | #CR points | rank span s(CR) mod D | #words | rank span(words) mod D | words in span s(CR)+D? |
|---|---|---|---|---|---|---|---|
| 1 | 2 | 2 | 6 | 2 | 3 | 2 | YES |
| 2 | 5 | 40 | 30 | 10 | 14 | 7 | YES |
| 3 | 9 | 240 | 90 | 30 | 39 | 17 | YES |
| 4 | 14 | 910 | 210 | 70 | 84 | 34 | YES |

Ranks are exact integer ranks (flint fmpz_mat, PROVED). Separately, every weight-3 word (any repeats/divergences) reduces via shuffle relations to the convergent distinct-letter words + products (exact, PROVED for s <= 4). So: **every weight-3 iterated-integral symbol with letters in {0,1,p_1,...,p_s} is, modulo products, a combination of Li_3-symbols at cross-ratios of {0,1,inf,p_1,...,p_s}, for s <= 4.** Combined with Goncharov's theorem that every weight-3 multiple polylogarithm reduces to Li_3's (and Rudenko's theorem that weight-3 MPLs are Li_3's at single cross-ratios + products), this proves the rule at the symbol level for every polygon with at most 4 far kink lines per edge - which covers all known examples (hexagon s = 1, octagons and decagon s = 2, 12-gon s = 4).

**Theorem candidate:** P'''(1) of any centrally symmetric convex polygon is a combination of Li_3 at the cross-ratios of the points where the kink lines cut the line of each edge, zeta(3), and products of lower weight. (Proved at the symbol level for s <= 4 far kink lines per edge; numerical/symbol-exact in all known closed forms.)

## 2. Task 2 -- the cobracket of the far kink lines

Framework: N3 = P'''(1) minus the local terms of prop:edgelaw, modulo products. delta(N3) = sum c_i (z_i ^ (1-z_i)) (x) z_i; computed from symbol-proved presentations. Basis-independence verified computationally: the generic-octagon delta from the OCT9d 44-term basis and from the independent cross-ratio presentations are identical (exact, PROVED); the 12-gon delta is invariant under the duplication presentation change z_4 = z_1^2 (exact, PROVED); every delta lands in ker(^3) = Lie^3(V) (exact check).

**Controls (all PROVED/pass):** delta{1/2}_3 = 0 (consistent with Li_3(1/2) = 7/8 zeta(3) + products); delta of the Landen triple {x}+{1-x}+{1-1/x} is 0 for 7 sample rationals; the hexagon's two arguments are cross-ratios; delta vanishes at r = 2 and r = 1/2 for any coefficients.

### 2(a). Generic rational octagon, edge by edge and total. PROVED; delta != 0 for every edge.
delta_e from the symbol-proved OCT9d coefficients (identical to the cross-ratio presentation):

- edge 1 (kappa 1/3, 1/9): 9 terms over letters {2,3,5}:
  (2^3): 5976/3125*2 - 6516/3125*3 + 2448/3125*5 ; (2^5): -6192/15625*2 - 5328/15625*3 + 6012/15625*5 ; (3^5): -17568/15625*2 + 5688/3125*3 - 2844/3125*5
- edge 2 (kappa 5/3, -5/9): 8 terms:
  (2^3): 78264/10985*2 - 11448/2197*3 + 144/10985*5 ; (2^5): -62136/10985*2 + 35928/10985*3 + 3168/10985*5 ; (3^5): 35784/10985*2 - 5256/2197*3
- edge 3 (kappa 1/4, 1/24): 6 terms:
  (2^5): (262400*2 + 75200*3 - 158200*5)/24389 ; (3^5): (75200*2 + 19200*3 - 43880*5)/24389
- edge 0 (kappa 3/2, -7/12): 24 terms over {2,3,5,7} (denominators 1681 = 41^2 and 68921 = 41^3), e.g. (2^3): (1800*2 - 5160*3 + 8820*5 - 1360*7)/68921, (2^7): (-26160*2 - 4180*3 - 4180*5 + 7000*7)/68921, ...  (full list in delta_oct9d.pkl and verify_pgen1.py output)

The total P'''(1) has the four edge contributions with Q-linearly independent scalar factors (sqrt41, 1, sqrt13, sqrt29), so the total delta vanishes iff each edge's delta vanishes: **the generic octagon's non-local zeta(3) is basis-dependent on every edge** (matching the paper's remark).

**Geometric formula.** For each edge, delta(N3)|_e = sum_{z in CR_e} c_z (z ^ (1-z)) (x) z with the explicit cross-ratio presentation coefficients (CR_presentations.pkl): a sum of {cross-ratio}_2 (x) cross-ratio over the edge's own marked-point cross-ratios. The wedge element is always in ker(^3), i.e. in Lie^3 of the edge's letter space (Witt dimensions 8, 8, 5, 10 for the four edges' CR spans).

### 2(b). Regular octagon. PROVED; delta != 0.
From the 6 Li_3 coefficients of the proved 34-term form (letters r2 = sqrt2, u = 1+sqrt2; wedge = r2^u only):
  delta(N3) = (r2^u)(x)r2 * (-5/2 + 15 sqrt2/8) + (r2^u)(x)u * (-11/4 + 2 sqrt2).
Invariant under the lattice of relations (Landen triples {sqrt2-1, 2-sqrt2, -sqrt2/2} and {sqrt2/2, 1-sqrt2/2, 1-sqrt2} have delta = 0; the duplication pair {sqrt2/2}+{-sqrt2/2} has delta = 0). The non-local zeta(3) (39 sqrt2 - 42)/4 * u is genuinely basis-dependent.

### 2(c). Decagon. PROVED; delta = 0.
Letters collapse to span{phi}: every z, 1-z is +-phi^k, so every wedge z^(1-z) = 0. The non-local part is in ker delta, hence canonical. With the golden reduction (NUMERICAL, two precisions: residual 1.5e-151 at 150 digits, 8.7e-252 at 250 digits)
  chi_3(1/phi) = (Li_3(1/phi)-Li_3(-1/phi))/2 = (25 sqrt5/48) L(chi_5,3) - (pi^2/12) log phi + (1/4) log^3 phi
(the golden analog of the paper's eq:chi3 for sqrt2; eq:chi3 itself re-verified to 4.9e-150 with correctly computed L(chi_8,3) = 0.958380454563094562...), and the sum Li_3(1/phi)+Li_3(-1/phi) = (zeta(3) - log(phi) pi^2/6)/5 + log^3(phi)/6 (DEC1, so alpha = 1/10 exact), the decagon's non-local part is
  **N3 = (1/64 + 31 sqrt5/320) * zeta(3) + (875/768 - 125 sqrt5/192) * L(chi_5,3)**
(NUMERICAL, verified against the direct single-valued evaluation of DEC1's presentation: difference 8e-41 at 120 digits).

### 2(d). 12-gon. PROVED; delta != 0.
From DOD5's 8 symbol-proved Li_3 coefficients (letters g2 = 1+sqrt3, g3 = sqrt3, u = 2+sqrt3):
  delta(N3) = (g2^u)(x)[(5/4+5sqrt3/4) g2 + (-1/4-sqrt3/4) g3 + (1/4-3sqrt3/4) u]
            + (g2^g3)(x)[(5-21sqrt3/4) g2 + (-37/8+29sqrt3/8) g3 + (-2+9sqrt3/4) u]
            + (g3^u)(x)[(7/4-13sqrt3/4) g2 + (-9/4+17sqrt3/8) g3 + (-17/16+25sqrt3/16) u]
i.e. 9 + 9 terms (rational and sqrt3 parts). Invariant under the duplication presentation z_4 = z_1^2 (exact). The paper's L(chi_12,3) coefficient 459/16 - 93sqrt3/16 is basis-dependent (delta != 0).

### 2(e). Hypothesis test: "the non-local part of the regular 2m-gon is delta-closed after summing over edges".
- m = 2, 3: trivial (fully local).
- m = 4 (octagon): **FAILS** (delta != 0, §2(b)).
- m = 5 (decagon): **HOLDS**, trivially: the golden letters collapse to span{phi}.
- m = 6 (12-gon): **FAILS** (delta != 0, §2(d)).
So the hypothesis is FALSE in general; the decagon is a special case where the whole non-local part is canonical.

## 3. Task 3 -- single-valued L_3 (basis-independence test)

Convention: for real x != 0,1, L_3(x) = Re Li_3(x) - log|x| Re Li_2(x) + (1/3) log^2|x| * Li_1(x) with Li_1(x) = -log|1-x| (Zagier's single-valued trilogarithm; equals Li_3(x) - log|x| Li_2(x) - (1/3) log^2|x| log(1-x) for x < 1). Verified: L_3(1/2) = 7/8 zeta(3) exactly; Landen L_3(x)+L_3(1-x)+L_3(1-1/x) = zeta(3) at x = 0.3, 0.7, -0.4, 2, -2, 2.5 (residuals < 1e-80); duplication L_3(x^2) = 4L_3(x)+4L_3(-x) at 0.6, -0.8 (residual 0.0); inversion L_3(x) = L_3(1/x) at x = 2,3,...,36, -2, -3/2, -7/5 (residuals < 1e-79); cross-checked against a series implementation.

For xi = a*zeta(3) + sum b_i Li_3(z_i) of the non-local part of each edge of the generic octagon, computed in the OCT9c and OCT9d bases:
- **edges 1, 2, 3: IDENTICAL values** (differences 5.7e-100, 2.6e-100, 1.7e-99 at 100 digits): L_3(xi_1) = 0.02367462581917732081212099120162979836462, L_3(xi_2) = 0.4874364362702959625368298252094181146908, L_3(xi_3) = -0.02102185222421880354518585875934411585906.
- **edge 0: the two bases differ by exactly (80/206763) zeta(3)** (pslq at 200 digits; 206763 = 3*41^3): L_3(xi_0)^OCT9c = 0.04876528948758063, L_3(xi_0)^OCT9d = 0.04830019393251146.
- All three edge-0 presentations (OCT9c 36-term, OCT9d 20-term, and our 25-term cross-ratio presentation) have the SAME total symbol (verified exactly, 0 difference) and pairwise L_3-differences that are rational zeta(3)-multiples (pslq at 200 digits): Sc-Sd = -(31119817701814913/654492976038000) zeta(3), S_CR-Sd = -(6228044/5169075) zeta(3). So they represent the same B_3-class; the (80/206763) zeta(3) shift between the two sessions' invisible-constant splittings is a rational zeta(3)-trade between the Li_3 basis and the invisible constant - exactly the basis-dependence the paper describes ("the rest of the zeta(3) depends on the basis"), and consistent with delta(N3) != 0 for edge 0: with a nontrivial cobracket, the zeta(3)-share of a presentation is not canonical.
**Explanation of the invariant:** L_3 descends to the Bloch group B_3 (it satisfies all trilogarithm functional equations); the B_3-class of an edge's non-local part is basis-independent (verified: identical symbols, rational-zeta(3) L_3-differences); the invariant real number is the single-valued evaluation of the class, defined up to rational zeta(3)-trades from different invisible-constant splittings; for edges 1-3 the two splittings agree exactly.

## 4. Evidence table

| check | label | residual / exactness |
|---|---|---|
| controls: delta{1/2}=0; Landen delta=0 (7 pts); hexagon CR; delta=0 at r=2,1/2 | PROVED | exact |
| generic octagon: TW3_e in span s(CR_e)+D (4 edges) | PROVED | exact ranks |
| CR-only presentations of 4 edges (25,18,18,14 terms) | PROVED | TW3 - sum in D, exact |
| delta per edge: OCT9d basis == CR basis | PROVED | exact, identical |
| regular octagon delta; decagon delta = 0; 12-gon delta | PROVED | exact |
| delta invariant under duplication (12-gon) | PROVED | exact |
| delta in ker(^3) (octagon, 12-gon) | PROVED | exact |
| alphabet proof s <= 4 (words in span s(CR)+D; shuffle reduction) | PROVED | exact flint ranks |
| 12-gon: 8 DOD5 points in CR-span (full 8-dim); 7/8 in CR, z4 = z1^2 | PROVED | exact |
| decagon pair (a,b) vs val(N3) | NUMERICAL | 8e-41 at 120 dps |
| golden chi_3 identity (25sqrt5/48) | NUMERICAL | 1.5e-151 / 8.7e-252 (150/250 dps) |
| eq:chi3 control | NUMERICAL | 4.9e-150 |
| D8/L(chi8,3) = 4sqrt2/3 | NUMERICAL | exact at 150, 250 dps |
| L_3(xi_e) OCT9c == OCT9d, edges 1-3 | NUMERICAL | 5.7e-100 ... 1.7e-99 |
| edge-0 L_3 trade = (80/206763) zeta(3); pairwise q rational | NUMERICAL | pslq at 200 dps |

## 5. Files
- `verify_pgen1.py` (< 5 min): controls, cross-ratio checks, delta computations (all polygons, per-edge and total), alphabet-proof ranks for s <= 4, task-3 numbers.
- `CR_presentations.pkl` (cross-ratio presentations of the 4 generic-octagon edges), `delta_oct9d.pkl`, `CR_regular.pkl`, `field_tables.pkl`, `task3_L3.pkl`, `qfvec.pkl`.
- `symbolproof_engine.py` (formal-letter proof engine), step_*.py (session log).
"""
open("result.md","w").write(result)
print("result.md written:", len(result), "chars")
