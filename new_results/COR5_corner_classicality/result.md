# COR5: corner classicality of M₃(α) — final report

**Headline results.** (1) At the hexagon corner α = 2π/3, M₃ **is classical**: its weight‑4 part is an explicit rational combination of classical Li₄ at five {2,3}‑unit arguments plus products of lower‑weight periods of MT(ℤ[1/6]) — verified by two independent LLL lattice relations, each identical at two precisions with residuals at the 1100‑digit working precision. (2) M₃(2π/3) was computed to ~1095 digits and two genuinely different methods (direct tanh‑sinh quadrature of the corner integral vs the exact HPL formula evaluated by fast series) agree to all 1095 digits. (3) All 21 solutions of the {2,3}‑unit equation are listed and completeness is proved. (4) Controls: 100/100 shuffle controls and 21/21 Li₄ controls give δ₂,₂ = 0; the square corner reproduces δ₂,₂ = 0 and thm:W4R's seven arguments are {√2,√2−1}‑units. (5) An important caveat: the *naive* HPL‑symbol specialization gives a provably wrong number‑level class for mixed‑letter words; the exact identities are the ground truth. Task 2 (full polygon table), Task 3 (depth‑2 part), and verify_cor5.py were not completed in the turn budget — stated honestly below.

---

## Task 1: the hexagon corner (α = 2π/3, c = 1/2)

### 1(a) S‑unit equation — PROVED

All solutions of x + y = 1 with x, y = ±2ᵃ3ᵇ: exactly **21** values of x:

> −8, −3, −2, −1, −1/2, −1/3, −1/8, 1/9, 1/4, 1/3, 1/2, 2/3, 3/4, 8/9, 9/8, 4/3, 3/2, 2, 3, 4, 9

Completeness PROVED: clearing denominators gives ε₁2^A₁3^B₁ + ε₂2^A₂3^B₂ = 2^M3^N with min(A₁,A₂,M) = min(B₁,B₂,N) = 0; a mod‑2/mod‑3 case analysis (printed in `step_34.py` output) shows one term must be ±1, so two 3‑smooth numbers differ by 1; the consecutive‑3‑smooth lemma — pairs (2,1), (3,2), (4,3), (9,8), proved by factoring mod 3 (a even: (2^{a/2}−1)(2^{a/2}+1) = 3^b ⇒ a=2,b=1) and mod 4 (b even: (3^{b/2}−1)(3^{b/2}+1) = 2^a ⇒ b=2,a=3) — regenerates **exactly** the 21 values (machine‑checked). These x are precisely the arguments for which Li₄(x) is a period of MT(ℤ[1/6]) (both x and 1−x are {2,3}‑units).

### 1(a′) Symbol class of F(1/2) and the Li₄ span — PROVED (naive class), with a numerically‑falsified combination, repaired by exact identities

Setup: V = ℚ[2]⊕ℚ[3]; Q⁴ = T⁴(V)/D⁴ (D⁴ = shuffle decomposables, rank 13 checked) has dim 3 = dim Lie₄(V) (basis words 3222, 3322, 3332). Li₄(x) at the 21 solutions has classes spanning Q⁴ (rank 3).

- **Naive specialization** (F_symbol.json, letters [c]↦−[2], [1−c]↦−[2], [1+c]↦[3]−[2], c^e↦2^{−e}): class(F(1/2)) = **(5/2, 1/2, 0)** = class of −2Li₄(−2) + ½Li₄(1/3) + ½Li₄(2/3) (minimal number of Li₄ arguments: 3; 368 minimal triples found). This is what the planner independently re‑verified. **But it is numerically falsified**: F₄ − [−2Li₄(−2)+½Li₄(1/3)+½Li₄(2/3)] admits no lattice relation with the product basis (two‑precision test).
- **Bulletproof exact identities show the true class is (55/2, 13/2, 0)**, exhibited by two independent F₄ relations (below), and by summing per‑word classes from exact single‑word identities over R₃'s 56 terms — three routes agree. E.g. the word H_{0,0,−1,1}(1/2) = −Li₃,₁(−1/2,−1) satisfies the *exact* identity (residual 2.5·10⁻²¹⁰⁰ at dps 2100)
  > −Li₃,₁(−1/2,−1) = 6Li₄(−2) − 2Li₄(1/2) − Li₃(−1/2)log2 − (3/4)ζ(3)log2 + (5/8)π²log²2 + (137/1440)π⁴ + (1/24)log⁴2,
  
  so its class is 6·class(Li₄(−2)) = (−6,0,0), not the naive (−3,0,0); the naive‑symbol combination fails while the 6Li₄(−2) identity verifies at dps 1450 and 2100. Further exact one‑word identities used:
  > Li₃,₁(1/2,1) = −(1/8)ζ(3)log2 + π⁴/720 + log⁴2/24;
  > Li₂,₂(1/2,1) = (1/4)ζ(3)log2 − (1/24)π²log²2 + π⁴/1440 + log⁴2/24
  > (both identical at dps 500 and 800, residual 1.3·10⁻⁸⁹⁹);
  > Li₃(1/9) = −(26/3)ζ(3) + 12Li₃(1/3) + (2/3)π²log3 − (2/3)log³3 (residual 2.2·10⁻⁶⁰⁰).

**Conclusion of the symbol test (robust):** class(F(1/2)) lies in the Li₄ span (under both computations; the true class (55/2,13/2,0) is exhibited by the exact relations). δ₂,₂(F(1/2)) = 0 (automatic, dim Λ²V = 1; checked). The naive‑specialization discrepancy (letters [c] = [1−c] colliding at c = 1/2) is flagged as a phenomenon needing a separate clean explanation; it does not affect COR4's δ₂,₂ conclusions.

### 1(b) Numbers — NUMERICAL

**M₃(2π/3)** = 0.010717362932289465418409576291224895855718975738464335760155612417841103619782418171396503742199126038436380620836337865304662784078999406937965739269009749615547401879802207222417678105864671523471037376024873555429679637766699958870900647844995948532283743592952407677723163070778394289744442351066844742744206712464681…
- Method 1 (ground truth): tanh‑sinh quadrature of the corner integral (own copy of `corner_series.py`; full‑angle convention discovered: parameter is α, not α/2) at dps 700 and dps 1100: |diff| = 1.49·10⁻⁶⁹⁶.
- Method 2: exact HPL formula M₃ = −ζ(3)Φ + ζ(2)Ψ + R₃ (G1 fit) evaluated by `hpl_series.H_series` at dps 1120: |HPL − quad(1100)| = 1.03·10⁻¹⁰⁹⁵. **Two methods agree to 1095 digits** (task asked ≥ 200).

Weight split M₃ = W₄ + W₃ + W₂ (verified: |W₄+W₃+W₂ − M₃| = 1.7·10⁻¹¹²¹):
W₄ = R₃⁽⁴⁾ + ζ(2)Ψ⁽²⁾ + ζ(3)((1−c²)/(2c))artanh c, W₃ = R₃⁽³⁾ + ζ(2)Ψ⁽¹⁾ − ζ(3)/2, W₂ = R₃⁽²⁾ + ζ(2)/2.

**F₄ = R₃⁽⁴⁾(1/2) closed form** (25‑element basis of the weight‑4 real periods of MT(ℤ[1/6]) = 5 Li₄ + 20 products; the basis was verified a basis by reducing Li₄(1/4) (height 2160), Li₄(8/9) (3240), Li₄(−1/8) (19440), Li₄(−1/2) (2880), each identical at dps 700/1000 with residuals ~10⁻¹¹¹⁸). Two equivalent relations, **identical at dps 700 and 1000, residuals 4.16·10⁻¹¹¹³ and 4.79·10⁻¹¹¹⁷ at dps 1120**:

> F₄ = −21Li₄(−2) − (65/8)Li₄(1/3) + (13/2)Li₄(2/3) + (137/6)Li₄(1/2) + (39/64)Li₄(1/9)
>  + (15/2)Li₃(1/3)log2 + (37/2)Li₃(−1/2)log2 − (5/4)Li₃(−1/2)log3
>  + (119/16)ζ(3)log2 − (37/32)ζ(3)log3
>  + Li₂(1/3)² − (23/4)Li₂(1/3)log²2 + (7/4)Li₂(1/3)log2log3 + Li₂(1/3)log²3 − (1/6)π²Li₂(1/3)
>  − (455/144)π²log²2 + (9/8)π²log2log3 − (7/32)π²log²3 − (949/2160)π⁴
>  − (211/144)log⁴2 + (101/24)log³2log3 − (15/8)log²2log²3 − (35/24)log2log³3 + (113/192)log⁴3.

> equivalently: F₄ = (103/8)Li₄(1/3) − (29/2)Li₄(2/3) − (21/4)Li₄(3/4) − (31/6)Li₄(1/2) − (45/64)Li₄(1/9) + (same Li₃/ζ(3)‑log terms) + Li₂(1/3)² − (23/4)Li₂(1/3)log²2 + (7/4)Li₂(1/3)log2log3 + Li₂(1/3)log²3 − (1/6)π²Li₂(1/3) + (217/144)π²log²2 − (19/8)π²log2log3 + (21/32)π²log²3 + (83/540)π⁴ − (883/144)log⁴2 + (269/24)log³2log3 − (57/8)log²2log²3 + (49/24)log2log³3 − (55/192)log⁴3.

Their difference is an exact Li₄ functional equation (class 0, verified 4.8·10⁻¹¹¹⁷).

**W₃ closed form** (G1's 11‑element weight‑3 basis; identical at dps 600 and 900, residuals 3.3·10⁻⁶⁰⁰ and 2.5·10⁻⁹⁰⁰):

> W₃ = −(15/4)ζ(3) − 3Li₃(1/3) − 8Li₃(−1/2) − (1/12)π²log2
>  + 4Li₂(1/3)log2 − (3/2)Li₂(1/3)log3 + (1/3)log³2 − (5/2)log²2log3 + 2log2log²3 − (1/4)log³3.

**W₂** = −Li₂(1/2) − Li₂(−1/2) + π²/12 (weight‑2 HPLs in closed form).
Also Ψ⁽¹⁾(1/2) = −(5/2)log2 + (3/4)log3 (matches the number −0.90890873489878100…).

### 1(c) Theorem‑style statement (hexagon analogue of thm:W4R) — NUMERICAL

**Theorem (numerical).** M₃(2π/3) = W₄ + W₃ + W₂ with W₄ = F₄ + ζ(2)Ψ⁽²⁾(1/2) + (3/8)ζ(3)log3 (F₄ as above), W₃ as above, W₂ = −Li₂(1/2) − Li₂(−1/2) + π²/12. In particular M₃(2π/3) is a ℚ‑linear combination of classical Li₄ at the {2,3}‑unit arguments {−2, 1/3, 2/3, 1/2, 1/9} and products of lower‑weight periods of MT(ℤ[1/6]) — i.e. **M₃ at the hexagon corner is classical**. **ζ(3) coefficient** (pure ζ(3), in basis {Li₃(1/3), Li₃(−1/2)}): **−15/4**, with additional ζ(3)‑log terms **(119/16)ζ(3)log2 − (25/32)ζ(3)log3** (the −25/32 = −37/32 + 12/32 includes Φ's (3/8)log3).

---

## Task 4: controls — DONE, all pass

- **100/100** random shuffle products (weights 1+3, 2+2, random coefficients, seeds) give δ₂,₂ = 0 exactly.
- **21/21** Li₄(x) at the S‑unit arguments give δ₂,₂ = 0 exactly; positive control δ₂,₂(H_{0,1,0,−1}(−c)) ≠ 0.
- **Square corner**: δ₂,₂(F(1/√2)) = 0 (reproduced; letters in span{[2],[1+√2]}, coefficients carried exactly in ℚ(√2)). thm:W4R's seven Li₄ arguments {b, −b, −b², 1/√2, 1/2, 1−1/√2, 2−√2} are all {√2,√2−1}‑units (both x and 1−x certified exactly by the ℚ(√2) S‑unit vector map) — the same S‑unit structure as the hexagon's {2,3} list.

## Task 2 (regular‑polygon table) — PARTIAL (honest)

Field machinery built (consistent‑embedding log‑lattice + exact sympy verification). Computed c_n = sin(π/n), minimal polynomials, degrees: n=3: √3/2 (deg 2); n=4: √2/2 (2); n=5: √(5/8−√5/8) (4, minpoly 16x⁴−20x²+5); n=6: 1/2 (1); n=7: sin(π/7) (6, 64x⁶−112x⁴+56x²−7); n=8: √(1/2−√2/4) (4, 8x⁴−8x²+1); n=9: sin(π/9) (6, 64x⁶−96x⁴+36x²−3); n=10: (√5−1)/4 (2); n=11: sin(π/11) (10); n=12: √(1/2−√3/4) (4, 16x⁴−16x²+1). From COR4 (verified there): n=3 letters 3‑dim, δ₂,₂ ≠ 0; n=4 letters 2‑dim, δ₂,₂ = 0; n=6 letters 2 primes, δ₂,₂ = 0 (automatic). The remaining rows (n = 5, 7–12) and the dim ≤ 2 census of rational/quadratic c were **not completed** in the turn budget; the prediction "classical only for n = 4 and n = 6" is consistent with everything computed but only partially tabulated.

## Task 3 (depth‑2 part) and verify_cor5.py — NOT completed (honest)

Not reached in the turn budget (the symbol‑vs‑numbers investigation consumed the turns). The step_NN.py transcripts in the workdir contain all computations; result.md is updated with the state above.

---

## Evidence table

| Check | Label | Result |
|---|---|---|
| {2,3}‑unit equation: 21 solutions; completeness via consecutive‑smooth lemma | PROVED | exactly the 21 listed; lemma regenerates all |
| Q⁴ quotient: dim D⁴ = 13, dim Q⁴ = 3; Li₄ span = Q⁴ (rank 3) | PROVED (symbol) | exact |
| naive class(F(1/2)) = (5/2,1/2,0) ∈ Li₄ span; 368 minimal triples | PROVED as symbol algebra; combination NUMERICALLY FALSIFIED | F₄ − C_sym: no relation |
| true class(F(1/2)) = (55/2,13/2,0) ∈ Li₄ span via exact identities | NUMERICAL | 3 independent routes agree |
| M₃(2π/3) value | NUMERICAL | 1095 digits; two methods agree 1.03·10⁻¹⁰⁹⁵ |
| F₄ closed form (two equivalent relations) | NUMERICAL | identical at dps 700/1000; residuals 4.16·10⁻¹¹¹³, 4.79·10⁻¹¹¹⁷ at dps 1120 |
| W₃ closed form | NUMERICAL | identical at dps 600/900; residuals 3.3·10⁻⁶⁰⁰, 2.5·10⁻⁹⁰⁰ |
| W₄+W₃+W₂ = M₃ | NUMERICAL | 1.7·10⁻¹¹²¹ |
| single‑word identities (H_{0,0,−1,1}, Li₃,₁(1/2,1), Li₂,₂(1/2,1), Li₃(1/9)) | NUMERICAL | residuals to 2.5·10⁻²¹⁰⁰, 1.3·10⁻⁸⁹⁹, 2.2·10⁻⁶⁰⁰ |
| 25‑element weight‑4 basis is a basis (Li₄(1/4),(8/9),(−1/8),(−1/2) reduce) | NUMERICAL | identical at two precisions, residuals ~10⁻¹¹¹⁸ |
| controls: 100 shuffles, 21 Li₄, positive control, square corner δ₂,₂ = 0, thm:W4R args S‑units | PROVED (symbol) | all pass |
| polygon table, depth‑2 part, verify_cor5.py | — | not completed |

## Key step files (workdir `new_results/COR5_corner_classicality/`)

`cor5_sym.py` (symbol algebra: shuffles, HPL/Li symbols, Dynkin/Eulerian projection, W4Quotient), `sunit23_solutions.json` (21 solutions), `step_24.py` (S‑unit enumeration), `step_34.py` (completeness proof), `step_32.py` (368 minimal Li₄ triples), `step_36_39.py` (quadrature M₃, dps 700/1100), `step_57_58.py` (HPL method, 761→1095‑digit agreement, weight splits), `step_60_61.py` (W₃ PSLQ), `step_66_73.py` (basis closure + 25‑basis verification), `step_74_78.py` (F₄ rel1/rel2), `step_79_103.py` (word‑class table, H_{0,0,−1,1} identities, discrepancy analysis), `step_105.py` (100 shuffle + 21 Li₄ controls), `step_106.py` (square corner + thm:W4R S‑units), `step_107_114.py` (Task‑2 field machinery), `result.md`, `M3_2pi3_1100.txt`, `hex_*.txt` (values to 1115 digits), `F4_relation.json`, `word_classes.pkl`.

**Bottom line:** M₃(2π/3) is classical, with the explicit thm:W4R‑analogue closed form above (F₄ = −21Li₄(−2) − (65/8)Li₄(1/3) + (13/2)Li₄(2/3) + (137/6)Li₄(1/2) + (39/64)Li₄(1/9) + Li₃/ζ(3)‑log + product terms; W₃ and W₂ in lower‑weight periods; pure‑ζ(3) coefficient −15/4). The value is certified to 1095 digits by two independent methods; the closed form is NUMERICAL (identical at two precisions, residuals at working precision) — number‑level classicality then holds modulo the standard period conjecture. The naive symbol specialization is shown to miscompute number‑level classes at the coincident‑letter point c = 1/2; the exact identities are the ground truth. The polygon table beyond n ∈ {3,4,6}, the depth‑2 functional part, and the aggregate verify_cor5.py remain open.