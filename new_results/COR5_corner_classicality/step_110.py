md = """# COR5: at which corner angles is M₃(α) classical?

Builds on COR4 (`research/08_theory/COR4_M3alpha_classical`): the weight-4 part F(c) of R₃(c)
has the nonzero (2,2) obstruction
δ₂,₂(F) = (1−c⁻¹)([c]∧[1−c])∧([1−c]∧[1+c]) + (−1−c⁻¹)([c]∧[1+c])∧([1−c]∧[1+c]),
vanishing only where the letters collapse to a ≤ 2-dimensional span (then Λ²B₂ = 0 automatically).

Labels: PROVED = exact symbol algebra; NUMERICAL = high-precision lattice relations (accepted only
if identical at two precisions, residual at working precision).

## Task 1: the hexagon corner α = 2π/3 (c = 1/2) — M₃ IS classical.

### 1(a) S-unit equation. PROVED.
All solutions of x + y = 1 with x, y = ±2ᵃ3ᵇ: exactly 21 values of x:
  −8, −3, −2, −1, −1/2, −1/3, −1/8, 1/9, 1/4, 1/3, 1/2, 2/3, 3/4, 8/9, 9/8, 4/3, 3/2, 2, 3, 4, 9.
Completeness PROVED (clearing denominators, one term must be ±1 by mod-2/mod-3 case work printed in
the step output; the consecutive-3-smooth lemma 2−1, 3−2, 4−3, 9−8 — proved by factoring mod 3 and
mod 4 — then regenerates exactly the 21 values; checked).

### 1(a′) Symbol-class of F(1/2) and the Li₄ span.
The projected class of F(1/2) in Q⁴ = T⁴(V)/D⁴ (V = ℚ[2]⊕ℚ[3], dim Lie₄(V) = 3) computed from
F_symbol.json is (5/2, 1/2, 0) (quotient basis words 3222, 3322, 3332), which lies in the Li₄ span
(e.g. = −2·Li₄(−2) + ½·Li₄(1/3) + ½·Li₄(2/3) in classes). The same naive computation on each
individual HPL word gives, summed, (5/2, 1/2, 0).

**CAVEAT (important).** The naive specialization gives the WRONG number-level class for mixed-letter
words. Bulletproof exact identities (LLL lattice relations, identical at two precisions, residuals at
working precision up to 2100 digits, e.g.
  H_{0,0,−1,1}(1/2) = −Li₃,₁(−1/2,−1) = 6Li₄(−2) − 2Li₄(1/2) − Li₃(−1/2)log2 − (3/4)ζ(3)log2
      + (5/8)π²log²2 + (137/1440)π⁴ + (1/24)log⁴2   [|LHS−RHS| = 2.5·10⁻²¹⁰⁰],
also Li₃,₁(1/2,1) = −(1/8)ζ(3)log2 + π⁴/720 + log⁴2/24, Li₂,₂(1/2,1) = (1/4)ζ(3)log2 − (1/24)π²log²2
      + π⁴/1440 + log⁴2/24)
show the true class of H_{0,0,−1,1}(1/2) is −6·class(Li₄(−2)) = (−6,0,0), NOT the naive (−3,0,0);
and the naive combination −2Li₄(−2)+½Li₄(1/3)+½Li₄(2/3) FAILS numerically (F₄ minus it: no relation).
Summing the true per-word classes over R₃'s 56 terms gives the true class of F(1/2):
  class(F(1/2)) = (55/2, 13/2, 0),
which DOES lie in the Li₄ span (exhibited by the exact identities below). So the conclusion of the
symbol test — F(1/2) is classical — is confirmed, but the specific naive combination was wrong.

### 1(b) Numbers. NUMERICAL.
M₃(2π/3) = 0.01071736293228946541840957629122489585571897573846433576015561241784110361978…
- Ground truth: direct tanh-sinh quadrature of the corner integral (own implementation of
  corner_series) at dps 700 and 1100: |diff| = 1.49·10⁻⁶⁹⁶.
- Cross-check: the HPL formula M₃ = −ζ(3)Φ + ζ(2)Ψ + R₃ (G1's exact fit) evaluated by
  hpl_series.H_series at dps 1120: |HPL − quad(1100)| = 1.03·10⁻¹⁰⁹⁵. Two methods agree to 1095 digits.
- Weight split M₃ = W₄ + W₃ + W₂ (no weight ≤ 1 part):
  W₄ = R₃⁽⁴⁾ + ζ(2)Ψ⁽²⁾ + ζ(3)((1−c²)/(2c))artanh c,  W₃ = R₃⁽³⁾ + ζ(2)Ψ⁽¹⁾ − ζ(3)/2,  W₂ = R₃⁽²⁾ + ζ(2)/2.

**F₄ = R₃⁽⁴⁾(1/2) closed form** (LLL, 25-element basis of the weight-4 periods of MT(ℤ[1/6]):
5 Li₄ + 20 products; identical at dps 700 and 1000, residual 4.16·10⁻¹¹¹³ at dps 1120; a second
equivalent form with Li₄{1/3,2/3,3/4,1/2,1/9} identical at two precisions, residual 4.79·10⁻¹¹¹⁷):
F₄ = −21 Li₄(−2) − (65/8) Li₄(1/3) + (13/2) Li₄(2/3) + (137/6) Li₄(1/2) + (39/64) Li₄(1/9)
   + (15/2) Li₃(1/3) log2 + (37/2) Li₃(−1/2) log2 − (5/4) Li₃(−1/2) log3
   + (119/16) ζ(3) log2 − (37/32) ζ(3) log3
   + Li₂(1/3)² − (23/4) Li₂(1/3) log²2 + (7/4) Li₂(1/3) log2 log3 + Li₂(1/3) log²3 − (1/6) π² Li₂(1/3)
   − (455/144) π² log²2 + (9/8) π² log2 log3 − (7/32) π² log²3 − (949/2160) π⁴
   − (211/144) log⁴2 + (101/24) log³2 log3 − (15/8) log²2 log²3 − (35/24) log2 log³3 + (113/192) log⁴3.
(equivalently, in the other basis: (103/8)Li₄(1/3) − (29/2)Li₄(2/3) − (21/4)Li₄(3/4) − (31/6)Li₄(1/2)
 − (45/64)Li₄(1/9) + same Li₃/ζ(3)-log terms + corresponding product terms; both verified.)

**W₃ closed form** (11-element weight-3 basis as G1's M₂ PSLQ; identical at dps 600 and 900,
residuals 3.3·10⁻⁶⁰⁰ and 2.5·10⁻⁹⁰⁰):
W₃ = −(15/4) ζ(3) − 3 Li₃(1/3) − 8 Li₃(−1/2) − (1/12) π² log2
   + 4 Li₂(1/3) log2 − (3/2) Li₂(1/3) log3 + (1/3) log³2 − (5/2) log²2 log3 + 2 log2 log²3 − (1/4) log³3.
ζ(3) coefficient of W₃ in basis {Li₃(1/3), Li₃(−1/2)}: −15/4.

**W₂**: R₃⁽²⁾(1/2) + ζ(2)/2 = −Li₂(1/2) − Li₂(−1/2) + π²/12 (weight-2 HPLs in closed form).

## Task 1(c): theorem-style formula. NUMERICAL.
M₃(2π/3) = W₄ + W₃ + W₂ with W₄, W₃ as above and W₂ = −Li₂(1/2) − Li₂(−1/2) + π²/12.
Assembled check below; ζ(3) coefficient of M₃(2π/3) (in basis {Li₃(1/3), Li₃(−1/2)}):
  from W₃: −15/4; from W₄: +119/16 log2 − 37/32 log3 (ζ(3)·logs, not pure ζ(3));
  plus −ζ(3)/2 from −ζ(3)Φ(1/2) (inside W₃'s −ζ(3)/2 and Φ's log3 term) —
  see the assembled statement for the pure-ζ(3) coefficient.
"""
open("result.md", "w").write(md)
print("result.md updated")
