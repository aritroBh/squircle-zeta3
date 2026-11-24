md = """# COR5: at which corner angles is M₃(α) classical?

Work in progress. Builds on COR4 (`research/08_theory/COR4_M3alpha_classical`), where the weight-4
part F(c) of R₃(c) was shown to have the nonzero (2,2) obstruction
δ₂,₂(F) = (1−c⁻¹)([c]∧[1−c])∧([1−c]∧[1+c]) + (−1−c⁻¹)([c]∧[1+c])∧([1−c]∧[1+c]),
vanishing only where the letters collapse to a ≤ 2-dimensional span (then Λ²B₂ = 0 automatically).

Labels: PROVED = exact symbol algebra (number-level consequences modulo the period conjecture);
NUMERICAL = high-precision evidence.

## Task 1(a): hexagon corner α = 2π/3 (c = 1/2) — symbol-level classicality. DONE, PROVED.

**S-unit equation.** All solutions of x + y = 1 with x, y = ±2ᵃ3ᵇ: exactly 21 values of x:
  −8, −3, −2, −1, −1/2, −1/3, −1/8, 1/9, 1/4, 1/3, 1/2, 2/3, 3/4, 8/9, 9/8, 4/3, 3/2, 2, 3, 4, 9.
Completeness PROVED: after clearing denominators the equation is
e₁2^A₁3^B₁ + e₂2^A₂3^B₂ = 2^M3^N with min(A₁,A₂,M) = min(B₁,B₂,N) = 0; a mod-2/mod-3 case check
(printed in step output) shows one term must be ±1, so two 3-smooth numbers differ by 1; the
consecutive-3-smooth lemma (2−1, 3−2, 4−3, 9−8; proved by factoring mod 3 and mod 4) then applies.
The lemma regenerates exactly the 21 values (checked).

**Projection modulo products.** V = ℚ[2]⊕ℚ[3]; Q⁴ = T⁴(V)/D⁴ (D⁴ = shuffle decomposables) has
dimension 3 (= dim Lie₄(V); quotient basis words 3222, 3322, 3332 over letters (2,3)).
The decomposables have rank 13 (checked); shuffle products are 0 in Q⁴ (exhaustive check); the
class of S(Li₄(x)) = −[1−x][x]³ is nonzero (control).

**The test.** S(F(1/2)) (56 weight-4 HPL terms of R₃, specialised: [c]↦−[2], [1−c]↦−[2],
[1+c]↦[3]−[2], c^e↦2^{−e}) has 24 words. Its class in Q⁴ is (5/2, 1/2, 0).
The classes of Li₄(x) at the 21 S-unit solutions span Q⁴ (rank 3), and the class of F(1/2) lies
in their span — minimally with **3 arguments** (368 minimal solutions). PROVED (exact rational
arithmetic on symbols; rebuilt independently from R₃ and from F_symbol.json, both agree).

Explicit minimal combinations for π(S(F(1/2))):
  = −2·Li₄(−2) + ½·Li₄(1/3) + ½·Li₄(2/3)
  = 2·Li₄(−1/2) + ½·Li₄(1/3) + ½·Li₄(2/3)
  = −½·Li₄(1/3) − 3/2·Li₄(2/3) − ½·Li₄(3/4)
  = ¼·Li₄(1/4) + ½·Li₄(1/3) + ½·Li₄(2/3)
(all verified exactly). δ₂,₂(F(1/2)) = 0 as required (automatic, dim Λ²V = 1; checked).

## Task 1(b): numerics — IN PROGRESS
"""
open("result.md", "w").write(md)
print("result.md written")
