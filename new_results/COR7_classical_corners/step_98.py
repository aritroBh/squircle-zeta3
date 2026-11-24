md = r"""# COR7: which corner angles make M₃ classical?

Work dir: `research/08_theory/COR7_classical_corners/`. Read-only inputs: COR6 (true number-level
symbols), COR5 (hexagon closed form; 25-element weight-4 basis of MT(ℤ[1/6])), G1
(M₃ = −ζ(3)Φ(c) + ζ(2)Ψ(c) + R₃(c) at every angle), the paper's thm:W4R (square-corner ℚ(√2) basis).

Labels: **PROVED** = exact symbol algebra / exact certified integer relations (number-level
classicality then holds modulo the standard period conjecture); **NUMERICAL** = lattice relations
identical at two precisions with residual at the working precision.

**Headline.**
1. The classical-at-the-symbol-level corners of M₃ are: c₀ = 1/2 (hexagon, known), 1/√2 (square,
   known), **1/3** (COR6's third point; here made fully classical with an explicit 1095-digit closed
   form), and four corners **new** to this work: **√2−1**, **(√2−1)²**, **1/φ (golden)**, **1/φ³ = √5−2**.
   Rational c₀: completeness **proved** ({1/3, 1/2} only). Quadratic c₀: exact census over the 15
   class-number-1 fields ℚ(√D), D ≤ 50, coefficients |a|,|b|,e ≤ 16 (a *census*, not a completeness
   proof): exactly the five quadratic points above.
2. c = 1/3: M₃ computed to 1095 digits by two independent methods (|quad − HPL| = 5.0·10⁻¹⁰⁹⁸);
   explicit closed forms for W₄ (COR5's 25-element basis), W₃ (COR5's 11-element basis), W₂;
   the Li₄ part's Lie class equals COR6's T4 class **exactly**; ζ(3) coefficient −26.
3. c = √2−1: M₃ to 1095 digits by two methods (5.0·10⁻¹⁰⁹⁸); closed forms with the ℚ(√2) weight-4
   basis of thm:W4R (7 Li₄ + 22 products) with ℚ(√2) coefficients; Li₄ part's Lie class = exact
   class(Ω(√2−1)); ζ(3) coefficient −5/4 − 63√2/32.
4. Every dim-V = 3 control tested (12 split-field points in ℚ(√17), ℚ(√41); 1/φ²; (√2−1)⁴;
   rational 1/5; (√2−1)³ at dim 4) has **nonzero** true δ₂,₂ — no accidental vanishing seen.

---

## Task 1: the c = 1/3 corner, α = 2 arccos(1/3) = 2.4619188346815493643 rad (141.0575587°)

### 1(a) The value — two independent methods, 1095 digits

- Method 1 (ground truth): own tanh-sinh quadrature of the corner integral's u³-coefficient at the
  FULL angle α (convention of G1's corner_series.py). **Bug found and fixed**: the u-series integrand
  loses all precision at quadrature nodes g ≲ 10^(−dps) (the product (1+g)⁻¹√(1−2cg+g²) − 1 cancels);
  COR5's copy stalls at ~17 digits at dps ≥ 100 (demonstrated). Fix: boost the internal precision by
  ⌈−log₁₀ g⌉+15 digits at nodes g < 0.1 (`cor7_lib.integrand_series_hp`). With the fix the quadrature
  converges to full precision (checked at dps 100–700 against the HPL value to ~1e-102 … 1e-698).
- Method 2: G1's exact HPL formula evaluated by the exact recursive series (`cor7_lib.H_series`,
  COR4's code) at dps 1150.

  M₃ = 0.0045712107180721471097362832226723427474247040964… (`c13_M3_hpl_1150.txt`, 1145 digits;
  `c13_M3_quad_1100.txt`, 1095 digits). **|quad(1100) − HPL(1150)| = 5.01·10⁻¹⁰⁹⁸.**

### 1(b) Weight split (COR5 convention): M₃ = W₄ + W₃ + W₂, |split − M₃| = 0 at dps 1150

W₄ = R₃⁽⁴⁾ + ζ(2)Ψ⁽²⁾ + ζ(3)·((1−c²)/(2c))artanh c, W₃ = R₃⁽³⁾ + ζ(2)Ψ⁽¹⁾ − ζ(3)/2, W₂ = R₃⁽²⁾ + ζ(2)/2.
At c = 1/3: (1−c²)/(2c) = 4/3, artanh(1/3) = (log 2)/2; Ψ⁽¹⁾(1/3) = (21/2)log2 − (15/2)log3 (1.8e-1121).

### 1(c) Closed forms (LLL; identical at two precisions; residual at working precision)

**W₂** = −Li₂(1/3) − Li₂(−1/3) + π²/12  (exact from the HPL formula; residual 9.2·10⁻¹¹²² at dps 1120).

**W₃** (COR5's 11-element weight-3 basis; identical at dps 600 and 900; residual 6.5·10⁻¹⁰⁴⁹ at dps 1050;
stored `c13_W3_relation.json`):

> W₃ = −26·ζ(3) − 87·Li₃(1/3) − 130·Li₃(−1/2) + (21/4)π²log2 − (43/12)π²log3
>      + 67·Li₂(1/3)log2 − 45·Li₂(1/3)log3 + (65/3)log³2 − 65·log²2log3 + (109/2)log2log²3 − (37/3)log³3.

**W₄** (COR5's 25-element weight-4 basis = 5 Li₄ + 20 products; identical at dps 700 and 1000;
residual 1.28·10⁻¹¹¹⁶ at dps 1120; height 15840; stored `c13_F4_relation.json`):

> W₄ = −(22/3)Li₄(−2) − (17/2)Li₄(1/3) + (20/3)Li₄(2/3) + (110/9)Li₄(1/2) + (25/48)Li₄(1/9)
>      − 12·Li₃(1/3)log2 + (28/3)Li₃(1/3)log3 − (16/3)Li₃(−1/2)log2 + (16/3)Li₃(−1/2)log3
>      + (14/3)ζ(3)log2 − (25/9)ζ(3)log3
>      − 4·Li₂(1/3)² + (13/3)Li₂(1/3)log²2 − 6·Li₂(1/3)log2log3 + (2/3)Li₂(1/3)log²3 + (11/18)π²Li₂(1/3)
>      − (31/27)π²log²2 + (31/72)π²log²3 − (271/1296)π⁴
>      + (46/27)log⁴2 − (62/9)log³2log3 + (61/6)log²2log²3 − (17/3)log2log³3 + (103/144)log⁴3.

Total: M₃ − (W₄+W₃+W₂ closed forms) = 1.2·10⁻¹¹¹⁹ (dps 1120).

### 1(d) Lie class = COR6's T4 class (PROVED symbol algebra)

class(Ω(1/3)) = class(−(11/81)Li₄(−8) − 2Li₄(−3) + (20/3)Li₄(2/3)) = **(14, 20/3, 8/3)** in Q⁴ =
T⁴(ℚ²)/D⁴; the Li₄ part of the NUMERICAL W₄ closed form has exactly the same coords (checked with
COR6's W4Quotient; the coefficient 20/3 of Li₄(2/3) is literally shared).

### 1(e) ζ(3) coefficient (basis-dependent)

In COR5's basis {Li₃(1/3), Li₃(−1/2)}: pure-ζ(3) coefficient **−26**; ζ(3)-log terms in W₄:
(14/3)ζ(3)log2 − (25/9)ζ(3)log3.

---

## Task 2: census of classical-at-the-symbol-level corners

δ₂,₂ ∈ Λ²B₂ vanishes automatically when the letter span V of {2, c₀, 1−c₀, 1+c₀} has dim ≤ 2
(then Λ²V is 1-dimensional and Λ²(Λ²V) = 0).

### 2(a) Rational c₀ — PROVED complete

**Theorem.** For c₀ = a/b ∈ (0,1), gcd(a,b) = 1: rank⟨2, c₀, 1−c₀, 1+c₀⟩ ≤ 2 ⟺ c₀ ∈ {1/3, 1/2}.

**Proof.** For each odd prime p write w_p = (v_p(c₀), v_p(1−c₀), v_p(1+c₀)). Since gcd(a,b) = 1 and
p is odd, p divides at most one of a, b−a, a+b; and p | b ⇒ p ∤ a(b−a)(a+b). Hence the p-column of
the valuation matrix W is −k(1,1,1) (if p | b) or k·e_i (if p | a, b−a, or a+b respectively),
k ≥ 1. rank of the group = 1 + rank(W) (the prime 2 is isolated), so rank ≤ 2 ⇔ rank(W) ≤ 1.
*Case A: some odd p | b.* Any further odd prime q contributes a column k'e_i, independent of
(1,1,1) — contradiction. So a, b−a, a+b are powers of 2: a = 2^α, b−a = 2^β, a+b = 2^γ. Then
2^β(2^{γ−β}−1) = 2a = 2^{α+1}, forcing 2^{γ−β}−1 = 1: γ = β+1, β = α+1, b = 3·2^α; gcd(a,b) = 2^α
gives α = 0: **c₀ = 1/3** (G = ⟨2,3⟩ ✓). *Case B: b = 2^s.* Every odd-prime column is a unit
vector; rank ≤ 1 forces all odd primes into the same row i. i=0: 2b = 2^β + 2^γ with β ≠ γ has an
odd factor — impossible. i=1: 2b = 2^γ + 2^β p^k has the odd factor (2^{γ−β}+p^k) — impossible.
i=2: 2b = 2^α + 2^β forces α = β, b = 2^{α+1}, a+b = 3·2^α = 2^γ p^k ⇒ p = 3, k = 1; gcd gives
α = 0: **c₀ = 1/2**. No odd primes at all is impossible (2b = 2^β+2^γ ⇒ β = γ ⇒ a = 0). ∎

Machine check: enumeration to denominator 5000 finds exactly {1/3, 1/2} (verify: to 300).

### 2(b) Quadratic c₀ — exact census (search, not a completeness proof)

Searched: all 15 real quadratic fields ℚ(√D) with class number 1, D ≤ 50 squarefree
(D = 2,3,5,6,7,11,13,17,19,21,23,29,33,37,41), c₀ = (a+b√D)/e, b ≠ 0, |a|,|b|,e ≤ 16,
gcd(|a|,|b|,e) = 1, 0 < c₀ < 1. For every candidate the letter-span rank was computed **exactly**:
generators (fundamental unit ε; one element per rational prime in the letter norms, certified by
exact NF multiplication — `cor7_pt.py`) and exponent vectors found by a norm-guided search and
verified exactly. (A log-lattice filter was tried first and found to be **vacuous** for S = {2} in
quadratic fields: v₂(N(x))·log2 = log|N(x)| = Σlog|σᵢx|, so the LV image lies in a 2-plane; the
exact machinery replaced it.) Result: **exactly 5 points with letter-span rank ≤ 2**:

| c₀ | field | α = 2 arccos c₀ (rad) | α (deg) | dim V | true δ₂,₂ | Lie class in Li₄(S-units) span |
|---|---|---|---|---|---|---|
| 1/√2 (square, n=4) | ℚ(√2) | 1.5707963267948966192 (π/2) | 90.0 | 2 | = 0 ✓ | yes (coords (0,20,0), (−26,0,155)) |
| √2−1 | ℚ(√2) | 2.2874354808048409875 | 131.0603990 | 2 | = 0 ✓ | yes ((8,−2,−6), (−3,1,8)) |
| (√2−1)² | ℚ(√2) | 2.7967406581640951061 | 160.2414361 | 2 | = 0 ✓ | yes ((84,−100,128), (−4,−16,48)) |
| 1/φ (golden) | ℚ(√5) | 1.8091137886047627283 | 103.6545847 | 2 | = 0 ✓ | yes ((111/2,0,−12), (−51/2,0,12)) |
| 1/φ³ = √5−2 | ℚ(√5) | 2.6649577299700610204 | 152.6908305 | 2 | = 0 ✓ | yes ((84,−21,6), (−18,3,0)) |

For each: δ₂,₂(Ω(c₀)) = 0 confirmed with the true-symbol code (exact; COR6's code for the three
ℚ(√2)/ℚ(√5) points with dim V = 2, cor7_pt for all five). Li₄ classes at the fields' {2}-unit
arguments (29 for ℚ(√2), 25 for ℚ(√5), exactly enumerated) span Q⁴-quotient = ℚ³, and **every
field component of class(Ω(c₀)) lies in their ℚ-span** (exact linear algebra). E.g. the golden
corner: class(Ω(1/φ)) = class of (7/6 − √5/18)Li₄(φ³) + (−4 + 4√5)Li₄(φ²/2) + (8 − 8√5)Li₄(φ/2).

### 2(c) dim V = 3 (and one dim 4) controls — no accidental zeros

All with exact true-symbol δ₂,₂: ten points in ℚ(√17) (√17−4, (√17−3)/2, (5−√17)/2, (√17−3)/4,
(√17−1)/4, (5−√17)/4, (√17−1)/8, (√17+1)/8, (7−√17)/16, (7+√17)/16), two in ℚ(√41) ((√41−5)/2,
(√41−5)/8): all letter-span rank 3, δ₂,₂ **nonzero** (3 components each). Also 1/φ² (rank 3, 2 comp),
1/(2φ) = sin(π/10) (rank 3, 3 comp), (√2−1)⁴ (rank 3, 3 comp), (√2−1)³ (rank 4, 12 comp), √3/2
(rank 3, 2 comp), and rationals 1/5 (3 comp), 2/5, 3/7 (COR6 table: 15, 12 comp): δ₂,₂ **nonzero**
everywhere. COR6's deg ≥ 4 rows (n = 5,7,8,9,11,12, cos(π/12)) likewise all nonzero.

### 2(d) The clean COR6 table (n = 3…12, π/6, 2/5, 3/7, 1/3) — naive vs TRUE δ₂,₂

| point | field deg | dim V | letter span | naive δ₂,₂ | TRUE δ₂,₂ | differ? |
|---|---|---|---|---|---|---|
| n=3 triangle, c=√3/2 | 2 | 3 | 3 | ≠0 | ≠0 (2 comp) | yes |
| n=4 square, c=1/√2 | 2 | 2 | 2 | 0 | **0** | no |
| n=5, c=sin(π/5) | 4 | 5 | 4 | ≠0 | ≠0 (36) | yes |
| n=6 hexagon, c=1/2 | 1 | 2 | 2 | 0 | **0** | no |
| n=7, c=sin(π/7) | 6 | 7 | 4 | ≠0 | ≠0 (6) | yes |
| n=8 octagon, c=sin(π/8) | 4 | 4 | 3 | ≠0 | ≠0 (12) | yes |
| n=9, c=sin(π/9) | 6 | 7 | 4 | ≠0 | ≠0 (6) | yes |
| n=10, c=sin(π/10)=1/(2φ) | 2 | 3 | 3 | ≠0 | ≠0 (3) | yes |
| n=11, c=sin(π/11) | 10 | 11 | 4 | ≠0 | ≠0 (6) | yes |
| n=12, c=sin(π/12) | 4 | 4 | 3 | ≠0 | ≠0 (3) | yes |
| α=π/6, c=cos(π/12) | 4 | 4 | 3 | ≠0 | ≠0 (8) | yes |
| c=1/3 | 1 | 2 | 2 | 0 | **0** | no |
| c=2/5 | 1 | 4 | 4 | ≠0 | ≠0 (15) | yes |
| c=3/7 | 1 | 4 | 4 | ≠0 | ≠0 (12) | yes |

Zero δ₂,₂ among all tested points ⟺ letter span ≤ 2 ⟺ {n=4, n=6, 1/3} in COR6's table; COR7 adds
{√2−1, (√2−1)², 1/φ, 1/φ³}. The naive (function-level) δ₂,₂ differs from the TRUE (number-level)
one at every point where both are nonzero (11 rows), and at the three zero rows it happens to agree.

---

## Task 3: the c₀ = √2−1 corner (α = 2.2874354808048409875 rad, 131.0603990°) — ℚ(√2) closed form

M₃ = 0.007189129721152490348077741879814010022297… (`b21_M3_hpl_1150.txt`, 1145 digits, HPL series;
`b21_M3_quad_1100.txt`, 1095 digits, boosted tanh-sinh quadrature). **|quad − HPL| = 5.01·10⁻¹⁰⁹⁸.**
Weights: W₄ = 1.21768419197083980701…, W₃ = −1.94318495181611313625…, W₂ = 0.73268988956642581959….
With b = √2−1, L = log(1+√2), λ = log2, and coefficients in ℚ(√2) (p + q√2):

**W₂** = −½·Li₂(b²) + π²/12 — exact (duplication Li₂(b)+Li₂(−b) = ½Li₂(b²), residual 0 at dps 1120).

**W₃** (paper's 12-element weight-3 basis × {1,√2}; identical at dps 700 and 1000; residual
8.05·10⁻¹¹¹⁸ at dps 1120; stored `b21_W3_relation.json`):

> W₃ = (−5/4 − 63√2/32)ζ(3) + (−11−13√2)Li₃(b) + (5+4√2)Li₃(−b) + (14+16√2)Li₃(1/√2)
>      + (−4−4√2)L·Li₂(b) + (4+4√2)λ·Li₂(b)
>      + (−11/12 − 49√2/48)π²L + (29/48 + 3√2/4)π²λ
>      + (−7/24 − √2/3)λ³ + (−7/4 − 2√2)λ²L + (11/4 + 3√2)λL² + (−1/3 − √2/4)L³.

**W₄** (thm:W4R's ℚ(√2) weight-4 basis: 7 Li₄ {b, −b, −b², 1/√2, 1/2, 1−1/√2, 2−√2} + 22 products
{L·(Li₃(b), Li₃(−b), Li₃(1/√2), ζ(3)), λ·(same), Li₂(b)², Li₂(b)·(L², λL, λ², π²), π⁴, π²·(L², λL, λ²),
L⁴, L³λ, L²λ², Lλ³, λ⁴} × {1,√2}; 59-dim lattice; identical at dps 700 and 1000; residual
1.35·10⁻¹¹¹⁶ at dps 1120; height 276480; stored `b21_W4_relation.json`):

> W₄ = (22 − 31√2/2)Li₄(b) + (24 − 27√2/2)Li₄(−b) + (−5/2 + 13√2/8)Li₄(−b²) − √2·Li₄(1/√2)
>      + (3/4)Li₄(1/2) + (−4+4√2)Li₄(1−1/√2) + (−2+3√2)Li₄(2−√2)
>      − 3L·Li₃(b) + L·Li₃(−b) − 2√2·L·Li₃(1/√2) + (53/16 + 2√2)L·ζ(3)
>      + 2λ·Li₃(b) + (−2+√2)λ·Li₃(−b) + (−2+2√2)λ·Li₃(1/√2) + (−25/16 − √2)λ·ζ(3)
>      + (2+√2)Li₂(b)² + (2+√2/2)Li₂(b)L² + (−2+√2/2)Li₂(b)λL + (1/4−√2/4)Li₂(b)λ² + (−1/2−√2/8)π²Li₂(b)
>      + (49/960 − 29√2/11520)π⁴ + (−5/8+√2/48)π²L² + (9/16−11√2/48)π²λL + (−3/16+5√2/32)π²λ²
>      + (1−√2/48)L⁴ + (−5/4+√2/4)L³λ + (5/8−√2/4)L²λ² + (−7/24+13√2/48)Lλ³ + (13/192−13√2/384)λ⁴.

**Lie-class check (PROVED symbol algebra):** the Li₄ part of this NUMERICAL closed form has
Q⁴-components (8, −2, −6) and (−3, 1, 8) — **exactly** class(Ω(√2−1)) from the true-symbol code.

Total: M₃ − (W₄+W₃+W₂ closed forms) = 5.5·10⁻¹¹²⁰ (dps 1120).
**ζ(3) coefficient** (basis {Li₃(b), Li₃(−b), Li₃(1/√2)}): pure ζ(3): **−5/4 − 63√2/32**;
ζ(3)-log terms in W₄: (53/16 + 2√2)L·ζ(3) + (−25/16 − √2)λ·ζ(3).

---

## Task 4: verify_cor7.py (runs in 7 s; 35/35 PASS)

`verify_cor7.py` checks: A1–A6 (c=1/3: split consistency, the three closed forms at dps 420 against
the stored 1145-digit values, quadrature at dps 60 vs stored, the exact T4 class equality, the −26
ζ(3) coefficient); B1–B4 (rational census to denom 300; exact rank + δ₂,₂ = 0 for the 5 quadratic
corners; Li₄-span = ℚ³ and class(Ω) membership for all 5; dim-3/4 controls nonzero); C1–C6 (√2−1:
the three closed forms, the W₂ duplication form, the Li₄-part Lie class match, quadrature vs
stored); D (recomputes the degree ≤ 2 rows of COR6's table with COR6's code and compares
true/naive/differ against `cor6_table.json`). Run: `TOTAL: 35/35 PASS, time 7s`.

## Evidence table

| Check | Label | Result |
|---|---|---|
| M₃(2 arccos(1/3)): quadrature (1100) vs HPL (1150) | NUMERICAL | 5.01·10⁻¹⁰⁹⁸ (1095 digits) |
| W₄ closed form (c=1/3), 25-element basis | NUMERICAL | identical dps 700/1000; residual 1.28·10⁻¹¹¹⁶ @1120 |
| W₃ closed form (c=1/3), 11-element basis | NUMERICAL | identical dps 600/900; residual 6.5·10⁻¹⁰⁴⁹ @1050 |
| W₂ = −Li₂(1/3) − Li₂(−1/3) + π²/12 | PROVED (from HPL formula) | 9.2·10⁻¹¹²² |
| M₃ = W₄+W₃+W₂ closed forms (c=1/3) | NUMERICAL | 1.2·10⁻¹¹¹⁹ |
| Li₄ part class = COR6 T4 = class(Ω(1/3)) | PROVED (symbol) | coords (14, 20/3, 8/3) exact |
| rational census: only {1/3, 1/2} | PROVED (prime-pattern case analysis) | enumerated to denom 5000 |
| 5 quadratic corners: letter span 2, true δ₂,₂ = 0 | PROVED (exact certified symbols) | all pass |
| Li₄(S-units) span = ℚ³ and ∋ class(Ω), all 5 corners | PROVED (exact linear algebra) | 29/25 args; all components in span |
| dim-3 controls (12 split-field, 1/φ², (√2−1)⁴, 1/5; (√2−1)³ dim 4) | PROVED (exact) | δ₂,₂ ≠ 0 everywhere |
| M₃(2 arccos(√2−1)): quadrature (1100) vs HPL (1150) | NUMERICAL | 5.01·10⁻¹⁰⁹⁸ |
| W₄ closed form (√2−1), thm:W4R basis × {1,√2} | NUMERICAL | identical dps 700/1000; residual 1.35·10⁻¹¹¹⁶ @1120 |
| W₃ closed form (√2−1), 12-element basis × {1,√2} | NUMERICAL | identical dps 700/1000; residual 8.05·10⁻¹¹¹⁸ @1120 |
| W₂ = −½Li₂(b²) + π²/12 | PROVED (duplication) | residual 0 @1120 |
| Li₄ part class = class(Ω(√2−1)) | PROVED (symbol) | (8,−2,−6), (−3,1,8) exact |
| M₃ = W₄+W₃+W₂ closed forms (√2−1) | NUMERICAL | 5.5·10⁻¹¹²⁰ |
| verify_cor7.py | — | 35/35 PASS, 7 s |
| quadrature integrand cancellation (COR5's copy stalls at ~17 digits at dps ≥ 100) | demonstrated; fixed in cor7_lib | before/after: 1e-17 → 1e-298+ |
| closed forms for the golden corner and (√2−1)² numbers | — | not attempted (symbol level only) |

## Key files (workdir)

`cor7_lib.py` (M₃ formula, boosted quadrature, LLL), `cor7_pt.py` (exact certified letter vectors for
quadratic fields), `verify_cor7.py`, `result.md`, values `c13_*_hpl_1150.txt`, `c13_M3_quad_1100.txt`,
`b21_*_hpl_1150.txt`, `b21_M3_quad_1100.txt`, relations `c13_F4_relation.json`, `c13_W3_relation.json`,
`b21_W3_relation.json`, `b21_W4_relation.json`, `golden_li4_combo.json`, `census_rank2.json`,
`lie_classes.pkl`. Key steps: cor7_lib (step_30–32), quadrature fix (step_43–48), c13 values
(step_49–50), c13 LLL (step_51–52), T4 (step_53), census (step_56–73), Li₄ spans (step_74–75),
b21 values (step_80), golden combo (step_84), b21 LLL (step_85–87), b21 checks (step_88–89),
verify_cor7.py (step_91–96).

**Bottom line.** The classical corners of M₃ at the symbol level are exactly: rational 1/2, 1/3
(completeness proved), and — in the searched region — quadratic 1/√2, √2−1, (√2−1)² (field ℚ(√2)),
1/φ, 1/φ³ (field ℚ(√5)). All have dim V = 2, true δ₂,₂ = 0, and Lie class in the Li₄ span at the
field's S-units. For 1/3 and √2−1 the number M₃ is classical explicitly (1095 digits, two methods;
closed forms identical at two precisions; Li₄ part's Lie class matching the exact symbol class).
"""
open('result.md', 'w').write(md)
print("result.md written:", len(md), "bytes")
