md = r"""# COR7: which corner angles make M₃ classical?

Work dir: `research/08_theory/COR7_classical_corners/`. Inputs (read-only): COR6 (true point symbols),
COR5 (hexagon closed form + 25-element weight-4 basis of MT(ℤ[1/6])), G1 (M₃ = −ζ(3)Φ + ζ(2)Ψ + R₃, all angles).

Labels: **PROVED** = exact symbol algebra (number-level classicality then modulo the period conjecture);
**NUMERICAL** = lattice relations identical at two precisions with residual at the working precision.

**Headline (status: Tasks 1 done; 2–4 in progress).**
1. At **c = 1/3** (α = 2 arccos(1/3) = 2.46191883468154936426…), M₃ **is classical**: explicit closed forms
   for the weight-4, weight-3, weight-2 parts in the COR5 bases of MT(ℤ[1/6]) (NUMERICAL, identical at
   dps 700/1000 (W₄) and 600/900 (W₃), residuals at working precision). M₃ computed to 1095 digits by two
   independent methods (own tanh-sinh quadrature with a precision-boosted integrand vs the exact HPL
   formula): |quad − HPL| = 5.01·10⁻¹⁰⁹⁸. The Li₄ part's Lie class equals COR6's T4 class **exactly**
   (PROVED symbol algebra): class = −(11/81)Li₄(−8) − 2Li₄(−3) + (20/3)Li₄(2/3), coords (14, 20/3, 8/3) ∈ Q⁴.

---

## Task 1: the c = 1/3 corner (α = 2 arccos(1/3))

### 1(a) The value (two independent methods, 1095 digits)

- Method 1 (ground truth): tanh-sinh quadrature of the corner integral
  M(u;α) = ∫₀¹[(1+g)^(−1−u)√(1+g²ᐧ²⁻²ᵘ − 2cos α·g¹⁻ᵘ) − 1]g^(u−1)dg, u³-coefficient, full angle
  α = 2 arccos(1/3). **Fix vs COR5's copy of corner_series.py**: the u-series integrand suffers
  catastrophic cancellation ((1+g)^(−1)√(1−2cg+g²) − 1) at nodes g ≲ 10^(−dps); we boost the internal
  precision by ⌈−log₁₀ g⌉+15 digits at nodes g < 0.1 (`cor7_lib.integrand_series_hp`). Without the boost
  the quadrature stalls at ~17 digits at dps ≥ 100 (verified). With it: |quad(1100) − HPL(1150)| = 5.01·10⁻¹⁰⁹⁸.
- Method 2: G1's exact HPL formula M₃ = −ζ(3)Φ(c) + ζ(2)Ψ(c) + R₃(c) evaluated by `cor7_lib.H_series`
  (COR4's exact recursive series) at dps 1150.

  M₃(2 arccos(1/3)) = 0.0045712107180721471097362832226723427474247040964…  (`c13_M3_hpl_1150.txt`, 1145 digits; `c13_M3_quad_1100.txt`, 1095 digits; agree to 5.01·10⁻¹⁰⁹⁸)

### 1(b) Weight split (COR5 convention; W₄+W₃+W₂ = M₃ to 1.7·10⁻¹¹²¹-style exactness of the split)

W₄ = R₃⁽⁴⁾ + ζ(2)Ψ⁽²⁾ + ζ(3)·((1−c²)/(2c))artanh c,  W₃ = R₃⁽³⁾ + ζ(2)Ψ⁽¹⁾ − ζ(3)/2,  W₂ = R₃⁽²⁾ + ζ(2)/2.

At c = 1/3:  W₄ = 1.29119734012229401049…, W₃ = −2.05191305933908006632…, W₂ = 0.76528692993485820294….

### 1(c) Closed forms (LLL; identical at two precisions; residual at working precision)

**W₂** = −Li₂(1/3) − Li₂(−1/3) + π²/12  (exact from the HPL formula; |W₂ − rhs| = 9.2·10⁻¹¹²²).

**W₃** (COR5's 11-element weight-3 basis; identical at dps 600 and 900; residual 6.5·10⁻¹⁰⁴⁹ at dps 1050):

> W₃ = −26·ζ(3) − 87·Li₃(1/3) − 130·Li₃(−1/2) + (21/4)π²log2 − (43/12)π²log3
>      + 67·Li₂(1/3)log2 − 45·Li₂(1/3)log3 + (65/3)log³2 − 65·log²2log3 + (109/2)log2log²3 − (37/3)log³3.

**W₄** (COR5's 25-element weight-4 basis = 5 Li₄ + 20 products; identical at dps 700 and 1000; residual
1.28·10⁻¹¹¹⁶ at dps 1120; full data `c13_F4_relation.json`):

> W₄ = −(22/3)Li₄(−2) − (17/2)Li₄(1/3) + (20/3)Li₄(2/3) + (110/9)Li₄(1/2) + (25/48)Li₄(1/9)
>      − 12·Li₃(1/3)log2 + (28/3)Li₃(1/3)log3 − (16/3)Li₃(−1/2)log2 + (16/3)Li₃(−1/2)log3
>      + (14/3)ζ(3)log2 − (25/9)ζ(3)log3
>      − 4·Li₂(1/3)² + (13/3)Li₂(1/3)log²2 − 6·Li₂(1/3)log2log3 + (2/3)Li₂(1/3)log²3 + (11/18)π²Li₂(1/3)
>      − (31/27)π²log²2 + (31/72)π²log²3 − (271/1296)π⁴
>      + (46/27)log⁴2 − (62/9)log³2log3 + (61/6)log²2log²3 − (17/3)log2log³3 + (103/144)log⁴3.

Total: M₃ − (W₄+W₃+W₂ closed forms) = 1.2·10⁻¹¹¹⁹.

Also Ψ⁽¹⁾(1/3) = (21/2)log2 − (15/2)log3 (residual 1.8·10⁻¹¹²¹); (1−c²)/(2c) = 4/3, artanh(1/3) = (log2)/2.

### 1(d) Lie class = COR6's T4 class (PROVED symbol algebra)

COR6 (verify_cor6.py, T4): class(Ω(1/3)) = class(−(11/81)Li₄(−8) − 2Li₄(−3) + (20/3)Li₄(2/3)).
The Li₄ part of the NUMERICAL W₄ closed form has Q⁴-coords (14, 20/3, 8/3) — **exactly equal** to the
T4 coords and to class(Ω(1/3)) computed from the true-symbol code. (The coefficient 20/3 of Li₄(2/3)
is even literally the same.)

### 1(e) ζ(3) coefficient (basis-dependent)

In COR5's basis {Li₃(1/3), Li₃(−1/2)}: pure-ζ(3) coefficient **−26**;
ζ(3)-log terms (14/3)ζ(3)log2 − (25/9)ζ(3)log3.

---

## Tasks 2–4: IN PROGRESS (census, ℚ(√2) corner, final report + verify_cor7.py)
"""
open('result.md', 'w').write(md)
print("result.md written,", len(md), "bytes")
