# COR7 — which corner angles make M₃ classical? (complete report)

Work dir `new_results/COR7_classical_corners/`. Deliverables: `result.md` (full report), `verify_cor7.py` (**35/35 PASS, 7 s**), `cor7_lib.py`, `cor7_pt.py`, values/relations (`c13_*`, `b21_*`, `*_relation.json`, `golden_li4_combo.json`, `census_rank2.json`), and step files step_01–97. Labels: **PROVED** = exact symbol algebra / exact certified integer relations (number-level classicality then modulo the period conjecture); **NUMERICAL** = lattice relations identical at two precisions, residual at working precision.

## Headline

The classical-at-the-symbol-level corners of M₃ are: c₀ = 1/2 (hexagon), 1/√2 (square), **1/3** (COR6's third point, now fully classical with a 1095-digit closed form), and four **new** corners: **√2−1**, **(√2−1)²**, **1/φ (golden)**, **1/φ³ = √5−2**. Rational c₀: completeness **proved** ({1/3, 1/2} only). Quadratic c₀: exact census (class-number-1 fields D ≤ 50, coefficients ≤ 16 — a census, not a proof of completeness): exactly those five quadratic points. All have dim V = 2, true δ₂,₂ = 0, and Lie class in the Li₄ span at the field's S-units. Every dim-V = 3 control tested has nonzero true δ₂,₂ — no accidental vanishing.

## Task 1: c = 1/3 (α = 2 arccos(1/3) = 2.4619188346815493643 rad = 141.0575587°)

- **M₃ = 0.0045712107180721471097362832226723427474247040964…** — 1095 digits by two independent methods: own tanh-sinh quadrature at the full angle (with a precision-boost fix: COR5's integrand u-series loses all precision at nodes g ≲ 10^(−dps) from the cancellation (1+g)⁻¹√(1−2cg+g²) − 1; demonstrated stalling at ~17 digits at dps ≥ 100; fixed by boosting internal precision at small nodes) vs G1's exact HPL formula by exact series at dps 1150. **|quad(1100) − HPL(1150)| = 5.01·10⁻¹⁰⁹⁸.**
- Weight split (COR5 convention): W₄ = R₃⁽⁴⁾ + ζ(2)Ψ⁽²⁾ + ζ(3)·((1−c²)/(2c))artanh c, W₃ = R₃⁽³⁾ + ζ(2)Ψ⁽¹⁾ − ζ(3)/2, W₂ = R₃⁽²⁾ + ζ(2)/2. At c=1/3: (1−c²)/(2c) = 4/3, artanh(1/3) = (log2)/2; Ψ⁽¹⁾(1/3) = (21/2)log2 − (15/2)log3 (resid 1.8·10⁻¹¹²¹).
- **W₂** = −Li₂(1/3) − Li₂(−1/3) + π²/12 (exact; 9.2·10⁻¹¹²²).
- **W₃** (COR5 11-element basis; identical dps 600/900; residual 6.5·10⁻¹⁰⁴⁹ @1050): W₃ = −26ζ(3) − 87Li₃(1/3) − 130Li₃(−1/2) + (21/4)π²log2 − (43/12)π²log3 + 67Li₂(1/3)log2 − 45Li₂(1/3)log3 + (65/3)log³2 − 65log²2log3 + (109/2)log2log²3 − (37/3)log³3.
- **W₄** (COR5 25-element basis; identical dps 700/1000; residual 1.28·10⁻¹¹¹⁶ @1120): W₄ = −(22/3)Li₄(−2) − (17/2)Li₄(1/3) + (20/3)Li₄(2/3) + (110/9)Li₄(1/2) + (25/48)Li₄(1/9) − 12Li₃(1/3)log2 + (28/3)Li₃(1/3)log3 − (16/3)Li₃(−1/2)log2 + (16/3)Li₃(−1/2)log3 + (14/3)ζ(3)log2 − (25/9)ζ(3)log3 − 4Li₂(1/3)² + (13/3)Li₂(1/3)log²2 − 6Li₂(1/3)log2log3 + (2/3)Li₂(1/3)log²3 + (11/18)π²Li₂(1/3) − (31/27)π²log²2 + (31/72)π²log²3 − (271/1296)π⁴ + (46/27)log⁴2 − (62/9)log³2log3 + (61/6)log²2log²3 − (17/3)log2log³3 + (103/144)log⁴3. Total M₃ − (W₄+W₃+W₂) = 1.2·10⁻¹¹¹⁹.
- **Lie class (PROVED):** Li₄ part's Q⁴ coords = **(14, 20/3, 8/3)** = COR6's T4 class −(11/81)Li₄(−8) − 2Li₄(−3) + (20/3)Li₄(2/3) = class(Ω(1/3)) (exact symbol algebra).
- **ζ(3) coefficient (basis-dependent):** pure −26 (basis {Li₃(1/3), Li₃(−1/2)}); ζ(3)-log terms (14/3)ζ(3)log2 − (25/9)ζ(3)log3.

## Task 2: census (letters {2, c₀, 1−c₀, 1+c₀}; δ₂,₂ = 0 automatic when dim V ≤ 2)

**Rational — PROVED complete.** For c₀ = a/b (gcd 1): rank ≤ 2 ⟺ c₀ ∈ {1/3, 1/2}. Proof: for odd p, the p-column of the valuation matrix W is −k(1,1,1) if p | b (p then divides none of a, b−a, a+b) and k·e_i if p | a (resp. b−a, a+b); group rank = 1 + rank(W). Case A (odd p | b): any other odd prime contradicts rank ≤ 1, so a, b−a, a+b are powers of 2; 2^β(2^{γ−β}−1) = 2a forces γ = β+1, β = α+1, b = 3·2^α, gcd ⇒ α = 0: c₀ = 1/3. Case B (b = 2^s): all odd primes must sit in one row; rows 0 and 1 give 2b or 2b/p^k an odd factor > 1 (impossible); row 2 gives 2b = 2^α+2^β ⇒ α = β, b = 2^{α+1}, a+b = 3·2^α ⇒ p = 3, α = 0: c₀ = 1/2. ∎ (Enumerated to denominator 5000: exactly {1/3, 1/2}.)

**Quadratic — exact census (search).** All 15 class-number-1 fields ℚ(√D), D ≤ 50; c₀ = (a+b√D)/e, |a|,|b|,e ≤ 16. A log-lattice filter proved **vacuous** for S = {2} (v₂(Nx)log2 = log|Nx| = Σ log|σᵢx|: image lies in a 2-plane) — replaced by exact certified generators + norm-guided, exactly verified exponent vectors (`cor7_pt.py`). Exactly **5 rank-2 points**:

| c₀ | field | α (rad) | α (deg) | dim V | true δ₂,₂ | Lie-class coords (in Li₄(S-unit) span, exact) |
|---|---|---|---|---|---|---|
| 1/√2 (square) | ℚ(√2) | 1.5707963267948966192 (π/2) | 90.0 | 2 | 0 ✓ | (0,20,0), (−26,0,155) |
| √2−1 | ℚ(√2) | 2.2874354808048409875 | 131.0603990 | 2 | 0 ✓ | (8,−2,−6), (−3,1,8) |
| (√2−1)² | ℚ(√2) | 2.7967406581640951061 | 160.2414361 | 2 | 0 ✓ | (84,−100,128), (−4,−16,48) |
| 1/φ (golden) | ℚ(√5) | 1.8091137886047627283 | 103.6545847 | 2 | 0 ✓ | (111/2,0,−12), (−51/2,0,12) |
| 1/φ³ = √5−2 | ℚ(√5) | 2.6649577299700610204 | 152.6908305 | 2 | 0 ✓ | (84,−21,6), (−18,3,0) |

Li₄ classes at the {2}-unit arguments (29 for ℚ(√2), 25 for ℚ(√5), exactly enumerated) span Q³, and **every field component of class(Ω) lies in their ℚ-span** (exact). Golden corner explicitly: class(Ω(1/φ)) = class of (7/6 − √5/18)Li₄(φ³) + (−4+4√5)Li₄(φ²/2) + (8−8√5)Li₄(φ/2).

**dim V = 3/4 controls (all exact):** 10 points in ℚ(√17), 2 in ℚ(√41) (all rank 3, δ₂,₂ ≠ 0, 3 components each); 1/φ² (rank 3, 2 comp), 1/(2φ) (rank 3, 3 comp), (√2−1)⁴ (rank 3), (√2−1)³ (rank 4), √3/2, rationals 1/5, 2/5, 3/7: δ₂,₂ ≠ 0 everywhere. **No accidental zeros.**

**Clean COR6 table** (deg = field degree; naive vs TRUE δ₂,₂; ncomp in parentheses): n=3 (2, dimV 3, ≠0/≠0, differ), n=4 (2, 2, 0/0), n=5 (4, 5, ≠0(36), differ), n=6 (1, 2, 0/0), n=7 (6, 7, ≠0(6), differ), n=8 (4, 4, ≠0(12), differ), n=9 (6, 7, ≠0(6), differ), n=10 (2, 3, ≠0(3), differ), n=11 (10, 11, ≠0(6), differ), n=12 (4, 4, ≠0(3), differ), α=π/6 (4, 4, ≠0(8), differ), 1/3 (1, 2, 0/0), 2/5 (1, 4, ≠0(15), differ), 3/7 (1, 4, ≠0(12), differ). Naive = true only at the three zero rows; differs at all 11 nonzero rows.

## Task 3: c₀ = √2−1 closed form (α = 2.2874354808048409875 rad)

**M₃ = 0.007189129721152490348077741879814010022297…** — 1095 digits by two methods: |quad(1100) − HPL(1150)| = 5.01·10⁻¹⁰⁹⁸ (quadrature 678 s). With b = √2−1, L = log(1+√2), λ = log2, coefficients in ℚ(√2):

- **W₂** = −½Li₂(b²) + π²/12 — exact (duplication Li₂(b)+Li₂(−b) = ½Li₂(b²); residual 0 @1120).
- **W₃** (12-element basis × {1,√2}; identical dps 700/1000; residual 8.05·10⁻¹¹¹⁸ @1120): W₃ = (−5/4 − 63√2/32)ζ(3) + (−11−13√2)Li₃(b) + (5+4√2)Li₃(−b) + (14+16√2)Li₃(1/√2) + (−4−4√2)L·Li₂(b) + (4+4√2)λ·Li₂(b) + (−11/12 − 49√2/48)π²L + (29/48 + 3√2/4)π²λ + (−7/24 − √2/3)λ³ + (−7/4 − 2√2)λ²L + (11/4 + 3√2)λL² + (−1/3 − √2/4)L³.
- **W₄** (thm:W4R basis: 7 Li₄ {b, −b, −b², 1/√2, 1/2, 1−1/√2, 2−√2} + 22 products, × {1,√2}; 59-dim lattice; identical dps 700/1000; residual 1.35·10⁻¹¹¹⁶ @1120; height 276480): W₄ = (22 − 31√2/2)Li₄(b) + (24 − 27√2/2)Li₄(−b) + (−5/2 + 13√2/8)Li₄(−b²) − √2·Li₄(1/√2) + (3/4)Li₄(1/2) + (−4+4√2)Li₄(1−1/√2) + (−2+3√2)Li₄(2−√2) − 3L·Li₃(b) + L·Li₃(−b) − 2√2·L·Li₃(1/√2) + (53/16 + 2√2)L·ζ(3) + 2λ·Li₃(b) + (−2+√2)λ·Li₃(−b) + (−2+2√2)λ·Li₃(1/√2) + (−25/16 − √2)λ·ζ(3) + (2+√2)Li₂(b)² + (2+√2/2)Li₂(b)L² + (−2+√2/2)Li₂(b)λL + (1/4−√2/4)Li₂(b)λ² + (−1/2−√2/8)π²Li₂(b) + (49/960 − 29√2/11520)π⁴ + (−5/8+√2/48)π²L² + (9/16−11√2/48)π²λL + (−3/16+5√2/32)π²λ² + (1−√2/48)L⁴ + (−5/4+√2/4)L³λ + (5/8−√2/4)L²λ² + (−7/24+13√2/48)Lλ³ + (13/192−13√2/384)λ⁴.
- **Lie class (PROVED):** Li₄ part's components (8,−2,−6) and (−3,1,8) = exact class(Ω(√2−1)). Total M₃ − (W₄+W₃+W₂) = 5.5·10⁻¹¹²⁰. **ζ(3) coefficient:** −5/4 − 63√2/32 (pure, basis {Li₃(b), Li₃(−b), Li₃(1/√2)}); ζ(3)-log terms (53/16 + 2√2)L·ζ(3) + (−25/16 − √2)λ·ζ(3).

## Task 4: verify_cor7.py (35/35 PASS, 7 s < 5 min)

Checks A1–A6 (c=1/3 closed forms at dps 420 vs stored 1145-digit values, quadrature at dps 60 vs stored, exact T4-class equality, ζ(3) = −26), B1–B4 (rational census to denom 300; exact rank + δ₂,₂ = 0 for the 5 quadratic corners; Li₄ span = Q³ and membership for all 5; dim-3/4 controls nonzero), C1–C6 (√2−1 closed forms, W₂ duplication, Li₄-part Lie class match, quadrature vs stored), D (recomputes COR6's degree ≤ 2 rows, comparing true/naive/differ with `cor6_table.json`). Run output: `TOTAL: 35/35 PASS, time 7s`.

## Evidence table (all numbers printed by code)

| Check | Label | Result |
|---|---|---|
| M₃(2 arccos(1/3)), 2 methods | NUMERICAL | 1095 digits; 5.01·10⁻¹⁰⁹⁸ |
| W₄ (c=1/3), 25-basis | NUMERICAL | identical dps 700/1000; 1.28·10⁻¹¹¹⁶ @1120 |
| W₃ (c=1/3), 11-basis | NUMERICAL | identical dps 600/900; 6.5·10⁻¹⁰⁴⁹ @1050 |
| W₂ (c=1/3) | PROVED from HPL formula | 9.2·10⁻¹¹²² |
| Li₄ part = COR6 T4 class | PROVED | coords (14, 20/3, 8/3) exact |
| rational census {1/3, 1/2} | PROVED (case analysis) | enumerated to denom 5000 |
| 5 quadratic corners: rank 2, δ₂,₂ = 0 | PROVED (exact certified) | all pass |
| Li₄(S-unit) span ⊇ class(Ω), 5 corners | PROVED (exact) | 29/25 args, all components |
| dim-3/4 controls (17 points) | PROVED | δ₂,₂ ≠ 0 everywhere |
| M₃(2 arccos(√2−1)), 2 methods | NUMERICAL | 1095 digits; 5.01·10⁻¹⁰⁹⁸ |
| W₄ (√2−1), thm:W4R basis | NUMERICAL | identical dps 700/1000; 1.35·10⁻¹¹¹⁶ @1120 |
| W₃ (√2−1) | NUMERICAL | identical dps 700/1000; 8.05·10⁻¹¹¹⁸ @1120 |
| W₂ (√2−1) duplication | PROVED | residual 0 @1120 |
| Li₄ part = class(Ω(√2−1)) | PROVED | (8,−2,−6), (−3,1,8) exact |
| quadrature cancellation bug (COR5's copy) | demonstrated, fixed | 1e-17 → 1e-298 at dps 300 |
| golden corner & (√2−1)² number-level closed forms | — | not attempted (symbol level only) |

**Key step files** (auto-saved in the workdir): cor7_lib + parse fix (step_30–33), formula validation vs COR5 (step_33, 35), quadrature diagnosis/fix (step_43–48), c13 values by two methods (step_49–50), c13 W₃/W₄ LLL (step_51–52), T4 class (step_53), assembly + Ψ⁽¹⁾ (step_54), rational census (step_56), LV-filter experiments and its vacuity (step_57–63), cor7_pt exact machinery (step_64–68), exact census → 5 corners (step_69), dim-3 controls (step_70–73), S-unit Li₄ args + Lie classes (step_74–75), thm:W4R basis from the paper (step_77–79), b21 values (step_80), angles + golden Li₄ combo (step_82–84), b21 W₃/W₄ LLL (step_85–87), b21 Lie class + assembly (step_88), b21 quadrature (step_89), relations stored (step_52, 90), verify_cor7.py written and run to 35/35 (step_91–96). Libraries: `cor7_lib.py`, `cor7_pt.py`; verify: `verify_cor7.py`; report: `result.md`.