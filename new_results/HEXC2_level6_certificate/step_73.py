result = '''# HEXC2 — exact level-6 certificate for the hexagon-corner closed form of M₃(2π/3)

**Status: PROVED (computer-assisted, exact).** `verify_hexc2.py` runs in ~35 s: **PART A PASS + PART B PASS**.

The three weight-homogeneous identities behind the closed form of M₃(2π/3),
- **ID₂ (weight 2):** −H₀₁ − H₀₋₁ + ζ(2)/2 = log²2/2 + Li₂(1/3) + ½log²(3/2)   (HPL words at ½)
- **ID₃ (weight 3):** R3_3 + ζ(2)Ψ1 − ζ(3)/2 = −(15/4)ζ(3) − 3Li₃(1/3) − 8Li₃(−½) + products
- **ID₄ (weight 4):** F4 = −21Li₄(−2) − (65/8)Li₄(1/3) + (13/2)Li₄(2/3) + (137/6)Li₄(1/2) + (39/64)Li₄(1/9) + products

are each an **exact rational combination of proved relations among level-6 coloured MZV values**.
This is the hexagon analogue of the square's Theorem thm:W4R certificate (A3h5): where the square used the
Cayley map onto 8th roots of unity, the hexagon uses a Möbius map onto 6th roots of unity.

## 1. The key fact (correcting HEXC1)

HEXC1 claimed no Möbius map sends {0, ½, 1, −1, ∞} into roots of unity. **Wrong.** Exact enumeration over
ℚ(ζ), ζ = e^{iπ/3} (`step_02.py`): there are exactly **12** such maps; all land in μ₆, in fact on μ₆∖{−1}.
The task map is the first:

> f(z) = (ζ²z + 1)/(ζ⁴z + 1):  0 ↦ 1, ½ ↦ ζ, 1 ↦ ζ², −1 ↦ ζ⁵, ∞ ↦ ζ⁴.   (exponents (0,1,2,5,4))

The 12 = 6 rotations × 2 orientations (`mobius_maps.json`, exact over ℚ(ζ), verified numerically in PART B).

## 2. The conversion (Task 1)

Under g = f⁻¹ the HPL letters {0, ±1} (dz/z, dz/(1−z), dz/(1+z)) pull back to dlog forms with singular
points f(0), f(±1), f(∞) ∈ μ₆. With (j) = ζʲdu/(1−ζʲu) = cmzv6 letter j (and dlog(u−ζᵏ) = −letter(−k)):

- 0 ↦ (2) − (0)      [dlog(u−1) − dlog(1−ζ²u)]
- 1 ↦ (4) − (2)      [−dlog(u−ζ²) + dlog(1−ζ²u)... = −dlog(g−1)]
- −1 ↦ (2) − (1)     [dlog(g+1)]

So every H_w(½) becomes an iterated integral from 1 to ζ along the image arc (a unit-circle arc; the
sector between it and the radial path 1→0→ζ contains none of the singular points {1,ζ²,ζ⁴,ζ⁵}, so the
two are homotopic — no branch ambiguity). Path composition through 0 and rescaling u = ζs,

  I(w; ζᶜ→ζᵉ) = Σₖ (−1)ⁿ⁻ᵏ reg[I(w[:k] shifted+e; 0→1)] ш reg[I(rev(w[k:]) shifted+c; 0→1)]

(first letter = outermost; composition puts the first part of the word on the FAR side), writes each as an
exact ℚ-combination of level-6 CMZV words. All words are convergent (no fit word ends in 0), so the T=0
regularization of the split factors is consistent and the composition theorem applies.

The closed-form side enters through cross-ratio frames: every needed argument is a cross-ratio of
{0, ½, 1, −1, ∞}: 1/3 (6 reps), −1/3 (2), 2/3 (6), ½ (8), −½ (6), **−2 (6)** — no inversion needed.
Li_n(z) = −G(0ⁿ⁻¹,1; z) = −∫_{p0}^{e}[dlog(x−p0)−dlog(x−pinf)]ⁿ⁻¹[dlog(x−p1)−dlog(x−pinf)], pulled back
letterwise. Only 1/9 is not a cross-ratio: Li₄(1/9) = 8(Li₄(1/3) + Li₄(−1/3)) is used as an input theorem
(absolute convergence: Li₄(x)+Li₄(−x) = 2Σx²ᵏ/(2k)⁴ = Li₄(x²)/8). Elementary atoms: log2 = −(3),
log3 = −(2)−(4), iπ = 3(1), ζ(2) = (−1,0), ζ(3) = (−1,−1,0), π² = −(iπ)², π⁴ = (iπ)⁴ = 1944·(1,1,1,1).

**Validation (PART B, 62 dps, independent of the proof):** all 38 distinct HPL-word conversions vs
`H_series` — worst residual 1.19e-65; all 9 Li atoms vs `mpmath.polylog` — worst 1.12e-65; 7 elementary
atoms; the Möbius images. Two independently written CMZV evaluators (mpmath log-series vs FLINT acb_poly)
agree to 1.95e-63 on random words. Branch and regularisation choices are validated by these checks.

## 3. The target vectors (Task 2)

Exact ℚ-vectors on admissible level-6 words (HPL side by §2; closed form by atoms + shuffle products;
c-e powers are rationals):

| target | definition | support | numeric value (62 dps) |
|---|---|---|---|
| T₂ | R3_2 + ζ(2)/2 − W2_closed | 26 | 1.34e-66 |
| T₃ | R3_3 + ζ(2)Ψ1 − ζ(3)/2 − W3_closed | 177 | 2.71e-65 |
| T₄ | F4 − F4_closed | 920 | 2.32e-64 |

## 4. The certificate (Task 3)

Relations (all proved theorems), built by HEXC1's `cmzv6.make_relations`: regularized double shuffle
(Ihara–Kaneko–Zagier), Zhao's regularized distribution d ∈ {2,3,6}, weight-1 branch identities
(1)+(5) = 0, (1)−(5)−2(2)+2(4) = 0, and shuffle lifts of lower-weight relations.

| weight | words | relation rows | rank (mod p) | quotient | Deligne N=6 | missing |
|---|---|---|---|---|---|---|
| 1 | 5 | 4 | 2 | 3 | 3 | 0 |
| 2 | 36 | 57 | 28 | 8 | 8 | 0 |
| 3 | 252 | 682 | 229 | 23 | 21 | 2 |
| 4 | 1764 | 7869 | 1695 | 69 | 55 | 14 |

(Deligne count checked: Hilbert series (1−x)⁻¹·(1−x)/(1−3x+x²) = 1/(1−3x+x²) = 1+3x+8x²+21x³+55x⁴+…
from 2 weight-1 generators, one generator in each weight ≥ 2, and powers of 2πi.)

For each weight: independent proved rows selected mod p (transpose RREF), pivot columns mod p, exact
rational solve y·S = T on the pivot block (FLINT `fmpz_mat.solve`), then **exact rational verification
y·S = T on ALL words** (FLINT `fmpq_mat`). Certificate files: `cert_w23.json`, `cert_w4.json`
(row indices prow + rational multipliers y).

| weight | certificate rows | target support | multiplier height | exact residual |
|---|---|---|---|---|
| 2 | 28 | 26 | ≤ 2 digits | 0 on all 36 words |
| 3 | 229 | 177 | ≤ 9 digits | 0 on all 252 words |
| 4 | 1695 | 920 | ≤ 142 digits (1692 nonzero) | 0 on all 1764 words |

**No frame double-shuffle rows were needed**: the targets lie in the span of the standard relations at all
three weights — as at the square ("M₃ needed none of the missing ones"), the 2 and 14 missing relations at
weights 3 and 4 play no role. As a self-containedness check, ζ(2) = π²/6 itself lies in the weight-2 span
(vector (−1,0)+3(1,1), checked mod p), so even Euler's identity comes out of the level-6 relations.

**π⁴ comes out exactly.** The closed form contributes −(949/2160)π⁴ (as 1944·(1,1,1,1), the shuffle
normalisation of (iπ)⁴ = 1944·val(1,1,1,1) = π⁴); the coefficient of (1,1,1,1) in T₄ is exactly
5148/5 (closed-form part alone: −8541/10), and after the exact reduction it is 0 — as is every one of
the 1764 word coefficients. No Borel argument and no denominator bound is needed anywhere.

## 5. What the proof uses

1. The G1 HPL fits (the definitions of R3_3, F4, Ψ1 at c = ½) and COR5's closed forms — inputs.
2. Theorems: regularised double shuffle and Zhao's regularised distribution for level-6 CMZV values;
   shuffle products of convergent regularised iterated integrals; path composition of regularised
   iterated integrals; substitution under the Möbius map; the one-line Li₄ distribution identity at x = 1/3.
3. Correctness of the exact conversion code (Möbius pullback, cross-ratio frames, path composition):
   audited and validated numerically at 62 dps (worst residual 1.19e-65 over all conversions), but not
   formally verified — same caveat as the square's certificate.
4. Exact rational arithmetic for the final identity (FLINT `fmpq_mat`); no numerics enter the proof.

## 6. Evidence table

| check | label | residual / digits |
|---|---|---|
| 12 Möbius maps {0,½,1,−1,∞} → μ₆ | PROVED (exact over ℚ(ζ)) | 12 = 6×2 enumerated exactly |
| task map images 0↦1, ½↦ζ, 1↦ζ², −1↦ζ⁵, ∞↦ζ⁴ | PROVED | 3.75e-66 numeric |
| 38 HPL-word conversions (62 dps) | validated | worst 1.19e-65 |
| 9 Li-atom conversions (62 dps) | validated | worst 1.12e-65 |
| 7 elementary atoms (log2, log3, iπ, ζ2, ζ3, π², π⁴) | validated | ≤ 3.04e-64 |
| T₂, T₃, T₄ numeric values | validated | 1.34e-66, 2.71e-65, 2.32e-64 |
| ζ(2) = π²/6 in weight-2 span | PROVED (in span) | mod-p membership |
| ID₂ exact certificate | **PROVED** | residual 0 on 36 words, 28 rows |
| ID₃ exact certificate | **PROVED** | residual 0 on 252 words, 229 rows |
| ID₄ exact certificate | **PROVED** | residual 0 on 1764 words, 1695 rows |
| quotient dims 3, 8, 23, 69 vs Deligne 3, 8, 21, 55 | checked | mod p = 1000000009 |
| two independent CMZV evaluators | cross-checked | 1.95e-63 |
| Li₄(1/9) distribution input | validated | 1.19e-66 |

## 7. Files

- `hexc2_common.py` — maps, letter pullbacks, path composition `path6`, conversions `hpl_to_cmzv`,
  `li_atom` (cross-ratio frames), elementary atoms, target builders `target2/3/4`.
- `fast_cmzv.py` — FLINT acb_poly CMZV evaluator (audit); `numeric_cmzv.py` — independent mpmath evaluator.
- `exact_linalg.py` — independent-row selection, pivot columns, exact solve, exact verify.
- `cert_w23.json`, `cert_w4.json` — the certificates (proved-row indices + rational multipliers).
- `verify_hexc2.py` — < 5 min: PART A exact certificate check, PART B 62-digit numeric audit.
- `mobius_maps.json` — the 12 maps (exponent images), enumerated exactly over ℚ(ζ).
- step files: `step_02` (map enumeration), `step_10` (cross-ratio table), `step_20…32` (conversion code +
  letter/Li validations), `step_37…40` (evaluator + 70-dps audits), `step_42…44` (lattice + membership mod p),
  `step_45…54` (exact certificates w2/w3/w4), `step_57` (ζ(2)=π²/6 + heights), `step_60…66` (fast evaluator),
  `step_69…` (verifier + result).

**Bottom line:** ID₂, ID₃ and ID₄ — hence the hexagon-corner closed form of M₃(2π/3) — are PROVED,
by an exact level-6 CMZV certificate through the Möbius map onto sixth roots of unity, with the π⁴
coefficient coming out exactly and no use of the missing 2/14 relations or of any denominator bound.
'''
open('result.md','w').write(result)
with open('progress.log','a') as f:
    f.write('step_71: result.md written; verify_hexc2.py PASS 34s (PART A exact 3/3, PART B numeric audit)\n')
print('result.md written', len(result))
