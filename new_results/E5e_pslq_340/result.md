# E5e — non-local ζ(5) of the b = 0 hexagon arcs: **RELATION FOUND AND CERTIFIED**

## Verdict
`X = P^(5)_arc/√5` lies in the ℚ-span of the 20-element weight ≤ 5 alternating-MZV basis with a **primitive integer relation, maxcoeff = 2,574,432,000**, returned **bit-identical by PSLQ at 250, 280, 310 (and 340) digits**, residual at 340 digits `5.69e-330 ≤ maxcoeff·10^(−330) = 2.57e-321` (margin ~10^8), and recovered **exactly** (as ± the vector) by an independent flint LLL. This resolves E5d's failure: the true maxcoeff 2.57×10^9 sat between E5d's excluded bound (≳2×10^8) and its 187-digit detection ceiling (~8×10^8).

## Data check
- `b0_p5_convB_350.json`: e1 (r=2), e2 (r=1/2), 337 decimals each; **agree through all first 336 decimals** (`e1−e2 = 4.0e-337` abs, 4.69e-340 rel), consistent with E5d's conv log (9.49e-339 / 6.43e-339); matches E5d's 187-digit prefix. ⇒ ~336 reliable digits.
- `X = −381.140137469830729527176573357265848015393511891093057029…`; basis recomputed at 360 dps from `../E5c_nonlocal_z5/mzv.py` (20 elements) → `consts360.json`.

## Evidence table (all numbers printed by code)

| step file | run | precision/tol | result | residual | maxcoeff |
|---|---|---|---|---|---|
| step_08.py | PSLQ [X,20-basis] | 250 d / 1e-225 | **R** | −8.7873e-244 | 2,574,432,000 |
| step_08.py | PSLQ | 280 d / 1e-255 | **same R** | −1.6175e-272 | 2,574,432,000 |
| step_08.py | PSLQ | 310 d / 1e-285 | **same R** | +2.5154e-302 | 2,574,432,000 |
| step_10.py | PSLQ | 340 d / 1e-315 | **same R** | +5.6841e-330 | 2,574,432,000 |
| step_09.py | 340-dps acceptance test | e1 data | PASS | 5.6841e-330 ≤ 2.574432e-321 | — |
| step_09.py | 340-dps, e2 leg (independent) | e2 data | PASS | 4.2853e-330 | — |
| step_10.py | **exact rational** residual (350-digit truncated decimals) | — | PASS | 5.6886e-330 | — |
| step_11.py | flint LLL δ=0.99, scale 10^300 | — | **row 0 = −R exactly** (norm 3.136e9; rows 1–5: 5.3e14–6.6e14, Dirichlet scale 1.9e14) | last coord 5.46e8 ≈ rounding | — |

`R` (for `[X, 1, log2, π², log2², ζ3, π²log2, log2³, π⁴, ζ3log2, π²log2², log2⁴, Li4(1/2), ζ5, ζ2ζ3, π⁴log2, ζ3log2², π²log2³, log2⁵, Li5(1/2), Li4(1/2)log2]`), gcd = 1:
`[7680000, 0, 1175040000, 49152000, 2574432000, −501811200, −111744000, 927640800, 6229440, 397036800, 45175200, 287114250, 381542400, 77518080, 52715520, −1240896, −1751040, −10005440, −2977497, −24330240, −56401920]`

Residuals scale like `maxcoeff·10^(−dps)` and sit ~10^19 below tolerance — true-relation signature (contrast E5d's spurious tolerance-level vectors).

## Exact closed form (re-checked at 340 digits: `X − closed = 7.37e-337`, `P5 − √5·closed = 1.66e-336` — at the data's own accuracy)

`P^(5)_arc = √5·X`,

```
X = −153·log2 − (32/5)·π² − (26817/80)·log2² + (3267/50)·ζ(3) + (291/20)·π²log2 − (386517/3200)·log2³
    − (6489/8000)·π⁴ − (20679/400)·ζ(3)log2 − (18823/3200)·π²log2² − (382819/10240)·log2⁴ − (1242/25)·Li4(1/2)
    − (20187/2000)·ζ(5) − (858/125)·ζ(2)ζ(3) + (6463/40000)·π⁴log2 + (57/250)·ζ(3)log2²
    + (31267/24000)·π²log2³ + (992499/2560000)·log2⁵ + (396/125)·Li5(1/2) + (918/125)·Li4(1/2)log2
```
(common denominator 7680000 = 2^12·3·5^4).

## ζ(5), ζ(2)ζ(3) coordinates in THIS basis vs local law

| constant | coord in X | local coord | non-local part |
|---|---|---|---|
| ζ(5) | **−20187/2000** = −10.0935 | −2079/250 = −8.316 | **N^(5) = −711√5/400 ≈ −3.9746108300058761854** (for P; −711/400 for X) |
| ζ(2)ζ(3) | **−858/125** = −6.864 | −477/250 = −3.816 | **N^(2,3) = −1239√5/250 ≈ −11.081952896488957735** (for P; −1239/250 for X) |

Both defects are nonzero and large. **Basis-dependence:** Li5(1/2) (+396/125), Li4(1/2) (−1242/25), Li4(1/2)log2 (+918/125) all occur with nonzero coefficients, so the ζ(5)/ζ(2)ζ(3) coordinates — and hence N^(5), N^(2,3) — change under a change of weight-5 basis (e.g. eliminating Li5(1/2) by its standard identity); they are coordinates in this basis, not invariants.

## Consistency with E5d
E5d proved maxcoeff ≳ 2×10^8 and a detection ceiling ~8×10^8 at 187 digits; the true maxcoeff 2.57×10^9 lies exactly in that gap, and the relation already locks at 250 digits, as E5d's scaling estimate predicted.

## Files (work dir `E5e_pslq_340/`)
`step_05.py` (data/agreement check) · `step_07.py` (360-dps constants → `consts360.json`) · `step_08.py` (PSLQ 250/280/310 → `pslq_runs.json`) · `step_09.py` (340-dps residuals e1/e2, gcd, coordinates → `relation.json`) · `step_10.py` (exact-rational residual; PSLQ @340) · `step_11.py` (flint LLL cross-check) · `step_12.py` (closed-form 340-digit re-check; N^(5), N^(2,3)) · `result.md` (this report).