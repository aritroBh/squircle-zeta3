# COR13C — PROVED: exact level-6 coloured-MZV certificate for COR7's closed form of M₃ at the corner c₀ = 1/3

**Verdict: PROVED (computer-assisted, exact), by the HEXC2 method.** The three weight-homogeneous identities ID₂, ID₃, ID₄ behind COR7's closed form of M₃(2 arccos(1/3)) are each an exact rational combination of proved relations among level-6 coloured MZV values, through Möbius maps onto sixth roots of unity. `verify_cor13c.py` runs in **28 s** (< 3 min): PART A (exact map enumeration + exact certificate, residual 0 on all words at all three weights) + PART B (62-digit numeric audit) both PASS; fresh re-run reproduced 27/27 PASS.

## Verdict table

| identity | statement (W-splits per COR5 convention, HPL fits from G1 theorem flat/M3angle) | verdict |
|---|---|---|
| **ID₂** | R₃⁽²⁾(1/3) + ζ(2)/2 = −Li₂(1/3) − Li₂(−1/3) + π²/12 | **PROVED** (exact; 28 proved rows; residual 0 on all 36 weight-2 words) |
| **ID₃** | R₃⁽³⁾(1/3) + ζ(2)Ψ⁽¹⁾(1/3) − ζ(3)/2 = −26ζ(3) − 87Li₃(1/3) − 130Li₃(−1/2) + (21/4)π²log2 − (43/12)π²log3 + 67Li₂(1/3)log2 − 45Li₂(1/3)log3 + (65/3)log³2 − 65log²2log3 + (109/2)log2log²3 − (37/3)log³3 | **PROVED** (229 rows; residual 0 on all 252 weight-3 words) |
| **ID₄** | R₃⁽⁴⁾(1/3) + ζ(2)Ψ⁽²⁾(1/3) + (2/3)ζ(3)log2 = −(22/3)Li₄(−2) − (17/2)Li₄(1/3) + (20/3)Li₄(2/3) + (110/9)Li₄(1/2) + (25/48)Li₄(1/9) − 12Li₃(1/3)log2 + (28/3)Li₃(1/3)log3 − (16/3)Li₃(−1/2)log2 + (16/3)Li₃(−1/2)log3 + (14/3)ζ(3)log2 − (25/9)ζ(3)log3 − 4Li₂(1/3)² + (13/3)Li₂(1/3)log²2 − 6Li₂(1/3)log2log3 + (2/3)Li₂(1/3)log²3 + (11/18)π²Li₂(1/3) − (31/27)π²log²2 + (31/72)π²log²3 − (271/1296)π⁴ + (46/27)log⁴2 − (62/9)log³2log3 + (61/6)log²2log²3 − (17/3)log2log³3 + (103/144)log⁴3 | **PROVED** (1695 rows; residual 0 on all 1764 weight-4 words) |

Bonus (as in HEXC2): **ζ(2) = π²/6 is itself exact** in the weight-2 span (`cert13_euler.json`, 28 rows, residual 0), and **π⁴ comes out exactly** (T₄'s coefficient of (1,1,1,1) is 2395/6, killed by the certificate — no Borel/denominator argument).

## Step 1 — Möbius maps: they exist. Count per level (EXACT over ℚ(ζ_N))

Points needed: the HPL words use letters {0,±1}, base point 0, endpoint 1/3, and ∞ enters the pullback formula — so the configuration is exactly **{0, 1/3, 1, −1, ∞}**. The Li-atom arguments {1/3, −1/3, 2/3, 1/2, −1/2, −2} are **all cross-ratios** of these five points (the cross-ratio set is the three S₃-orbits {−2, 3, −1/2, 3/2, 1/3, 2/3}, {−3, 4, −1/3, 4/3, 1/4, 3/4}, {−1, 2, 1/2}); Li₄(1/9) enters only via the input theorem Li₄(1/9) = 8(Li₄(1/3) + Li₄(−1/3)). No sixth point (in particular not −1/3 as a *point*) is needed.

| level N | injective maps {0,1/3,1,−1,∞} → μ_N | extra maps into μ_N ∪ {0,∞} |
|---|---|---|
| 6 | **12** | 0 |
| 12 | **24** | 0 |
| 18 | **36** | 0 |
| 24 | **48** | 0 |
| 30 | **60** | 0 |
| 36 | **72** | 0 |

Exactly **2N** maps at each level, enumerated EXACTLY over ℚ(ζ_N) two independent ways: (i) direct Fraction arithmetic over ℚ(ζ₆) (HEXC2 step_18 style, N=6: 12 maps, matches stored `mobius_maps_13.json`); (ii) mod-p filter with sharp primes ((p−1)/N ≥ 50, a complete superset of candidates) + exact verification in ℚ(ζ_N) (`cyclo.py`). Structure (exact): all maps are shifts (multiplication of all images by ζᵏ) of two conjugate base maps with exponent tuples (for (0, 1/3, 1, −1, ∞)) **(0, N/6, N/3, 2N/3, N/2)** and **(0, 5N/6, 2N/3, N/3, N/2)**. At N=6 each map omits exactly one root, 2 maps per omitted root — structurally identical to the hexagon corner. Task map used: (0,1,2,4,3), i.e. 0↦1, 1/3↦ζ, 1↦ζ², −1↦ζ⁴, ∞↦ζ³=−1.

## Step 2 — Conversion (validated numerically ≥ 60 digits, two independent evaluators)

HPL letters pull back as 0 ↦ (3)−(0), 1 ↦ (4)−(3), −1 ↦ (3)−(2); iterated integrals from 1 to ζ become level-6 CMZV words by shuffle-regularized path composition through 0 (HEXC2 `path6`, first letter outermost); Li atoms via cross-ratio frames; atoms log2 = −(3), log3 = −(2)−(4), iπ = 3·(1), ζ(2) = (−1,0), ζ(3) = (−1,−1,0). All 38 HPL words at 1/3, all 9 Li atoms, 7 elementary atoms validated at 70 dps (mpmath evaluator) and 62 dps (acb_poly evaluator); targets T₂/T₃/T₄ (supports 26/170/1083) evaluate to 2.65e-66 / 4.38e-64 / 5.81e-65.

## Step 3 — Certificate (exact)

Targets T₂/T₃/T₄ lie in the span of cmzv6's proved relations (regularized double shuffle [Ihara–Kaneko–Zagier], Zhao distribution d ∈ {2,3,6}, weight-1 branch rows, lifts): membership first checked mod p = 1000000009 and 1000000097, then solved exactly over ℚ and verified exactly with FLINT fmpq_mat on **every** admissible word. Quotient dims are 3, 8, 23, 69 (same as HEXC2; the 2 and 14 relations missing vs Deligne's 3,8,21,55 are not needed). **Certificate sizes:** weight 2: 28 rows, 25 nonzero multipliers, height ≤ 2 digits; weight 3: 229 rows, height ≤ 10 digits; weight 4: 1695 rows, 1692 nonzero multipliers, height ≤ 141 digits (median 139).

## Evidence table (all numbers printed by my code)

| check | label | residual / digits |
|---|---|---|
| COR7's ID₂/ID₃/ID₄ vs stored 1150-digit HPL values | NUMERICAL (input) | 7.64e-152 / 5.5e-150 / 5.5e-150 @150dps |
| map counts N=6,12,18,24,30,36 | **PROVED** (exact over ℚ(ζ_N)) | 12, 24, 36, 48, 60, 72 = 2N; 0 maps need 0/∞ |
| 12 Möbius maps, images on {1/3, −1} | PROVED exact; numeric | 6.73e-71 @70dps; 3.2e-66 @62dps |
| 38 HPL-word conversions at 1/3 | validated | worst 1.23e-70 @70dps; 2.77e-66 @62dps |
| 9 Li-atom conversions | validated | worst 1.86e-70 @70dps (Li₄(−2)); 2.14e-65 @62dps |
| 7 elementary atoms | validated | ≤ 1.31e-68 @70dps; ≤ 3.04e-64 @62dps |
| T₂, T₃, T₄ numeric values | validated | 2.65e-66, 4.38e-64, 5.81e-65 @62dps |
| **ID₂ / ID₃ / ID₄ exact certificates** | **PROVED** | residual 0 on 36/252/1764 words; 28/229/1695 proved rows |
| ζ(2)=π²/6 in weight-2 span | **PROVED** | residual 0, exact |
| π⁴ exact | **PROVED** | T₄ coeff of (1,1,1,1) = 2395/6 → killed exactly |
| quotient dims 3,8,23,69 | checked | mod p ∈ {1000000009, 1000000097} |
| two independent CMZV evaluators (mpmath vs acb_poly) | cross-checked | 1.9e-63 (12 random words) |
| Li₄(1/9) distribution (input theorem) | validated | 1.13e-71 |
| `verify_cor13c.py` | fresh run | **28 s, 27 PASS, 0 FAIL** |

## Files (new_results/COR13C_third_corner_certificate)

- `verify_cor13c.py` (28 s < 3 min), `cert13_w23.json` (9,868 B), `cert13_w4.json` (473,712 B), `cert13_euler.json`, `mobius_maps_13.json`, `map_counts.json`, `progress.log`
- `cor13c_common.py` (maps/pullbacks/cross-ratio table/li_atom/atoms/targets), `cyclo.py` (exact ℚ(ζ_N)), `exact_linalg.py`, `fast_cmzv.py` + `numeric_cmzv.py` (two independent evaluators), `w4_pivots_13.pkl`, `w4_y_13.pkl` (checkpoints)
- Producing steps: `step_21.py` (12 maps exact over ℚ(ζ₆)), `step_22.py` (exact counts 2N for N ≤ 36, mod-p filter + exact verification), `step_23.py` (conversion module; target supports 26/170/1083), `step_24.py` (70-dps audits, 16 min), `step_25.py` (targets ≈ 0, evaluator cross-check), `step_26.py`/`step_27.py` (membership mod p, two primes at weight 4), `step_28.py` (exact certs w2/w3), `step_29.py`/`step_30.py` (weight-4 solve + exact verification on 1764 columns), `step_31.py` (Euler cert + π⁴ coefficient 2395/6), `step_32.py` (verifier), `step_34.py` (report data), `step_35.py` (fresh run 28 s, 27/27 PASS).

**Caveat (identical to HEXC2's):** the conversion code (Möbius pullback, cross-ratio frames, path composition) is validated numerically to ≥ 60 digits on every piece but not formally verified; the final identities themselves are checked with exact rational arithmetic only, against relation rows that are proved theorems (RDS, Zhao distribution, branch, lifts).