# HEXC2 — PROVED: exact level-6 certificate for the hexagon-corner closed form of M₃(2π/3)

**Verdict: PROVED (computer-assisted, exact).** The three weight-homogeneous identities ID₂, ID₃, ID₄ behind COR5's closed form of M₃(2π/3) are each an exact rational combination of proved relations among level-6 coloured MZV values, through the Möbius map onto sixth roots of unity. `verify_hexc2.py` runs in **22 s** (< 5 min): PART A (exact certificate, residual 0 on all words at all three weights) + PART B (62-digit numeric audit) both PASS.

## The statements proved

- **ID₂ (wt 2):** −H₀,₁(½) − H₀,₋₁(½) + ζ(2)/2 = log²2/2 + Li₂(1/3) + ½log²(3/2)
- **ID₃ (wt 3):** R3_3 + ζ(2)Ψ1 − ζ(3)/2 = −(15/4)ζ(3) − 3Li₃(1/3) − 8Li₃(−½) + products
- **ID₄ (wt 4):** F4 = −21Li₄(−2) − (65/8)Li₄(1/3) + (13/2)Li₄(2/3) + (137/6)Li₄(1/2) + (39/64)Li₄(1/9) + products

## Key facts and method

1. **The key fact (correcting HEXC1):** exactly **12** Möbius maps send {0,½,1,−1,∞} into μ₆ (enumerated exactly over ℚ(ζ), step_18; each map omits exactly one root, 2 maps per root). The task map f(z) = (ζ²z+1)/(ζ⁴z+1) sends 0↦1, ½↦ζ, 1↦ζ², −1↦ζ⁵, ∞↦ζ⁴.
2. **Conversion:** HPL letters pull back as 0↦(2)−(0), 1↦(4)−(2), −1↦(2)−(1); iterated integrals from 1 to ζ become level-6 CMZV words by path composition through 0 (with the correct first-letter-outermost composition order — an initial sign/order bug was caught and fixed via numeric checks). Closed-form arguments {1/3, −1/3, 2/3, ½, −½, **−2**} are all cross-ratios of the configuration (step_19); only 1/9 is not — covered by the one-line distribution identity Li₄(1/9) = 8(Li₄(1/3)+Li₄(−1/3)).
3. **Certificate:** targets T₂/T₃/T₄ (supports 26/177/920) lie in the span of cmzv6's proved relations (RDS, Zhao distribution d∈{2,3,6}, weight-1 branch, lifts) at all weights. **No frame double-shuffle rows were needed** — the 2 and 14 relations missing vs Deligne's 3,8,21,55 (quotients 3,8,23,69 confirmed mod p) play no role, exactly as at the square. Multipliers from exact solves (weight 4: 1695×1695, heights ≤ 142 digits), verified exactly on all words with FLINT fmpq_mat. **π⁴ comes out exactly** (coefficient of (1,1,1,1) in T₄ = 5148/5 → 0 after reduction; no Borel/denominator argument). Bonus: ζ(2) = π²/6 is itself exact in the weight-2 span (`cert_euler.json`).

## Evidence table (all numbers printed by my code)

| check | label | residual / digits |
|---|---|---|
| 12 Möbius maps, task-map images | PROVED (exact over ℚ(ζ)) | 3.75e-66 numeric |
| 38 HPL-word conversions | validated | worst 1.43e-70 @70dps; 1.19e-65 @62dps |
| 9 Li-atom conversions | validated | worst 2.20e-70 @70dps; 1.12e-65 @62dps |
| 7 elementary atoms | validated | ≤ 3.04e-64 |
| T₂, T₃, T₄ numeric values | validated | 1.34e-66, 2.71e-65, 2.32e-64 |
| ID₂ / ID₃ / ID₄ exact certificates | **PROVED** | residual 0 on 36/252/1764 words; 28/229/1695 proved rows |
| ζ(2)=π²/6 in weight-2 span | **PROVED** | residual 0, exact |
| quotient dims 3,8,23,69 vs Deligne 3,8,21,55 | checked | mod p = 1000000009 |
| two independent CMZV evaluators (mpmath vs acb_poly) | cross-checked | 1.95e-63 |

## Files (new_results/HEXC2_level6_certificate)

- `result.md` (full write-up), `verify_hexc2.py` (22 s), `cert_w23.json`, `cert_w4.json`, `cert_euler.json`, `mobius_maps.json`, `progress.log`
- `hexc2_common.py` (maps/pullbacks/path6/hpl_to_cmzv/li_atom/atoms/targets), `fast_cmzv.py` + `numeric_cmzv.py` (two independent evaluators), `exact_linalg.py`
- Producing steps: `step_18.py` (12 maps over ℚ(ζ)), `step_19.py` (cross-ratio table), `step_21…33` (conversion module + letter/Li sign fixes), `step_34…36` (mpmath evaluator + composition-order fix), `step_37.py`/`step_39.py` (70-dps audits: 38 HPL words worst 1.43e-70; 9 Li atoms worst 2.20e-70), `step_40.py` (targets), `step_43.py`/`step_45.py` (lattice, membership mod p), `step_53.py` (exact certs w2/w3), `step_55.py`/`step_56.py` (weight-4 solve + exact verification on 1764 columns), `step_57.py` (ζ(2)=π²/6, heights 142 digits, cert_w4.json), `step_64…67` (fast evaluator + audits), `step_71/72` (verifier + fresh run), `step_78` (Euler exact cert), `step_80` (final fresh PASS 22 s).

**Caveat (same as the square's certificate):** the exact conversion code (Möbius pullback, cross-ratio frames, path composition) is validated numerically to ≥ 60 digits on every conversion but not formally verified; the final identity itself is checked with exact rational arithmetic only.