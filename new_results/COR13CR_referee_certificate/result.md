# COR13CR — adversarial referee report on COR13C (exact level-6 certificate for the closed form of M₃ at c₀ = 1/3)

**Verdict: PROVED (computer-assisted, exact).** Every claim in COR13C is confirmed end-to-end with fully independent code: the three weight-homogeneous identities ID₂, ID₃, ID₄ behind COR7's closed form of M₃(2 arccos(1/3)) are exact rational combinations of proved relations among level-6 coloured-MZV values (regularized double shuffle, Zhao distribution d ∈ {2,3,6}, weight-1 branch identities, shuffle lifts), through the Möbius map 0 ↦ 1, 1/3 ↦ ζ, 1 ↦ ζ², −1 ↦ ζ⁴, ∞ ↦ ζ³ = −1 (ζ = e^{iπ/3}); residual **exactly 0 on every word** (36/252/1764 words at weights 2/3/4), certificate sizes 28/229/1695 rows, multiplier heights ≤ 141 digits; no relations beyond the standard proved set enter (quotient dims 3, 8, 23, 69 vs Deligne's 3, 8, 21, 55). **No errors and no display typos were found** — every displayed formula in COR13C's result.md matches the certified JSON data exactly. The referee adds two strengthenings (as HEXC2R did for HEXC2): the 36 pulled-back 1-form identities are verified **exactly over ℚ(ζ)** (upgrading the conversion's analytic input from numeric to exact), and frame independence of the targets is **proved** (frame differences lie in the proved relation span). The remaining numeric component is the conversion's path-composition/regularization machinery, audited to ≥ 65 digits on every piece with an independent evaluator, with a demonstrated discriminating wrong-order test.

## What was checked, task by task

### Task 1 — certificate algebra: CONFIRMED (exact, three backends + numeric row validity)
- My relation generators (HEXC2R's reflib, written from the definitions) reproduce cmzv6's proved row sets **exactly as sets** at weights 1–4: 4/57/682/7869 primitive rows [step_38].
- Rebuilt rows and targets, re-ran the certificates: Σ yᵢ·row(prowᵢ) − T = **0 on every word** at weights 2, 3, 4 and for the Euler bonus, with (a) FLINT fmpq_mat (all), (b) an independent pure-Fraction backend (w2, w3, Euler), (c) my own modular arithmetic mod 5 primes (w4). Target supports 26/170/1083 as stated; T₄ coefficient of (1,1,1,1) = 2395/6 (π⁴ comes out exactly) [step_39].
- Multiplier statistics reproduced exactly: w2: 25/28 nonzero, height ≤ 2 digits; w3: 229/229 nonzero, height ≤ 10 digits; w4: 1692/1695 nonzero, height ≤ 141 digits, median 139; Euler: 22/28 nonzero, height ≤ 2 [step_39, step_40].
- **Row validity**: all 1980 certificate rows (28+229+1695, plus the 28 Euler rows from the same weight-2 pool) plus stratified unused samples vanish at 65 dps with my own truncated-iterated-sum evaluator: worst 3.92e-65 (w2), 1.9e-64 (w3), 1.37e-63 (w4), classified by family — w2: 20 DS + 5 dist₂ + 1 dist₃ + 2 lift₁; w3: 203 DS + 11 dist₂ + 2 dist₃ + 13 lift₁; w4: 1612 DS + 26 dist₂ + 5 dist₃ + 50 lift₁ + 2 lift₂; **no frame rows** [step_43, step_59]. Evaluator cross-validated against COR13C's two independent evaluators: worst 9.3e-63 / 7.29e-63 on 150 random words [step_41].
- **Tamper test**: y₀ → y₀+1 gives 2 bad words; y₁₀₀ → 2y₁₀₀ gives 8 bad words — the verification genuinely discriminates [step_56].

### Task 2 — conversion: CONFIRMED (exact equality of rational vectors; numeric to 70 digits; analytic input now exact)
- My own conversion (my 12-map enumeration, my pullback, reflib's Chen path composition through 0 + shuffle regularization, my cross-ratio frames from {0, 1/3, 1, −1, ∞}) produces rational vectors **exactly identical** to cor13c_common's on all 38 HPL words and all 9 Li atoms (with matched frames) [step_44]. Pullbacks at the task map: 0 ↦ (3)−(0), 1 ↦ (4)−(3), −1 ↦ (3)−(2), as claimed.
- Numeric audit at 70 dps with my evaluator: 38 HPL conversions worst 6.13e-71; 9 Li atoms worst 1.21e-70 (Li₄(2/3)); 7 elementary atoms worst 1.18e-68. **Wrong composition order fails**: 36/38 words off (worst 0.905) — the audit discriminates the convention [step_47].
- Targets: mine == theirs exactly at weights 2, 3 (frames coincide) and at weight 4 with matched frames [step_44, step_46].
- **Referee strengthening 1**: all 36 pulled-back 1-form identities (12 maps × 3 HPL letters) verified **exactly over ℚ(ζ)** as rational-function identities, cross-multiplied in my own ℚ(ζ) polynomial arithmetic [step_49].
- **Referee strengthening 2 — frame independence**: with my own (different) rep-0 frames for z = 2/3, 1/2, −2, my targets vanish numerically (1.2e-65 / 3.55e-64 / 4.98e-64 at 65 dps) and the weight-4 difference (662 words) lies **exactly in the proved span** (rank 1695 = 1695 mod 3 primes; weights 2,3 identical) [step_45, step_46].
- Map independence spot check: target₂ vanishes under maps 3, 7, 11 (≤ 6.9e-61); target₃ under map 5 (4.55e-59) [step_56].

### Task 3 — the identities and COR7's closed form: CONFIRMED (250 digits for IDs; 1097+ digits for M₃)
- ID₂, ID₃, ID₄ recomputed at **250 dps** with my own HPL values at 1/3 (reflib truncated sums) and my own Li values (direct series; Jonquière inversion for Li₄(−2), cross-checked against mpmath to 1.22e-250): residuals **2.62e-251 / 7.02e-249 / 3.0e-249**. Auxiliaries: Ψ⁽⁰⁾(1/3) = 1/2 exactly (0.0), |Ψ⁽¹⁾(1/3) − ((21/2)log2 − (15/2)log3)| = 5.24e-251 [step_50].
- **M₃ assembly at 1150 dps** (my Li values, my HPL fit evaluation): |M₃(closed form) − c13_M3_quad_1100.txt| = 9.63e-1098 (**1097 digits**, the reference's own accuracy); |closed − c13_M3_hpl_1150| = 2.31e-1148 (**1147 digits**); |my G1-fit assembly − closed form| = 1.52e-1148; components W₂/W₃/W₄ vs stored references ≤ 1.43e-1145. M₃ = 0.00457121071807214710973628322267234274742470409643267589503504… [step_51].
- Li₄(1/9) = 8(Li₄(1/3)+Li₄(−1/3)) (the distribution input): 6.06e-211 at 210 dps with my own series [step_37].

### Task 4 — frame independence and the map count: CONFIRMED (exact)
- **N = 6**: my own exact enumeration over ℚ(ζ₆) (my own pair arithmetic, z² = z−1): exactly **12** maps, identical to mobius_maps_13.json; task map images 0↦1, 1/3↦ζ, 1↦ζ², −1↦ζ⁴, ∞↦ζ³ verified directly; each map omits exactly one root, 2 maps per root [step_33].
- **2N counts** for N = 6, 12, 18, 24, 30, 36: **12, 24, 36, 48, 60, 72** — my own mod-p filter (sharp primes, (p−1)/N ≥ 50: p = 307, 601, 919, 1201, 1531, 1801) + exact verification in my own ℚ(ζ_N) (my Cyclo class); **0** maps need 0 or ∞; exponent tuples match map_counts.json exactly at every level; all maps are shifts of the two conjugate base maps (0, N/6, N/3, 2N/3, N/2) and (0, 5N/6, 2N/3, N/3, N/2) [step_34, step_36].
- Cross-ratios of {0, 1/3, 1, −1, ∞}: 15 values forming the three stated S₃-orbits; all six Li arguments {1/3, −1/3, 2/3, 1/2, −1/2, −2} present (4/2/4/12/4/4 frames); 1/9 absent [step_37].
- Quotient dims **3, 8, 23, 69** (ranks 2/28/229/1695) vs Deligne 3, 8, 21, 55 — the 2+14 missing relations are not used; my own GF(p) elimination, two primes [step_52, step_53].

### Task 5 — verdict
**PROVED (computer-assisted, exact)** — no repairs needed. Exact (proof-level): the 12 maps and 2N counts (ℚ(ζ_N) arithmetic), the cross-ratio table, the relation row sets, the certificates Σ y·rows = T with residual exactly 0 on every word (three arithmetic backends), the Euler bonus ζ(2) = π²/6, π⁴ coming out exactly, quotient dims, the 36 pulled-back 1-form identities (referee strengthening), frame independence within the proved span (referee strengthening). Numerical only (as declared in COR13C's own caveat, and audited here with an independent evaluator at 65–70 dps on every piece, wrong-order discrimination demonstrated): the conversion's path-composition/regularization machinery; the CMZV/HPL evaluators; the ID₂/₃/₄ residuals at 250 dps and the M₃ value at 1150 dps (inputs from G1/COR7 re-confirmed to 250/1147 digits). COR13C's own verifier re-ran fresh in my hands: 25 s, 27/27 PASS, rc = 0 [step_54]. My own aggregate verifier `verify_cor13cr.py` (580 s): ALL PASS [step_57, step_58]. Displayed ID₃/ID₄ formulas match the certified JSONs exactly — no typos (unlike HEXC2, which had one) [step_60].

## The sentence the paper can carry

> **Theorem (third corner c₀ = 1/3, computer-assisted).** M₃(2 arccos(1/3)) = W₄ + W₃ + W₂ with W₂ = −Li₂(1/3) − Li₂(−1/3) + π²/12, W₃ = −26ζ(3) − 87Li₃(1/3) − 130Li₃(−1/2) + (21/4)π²log 2 − (43/12)π²log 3 + 67Li₂(1/3)log 2 − 45Li₂(1/3)log 3 + (65/3)log³2 − 65log²2 log 3 + (109/2)log 2 log²3 − (37/3)log³3, and W₄ = F₄ with F₄ = −(22/3)Li₄(−2) − (17/2)Li₄(1/3) + (20/3)Li₄(2/3) + (110/9)Li₄(1/2) + (25/48)Li₄(1/9) plus products of lower-weight periods over ℚ as displayed in COR7/COR13C; in particular M₃ at c₀ = 1/3 is classical. Each weight-homogeneous identity is an exact rational combination of proved relations among level-6 coloured MZV values (regularized double shuffle, Zhao's distribution relations of orders 2, 3, 6, the weight-one branch identities, and shuffle lifts), through the Möbius map sending (0, 1/3, 1, −1, ∞) to (1, ζ, ζ², ζ⁴, ζ³); the residual is exactly 0 on every word (36/252/1764 words of weights 2/3/4), the multipliers come from exact rational solves (heights ≤ 141 digits), π⁴ and ζ(2) = π²/6 come out exactly, and no relations beyond the standard proved ones enter (the 2 weight-three and 14 weight-four relations missing from the standard set play no role). The conversion's analytic input — the 36 pulled-back differential forms — is verified exactly over ℚ(ζ), and the targets are frame-independent modulo the proved relations; the remaining conversion machinery is audited numerically to ≥ 65 digits on every piece with independent evaluators.

## Evidence table (all numbers printed by my code)

| check | step | label | residual / result |
|---|---|---|---|
| relation row sets: reflib == cmzv6, weights 1–4 | step_38 | **exact** | 4/57/682/7869 rows, sets equal |
| certs Σ yᵢrows − T = 0 on every word, w2/w3/w4 + Euler | step_39 | **exact** | residual 0 on 36/252/1764/36 words; fmpq + pure-Fraction + mod 5 primes |
| multiplier stats | step_39/40 | exact | w2 25/28 nz ≤2 dig; w3 229/229 ≤10; w4 1692/1695 ≤141 (med 139); euler 22/28 ≤2 |
| T₄ coeff of (1,1,1,1) (π⁴ exact) | step_39 | exact | 2395/6, killed by certificate |
| tamper test (perturbed certs fail) | step_56 | exact | 2 and 8 words off |
| my evaluator vs fast_cmzv / numeric_cmzv | step_41 | cross-check | worst 9.3e-63 / 7.29e-63 (150 words, 65 dps) |
| all 1980 cert rows + unused vanish numerically | step_43 | numeric 65 dps | w2 3.92e-65; w3 1.9e-64; w4 1.37e-63 |
| family usage (no frame rows) | step_43/59 | exact count | w2 20+5+1+2; w3 203+11+2+13; w4 1612+26+5+50+2 |
| 12 maps exact == mobius_maps_13.json; task-map images; 2 per omitted root | step_33 | **exact over ℚ(ζ₆)** | identical |
| 2N counts N = 6,12,18,24,30,36; 0 maps need 0/∞ | step_34/36 | **exact** | 12/24/36/48/60/72 == map_counts.json; two shift orbits |
| cross-ratios: 6 Li args present, 1/9 absent, 3 orbits | step_37 | exact rationals | 15 values, orbits as stated |
| Li₄(1/9) = 8(Li₄(1/3)+Li₄(−1/3)) | step_37 | numeric | 6.06e-211 @210 dps |
| my conversion == cor13c (38 HPL + 9 Li + targets) | step_44/46 | **exact** | rational vectors identical |
| conversion audit: HPL / Li / atoms | step_47 | numeric 70 dps | 6.13e-71 / 1.21e-70 / 1.18e-68 |
| wrong composition order fails | step_47 | numeric | 36/38 words off (worst 0.905) |
| 36 pulled-back 1-form identities | step_49 | **exact over ℚ(ζ)** | 36/36 hold (12 maps × 3 letters) |
| frame independence | step_46 | numeric + **exact** | my-frame targets ≤ 4.98e-64 @65 dps; w4 diff in proved span mod 3 primes |
| map independence | step_56 | numeric 60 dps | target₂ vanishes under maps 3,7,11; target₃ under map 5 |
| ID₂ / ID₃ / ID₄ residuals | step_50 | numeric 250 dps | 2.62e-251 / 7.02e-249 / 3.0e-249 |
| M₃ closed form vs quad(1100) / hpl(1150); my assembly | step_51 | numeric 1150 dps | 1097 / 1147 / 1147 digits; components ≤ 1.43e-1145 |
| quotient dims 3, 8, 23, 69 vs Deligne 3, 8, 21, 55 | step_52/53 | exact mod 2 primes | ranks 2/28/229/1695; missing 0,0,2,14 |
| verify_cor13c.py fresh run | step_54 | — | 25 s, 27/27 PASS, rc 0 |
| verify_cor13cr.py (my aggregate verifier) | step_57/58 | — | 580 s, ALL PASS, rc 0 |
| displayed ID₃/ID₄ == certified JSONs | step_60 | exact | no display typos |

## Files (new_results/COR13CR_referee_certificate)

- `reflib13.py` — the referee's c₀ = 1/3-specific independent library (my 12-map enumeration over ℚ(ζ₆), HPL pullbacks at 1/3, cross-ratio frames from {0, 1/3, 1, −1, ∞}, Li atoms, elementary atoms, targets); reuses HEXC2R's `reflib.py` for the generic machinery [written in step_32].
- `mymaps.py` — my own cyclotomic arithmetic ℚ(ζ_N) + mod-p map counting [step_36].
- `verify_cor13cr.py` — the 580 s aggregate verifier (PART A exact; PART B numeric audit) [step_57].
- Producing steps: step_33 (12 maps exact), step_34/36 (2N counts, N ≤ 36), step_37 (cross-ratios, orbits, Li₄(1/9)), step_38 (row-set equality), step_39/40 (certificates exact, heights), step_41 (evaluator cross-check), step_43 (row validity 65 dps, families), step_44–46 (conversion exact + frame independence), step_47 (70-dps audit + wrong order), step_49 (1-forms exact), step_50 (IDs 250 dps), step_51 (M₃ 1150 dps), step_52/53 (dims), step_54 (their verifier fresh), step_56 (tamper + map independence), step_58 (my verifier fresh), step_59 (Euler/w2 families), step_60 (display check). One line per finished step in `progress.log`.