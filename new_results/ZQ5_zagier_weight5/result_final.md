# ZQ5: explicit weight-5 formulas for L(χ_D,5) — complete result

**Everything below is NUMERICAL; nothing is a proof.** All relations are integer relations found by LLL and verified to the stated precision with exact symbol/cobracket certificates.

## Headline results

| D | field | S | formula? | terms | max num | lcm den | residual (fresh) | symbol 0? | cobracket 0? |
|---|-------|---|----------|-------|---------|---------|------------------|-----------|--------------|
| 8 (control) | Q(√2) | {2} | **YES** | 24 (11 Li5 + ζ5 + 12 prod) | 1955 | 1916928 = 2¹⁴·3²·13 | −5.96e-1216 @1200d | YES | YES |
| 5 | Q(√5) | {5,2} | **YES** | 7 (4 Li5 + 0 ζ5 + 3 prod) | 1728 | 103125 = 3·5⁵·11 | +3.39e-1416 @1400d | YES | YES |
| 12 | Q(√3) | {3,2} | **YES** | 88 (53 Li5 + ζ5 + 34 prod) | 913 | 15552 = 2⁶·3⁵ | −3.00e-1416 @1400d | YES | YES |
| 13 | Q(√13) | {13,2} | NO | — | — | — | — | — | — |
| 21 | Q(√21) | {3,7,2} | NO | — | — | — | — | — | — |
| 29 | Q(√29) | {29,2} | NO | — | — | — | — | — | — |
| 17 | Q(√17) | {17,2} (2 split) | NO | — | — | — | — | — | — |

S = ramified primes plus 2 throughout, with letters: D=8: log√2, log(1+√2); D=5: log φ, log √5, log 2; D=12: log(2+√3), log √3, log(1+√3); D=13: log((3+√13)/2), log √13, log 2; D=21: log((5+√21)/2), log √21, log((3+√21)/2), log 2; D=29: log((5+√29)/2), log √29, log 2; D=17: log(4+√17), log √17, log((3+√17)/2), log((3−√17)/2).

## The formulas

**D = 5** (ρ = 1/φ = (√5−1)/2), L(χ5,5) = √5·T with
> T = (1728/6875)[Li5(ρ) − Li5(−ρ)] − (256/6875)[Li5(ρ³) − Li5(−ρ³)] + (648/34375)·log⁵φ − (48/6875)·π²·log³φ + (64/20625)·π⁴·log φ.

Equivalent (unladdered) form: T = (1728/6875)Li5(1/φ) − (1728/6875)Li5(−1/φ) − (512/6875)Li5(φ⁻³) + (16/6875)Li5(φ⁻⁶) + same 3 products, using Li5(ρ⁶) = 16[Li5(ρ³) + Li5(−ρ³)]. Note: no ζ(5) term, no log 2 / log √5 products although 2 ∈ S. Factorisations: v(ρ) = [−1,0,0], v(1−ρ) = [−2,0,0]; v(−ρ) = [−1,0,0], v(1+ρ) = [1,0,0]; v(ρ³) = [−3,0,0], v(1−ρ³) = [−2,0,1]; v(ρ⁶) = [−6,0,0], v(1−ρ⁶) = [−3,0,2] in letters (φ, √5, 2). Ladder form verified fresh: residual 3.9e-1316 @1300d; orig−alt = 1.8e-1317.

**D = 8** (control): T = L(χ8,5)/√2 = −(645/16384)ζ(5) − (75/104)Li5(−1/√2) − (35/208)Li5(1−√2) + (25/104)Li5(2−√2) + (275/3328)Li5((1−√2)/2) − (225/1664)Li5(1/2) + (665/6656)Li5(−3+2√2) − (375/3328)Li5(−2+2√2) + (25/416)Li5((2−√2)/4) − (385/106496)Li5(17−12√2) − (15/3328)Li5((4−3√2)/8) + (5/3328)Li5(−16+12√2) + (1255/79872)λ₂⁵ − (925/39936)λ₂⁴L + (−25/9984)(λ₂³L² − λ₂²L³) + (25/39936)λ₂L⁴ + (31/26624)L⁵ + (1175/119808)λ₂³π² − (175/19968)λ₂²Lπ² + (25/19968)λ₂L²π² − (55/39936)L³π² + (1955/479232)λ₂π⁴ + (175/159744)Lπ⁴, with λ₂ = log√2, L = log(1+√2). Agrees with M4i's `L_relation.json` numerically (M4i−mine = −1.24e-1215; M4i's own formula rechecked with my code to −1.0e-845 @830d; planner independently verified mine to −7.0e-1331 @1300d). M4i's 13 arguments are all in my 16-element exceptional list and their coefficient vector lies in my saturated kernel (rank-4 symbol map, 12-dim kernel).

**D = 12**: T = L(χ12,5)/√3 = (253/7776)ζ(5) + Σ (53 Li5 terms) + Σ (34 products). Li5 terms (x : coefficient): (−1+2√3/3): 77/5184, (1−2√3/3): −19/216, (1−√3/3): −77/2592, (5/2−3√3/2): 29/2592, (2−√3): 253/7776, (−2+√3): −17/1944, (−1+√3): 469/3888, (1−√3): −913/7776, (−3+2√3): −11/864, (3−2√3): −59/864, (1/2−√3/6): −7/162, (√3/3): −19/108, (−√3/3): −10/27, (−1/2+√3/2): 331/1296, (1/2−√3/2): 13/648, (3/2−√3/2): 5/432, (1/2+√3/6): 23/2592, (7−4√3): −5/288, (−7+4√3): 2/243, (−5+3√3): −7/96, (4−2√3): 137/3888, (−6+4√3): −13/432, (2/3−√3/3): 1/288, (1/3−√3/3): 37/1296, (2/3): −1/162, (−1/3): −5/144, (1/3+√3/3): −103/2592, (1/2−√3/3): 47/648, (1/2): −509/3888, (−1/2): −499/7776, (√3/2): 25/5184, (−√3/2): 47/2592, (−26+15√3): 11/3888, (10−6√3): 113/7776, (2/3−2√3/9): −31/7776, (2−2√3/3): 103/2592, (1/2−5√3/18): −1/324, (1/3−2√3/9): 41/7776, (−5/4+3√3/4): −1/108, (1/3+2√3/9): 103/15552, (1/4−√3/4): −103/3888, (3/4−√3/4): −35/1296, (9/4−3√3/4): 23/7776, (1/4+√3/4): 19/1944, (1/2−√3/4): 7/432, (1/2+√3/4): −1/864, (−20+12√3): −37/7776, (−4/3+4√3/3): −7/2592, (5/8−3√3/8): 23/7776, (−96+56√3): 1/1296, (8/9): 1/2592, (1/2−7√3/24): −5/2592, (−1/8): 5/7776. Products (e = log(2+√3), s = log√3, u = log(1+√3)): e⁵: 1/288, e⁴s: −173/15552, e⁴u: −61/15552, e³s²: 91/7776, e³su: 47/2592, e³u²: −65/15552, e²s³: −17/1944, e²s²u: −53/5184, e²su²: −1/2592, e²u³: 31/3888, es⁴: 295/15552, es³u: −23/1296, es²u²: 1/72, esu³: −5/324, eu⁴: −101/15552, s⁵: 7/1296, s⁴u: −59/1944, s³u²: 617/15552, s²u³: −371/15552, su⁴: 95/7776, u⁵: 7/1944, e³π²: −53/15552, e²sπ²: 103/7776, e²uπ²: −1/2592, es²π²: −17/1944, esuπ²: −77/2592, eu²π²: 73/5184, s³π²: −7/1296, s²uπ²: 313/15552, su²π²: 73/7776, u³π²: −79/7776, eπ⁴: 13/5184, sπ⁴: 7/7776, uπ⁴: −11/7776. (Alternative representative: 86 terms, max num 1907.)

## Method and certificates

Pipeline (generic in D): exact Q(√D) arithmetic (Fraction-based) → S-unit factorisation (norm valuations + unit-power recovery, tested on 20–30 random S-units per D, plus rejection cases) → exceptional S-units with |x| < 1 (B = 8; completeness checked at B = 10 for D=21 and B = 14 for D=13) → exact weight-5 symbols −(1−x)⊗x^⊗4 → saturated integer kernel (HNF/saturation with SNF certificates = all 1s) of "symbol ∈ Sym⁵" → numerics at 1150 digits (Li5 direct series, Jonquière expansion for 0.9 < x < 1, cross-checked vs polylog; L via Hurwitz ζ, cross-checked at P and P+150 and against M4i's 1650-digit L(χ8,5) to 1007 digits) → LLL relation search (flint, gram='exact', planted-relation test passed) at 520 and 830 digits → minimisation over the exact-zero freedom (duplication relations + 2–27-dim zero lattices) → verification.

Per-D verification (a)–(e):
- **D=8**: (a) exact relation at 520/830 differ only inside the 5-dim exact-zero span (span test: True); residuals 4.9e-1161/1.6e-1161; (b) exact symbol ZERO; (c) (3,2)-cobracket ZERO (full & cut, total & Li5-part & products-part; each Li5 symbol individually kills δ); (d) all 11 used args have x, 1−x S-units (table in step_21 output); (e) 24 terms, max num 1955, lcm 2¹⁴·3²·13.
- **D=5**: (a) two-precision in zero span: True, a0 = 34375 = 5⁵·11 both; (b) symbol ZERO; (c) cobracket ZERO (both modes, total and Li5-part); (d) 4 args factorised (above); (e) 7 terms, max num 1728, lcm 3·5⁵·11. Planner independently verified to −1.8e-1530 @1500d.
- **D=12**: (a) two-precision in zero span: True, a0 = −15552 both; (b) symbol ZERO; (c) cobracket ZERO; (d) 53 args factorised (step_41 output); (e) 88 terms, max num 913, lcm 2⁶·3⁵.

**Unit-only test (task 2)**: searched [T, ζ(5), Li5(±ε⁻ᵏ) k = 1..6, log⁵ε, π²log³ε, π⁴logε] at 520 and 830 digits. D=5: YES (the formula above is unit-only). D=8: **NO** exact relation at either precision (0 found). D=12: **NO** exact relation at either precision (0 found). (Detection range ~1e⁵⁰ coefficients at dim 17, 830 digits.)

**Negative results (task 3)**: D=13: only 4 exceptional S-units exist (complete at B = 14); kernel dim 1; no exact relation at 520/830 (both targets L/√13 and L). D=29: same situation (4 candidates, kernel dim 1, no relation). D=21: 78 exceptional S-units (complete at B = 10), kernel dim 17; no exact relation at 520/830, target L/√21 or L (values certified good to ~1e-1150; mpmath's polylog, not our series, was the imprecise one near x = 1). D=17: 147 exceptional S-units, kernel dim 63; no exact relation at 520/830, target L/√17 or L. Whether these fields have formulas with much larger coefficients (> ~1e⁴⁵), enlarged S, or none at all remains open.

## Files and rerun commands

Work dir `new_results/ZQ5_zagier_weight5`, PY = `<project>/.venv/bin/python`:
- `w5lib.py` — exact Q(√D) arithmetic, S-unit factorisation, exceptional-unit enumeration, characters, Hurwitz L, Li5 series/Jonquière, symbols, shuffle products, saturation (HNF/SNF), (3,2)-cobracket.
- `pipe.py` — pipeline: phase_enum / phase_kernel / phase_numerics / phase_lll, products, minimisation helpers.
- `out/d{D}/`: `candidates.json`, `kernel.json`, `vals_1150.pkl`, `rel_min*.json`, `atom_lll.pkl` (D=5, 12, 21), plus `kernel_b10.json`, `vals_b10.pkl` (D=21).
- `steps/step_NN.py` — the steps that produced every printed number (key ones: step_06_tests.py, step_08_d8_enum.py, step_16_d8_min.py, step_21_d8_verify.py, step_26_shufflefix rerun, step_28–33 (D=5), step_35–41 (D=12), step_42_unitonly.py, step_43–45 (D=13), step_46 (D=29), step_48–52 (D=21), step_53–58 (D=17), step_59_ladder.py).

Rerun examples:
```
nice $PY steps/step_21_d8_verify.py    # D=8 control (a)-(e) + M4i agreement + fresh 1200d
nice $PY steps/step_33_d5_verify.py    # D=5 (b)-(e) + fresh 1400d
nice $PY steps/step_41_d12_verify.py   # D=12 (a)-(e) + fresh 1400d
nice $PY steps/step_42_unitonly.py     # unit-only tests D=8/D=12
nice $PY steps/step_59_ladder.py       # D=5 ladder form + final formula printouts
nice $PY steps/step_44_d13_lll.py      # D=13 negative result
nice $PY steps/step_51_d21_b10.py      # D=21 negative result (B=10)
nice $PY steps/step_57_d17_lll830.py   # D=17 negative result
```

Notable implementation pitfalls encountered (for reuse): `mp.mp.dps` (not `mp.dps`); unary `+v` in mpmath rounds to current dps (precision trap in our early code); `str(mpf)`/`mp.mpf(str)` at ambient dps truncated the value cache to 17 digits (fixed by storing `mp.nstr(v, P+25)` and parsing at high dps); shuffle-product symbols need multiplicity factors ∏ multᵢ!; π²-products have zero symbol (not identity); letters that are negative field elements need log|·| (mp.log(abs(...))); flint LLL needs gram='exact' at these precisions; spurious LLL rows are separated from exact relations only by a ~1e-1080 residual cut.

**Conclusion**: weight-5 Zagier-type formulas for L(χ_D,5) exist numerically for D = 8 (control, matching M4i), D = 5 (a striking 7-term unit-power formula in φ), and D = 12 (88 terms); none were found for D = 13, 21, 29, 17 with S = ramified ∪ {2} and candidates up to B = 8–14. The D = 5 formula is the analogue of Zagier's 1991 weight-3 formula for Q(√5); the D = 5 and D = 12 formulas appear to be the first weight-5 real-quadratic formulas beyond Q(√2).