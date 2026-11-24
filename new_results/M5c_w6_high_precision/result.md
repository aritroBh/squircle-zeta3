# M5c — W6 − N6 at 3000–3300 digits with reserve Li₆ extras

**Verdict (NUMERICAL, nothing here is a proof): W6 is NOT closed.** At 3500‑digit values, the 179‑column kernel LLL at **3000** and **3300** digits found no target relation, and neither did the enlarged lattices with the reserve Li₆ arguments **(4−3√2)/8**, **1−b⁴**, **2b** added cumulatively (181, 183, 185 columns at 3300 digits). In every span tried, all reduced rows are spurious with coefficients ~10¹⁹ and residuals ~1e‑3283, so any genuine relation in these spans has coefficients ≳ 10¹⁹ — far above the weight‑5 pattern (max numerator 3.6·10⁷), strongly supporting the "missing symbol‑invisible directions" explanation. The remaining reserve arguments b⁴, b², −1/√2 were **provably redundant** (exact duplication combos of base atoms, verified to ~1e‑3610) and were skipped. The pipeline itself was validated by a full‑lattice planted‑relation control (PASS).

## 1. Precision (all targets met or exceeded)

| quantity | precision | control | result |
|---|---|---|---|
| W6 @11500 bits | radius 5.0e‑3459, tail 1.0e‑3474 (3431 digits) | vs quad5_100 | 4.3e‑101 ✓ |
| W6 @11800 bits | radius 2.6e‑3549, tail 5.7e‑3565 (3522 digits) | agreement 11500/11800 | **2.3e‑3432** ✓ |
| W6 @11800 bits | | vs W6_400.txt / W6_5400.json | 3.7e‑401 / 3.5e‑1651 ✓ |
| L(χ8,5) | 3600 digits, Hurwitz‑ζ @12600 bits, rad 5.3e‑3793 | vs Lchi8_5_2500.json | 3.4e‑2511 ✓ |
| 45 leaves + 149 atoms | 3500 digits (ctx 12092 bits) | two‑precision (12700 vs 12092 bits) | all ≤ 9e‑3559 after fix ✓ |
| Li₆ duplication | 1/√2 pair vs 1/2; b pair vs b² | | 4.4e‑3611 / 2.8e‑3611 ✓ |
| Li₃,₂(b,1) | mpl_general vs M4c `w5_search.mpl` | | 9.5e‑3542 ✓ |
| N6 atoms / N6 sum | vs M5a `n6_values.json` | | ≤4.6e‑77 / 4.3e‑75 ✓ |

**Bug found and fixed (affects M5b retroactively).** `w6_values.py::mpl_general` chose the truncation N from the rate |xy|, but for |y|≤1 the outer series decays like |x|ⁿ (H(n)→Li_t(y)=O(1)). Li₄,₂(b,√2/2) — one of the four N6 atoms — was therefore accurate to only ~1e‑2582 in the fresh 3500‑digit cache (old−new = 2.3e‑2582) and to only **~1e‑1751 in M5b's 2400‑digit cache, which invalidated M5b's own 1900‑digit LLL** (its 1300/1600‑digit runs were unaffected). Recomputed with rate |x|·max(1,|y|) and a rigorous ball tail: radius 3.7e‑3639; two‑N agreement 1.8e‑3539; cache patched (the other 7 mpl leaves were already fine, old−new ≤1.1e‑3584; no basis atom contains this leaf, only the N6 target does).

## 2. LLL ladder (179+ columns, gram='exact', 3500‑digit values, W6 = w6_11800)

| lattice | cols | time | best rows (spurious) | genuine relation |
|---|---|---|---|---|
| base @3000 digits | 179 | 3250 s | maxcoef 4.9e17–8.2e17, residual ~1e‑2984 | none |
| base @3300 digits | 179 | 3916 s | maxcoef 2.3e19–3.0e19, residuals 3.0e‑3282…3.4e‑3283 | none |
| +(4−3√2)/8 @3300 | 181 | 3200 s | maxcoef 1.5e19–2.5e19, residual ~1e‑3283 | none |
| +1−b⁴ @3300 | 183 | 5770 s | maxcoef 9.7e18–1.3e19, residual ~1e‑3283 | none |
| +2b @3300 | 185 | 4909 s | maxcoef 6.6e18–9.9e18, residual ~1e‑3284 | none |

- Planted‑relation controls: 3‑column (PASS, every run) and **full 179‑column lattice @1300 digits: planted 6‑term relation recovered exactly, residual 2.5e‑1538 (1128 s)** — so a genuine relation with coefficients ≲1e19 WOULD have been found.
- det^{1/n} scale: ~10^17.3 @3000, ~10^18.6 @3300; genuine coefficients exceed ~10^19 in all five spans.
- Extras b⁴(13), b²(7), −1/√2(6) skipped after verifying exact duplication: Li₆(b²)=32(Li₆(b)+Li₆(−b)) residual 8.5e‑3611; Li₆(b⁴)=1024‑base‑combo residual 2.8e‑3608; Li₆(−1/√2)=Li₆(1/2)/32−Li₆(√2/2) residual 5.1e‑3611. A {14,15,13} run was started, then killed when 13 was proved redundant (span ≡ {14,15}).
- Symbol stages for all cumulative extra sets: rank stays 60, target in span, saturated kernel dims 89 (base), 90, 91, 92, 93, 94, 95 (SNF invariants all 1).

## 3. How far it got / what remains

Values are solid to 3500 digits and the LLL machinery is validated, so the obstruction is structural, not numerical: the current atom set (products + Li₆ at the 10 base args + the three genuinely‑new pool extras) does not contain the target's zero‑symbol constants at any findable height. Next steps (not run): (a) 3600–4200‑digit LLL (resolves ~1e20–1e22; needs W6 ≈13400 bits, leaves ≈3950 digits — the value scripts support this directly); (b) new kernel columns from zero‑symbol depth‑2 constants (Li₃,₃/Li₄,₂ combos) or other primitives — analysis had just begun (w4_symbol.py inspection) when the turn budget ended.

## 4. Files (work dir `new_results/M5c_w6_high_precision/`)

- Scripts written this session: `w6c_lll.py` (driver), `w6c_values_extra.py`, `w6c_control.py`, `w6c_assemble.py`, `w6c_verify.py`, `w6c_indep_check.py`, `launch_*.sh`; patched copy `w6_values.py` (Lchi8_5_3600); M5b scripts copied unmodified otherwise. Tool‑call transcripts: `step_01.py`…`step_~70.py` (env checks, launches, mpl fix, duplication checks, monitors).
- Data: `w6_11500.json`, `w6_11800.json` (+ `_partial`), `Lchi8_5_3600.json`, `w6_leaves_3500.pkl` (Li42 patched), `w6_atomvals_3500.pkl`, `mpl_recheck.json`, `w6x*_symbol_solve.json` (7 sets), `w6_leaves_extra_w6x*.pkl`, `w6_atomvals_3500_x*.pkl`, `w6c_lll{,_w6x14,_w6x14_15,_w6x14_15_9}_{3000,3300}.json`, `progress.log`.

Rerun commands (`PY='<project>/.venv/bin/python'`, all detached with `nice -n 10`):
```
$PY -B m5b_precision.py w6 11500 11800      # W6 quadrature (262s/553s)
$PY -B w6_values.py 3500                    # leaves+atoms (99s); then re-patch Li42 via mpl_recheck.json
$PY -B w6c_lll.py 3000            # base LLL; also 3300, and '... 3300 w6x14_15_9'
$PY -B w6_solve.py 14 15 9        # symbol stage for extras
$PY -B w6c_values_extra.py 14 15 9
$PY -B w6c_control.py 1300        # full-lattice planted control
```
All claims above are NUMERICAL evidence printed by these runs; none constitute a proof.