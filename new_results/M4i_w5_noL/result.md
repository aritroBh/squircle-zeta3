# M4i — W5 without L(χ8,5), and the explicit L-relation (Kimi K3)

**Work dir:** `new_results/M4i_w5_noL` · `PY='<project>/.venv/bin/python'` · M4f copied (34 files), never edited in place.

## Verdict

**L(χ8,5) is droppable (numerically, with exact symbol and exact (3,2)-cobracket certificates).**
1. `drop_L.py` runs as-is: **`target in span(kv, u0..u6) WITHOUT L: True`** (tbar coords on L, s2L = −35426/315, 2152/15).
2. The explicit relation below puts L(χ8,5) in the Q(√2)-span of Li₅ at the M4e arguments + ζ(5) + products — in fact **all 26 coefficients are pure √2**, i.e. **L(χ8,5)/√2 is a plain Q-combination of the classical atoms**. Verified: exact symbol 0, exact cobracket 0, and fresh 1620-digit value agreeing to **2.8e-1741** with an independent Hurwitz-ζ evaluation of L.
3. `R5_noL.json` (no L column anywhere) closes W5 = N5 + R5 with the same fresh-value residual as R5_best: **−2.1052657e-1596**.

## The L relation (`L_relation.json`, reductions Li₅(b²)→16[Li₅(b)+Li₅(−b)], Li₅(−1)→−15ζ(5)/16, Li₅(−1/√2)→Li₅(1/2)/16−Li₅(1/√2) applied)

**26 terms · max |numerator| = 4815 · lcm denominator = 638,976 = 2¹⁴·3·13** (13 Li₅ of the 14 arguments — Li₅(1−b⁴) unused — plus ζ(5) and 12 products). Compact u-form found by exact rref in the 14-dep lattice:
√2·L = (255/1664)u₀ − (395/1664)u₁ + (75/416)u₂ + (25/416)u₃ − (15/26624)u₄ − (5/416)u₅ + (21-term kv combination), then expanded to atoms and minimised over the 12-dim exact-zero freedom (LLL gram=exact + rational descent + restarts; raw rref point was already optimal).

10 largest coefficients (term: [p, q], coefficient = p+q√2; all p=0):
| |num| | term | coefficient |
|---|---|---|
| 4815 | z5 | (4815/106496)√2 |
| 915 | Li₅(1/2) | −(915/3328)√2 |
| 395 | Li₅(−2+2√2) | −(395/3328)√2 |
| 335 | λ·π²·π² | (335/159744)√2 |
| 265 | λ⁵ | −(265/212992)√2 |
| 255 | Li₅(1/2−√2/2) | (255/3328)√2 |
| 185 | Li₅(−1+√2) | −(185/208)√2 |
| 145 | Li₅(−3+2√2) | (145/3328)√2 |
| 125 | λ·L⁴ | −(125/26624)√2 |
| 115 | λ³·L·π²… (tie: λ²·L², L²·π²·…) | ∓(115/13312)√2 |

Full term list (λ=log2, L=log(1+√2)): λ⁵ −265/212992; λ⁴L 25/6656; λ³L² −115/13312; λ²L³ 105/13312; λL⁴ −125/26624; L⁵ 67/26624; λ³π² 95/26624; λ²Lπ² −105/13312; λL²π² 85/13312; L³π² −115/39936; λπ²π² 335/159744; Lπ²π² 79/159744; Li₅(−1+√2) −185/208; Li₅(1−√2) −55/52; Li₅(−3+2√2) 145/3328; Li₅(√2/2) 75/104; Li₅(1/2) −915/3328; Li₅(1−√2/2) −5/208; Li₅(2−√2) 45/208; ζ(5) 4815/106496; Li₅(1/2−√2/2) 255/3328; Li₅(−2+2√2) −395/3328; Li₅(1/2−√2/4) 75/832; Li₅(1/2+√2/4) 25/832; Li₅(17−12√2) −15/53248; Li₅(1/2−3√2/8) −5/832 — each times √2.

## R5_noL vs R5_best

| form | terms | Li₅ args | products | max \|num\| | lcm den | L(χ8,5) |
|---|---|---|---|---|---|---|
| R5_best | 70 | 14 | 54 (+ζ5) | 80,354 | 11,612,160 | (−80354/315, 22520/3) |
| **R5_noL** | **69** | **14 (same set)** | **54 (+ζ5)** | **2,487,169** | **301,916,160** | **none** |

Optimisation path (L-free affine space: 84 columns, nullity 12): rref particular 166,920,115 → substituted-R5_best seed 10,469,591 → descent 2,487,169; reweighted LLL-CVP rounds + ~460 shake restarts gave no further improvement (accepted). The ~31× numerator growth vs R5_best is the price of eliminating L (L-coeff of R5_best · L-relation numerator 4815 ≈ 1.2M lower scale).

## Verification outputs (all printed by my code)

| check | result | precision / digits |
|---|---|---|
| drop_L exact span test | **True** (no L columns) | exact Q (flint rref) |
| L-rel: (a) exact full point symbol of RHS = 0 | **PASS** | exact |
| L-rel: (b) all 13 Li₅ reduce to 0 in (3,2) cobracket | **PASS** | exact |
| L-rel: smoke (cached S0 + fresh extras) | −1.1367e-1310 | 1320 dps (floor) |
| L-rel: (c0) Hurwitz 8⁻⁵[ζ(5,1/8)+ζ(5,7/8)−ζ(5,3/8)−ζ(5,5/8)] − Lchi8_5.json | 4.2434192e-1656 | fresh 1620 (file 1650 d) |
| L-rel: (c) RHS − L[Hurwitz] | **2.7833324e-1741** | fresh 1620 (dps 1740, arb 1680/1780 d) |
| R5_noL: exact symbol == target_symbol | **PASS** | exact |
| R5_noL: every Li₅ (incl. b⁴, 1−b⁴) cobracket-0; reduce(N5−target)=0 | **PASS / PASS** | exact |
| R5_noL: fresh R5 + N5 − W5_5400 | **−2.1052657e-1596** (identical to R5_best) | fresh 1620 |
| Cross: Z = R5_noL − R5_best − c_L·(Lrel−L), c_L=(−80354/315, 22520/3): exact symbol 0; value | **PASS; 7.5150e-1740** (54 slots, max num 3743) | exact; fresh 1620 |
| deps 1200+1500 combined rank | 14 (same span); L/s2L projection rank 2 each | exact |

## Files and rerun commands (cwd = `new_results/M4i_w5_noL`)

Artifacts: `L_relation.json`, `R5_noL.json` (R5_best format `{'terms': {name:[p,q]}}`), `lrel_state.pkl`, `state_noL.pkl`, `progress.log`; scripts `lrel_build.py, lrel_smoke.py, lrel_verify.py, nol_build.py, nol_optimize2.py, nol_export.py, nol_cross.py, deps_consistency.py` (+ copied M4f scripts; harness-saved `step_01.py…step_35.py`: step_17=drop_L run, step_19=lrel_build run, step_21/22=smoke, step_23=lrel_verify, step_25/26=nol_build, step_29/30=nol_optimize2, step_31=export, step_32=verify_json, step_33=nol_cross, step_34=deps_consistency).

```
nice $PY drop_L.py                  # L-free span test: True
nice $PY lrel_build.py 420          # build+minimise L relation -> L_relation.json
nice $PY lrel_smoke.py              # 1200-digit smoke: -1.14e-1310
nice $PY lrel_verify.py             # symbol=0, cobracket=0, Hurwitz cross-check, 2.78e-1741
nice $PY nol_build.py 60            # L-free space + first optimisation -> state_noL.pkl
nice $PY nol_optimize2.py 750       # reweighted CVP + restarts (resumable)
nice $PY nol_export.py              # -> R5_noL.json
nice $PY verify_json.py R5_noL.json # symbol, cobracket, fresh 1620d: -2.1052657e-1596
nice $PY nol_cross.py               # R5_noL-R5_best = c_L(Lrel-L)+exact-zero, 7.5e-1740
nice $PY deps_consistency.py        # 1200/1500 dep spans agree (rank 14)
```
---

## Planner verification (W4 session, 2026-10-27)
Both claims are checked with planner code that shares nothing with the M4c/M4e/M4f/M4i scripts: direct power series for
every Li_k and Li_{3,2}, N5 hard-coded, and L(χ8,5) from mpmath's Hurwitz ζ.
- **The L relation** (`../M4e_w5_lchi/planner/planner_check_Lrel.py L_relation.json 1700`):
  - 26 terms; every coefficient is purely a multiple of √2; max |numerator| 4815; lcm denominator 638976 = 2^14·3·13.
  - RHS − L(χ8,5) = 3.6e-1740 at 1700 digits, the working-precision limit.
  - So L(χ8,5)/√2 is a Q-combination of Li5 at 13 S-units of Q(√2), ζ(5), and λ, L, π² products. It is numerical; the
    exact symbol and cobracket are 0 (Kimi).
- **R5_noL** (`../M4e_w5_lchi/planner/planner_eval_terms.py R5_noL.json 1600`):
  - R5 + N5 − W5 = −2.1e-1596 against M4c's W5, and 2.8e-601 against M4a's independent W5.
  - It has 73 terms including N5; max numerator 2,487,169; lcm denominator 301,916,160.
- **Status:** L(χ8,5) is droppable (numerically). The paper keeps R5_best (max numerator 80,354), which is far smaller.
  The L relation is a candidate remark: an explicit weight-5 Zagier-type formula for ζ_{Q(√2)}(5)/ζ(5) = L(χ8,5) in
  real Li5 values.
