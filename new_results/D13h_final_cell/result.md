# D13h complete: certified decay theorem for ALL N ≥ 100

## THEOREM — PROVED (arb) †

For every integer N ≥ 100,

**|M_k^main + (√π/N)·e^{−√(πN)}·cos√(πN)| ≤ C·N^{−3/2}·e^{−√(πN)}** with **C = 617117.6**

† every step an arb enclosure (python-flint, ctx.dps = 70); the single classical citation is Binet's log-Γ remainder formula (see §4), itself numerically verified on the ray arg z = π/4 (worst ratio 0.7071 = 1/√2).

C = max of three certified suprema (each an arb upper bound):

| range | cover | sup C | attained at |
|---|---|---|---|
| [100, 400] | 29 ratio-1.05 cells, `certified_y38` (J=24, nsub=4) | **194.3809** | [373.3456, 392.0129] |
| [400, 10¹⁰] | 43 exactly-contiguous cells (planner's sweep; used, not redone) | **617117.5536** | [6.6344e9, 9.9516e9] |
| [10¹⁰, ∞) | ONE analytic cell y ∈ [0, y_J] (this work) | **49.5626** | — |

## Task 1 — FINAL CELL: C = 49.56260 for all N ≥ 10¹⁰ (PROVED/arb)

`step_10.py` (J=48 → C = 53.61569), `step_10_J96.py` (C = 50.85487), `step_10_J192.py` (**C = 49.56260**); J only refines, every J certifies. Every piece of step_38's `cell()` is written as y³·(regular factor) and every regular factor is evaluated on the box y ∈ [0, y_J], σ ∈ [0,1], ρ = 1+y²σ, t = ρy ∈ [0, 1.25334e-5] with **no division by y**:
- Factored: Φ3b = y²Ph3reg, Φ4b = y³Ph4reg (πρ⁴/3πρ⁵ rational forms), G1 = yG1reg, G2 = y²G2reg, c2 = yc2reg (box [1.1107 ± 6.0e-5]), om = y·omreg ([1.5708 ± 5.9e-5]).
- g0 via h₁-Binet: |g0| ≤ t^{5/2}(√π/4)2^{3/2}e^{ReE}, ReE ≤ (t/24)(−√2 + 459t²/(1−3T)) ≤ 0 (asserted).
- v-moments: v = u/x0 = √y·z, Gaussian e^{−K2reg z²} exact (K2reg = c2reg/ρ²); zmax_lo = 3.2886 > zstar = 3.
- Tails: d = w/√(2c2) (< x0/3 asserted), c2d² = w²/2 exact, e^{−c2d²} = y√(2/π)e^{−c_w/2} exact; curvature-drift corrections via Lip_η = 4.9953 and √y·w³ ≤ √y_J·w_J³; lamR/L ≥ √y·w_J·lamreg with lamRreg = 1.39009 > 0, lamLreg = 1.43862 > 0 (asserted).
- TLf: MmL_up = −0.249982 < 0 → TLf·N^{3/2} = 2.8e-19229; TL0·N^{3/2} = 9.9e-2666376435 (monotonicity asserted).
- D: z-factorisation z = e^{iπ/4}ρζ̃₃e^E e^{mS} e^{om²invp2/2} e^{iφ} verified **exact** (ratio 1 ± 1.3e-57); Dreg = 11.44228.

**Validation:** W1's C-part 15.3849 matches step_38's W1·N^{1.5} at N=10⁹ exactly; W2: true 7.70 ≤ ours 9.30 ≤ step_38's 11.15; W3: true 28.85 ≤ ours 33.11 ≤ step_38's 46.81; pointwise trend C(N): 155.5 (400) → 83.3 (1e5) → 57.5 (1e7) → 51.2 (1e9) → limit ~49 ✓.

## Task 2 — [100,400] tightening (PROVED/arb)

29 cells ratio ≤ 1.05 (`step_18.py`): C rises 125.69 → **sup 194.3809** < 200, so step_34's midpoint+reach W3 port was **not needed** (conditional, skipped per instructions). Full table printed / `tight_100_400.json`.

## Task 3 — sweep verification (not redone)

43 cells **exactly** contiguous [400, 10¹⁰]; provenance spot-checked by re-running `certified_y38` on cells 0, 20, 41 — exact match to all digits (`step_19.py`). Archived: `sweep_cells_400_1e10.json`.

## Evidence table (arb enclosures, ctx.dps = 70)

| quantity | value | radius | digits |
|---|---|---|---|
| final cell C (J=192) | 49.56259792396024857830698314833209853198249410057901442521712151382263 | 2.0e-69 | 68 |
| final cell C (J=96) | 50.85487159316713887613912557631075219361998757785207777044262300667321 | 4.5e-69 | 68 |
| final cell C (J=48) | 53.61569012880916513058358542384421037118105513550664683468919612081652 | 2.6e-69 | 68 |
| sup [100,400] | 194.3808878592793 (arb upper) | — | — |
| sup [400,10¹⁰] | 617117.5535603258 (arb upper) | — | — |
| z-factorisation ratio | 1.0 | 1.3e-57 | 56 |
| μ-bound worst ratio | 0.7071067811865476 (= 1/√2) | NUMERICAL | — |

## Rerun commands

```
PY='<home>)/.venv/bin/python'
cd research/08_theory/D13h_final_cell
$PY step_38.py        # D13g base: point evals + 4 cells [100,400]           (<1 s)
$PY step_10.py        # FINAL CELL J=48 : C = 53.61569, all N >= 1e10        (<1 s)
$PY step_10_J96.py    # FINAL CELL J=96 : C = 50.85487                       (~1 s)
$PY step_10_J192.py   # FINAL CELL J=192: C = 49.56260                       (~2 s)
$PY step_18.py        # [100,400] ratio-1.05: 29 cells, sup 194.3809         (~1 s)
$PY step_19.py        # spot-check planner sweep cells 0,20,41 (exact match)
```

**Files:** step_38.py (D13g base), step_08.py/step_09.py (μ/logcf/z-identity checks), step_10.py, step_10_J96.py, step_10_J192.py (final cell), step_16.py/step_17.py (W2/W3 validation vs true integrals), step_18.py (tightening), step_19.py (sweep spot-check), step_20.py (junction + refinement), final_cell{,_J96,_J192}.json, tight_100_400.json, sweep_cells_400_1e10.json, result.md, progress.log.

**Labels:** THEOREM = PROVED (arb) modulo classical Binet formula; the μ-inequality |μ(z)| ≤ 1/(12|z|cos(arg z)) and all remainder/monotonicity bounds proved in-derivation and machine-checked. NUMERICAL (support only): pointwise C(N) trend → ~49. Optional future work (not done, not needed): ratio-1.05 cover of [400,10¹⁰] (~330 cells) would lower the middle sup to O(10²–10³); the theorem stands with **C = 617117.6**.