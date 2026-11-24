**Bug found and fixed; full 14-row agreement restored.** The symbol map is linear over the field, so the weight-0 prefactors c₀^e act as scalars on the symbol (they do not drop out). My earlier disagreement was this c^e-scaling error, not the symbol recursion (which passed every control: differential formula to 10⁻¹⁵, shuffle identities, five-term relation, Liₙ symbols, function-level δ₂,₂ matching the paper's Prop. M3angle exactly).

Per-row verdicts (true number-level δ₂,₂ of Ω(c₀), constant letters included, c₀^e scaling applied; `verify_cor6x.py`, 141 s):

| point | deg | dim U | my δ₂,₂ | COR6 | verdict |
|---|---|---|---|---|---|
| n=3 (sin π/3) | 2 | 3 | NONZERO | NONZERO | AGREE |
| n=4 (1/√2, zero control) | 2 | 2 | ZERO | ZERO | AGREE |
| n=5 (sin π/5) | 4 | 4 | NONZERO | NONZERO | AGREE |
| n=6 (1/2) | 1 | 2 | ZERO | ZERO | AGREE |
| n=7 (sin π/7) | 6 | 4 | NONZERO | NONZERO | AGREE |
| n=8 (sin π/8) | 4 | 3 | NONZERO | NONZERO | AGREE |
| n=9 (sin π/9) | 6 | 4 | NONZERO | NONZERO | AGREE |
| n=10 (sin π/10) | 2 | 3 | NONZERO | NONZERO | AGREE |
| n=11 (sin π/11) | 10 | 4 | NONZERO | NONZERO | AGREE |
| n=12 (sin π/12) | 4 | 3 | NONZERO | NONZERO | AGREE |
| α=π/6 (cos π/12) | 4 | 3 | NONZERO | NONZERO | AGREE |
| c₀=1/3 | 1 | 2 | ZERO | ZERO | AGREE |
| c₀=2/5 | 1 | 4 | NONZERO | NONZERO | AGREE — values match exactly |
| c₀=3/7 | 1 | 4 | NONZERO | NONZERO | AGREE — values match exactly |

Controls: δ₂,₂(Li₄(x)) = 0, δ₂,₂(Li₂²) = 0, square zero control passes. **Task 3 resolved:** with the fix, the paper's Lie classes hold exactly mod shuffles — Ω(1/2) ≡ −21Li₄(−2) + (13/2)Li₄(1/3) + (13/2)Li₄(2/3) (decomposable difference, exact), and Ω(1/3) ≡ −(11/81)Li₄(−8) − 2Li₄(−3) + (20/3)Li₄(2/3) (decomposable difference, exact); the Ω(1/2) number-level closed form (five Li₄'s + COR5's products) confirmed with own code to 1.5·10⁻⁷². My earlier Task-3 DISAGREE is retracted (cause: the c^e prefactors were wrongly dropped instead of applied as field scalars). Deliverables in `new_results/COR6X_independent_corner_symbols/`: `cor6x_core.py`, `hpl_num.py`, `lll_rel.py`, `step_02`–`step_14.py`, `verify_cor6x.py`, updated `result.md`, `progress.log`.