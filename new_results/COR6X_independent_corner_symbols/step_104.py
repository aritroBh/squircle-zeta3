result = r'''# COR6X: independent referee report — true point symbols of Ω(c₀) at polygon corners (FINAL)

**Task.** Independently re-implement the true number-level symbol of Ω(c₀) (weight-4 part of R₃(c), the 56
weight-4 terms of `hpl_fit_k3_(0,0,0).json`), its δ₂,₂ projection, and the corner table; compare with COR6.
Own code only (`cor6x_core.py`), no COR6 imports.

## Machinery (cor6x_core.py), built from the definitions and validated

- **Regularized symbol of the number** G(a₁,…,a_n;x) = I(0;a₁,…,a_n;x): Goncharov recursion
  S(I) = Σᵢ S(I(omit aᵢ)) ⊗ ([prevᵢ−aᵢ] − [nextᵢ−aᵢ]) with the cyclic convention prev₁ = endpoint,
  next_n = start point; regularization aᵢ ↦ aᵢ+εⁱ, collision letters (ε) killed at the end.
  *Derived* from the total differential of iterated integrals (re-derived here by IBP) and checked
  **numerically** at weights 2–3 against direct quadrature of ∂I/∂(each parameter): agreement to 10⁻¹⁵.
- **H ↔ G:** H_w(x) = (−1)^{#ones} G(w;x), checked by direct quadrature on all 24 weight-4 word shapes to 10⁻²⁶.
- **Controls (all pass):** S(Li_n(x)) = −[1−x]⊗[x]…[x]; H_{0,0,0,1}(x) = Li₄(x); δ₂,₂ kills Li₄(x), all
  shuffle products; δ₂,₂ ≠ 0 on Li_{2,2}(x,y) and on Li_{3,1}(x,y); the five-term relation for Li₂ holds;
  15/15 shuffle identities S(u)шS(v) = S(uшv); symbol integrability holds for every word.
- δ₂,₂(w₁w₂w₃w₄) = (w₁∧w₂)∧(w₃∧w₄) ∈ Λ²(Λ²V), letters expanded in a prime/S-unit basis.
- Function level: δ₂,₂(Ω)(c) = (1−1/c)([c]∧[1−c])∧([1−c]∧[1+c]) − (1+1/c)([c]∧[1+c])∧([1−c]∧[1+c]),
  **matching the paper's Prop. prop:M3angle exactly**.
- Value level (hpl_num.py, independent series code): M₃(2π/3) via R₃+Ψ+Φ matches the stored series to
  2·10⁻⁴⁸; the weight-4 accounting at c=1/√2 reproduces W₄ = 0.79578505983224081132… to all digits.

## The c^e bug and its fix

The 56 terms are coeff·c^e·H_w(c). A first version of this project dropped the c₀^e prefactors from the
symbol ("weight 0 drops out") — **wrong**: the symbol map is **linear over the field**, so weight-0
constants c₀^e act as *scalars*: S(c₀^e·X) = c₀^e·S(X), they do not drop out. With the correct scaling
every check now passes. (COR6's pipeline had the scaling right; the paper's δ₂,₂ zero/nonzero table is
unchanged by the fix, but the δ₂,₂ *values* at rational points change and then match COR6 exactly.)

## The corner table (true point symbols, constant letters included, c₀^e scaling applied)

δ₂,₂ of Ω(c₀), letter group of {2, c₀, 1−c₀, 1+c₀} presented exactly (every relation verified as an exact
field identity; span dim certified by independence arguments + LLL):

| point | field (deg) | relations (verified) | dim U | δ₂,₂ | COR6 | verdict |
|---|---|---|---|---|---|---|
| n=3, c₀=√3/2 | ℚ(√3) (2) | 2²(1−c₀²)=1 | 3 | NONZERO (3 comp.) | NONZERO, lspan 3 | **AGREE** |
| n=4, c₀=1/√2 | ℚ(√2) (2) | 2(1−c₀²)=1, 2c₀²=1 | 2 | **ZERO** | ZERO, lspan 2 | **AGREE** |
| n=5, c₀=sin π/5 | ℚ(ζ₂₀)⁺ (4) | none | 4 | NONZERO (6 comp.) | NONZERO, lspan 4 | **AGREE** |
| n=6, c₀=1/2 | ℚ (1) | 2c₀=1, 2(1−c₀)=1 | 2 | **ZERO** | ZERO, lspan 2 | **AGREE** |
| n=7, c₀=sin π/7 | ℚ(ζ₂₈)⁺ (6) | none | 4 | NONZERO (6 comp.) | NONZERO, lspan 4 | **AGREE** |
| n=8, c₀=sin π/8 | ℚ(ζ₁₆)⁺ (4) | 2³c₀²(1−c₀²)=1 | 3 | NONZERO (3 comp.) | NONZERO, lspan 3 | **AGREE** |
| n=9, c₀=sin π/9 | ℚ(ζ₃₆)⁺ (6) | none | 4 | NONZERO (6 comp.) | NONZERO, lspan 4 | **AGREE** |
| n=10, c₀=sin π/10 | ℚ(√5) (2) | 2³c₀²(1+c₀)=1 | 3 | NONZERO (3 comp.) | NONZERO, lspan 3 | **AGREE** |
| n=11, c₀=sin π/11 | ℚ(ζ₄₄)⁺ (10) | none | 4 | NONZERO (6 comp.) | NONZERO, lspan 4 | **AGREE** |
| n=12, c₀=sin π/12 | ℚ(√2,√3) (4) | 2⁴c₀²(1−c₀²)=1 | 3 | NONZERO (3 comp.) | NONZERO, lspan 3 | **AGREE** |
| α=π/6, c₀=cos π/12 | ℚ(√2,√3) (4) | 2⁴c₀²(1−c₀²)=1 | 3 | NONZERO (3 comp.) | NONZERO, lspan 3 | **AGREE** |
| c₀=1/3 | ℚ (1) | (letters ⊂ span{[2],[3]}) | 2 | **ZERO** | ZERO, lspan 2 | **AGREE** |
| c₀=2/5 | ℚ (1) | none | 4 | NONZERO (15 comp., **values match COR6 exactly**) | NONZERO, lspan 4 | **AGREE** |
| c₀=3/7 | ℚ (1) | none | 4 | NONZERO (12 comp., **values match COR6 exactly**) | NONZERO, lspan 4 | **AGREE** |

All 14 rows: field degree, letter-span dimension, and the true δ₂,₂ zero/nonzero verdict **AGREE** with
`cor6_table.json`; at c₀ = 2/5 and 3/7 the full δ₂,₂ values match component-for-component. The naive
(function-substituted) δ₂,₂ also agrees on every row (ZERO exactly at n=4, n=6, c₀=1/3; NONZERO elsewhere).

## Task 3 (Lie classes) — RESOLVED

With the correct c₀^e scaling, the paper's Lie classes hold **exactly** (modulo shuffle products, exact
ℚ-linear algebra):
- Ω(1/2) ≡ −21Li₄(−2) + (13/2)Li₄(1/3) + (13/2)Li₄(2/3) (mod shuffles): **TRUE (decomposable difference)**.
- Ω(1/3) ≡ −(11/81)Li₄(−8) − 2Li₄(−3) + (20/3)Li₄(2/3) (mod shuffles): **TRUE (decomposable difference)**.
The number-level identity Ω(1/2) = −21Li₄(−2) − (65/8)Li₄(1/3) + (13/2)Li₄(2/3) + (137/6)Li₄(1/2) +
(39/64)Li₄(1/9) + (COR5's product terms) was also confirmed numerically with own code to 1.5·10⁻⁷².
My earlier "DISAGREE" on these two entries was caused by the c^e bug above and is hereby retracted.

## Files
- `cor6x_core.py` — symbol machinery; `hpl_num.py` — independent HPL evaluator; `lll_rel.py` — flint LLL.
- `step_02.py` (controls), `step_03.py` (HPL/quadrature conventions), `step_04.py` (rational points),
  `step_05a–g.py` (differential formula + number-level checks), `step_06a–e.py` (evaluator + W₄/M₃ checks),
  `step_07–09.py` (cyclotomic table, drop-ce version), `step_10–13b.py` (bug hunt), `step_14.py` (fixed table),
  `verify_cor6x.py` (final 14-row comparison + controls + Lie-class confirmations, 141 s).
'''
open("<project>/research/08_theory/COR6X_independent_corner_symbols/result.md","w").write(result)
open("<project>/research/08_theory/COR6X_independent_corner_symbols/progress.log","a").write(
"verify_cor6x.py (fixed c^e scaling, 141s): all 14 rows AGREE with cor6_table.json; values match at 2/5,3/7; Lie classes at 1/2 and 1/3 confirmed; Task-3 DISAGREE retracted\n")
print("done")
