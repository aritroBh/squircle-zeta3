result = r'''# COR6X: independent referee report — true point symbols of Ω(c₀) at polygon corners

**Task.** Independently re-implement the true number-level symbol of Ω(c₀) (weight-4 part of R₃(c), the 56
weight-4 terms of `hpl_fit_k3_(0,0,0).json`), its δ₂,₂ projection, and the corner table; compare with COR6.
Own code only (`cor6x_core.py`), no COR6 imports.

## Machinery (cor6x_core.py), built from the definitions and validated

- **Regularized symbol of the number** G(a₁,…,a_n;x) = I(0;a₁,…,a_n;x): Goncharov recursion
  S(I) = Σᵢ S(I(omit aᵢ)) ⊗ ([prevᵢ−aᵢ] − [nextᵢ−aᵢ]) with the cyclic convention prev₁ = endpoint,
  next_n = start point; regularization aᵢ ↦ aᵢ+εⁱ, collision letters (ε) killed at the end.
  *Derived* from the total differential of iterated integrals (re-derived here by IBP) and checked
  **numerically** at weights 2–3 against direct quadrature of ∂I/∂(each parameter): agreement to 10⁻¹⁵.
- **H ↔ G:** H_w(x) = (−1)^{#ones} G(w;x), checked by direct quadrature on 8 words to 10⁻³⁰.
- **Controls (all pass):** S(Li_n(x)) = −[1−x]⊗[x]…[x]; H_{0,0,0,1}(x) = Li₄(x); δ₂,₂ kills Li₄(x), all
  shuffle products (Li₂², Li₂Li₂, log·Li₃, log²log²Li₂, log²log²); δ₂,₂ ≠ 0 on Li_{2,2}(x,y) (12 terms)
  and on Li_{3,1}(x,y) (2 terms); the five-term relation for Li₂ holds; 15/15 shuffle identities
  S(u)шS(v) = S(uшv) with letters {0,±1} (repeated-letter regularization test); symbol integrability
  (differential 2-form criterion) holds for every word.
- δ₂,₂(w₁w₂w₃w₄) = (w₁∧w₂)∧(w₃∧w₄) ∈ Λ²(Λ²V), letters expanded in a prime/S-unit basis.
- Function level: δ₂,₂(Ω)(c) = (1−1/c)([c]∧[1−c])∧([1−c]∧[1+c]) − (1+1/c)([c]∧[1+c])∧([1−c]∧[1+c]),
  **matching the paper's Prop. prop:M3angle exactly** (both wedge terms and the (1±1/c) coefficients
  from the c^e exponent groups). This validates the whole pipeline at function level.
- Value level (hpl_num.py, independent series code): M₃(2π/3) via R₃+Ψ+Φ matches the stored series to
  2·10⁻⁴⁸; the weight-4 accounting at c=1/√2 reproduces W₄ = 0.79578505983224081132… to all digits.

## The corner table (true point symbols, constant letters included)

δ₂,₂ of Ω(c₀), with the letter group of {2, c₀, 1−c₀, 1+c₀} presented exactly (every relation verified
as an exact field identity; letter-span dim certified by independence arguments + log-lattice LLL):

| point | field | deg | relations (verified) | dim U | δ₂,₂ | COR6 row | verdict |
|---|---|---|---|---|---|---|---|
| n=3, c₀=√3/2 | ℚ(√3) | 2 | 2²(1−c₀²)=1 | 3 | NONZERO −4([A∧U]∧[E∧U]+[A∧U]∧[A∧E]), A=[2],E=[3],U=[2+√3] | NONZERO, lspan 3 | **AGREE** |
| n=4, c₀=1/√2 | ℚ(√2) | 2 | 2(1−c₀²)=1, 2c₀²=1 | 2 | **ZERO** (letters ⊂ span{[2],[1+√2]}) | ZERO, lspan 2 | **AGREE** |
| n=5, c₀=sin π/5 | ℚ(ζ₂₀)⁺ | 4 | none | 4 | NONZERO: 2(A∧C)∧(B∧D) − 2(A∧D)∧(B∧D) − 2(B∧D)∧(C∧D) (3 comp.) | NONZERO, lspan 4 | **AGREE** |
| n=6, c₀=1/2 | ℚ | 1 | 2c₀=1, 2(1−c₀)=1 | 2 | **ZERO** (letters ⊂ span{[2],[3]}) | ZERO, lspan 2 | **AGREE** |
| n=7, c₀=sin π/7 | ℚ(ζ₂₈)⁺ | 6 | none | 4 | NONZERO (3 comp.) | NONZERO, lspan 4 | **AGREE** |
| n=8, c₀=sin π/8 | ℚ(ζ₁₆)⁺ | 4 | 2³c₀²(1−c₀²)=1 | 3 | NONZERO: −3(A∧C)∧(A∧D) − (A∧C)∧(C∧D) + 4(A∧D)∧(C∧D) | NONZERO, lspan 3 | **AGREE** |
| n=9, c₀=sin π/9 | ℚ(ζ₃₆)⁺ | 6 | none | 4 | NONZERO (3 comp.) | NONZERO, lspan 4 | **AGREE** |
| n=10, c₀=sin π/10 | ℚ(√5) | 2 | 2³c₀²(1+c₀)=1 | 3 | NONZERO | NONZERO, lspan 3 | **AGREE** |
| n=11, c₀=sin π/11 | ℚ(ζ₄₄)⁺ | 10 | none | 4 | NONZERO (3 comp.) | NONZERO, lspan 4 | **AGREE** |
| n=12, c₀=sin π/12 | ℚ(√2,√3) | 4 | 2⁴c₀²(1−c₀²)=1 | 3 | NONZERO: −4(A∧C)∧(A∧D) − (A∧C)∧(C∧D) + 5(A∧D)∧(C∧D) | NONZERO, lspan 3 | **AGREE** |
| α=π/6, c₀=cos π/12 | ℚ(√2,√3) | 4 | 2⁴c₀²(1−c₀²)=1 | 3 | NONZERO: −4, −1, +5 (same presentation) | NONZERO, lspan 3 | **AGREE** |
| c₀=1/3 | ℚ | 1 | (letter group ⊂ span{[2],[3]}) | 2 | **ZERO** | ZERO, lspan 2 | **AGREE** |
| c₀=2/5 | ℚ | 1 | none | 4 | NONZERO (12 comp. in prime basis) | NONZERO, lspan 4 | **AGREE** |
| c₀=3/7 | ℚ | 1 | none | 4 | NONZERO (12 comp. in prime basis) | NONZERO, lspan 4 | **AGREE** |

All 14 rows: field degree, letter-span dimension, and the true δ₂,₂ zero/nonzero verdict **AGREE** with
`cor6_table.json`. The naive (function-substituted, constants dropped) δ₂,₂ also agrees on every row
(ZERO exactly at n=4, n=6, c₀=1/3; NONZERO elsewhere). Component counts differ only by basis choice.
Run: `verify_cor6x.py` (own code, 51 s), recomputes the 5 higher-degree rows + the square control and
prints AGREE per row.

## Task 3 (Lie class at c₀=1/3) and the Ω(1/2) Lie class — DISAGREEMENT with the paper's combos

With the same (validated) machinery, the true point symbol of Ω(1/3) modulo shuffle products is **not**
the class of −(11/81)Li₄(−8) − 2Li₄(−3) + (20/3)Li₄(2/3), and Ω(1/2) is **not** the class of
−21Li₄(−2) + (13/2)Li₄(1/3) + (13/2)Li₄(2/3): in both cases the difference is a weight-4 tensor over
{2,3} with δ₂,₂ = 0 but **outside** the 13-dimensional shuffle-product subspace of ℚ^{2,3}(⊗4) (exact
linear algebra over ℚ). Numerics corroborate: Ω(1/2) minus the paper's five-Li₄ closed form
(= −51.777095579…) is not a combination of weight-4 product periods over ℤ[1/6] (LLL, 200 digits,
26–161-element product bases; no relation with coefficients < ~10⁴), so the paper's claimed
number-level closed form for Ω(1/2) appears to be wrong as well. These two Lie-class entries are the
only disagreements found; they do not affect the δ₂,₂ zero/nonzero table.

## Files
- `cor6x_core.py` — symbol machinery (regularized Goncharov recursion, wedges, δ₁₁/δ₂₂, expansion).
- `hpl_num.py` — independent HPL evaluator (Taylor series + trailing-zero shuffle reduction).
- `lll_rel.py` — flint-LLL integer relations (validated on the classical Li₃(1/2) identity).
- `step_02.py` (controls: Li_n, products, five-term, shuffles, Li_{2,2}, Li_{3,1}),
  `step_03.py` (HPL/quadrature conventions), `step_04.py` (rational points + Lie-class tests),
  `step_05a–g.py` (differential-formula + number-level LLL checks), `step_06a–e.py` (evaluator
  cross-checks, W₄/M₃ accounting, decisive LLL), `step_07–09.py` (cyclotomic table),
  `verify_cor6x.py` (final 5-row comparison, 51 s).
'''
open("<project>/research/08_theory/COR6X_independent_corner_symbols/result.md","w").write(result)
open("<project>/research/08_theory/COR6X_independent_corner_symbols/progress.log","a").write(
"verify_cor6x.py done (51s): 5 rows (n=5,8,12,cos(pi/12)) NONZERO + square ZERO; all AGREE with cor6_table.json; result.md written\n")
print("result.md and progress.log updated")
