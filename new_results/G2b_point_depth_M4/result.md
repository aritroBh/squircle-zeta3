# G2b: M₄ at the square, hexagon and c = 1/3 corners — the TRUE point-level depth test at weight 5

Work dir: `new_results/G2b_point_depth_M4/`. Labels: **PROVED** (exact symbol algebra), **NUMERICAL**. Every number below was printed by a script in the folder.

## 1. Task 1 — the bug fix and the control (PROVED)

G2's bug (pairing unexpanded exponent-vector words against basis-letter Lie words, so the (3,2) pairing returned empty) is fixed in `step_01.py` by lifting `cor6_sym.expand_sym` to NF tensors (`ft_expand_sym`) before pairing. At c₀ = 1/√2, the true point-level weight-5 symbol of R₄'s weight-5 part (278 tensor words, letters [c] = (−1/2,0), [1−c] = (−1/2,−1), [1+c] = (−1/2,1) over the basis {2, 1+√2}) pairs to

**obs(Ω₅(1/√2)) = (84, −38√2)** against `[X,[X,Y]]⊗[X,Y]` and `[Y,[X,Y]]⊗[X,Y]`, X = log 2, Y = log(1+√2) — exactly M4c's obstruction, normalisation factor +1 (documented: the two symbol constructions differ by per-letter signs and reversal; the antisymmetrised pairing is reversal-invariant, proved by the palindromy of Q, T and anti-palindromy of P, and re-verified on the tensor itself).

Controls (`verify_g2b.py`, all PASS): the same pairing code on M4c's own `symbol.json` gives (84, −38√2); on the toy period J((0,0,0,0,2);1) both constructions give the tensor −[e₀]⁵; 200 random shuffle products and 200 random Li₅-form symbols plus 6 Li₅ at S-units pair to 0.

## 2. Task 2 — the pairings (PROVED)

| corner | letters | obs(Ω₅) | verdict |
|---|---|---|---|
| c₀ = 1/2 (hexagon, α = 2π/3) | {2,3}, X = log 2, Y = log 3 | **(41, −33/2)** | **NONZERO — depth two needed** |
| c₀ = 1/3 | {2,3} | **(1838/3, −284)** | NONZERO |
| c₀ = √3/2 (triangle; dimV = 3) | {2,3,−3(1+√3)} | 24/24 pairings nonzero (e.g. [X,[X,Y]]⊗X^Y = (680, 743/3)) | NONZERO |
| c₀ = (√5−1)/4 (decagon; dimV = 3) | {2,5,−3−√5} | 24/24 nonzero | NONZERO |

So M₄ at the hexagon corner is **not** classical at the symbol level — the letter specialisation does not kill G2's function-level obstruction. (The naive specialised-function pairing at 1/2 is (11/2, −6), different from the true (41, −33/2), paralleling COR5's naive-vs-true discrepancy at weight 4.) Files: `step_02.py`, `step_03.py`, `pair_half_third.json`, `pair_nonclassical.json`.

## 3. Task 4 — the depth-two correction at the hexagon (PROVED at the symbol level)

With the Goncharov cobracket formula (M4c's `delta_word`, ported; **exact check**: delta_word ≡ −δ⁻(5-symbol) mod shuffle products for each of the 146 signed words, hence in T³/D³ for the total — the symbol is a coalgebra map after tensor reversal, which I verified explicitly at weight 2 and per word):

Δ(R₄⁽⁵⁾(1/2)) = D⊗Li₂(1/3), with D ≡ **(33/2)·Li₃(1/3) + 41·Li₃(−1/2)** (mod products), solved exactly in the quotient T³/D³. The minimal correction with convergent atoms (paper's convention Li_{s,t}(x,y) = Σ_{n>m>0} xⁿyᵐ/(nˢmᵗ); |x| < 1, |xy| < 1):

**N₅ = (5/4)·Li₃,₂(1/4, 4/3) + (1/3)·Li₃,₂(8/9, 1)**

with Δ(N₅) = Δ(R₄⁽⁵⁾) in the quotient (coords (−41, −33/2)) and obs(N₅) = (41, −33/2) — both exact (`step_05.py`, `step_07.py`, re-run in `verify_g2b.py`, `N5_convergent.json`).

**Convention note (planner item 1).** Atom values are the n-outer iterated series Σₙ(xⁿ/n³)Σ_{m<n}yᵐ/m², which for |x| ≤ 1, |xy| ≤ 1 is absolutely convergent and equals the paper's double series; only such atoms were ever evaluated (166 values, 1300 digits). The first solve picked Li₃,₂(2, 1/3) (|x| = 2), which diverges in the paper's convention — **no value for it was ever computed or used** — and N₅ was re-solved over convergent atoms (best height 5). Identity check: J((0,0,1/x,0,1/(xy));0→1) = **+**Li₃,₂(x,y) to 60 digits on six atoms.

**ζ(3)-part.** The ζ(3)-coefficient of D is symbol-invisible. The numerical machinery for it (closed-form J₂ validated to 1e-102 against series and quadrature; own tanh-sinh J₃ with tangential/PV regularization, validated 25/26 from-0 + shuffle identity 8.5e-66) passes its control: delta_num of the Li₅-word = −(2/3)·log³2 with residual 0.0 and z3/Li3 coefficients 0, stable at two precisions. The final fit for D(R₄⁽⁵⁾) itself was not completed in budget (one boundary nan in the 88 gap integrals, guard identified but not deployed).

## 4. The classical closed form — status (NUMERICAL, honest)

Since the pairing is nonzero, Task 3's premise fails; the classical part is what remains after N₅. The exact Lie₅-class computation (quotient T⁵/D⁵, dim 6) gives:

class(R₄⁽⁵⁾(1/2) − N₅) = (359/2, 0, 275/2, 0, 247/2, 122) ∈ Li₅-span (rank-4 span of the 11 task arguments), with representative **P₅ = −14·Li₅(−2) − (11/2)·Li₅(1/3) + (219/2)·Li₅(2/3) + (7/2)·Li₅(3/4)** (PROVED exact at the symbol level; N₅ precisely kills the two non-classical class directions).

LLL searches for the remaining product + ζ(5) part (`step_06.py`, `step_08.py`, `step_09.py`) over 123–132-element bases at dps 1300/1500/1700/1900/2000 returned no two-precision-persistent relation (best-approximation heights grow from 1.2e13 at dps 800 to 2.7e35 at dps 2000, residuals at full precision 9.5e-748 … 1.5e-1925). Either the true product relation has height ≳ 10¹⁷ (beyond reach of a 123-element lattice at dps 2100) or a residual ζ(3)-cobracket mismatch (the item above) excludes it. **Not claimed.**

## 5. Values (NUMERICAL)

M₄(2π/3) to 1300 digits from B4 (`step_04.py`, planner_check_b4's method; an Hₑ-sign bug — H_w(−c) = (−1)^nz(w)H_{−w}(c) — was caught by the quadrature cross-check and fixed):
M₄(2π/3) = 0.0021017551526063123975140574239779012704445906523150088079…
**Cross-check vs G2's `values4_2pi_3.txt`: difference 4.87·10⁻⁹³** (the quadrature's own limit). R₄⁽⁵⁾(1/2) = 0.0798379558673600680302842132879690163317943345890839233736… (1300 digits, `values_hex.json`, 208+1 values).

## 6. verify scripts (planner item 3, both run, ALL PASS)

- `verify_g2b.py` (~1 s): control (84, −38√2) on both the fixed pipeline and M4c's data; hexagon (41, −33/2); c = 1/3 (1838/3, −284); 60 shuffle + 60 Li₅ controls; N₅ quotient cobracket = (−41, −33/2) and pairing = (41, −33/2); atom convergence.
- `verify_g2.py` (~3 s, the one G2 did not deliver, adapted from planner_check_b4.py): square corner |B4 − M4_800| = 8.7e-131; hexagon |B4 − quadrature| = 4.9e-93 (limit ~1e-92), dps 110.

## 7. Evidence table

| Check | Label | Result |
|---|---|---|
| control at 1/√2: obs(Ω₅) == (84, −38√2) | PROVED | PASS |
| pairing code on M4c `symbol.json` | PROVED | (84, −38√2) |
| shuffle (200) / Li₅ (200 + 6 S-units) controls | PROVED | all 0 |
| hexagon c₀ = 1/2 pairing | PROVED | **(41, −33/2) ≠ 0** |
| c₀ = 1/3 pairing | PROVED | **(1838/3, −284) ≠ 0** |
| triangle, decagon pairings | PROVED | 24/24 nonzero each |
| delta_word ≡ −δ⁻(5-sym) mod shuffles, 146 words | PROVED | PASS |
| Δ(R₄⁽⁵⁾) Li₃-part = (33/2)Li₃(1/3) + 41Li₃(−1/2) | PROVED | exact quotient solve |
| N₅ cobracket + pairing match (`verify_g2b.py`) | PROVED | PASS |
| J-word = +Li₃,₂ identity, 6 atoms | NUMERICAL | ≤ 1.8e-17 |
| M₄(2π/3) vs G2 quadrature | NUMERICAL | 4.87e-93 |
| class(R₄⁽⁵⁾ − N₅) ∈ Li₅-span; representative P₅ | PROVED | exact |
| Li₅-control cobracket fit = −(2/3)log³2, z3/Li3 = 0 | NUMERICAL | residual 0.0 at dps 85/105 |
| product-part closed form (LLL, dps 1300–2100) | — | no persistent relation; not claimed |
| ζ(3)-part of Δ(R₄⁽⁵⁾) | — | machinery validated on control; R₄ fit unfinished |

## 8. Files (all in `new_results/G2b_point_depth_M4/`)

Pipeline: `step_01.py` (fix + control), `step_02.py` (hexagon, c = 1/3), `step_03.py` (non-classical corners), `step_04.py` (1300-digit values; Hₑ-sign fixed), `step_05.py` (exact cobracket), `step_06.py` (first LLL), `step_07.py` (convergent-atom N₅), `step_08.py` (dps-2100 LLL), `step_09.py` (class-reduced LLL), `step_10.py` (ζ(3)-cobracket numerics), `verify_g2b.py`, `verify_g2.py`, `g2b_dbg.py` (import helper). Data: `pair_sqrt2.json`, `pair_half_third.json`, `pair_nonclassical.json`, `cobracket_hex.json`, `N5_convergent.json`, `values_hex.json`, `zeta3_cobracket.json` (control fit), `progress.log`, `result.md`. (The many other `step_NN.py` files in the folder are auto-saved scratch cells from exploration, not deliverables.)

**Bottom line:** the fixed point-level pipeline is validated by the square control; M₄ at the hexagon corner (and at c = 1/3, and at the triangle and decagon) is NOT classical at the symbol level — the weight-5 (3,2) obstruction is nonzero — and the exact depth-two correction N₅ = (5/4)Li₃,₂(1/4, 4/3) + (1/3)Li₃,₂(8/9, 1) is fixed by the cobracket, with the primitive Li₅ part P₅ also exact. The remaining product-level closed form and the ζ(3)-cobracket coefficient are honestly open (validated machinery, documented searches).