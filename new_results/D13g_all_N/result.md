# D13g — Certified decay theorem: box machinery extended (regular TLf + monotone Ĝ); all‑N completion NOT finished this session

## Honest status

The target theorem for **all N ≥ N₀ is not yet proved**. What is proved this session: a stronger certified box assembly (`step_38.py`) with the TLf interval‑dependency **eliminated** and the Ĝ eta‑box inflation **eliminated**, giving certified cell constants on all of **[100,400]** (sup C = **930.4**) and on scattered ratio‑1.5 cells up to N = 1.5×10¹⁰. The contiguous dyadic cover [400, 10¹⁰], the sup‑C < 200 tightening of [100,400], and the final analytic cell N ≥ 10¹⁰ remain to be executed (design fully worked out below).

## 1. New fixes this session (each verified by code)

**Fix A — regular TLf (PROVED identity).** The D13f §4 interval dependency in `TLf` came from `ReX2 = −a·x2 + θ − (N/2)log(1+vf)` (cancellation of two O(θ) boxes). Derived and implemented the regular factored form, exact in (y, ρ, η₂), x2 = η₂/(ρy), N = π/(2y²), θ = π/(√2·y):

> ReX2 = θ·(M1 − L2), M1 = 1 − η₂/(2ρ), L2 = (√2/4)·ureg·logcf(u), ureg = √2ρ/η₂ + ρ²y/η₂², u = y·ureg, logcf(u) = log(1+u)/u (decreasing ⇒ tight box via endpoint point‑evals).

Verified against the old expression: identical values, residuals ≤ 9.2e‑59 (N=100), ≤ 1.6e‑58 (N=400), ≤ 1.6e‑55 (N=10⁶) at ctx.dps=60. Effect on the D13f cells (TLf·N^1.5): 786.5→**150.1**, 1232.7→**113.8**, 2192.6→**78.0**, 16.9→**3.6**.

**Fix B — monotone Ĝ box (PROVED monotone, D13f).** `Gh(η,ρ,y)` is decreasing in η; naive interval evaluation over the wide η‑boxes (ηR ⊃ [1.08, 1.33]) inflated `Ĝ.upper` (0.58→0.81) and killed `lamR`'s sign for ratio‑1.5 cells at N ≥ 10⁶. Implemented `Ghbox` via endpoint point‑evals (asserts monotonicity each call). lamR now passes at N = 10⁶, 10⁸, 10¹⁰.

**Pitfall documented (NUMERICAL/diagnostic):** this environment's python‑flint 0.9.0 renders wide arb boxes misleadingly via `str` (e.g. `[1.0826,1.3274]` prints as `[1e+0 ± 0.328]`); `.mid()`/`.rad()` confirm the arithmetic is exact. Do not trust `str()` for wide boxes.

## 2. Certified cell constants (PROVED — arb enclosures, all sign asserts passed: om>0, c2>0, c2R/L>0, lamR/L/2>0, den>0, ex0<0, xl>0, M1−L2<0, Ĝ‑mono)

| N‑box | C (certified upper) | | N‑box | C |
|---|---|---|---|---|
| [100,150] | **679.5** | | [400,600] | **1270.0** |
| [150,225] | **754.2** | | [600,900] | **1701.4** |
| [225,340] | **930.4** | | [10⁶,1.5·10⁶] | **11332.6** |
| [340,400] | **340.5** | | [10⁸,1.5·10⁸] | **79184.0** |
| | | | [10¹⁰,1.5·10¹⁰] | **756712.0** |

Piece breakdown [100,150]·N^1.5: W1 78.7, W2 18.3, W3 478.0, TM 62.9, TR 63.0, TLm 65.1, TLf 150.1, TL0 0.0, D 129.6. W3 dominates everywhere (crude full‑segment Taylor V + q‑factor), and grows ~N^1/4·polylog on fixed‑ratio cells — this, not TLf, is now the tightening target.

**Theorem proved this session (PROVED):** for all N ∈ [100,400] (k = N−1 ∈ [99,399]),
> |M_k^main + (√π/N)·e^{−√(πN)}·cos√(πN)| ≤ **930.4** · N^{−3/2}·e^{−√(πN)}.

**Truth comparison (from D13f_sweep/mk_deg9.tsv, arb enclosures of M_k^main):** true C(N) = 0.594 (N=100), 1.530 (150), 1.597 (200), 1.127 (300), 1.065 (400); **max over [100,400] = 1.64712094341739875… at N = 114** (±4.1e‑60). Certified 930.4 ≫ true 1.65 ✓ consistent.

## 3. Evidence table

| check | result | residual/margin | digits |
|---|---|---|---|
| ReX2 regular form vs old (N=100,400,10⁶) | exact identity | ≤ 1.6e‑55 | 60 |
| Ĝ monotone endpoint box (assert in `Ghbox`) | passed every cell | — | 70 |
| All step_38 sign asserts (9 cells run) | passed | — | 70 |
| TLf·N^1.5 on [100,400]: before→after | 786/1233/2193/17 → 150/114/78/3.6 | 5–28× tighter | 70 |
| sup C on [100,400]: before→after | 2281.28 → **930.4** | 2.45× tighter | 70 |
| lamR at ratio 1.5, N=10⁶/10⁸/10¹⁰ | was FAIL → now passes | — | 70 |
| true max C on [100,400] | 1.6471209… @ N=114 | ±4.1e‑60 | 60 |

## 4. Remaining work (design complete, unexecuted — CONJECTURE only that these suffice)

1. **Dyadic sweep [400,10¹⁰]:** cells run in <0.05 s; use an adaptive driver: try ratio 1.5, halve the log‑ratio on assert failure or C > 3000. Expected ~10³ subcells (binding constraints: lamR needs √r−1 ≲ d/x0 ~ N^{−1/4}√log N; M1−L2 < 0 needs ε_y ≲ 0.1).
2. **[100,400] → sup C < 200:** rerun step_38 on ratio‑1.05 N‑cells (~29 cells) — TLf already fixed; W3's q‑factor exponent (c2up−c2side)u² shrinks ~5×; if still > 200, port step_34's midpoint+reach W2/W3 (sharp values: C(100)=42.4247, C(400)=46.5448 show headroom).
3. **Final analytic cell y∈[0,y_J]:** every piece is exactly y³·(regular factor); implement factored forms: Phi3b = πρ⁴y²(ρ²y²+3ηρy+3η²)/(η³(η²+√2ηρy+ρ²y²)^{3/2}), Phi4b similarly ~y³; G1 = 6ρ²y²/η²+(1.5√2+1)ρy/η, G2 = 6√2ρ³y³/η³+4ρ²y²/η²; g0 via the h₁ Binet bound |g| ≤ (√π/4)·2^{3/2}(ρy)^{5/2}e^{ReE}, ReE ≤ −3√2ρy/8 + 9ρ²y²/2 + ρy/6; W‑integrals by v = u/x0 moments with K2 = c2reg/(ρ²y); tails via w² = log N + c_w with e^{−c2R d²} = (√(2/π)·y)^{c2R/c2}e^{−c_w·c2R/2c2} (bounded box power); TLf/TL0 by manual monotonicity of y^{−3}e^{K/y} and e^{ex0·N}N^{3/2}/(N−1) (decreasing for y_J < |K|/3).

## 5. Files (in `research/08_theory/D13g_all_N/`)

- **step_38.py** — the box assembly to use: step_37 + Fix A (regular TLf, tight logcf) + Fix B (`Ghbox`); `certified_y38(y1,y2,J=24,nsub,c_w)`; `__main__` reproduces the four [100,400] cells.
- step_37.py, step_34.py — copied references (D13f).
- Tool‑call scripts step_01.py…step_20.py: step_11 (truth constants), step_13 (ReX2 identity check), step_14–16/19 (step_38 construction + cell runs), step_17–18 (flint str diagnostic), step_20 (final sweep probes).

**Rerun:** `cd research/08_theory/D13g_all_N && PY step_38.py` (the four cells, < 1 s); dyadic probes: `certified_y38(sqrt(pi/(2·n2)), sqrt(pi/(2·n1)), nsub=2)`.

**Bottom line:** the TLf and Ĝ box pathologies are fixed and verified; the certified sup on [100,400] is 930.4 (was 2281.28); ratio‑1.5 cells now certify at every scale up to 1.5×10¹⁰ (C ≤ 756712, weak). The contiguous dyadic cover, the < 200 tightening, and the degenerate final cell — hence the all‑N ≥ N₀ theorem — are designed but **not executed/proved this session** (turn limit hit during the flint‑str diagnostic).