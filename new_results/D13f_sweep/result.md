# D13f — Certified decay theorem: the arb box-gap is CLOSED; certified C on [100,400] (cells) + sharp points; dyadic/final-cell runs pending

## Bottom line

**(1) The D13e gap is fixed and verified.** The blow-up was isolated this session to a genuine **arb/FLINT artifact**: complex `digamma` over a box of input-width 0.035 returns an enclosure of width **±103** (≈2×10⁴× overestimate; `lgamma` similar, ±1.07). This poisoned `gapF`/`logg` (→ W1/W2/W3) and `gF` (→ g0m, q-factor, D). In `step_37.py` all complex digamma/lgamma/gamma **box** evaluations are eliminated and `λR` keeps its sign via the regular Ĝ-form. Every sign assertion (`om>0`, `c2>0`, `c2R/c2L>0`, `lamR/lamL/lam2>0`, `den>0`, `ex0<0`, `xl>0`) **passed on every cell**, so every printed C below is a valid arb enclosure.

**(2) Certified this session (PROVED, arb enclosures):**
- Sharp single-N (step_34, reproduced): **C(100) = 42.4247, C(400) = 46.5448**.
- Box assembly (step_37): points C(100) = 99.99, C(400) = 152.68; **cells: [100,150] → 1082.53, [150,225] → 1470.06, [225,340] → 2281.28, [340,400] → 345.68**.
- Hence **for all N ∈ [100,400]: |M_k^main + (√π/N)e^{−√(πN)}cos√(πN)| ≤ 2281.28·N^{−3/2}e^{−√(πN)}** (sup of the 4 cell constants). The cell constants are weak (W3 Taylor-crudeness + TLf interval-dependency, see §4) but **certified**.

**(3) NOT finished (turn limit):** the dyadic sweep 400→10¹⁰ and the final analytic cell y*∈[0,y_J], so "for all N ≥ N₀" remains open; the mk_deg9.tsv table check was not run (prior D13d numerics: truth ≈ 0.60–1.60·sc ≪ 2281·sc ✓ consistent).

## 1. The fixes (each derived here + verified at 50 digits)

- **γ real bounds** (replaces all box digamma/lgamma): from ψ(w)−ψ(w+3/2) = −(3/2)Σ_{n≥0} 1/((n+w)(n+w+3/2)) and |n+w| ≥ n+t, w = t(1+i), t = x/(2√2): **|γ′(x)| ≤ (3/4)ψ′(t)+1/x ≤ (3/4)(1/t²+1/t)+1/x**; writing 1/(n+w)²−1/(n+w+3/2)² = ∫₀^{3/2}2/(n+w+r)³dr: **|γ″(x)| ≤ (3/8)|ψ″(t)|+1/x² ≤ (3/8)(2/t³+1/t²)+1/x²**. Verified: x=7.9: |γ′|=0.308 ≤ 0.491; |γ″|=0.038 ≤ 0.0985.
- **λR regular form**: −ReΦ′(ηx₀,N) = **a − (πρ/2)Ĝ(η)**, Ĝ(η) = ρ(η+t√2)/(√2η(η²+√2ηt+t²)), t = ρy — pure algebraic substitution of N = π/(2y²), x = η/(ρy) (no cancellation, all terms positive; Ĝ decreasing). Verified: residual ≤ 7×10⁻¹⁸ over N=100…10⁸, η∈[0.5,2].
- **curv2** = |ReΦ″(ηx₀)|/2 (regular): verified ratio = 1.0 exactly.
- **|g(x)| strictly decreasing on x>0** (new, used for the W3 q-factor): Re γ′(x) = −(3/4)Σ Re[b/((n+w)(n+w+3/2))] − 1/x < 0 since each Re[…] = (A−t²+tB)/(√2D) > 0 (A = (n+t)(n+t+3/2) > t², B = 2(n+t)+3/2). Numerically confirmed (0.816→0.0679→0.00668→0.000122 for x=1,3,7.9,40).
- **|g(x)| ≤ Hr(x/√2)/x** via the beta integral Γ(w)/Γ(w+3/2) = (2/√π)∫₀¹t^{w−1}(1−t)^{1/2}dt (ratio 0.60–0.69, tight).
- **|Φ‴| ≤ Phi3b, |Φ⁗| ≤ Phi4b** (proved via (1+s)^m − s^m = Σ\binom{m}{k}s^k and |1+s|² = 1+√2x+x²; both monotone decreasing in x — used for box Taylor bounds M3/M4).
- **g0 complex enclosure by propagation**: g0box = g(x_ref)·exp(ball(r_g)), r_g = G1box·hx (x_ref = box midpoint) — used for D and (with min against the beta bound) for P0t.
- **h₁ two-sided bound derived for the final cell** (not yet implemented): Γ(z)/Γ(z+3/2) = z^{−3/2}e^{E(z)}, E(z) ∈ (−3/(8z) − 9/(8z³), −3/(8z) + 9/(8z²) + 1/(12z)) via Binet 0 < μ(z) < 1/(12z); regular in (y,σ) since 1/z = 2√2ρy.

## 2. step_37 results (arb intervals, ctx.dps=70; pieces normalized by N^{−3/2}e^{−θ})

| box | W1 | W2 | W3 | TM | TR | TLm | TLf | D | **C (certified)** |
|---|---|---|---|---|---|---|---|---|---|
| N≈100 (pt) | 16.16 | 5.79 | 91.27 | 12.14 | 8.13 | 16.89 | 5.71 | 2.69 | **99.99** |
| N≈400 (pt) | 15.50 | 11.08 | 190.29 | 4.04 | 7.31 | 9.16 | 0.47 | 1.99 | **152.68** |
| [100,150] | 78.7 | 18.3 | 478.0 | 62.9 | 63.0 | 65.1 | 786.5 | 129.6 | **1082.53** |
| [150,225] | 73.8 | 22.9 | 637.0 | 44.1 | 59.2 | 59.2 | 1232.7 | 143.5 | **1470.06** |
| [225,340] | 72.4 | 29.5 | 918.2 | 32.5 | 61.9 | 53.8 | 2192.6 | 168.2 | **2281.28** |
| [340,400] | 27.1 | 16.7 | 376.5 | 8.7 | 14.3 | 18.1 | 16.9 | 43.2 | **345.68** |

Certified interval uppers: C(100 pt) = [99.98994796256576624268362286701867512610352200355638096612116249128805 ± 1.91e-69]; C(400 pt) = [152.6802475086496487455891799400063328271125906692597037540867665109418 ± 4.90e-68]; cells: 1082.53051… ± 2.98e-67, 1470.05952… ± 4.22e-67, 2281.28071… ± 3.97e-67, 345.68338… ± 1.83e-68. (step_34 sharp points: [42.42472910… ± 8.14e-64], [46.54477870… ± 5.57e-63], 62 digits.)

## 3. Evidence table

| check | result | residual/margin | digits |
|---|---|---|---|
| arb digamma over box (width 0.035) | output width ±103 (artifact; documents the gap) | — | 70 |
| λR = a−(πρ/2)Ĝ(η) vs −ReΦ′ (N=100…10⁸, η=0.5…2) | exact identity | ≤ 7e-18 | 50 |
| curv2 vs \|ReΦ″\|/2 (N=100, 10⁴) | ratio = 1.0 | < 1e-6 | 50 |
| γ′ bound at x=7.9 | 0.308 ≤ 0.491 | margin 1.6× | 50 |
| γ″ bound at x=7.9 | 0.038 ≤ 0.0985 | margin 2.6× | 50 |
| beta bound \|g\| ≤ Hr(x/√2)/x (x=1…40) | holds | ratio 0.60–0.69 | 50 |
| \|g\| decreasing (x=1,3,7.9,40) | 0.816/0.0679/0.00668/0.000122 | monotone | 50 |
| step_37 sign asserts (all cells) | all pass | — | 70 |
| step_34 C(100), C(400) | 42.4247, 46.5448 | arb width ~1e-63 | 62 |

## 4. Status vs task list

| item | status |
|---|---|
| Fix box digamma/lgamma blow-up (γ real bounds) | **DONE, PROVED** (W1/W2/W3 no longer use complex box evals) |
| λR sign via regular Ĝ form | **DONE, PROVED** (asserts pass on all cells) |
| Cells on [100,400] (ratio 1.2–1.5) | **DONE, certified** (sup C = 2281.28) |
| Dyadic cells to N_J ≈ 10¹⁰ | NOT run (turn limit) — code path identical to cells |
| Final analytic cell y*∈[0,y_J] | NOT run; h₁ bound derived, design fixed (v-moments full-line; tails ~ N^{−r/2}; needs the regular degenerate-box forms) |
| mk_deg9.tsv check (N₀≤N≤400) | NOT run; prior NUMERICAL (D13d, 50 dps): truth 0.60–1.60·sc ≪ 2281·sc, consistent |

**Why the cells are weak (fixable, not a correctness issue):** W3 uses crude full-segment Taylor V-bounds (V ≤ M3u³/6+G1u vs step_34's midpoint+reach) and a Lipschitz q-factor — inflating W3 ≈ 8×; TLf's `ReX2 = −ax2 + θ − (N/2)log(1+vf)` suffers interval dependency between `θ` and `(N/2)log(1+vf)` over wide N-boxes (huge cancellation) — both curable by finer N/σ subdivision and the midpoint-reach technique ported from step_34.

## 5. Files (in `research/08_theory/D13f_sweep/`)

- **step_37.py** — the fixed box assembly (main new result); output `cert_boxes37.json`.
- step_34.py (re-verified: C(100)=42.4247, C(400)=46.5448, 62-digit arb intervals in `cert_single_N.json`), step_35.py (truths cross-check), step_36.py (old, gap).

**Summary of what is PROVED today:** the single-N inequality at N=100 (C=42.4247) and N=400 (C=46.5448); the box-certified inequality with C ≤ 2281.28 for all N ∈ [100,400]; and all supporting regular-form identities/bounds (λR-Ĝ, curv2, γ′/γ″ bounds, Phi3b/Phi4b, beta bound, |g| decreasing) — every number an arb enclosure printed above. The all-N≥N₀ completion (dyadic + final analytic cell) is architected with the arb-gap closed, and remains to be executed.