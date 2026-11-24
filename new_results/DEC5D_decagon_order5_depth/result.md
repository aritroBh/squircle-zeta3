# DEC5D — the regular decagon's P⁽⁵⁾(1): does it need a depth-two polylogarithm?

## VERDICT: **NO** — unlike the octagon, the decagon needs no depth-two polylogarithm at order five.

**PROVED at the symbol level** (exact constant-retaining Goncharov symbol + cobracket obstruction + exact classical-span membership); **number level modulo the period conjecture**. In sharp contrast, the same general code reproduces the octagon's nonzero obstruction (−1485+1035√2, 900−540√2).

---

## Task 1 — Value — NUMERICAL

`dec_engine.py` (Route-A d-series of `polygon_engine.py`, one wedge (−π/10, π/10), ×10 symmetry), two independent runs at dps 120 and dps 170 (final tanh-sinh level diffs 1.08e-117 and 1.85e-167):

- **P⁽⁵⁾(1) = −10.1644521920838037189075927518683458765038789365773310500737349139543150766489835168246984814886990478204766055284955…** (two-run agreement 1.6e-114; P…P⁽⁴⁾ agree to 116–126 digits)
- **P(1)** = 2.008114158862272798697976726327421423305944607218666812… = **10·L** (perimeter), L = √((25−11√5)/10); residual **2.9e-177**
- **P‴(1)** = −4.776200180882238416632022453444885836970773916879016004… = **60·L·S** with the paper/DEC1 15-term ℚ(√5) form; residual **6.6e-170**
- P⁽⁵⁾(1)/L = −50.61690415968495829801181983715576285946644170041575072005035186…

## Task 2 — Exact structure — PROVED algebra, NUMERICAL checks

Edge (vertical, x=−h): V0 = (−h, R sin(π/10)), h = 1/(2φ), R = 1/(2φ²sin(π/5)). `symeng.edge_params` gives j0=2, j1=3, nonloc={0,1,4}; exact ℚ(√5) data:
- local: a₂ = t/(2φ²), a₃ = (1−t)/(2φ²) (α=β=(3−√5)/4);
- non-local: a₁ = (1+t/φ)/(2φ²) (**κ₁ = 1/φ**), a₄ = (1−t/φ²)/(2φ) (**κ₄ = −1/φ²**), a₀ = 1/(2φ) **constant** (κ₀ = 0);
- s₀ = −L/2, h²/L² = (5+2√5)/4, E5 = L·(ℚ(√5)-coefficient polynomial).
Checks: substituted E5/L vs `symeng.E5_sym`: **2.3e-50** at 5 points; normalization P⁽⁵⁾ = 1200∫E5 (**2.7e-41**).
Expansion: **461 monomials** (a,b,c,d,e,f) in (log t, log(1−t), w₁=log(1+t/φ), w₄=log(1−t/φ²), log 2, log φ), degrees {1:6, 2:21, 3:56, 4:126, 5:252} (octagon's exact structure), coefficients L·(p+q√5)tⁿ. Assembly from 756 numerically integrated monomials matches the engine to **1.6e-101**.

## Task 3 — Letters and V — PROVED

w₁ = G(−φ; t), w₄ = G(φ²; t): letters **{0, 1, −φ, φ²}** = {0, 1, −1/κ₁, −1/κ₄}. Float-free S-unit factorisation (`logvec5`): every difference is ±φ^k (1+φ = φ², φ+φ² = φ³, φ²−1 = φ — the golden orbit closes on itself); integrand constants log(2φ²), log(2φ) → log 2, log φ; the log√5 direction never appears. **V = ⟨log 2, log φ⟩_ℚ, dim V = 2** (not 3).

## Task 4 — Depth test — PROVED at symbol level

Exact G-rewrite of all 756 integrals (1024 words, weights ≤ 5, no leading-1 words; numeric match ≥ 94–102 digits). Full number-level symbol with every (a_i−a_j) entry; Hall duality invertible (d=2: det 2; d=3: det 8 — the pairings separate Lie₃(V)⊗Lie₂(V)).

**Controls (all PASS, exact):** Li₅ at 14 S-units (incl. √5, √5/5, 2/φ — d=3 machinery, 24 pairings) → 0; 300+300 random shuffle products → 0; M4c W₅ → **(84, −38√2)**; octagon blocks and total → local (0,0), single-p = single-q = (−53955/8+4770√2, −10215/2+14445√2/4), mixed (48015/4−8505√2, 11115−15525√2/2), total **(−1485+1035√2, 900−540√2)** — exact match of OCT5D2.

**Decagon obstruction (X = log 2, Y = log φ):**

| block | d=3 (24 pairings) | d=2 M4c pairings | (4,1) component | classical span |
|---|---|---|---|---|
| local | 0 | (0, 0) | 0 | yes — pure products |
| single-1 | 0 | (0, 0) | 0 | yes — pure products |
| single-4 | 0 | (0, 0) | 0 | yes — pure products |
| mixed | 0 | (0, 0) | 0 | yes — pure products |
| **total W5/L** | **0** | **(0, 0)** | **0** | **yes — pure products** |

Every block and the total lie **exactly** in the classical span (shuffle products + Li₅; rank 30 of 32 — the 2 missing directions are precisely the depth-two space Lie₃⊗Lie₂), fits verified by exact substitution with residual 0; each solution uses **only shuffle generators** (the weight-5 symbol is a pure product — no Li₅, no depth two). Control: the octagon's total is NOT in its classical span.

## Task 5 — Report

- **No block carries an obstruction**; the whole weight-5 symbol is classical — opposite of the octagon (single-p, single-q, mixed blocks each carry depth two there, totalling (−1485+1035√2, 900−540√2)).
- Structural reason: the decagon's letter differences span a ONE-dimensional log group (φ^k only), mirroring DEC1's weight-3 result (Li₂,₁(−1/φ,−1/φ), Li₂,₁(1/φ²,−φ) reduce).
- **Bonus (NUMERICAL, 239-digit verified):** exact local-block decomposition (34 terms, residual 1.25e-239) with ζ(5)-coefficient **(73125−30825√5)/64** and ζ(2)ζ(3)-coefficient **(−64875+32775√5)/64** — exactly the paper's `prop:order5` values C₅⁽⁵⁾/L, C₅⁽²·³⁾/L for m=5 (a new order-five check of the local law).
- **Full W5/L classical closure: attempted, INCONCLUSIVE (honest).** LLL at 240 and 420 digits over a 32-atom classical basis produced only threshold candidates (maxc 1620/1935/2517/1838858 across runs); the fixed-vector two-precision test **rejects** them (the 240-digit candidate verifies at 1.81e-189 at *both* precisions instead of dropping to ~1e-360). The symbol-level proof stands; the explicit number-level fit needs ~600+ digits.

## Evidence table (all numbers printed by the named steps)

| # | Check | Residual/agreement | Digits | Producer |
|---|---|---:|---:|---|
| 1 | P⁽⁵⁾(1), dps 120 vs 170 | 1.6e-114 | 113+ | dec_engine.py (steps ~20–21) |
| 2 | P(1) = 10L | 2.9e-177 | 177 | step ~22 |
| 3 | P‴(1) = 60LS | 6.6e-170 | 170 | step ~22 |
| 4 | E5/L vs symeng, 5 pts | 2.3e-50 | 50 | step ~25 |
| 5 | assembly vs engine | 1.6e-101 | 101 | step ~28 |
| 6 | rewrite of 756 integrals vs table | 2.0e-102 | 102 | step ~51 |
| 7 | P5/L from rewrite vs engine | 2.2e-114 | 114 | step ~52 |
| 8 | G-table self-checks | 1e-116…1e-122; 1.5e-357 | 120/420 | gengine_dec.py (steps ~50, 92) |
| 9 | controls: Li₅, shuffles, M4c, octagon | exact | — | steps ~53–58 |
| 10 | **decagon obstruction: all vanish** | exact | — | steps ~59–61 |
| 11 | classical-span fits, all blocks+total | residual exactly 0 | — | steps ~65–67 |
| 12 | local block exact; ζ(5), ζ(2)ζ(3) = prop:order5 | 1.25e-239; exact | 239 | steps ~79–83 |
| 13 | W5/L closure 240/420 digits | rejected (1.81e-189 both) | 420 | steps ~87–97 |
| 14 | **verify_dec5d.py: 21 PASS, exit 0** | 186.5 s | — | verify_dec5d.py |

## Deliverables in the work dir

- **`verify_dec5d.py`** — 21 checks, runs 186.5 s (< 5 min), prints 21 PASS lines, exit 0 (fixed the earlier `P'''`-quote syntax bug).
- **`dec_engine.py`** (value engine), **`dec5d_lib.py`** (QS5 arithmetic, float-free logvec5/logvec2, generic symbol machinery + Hall obstruction + rewrite), **`gengine_dec.py`** (G-tables), **`result.md`** (full report).
- Artifacts: `E5_dec_subs.pkl`, `parsed_mons_dec.pkl`, `itab_dec110.pkl`, `rewrite_dec_exact.pkl`, `dec_block_symbols.pkl`, `dec_classical_fits.pkl`, `local_block_exact.pkl`, `w5L_420.pkl`, `p5_gold.pkl`, transcript `step_01.py…step_100.py` (key producers: ~19–22 value, ~23–28 exact structure+assembly, ~30–46 letters/V/duality, ~50–52 rewrite, ~53–58 controls, ~59–61 obstruction, ~65–68 span fits, ~79–83 local block, ~87–97 closure attempts, ~98–100 verifier+result).

**Labels:** PROVED — edge structure, letters/V, the vanishing of the (3,2) and (4,1) obstructions, classical-span membership, all controls (exact arithmetic). NUMERICAL — P⁽⁵⁾(1) (113+ digits, two levels), P(1)=10L, P‴(1)=60LS, assembly matches, the local-block closed form, and the (rejected) W5/L closure candidates. The headline statement "no depth two needed" is PROVED at the symbol level; at the number level it holds modulo the period conjecture.