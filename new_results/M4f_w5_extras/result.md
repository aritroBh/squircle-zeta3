# M4f — smaller coefficients for W5 by choosing the extra Li₅ arguments well

## Verdict

**Yes — a far better choice exists, and the picture is structurally different from what M4e assumed.** With the full pool of exceptional S‑units kept in play, the symbol‑invisible quotient is only **3‑dimensional**, not 5: the "five extras minimum" was an artifact of M4e dropping the ζ(5) and Li₅(b)‑family kernel columns. **Three extras suffice to close W5**, and using the extra freedom of the full 7‑element pool the optimized closed form has **max numerator 80,354 vs M4e's 3,061,997 — a 38× reduction** (the brief's keep‑M4e threshold was 10×), at 71 classical terms (+ the 4 exact Li₃,₂ terms of N5). The lcm denominator grows 1.5× (11,612,160 vs 7,741,440); a Pareto variant keeps M4e's exact denominator with max numerator 965,201 (3.2× reduction).

⚠️ **Session ended early (hard stop):** the final 1500/1590‑digit re‑verification, the drop‑L test, the LaTeX file, and the M4 assembly check were **not run**. Everything claimed below is backed by printed output as itemized in §5; open items have exact rerun recipes in §7.

## 1. The pool (step 1 — `pool.py`, printed)

Enumerating x = ±(√2)ᵐ(1+√2)ⁿ, |m|≤12, |n|≤8, x∈(−1,1)\{0}, with 1−x also an S‑unit: **exactly 16 exceptional units**. All 16 have Li₅‑symbol in the S0 span (exact rref, `target-in-span=True` for each). Nine are S0 members (u=0); the **7 extras candidates**, ordered by height h = max over {x, 1−x} of |m|+|n|:

| z | name | h | particular max\|num\|, lcm |
|---|---|---|---|
| (1−√2... ) 1/2−√2/2 | −b/2 | 3 | 24, 960 |
| −2+2√2 | 2b | 3 | 24, 960 |
| 1/2−√2/4 | (2−√2)/4 | 4 | 144, 3840 |
| 1/2+√2/4 | (2+√2)/4 | 4 | 144, 3840 |
| 17−12√2 | b⁴ | 7 | 1280, 5 |
| 1/2−3√2/8 | (4−3√2)/8 | 7 | 2100, 3840 |
| −16+12√2 | 1−b⁴ | 7 | 2100, 160 |

## 2. Overcomplete search (step 2 — `overcomplete.py 1200`, printed)

Columns: [target, kv(35), √2kv(35), L(χ8,5), √2L, u_z(7), √2u_z(7)] = 87; planted‑relation control **PASS**. Automatic span‑preserving dependency dropping found **14 exact value relations** (residuals 2e‑1310 … 7.2e‑1306), leaving 72 independent columns:

- `kv30 = 0` (duplication {±1/√2}→1/2); `16·kv16 + 15·kv17 = 0` (Li₅(−1) = −15ζ(5)/16, kv17 = ζ5); `kv18 − kv20 − kv22 − 2·kv27 = 0` ({±b}→b² duplication);
- **u(2b) = −u(−b/2) + kv‑combo** (odd‑weight inversion, exact to 1.4e‑1306);
- **u(b⁴) = pure kv‑combo** (duplication; height ~576/5, cf. M4e's 571);
- **u((4−3√2)/8) = (65/2)u₀ + 15u₂ + 5u₃ − (416/5)√2·L + kv‑combo**;
- **u(1−b⁴) = −(65/2)u₀ + 5u₂ + 15u₃ + (416/5)√2·L + kv‑combo** (so u₅+u₆ = 20(u₂+u₃)+kv).

Shortest target relation (LLL, printed): **W5‑coefficient c₀ = −5,806,080 = 2¹¹·3⁴·5·7**, maxcoef 106,522,250,688, residual 6.076e‑1289 — and it uses **only three extras**: u₀=−b/2, u₂=(2−√2)/4, u₃=(2+√2)/4, with coefficients **304, 272⅔, 326**. The target has full support on all 3 quotient directions ⇒ **3 extras is the minimum** (2 can span only a 2‑dim subspace).

## 3. Subset search (step 3 — `subset_exact.py`, `optimize*.py`, printed)

The 14 deps give exact coordinates (max coefficient 159,744 — findable height) of every column in the 72‑column basis, so **every** subset was screened by exact Q(√2) linear algebra: **76 of the 99 subsets with |S|≥3 close**; closing 3‑subsets are exactly the quotient bases (C(6,3)=20 minus dependents). For each, the 2|S|‑dim affine solution space was optimized by cascaded D‑scaled integer‑LLL/Babai + rational descent (score = max\|num\|, then lcm, then terms, after Li₅(b²)→16[Li₅(b)+Li₅(−b)], Li₅(−1)→−15ζ(5)/16).

**Comparison table** (classical terms; add 4 for N5's Li₃,₂):

| form | extras | max\|num\| | lcm den | terms |
|---|---|---|---|---|
| **best (full 14‑dim space optimum, `state_opt2.pkl`)** | 7 | **80,354** | 11,612,160 | 71 |
| earlier representative of same space (`state_opt.pkl`) | 7 | 133,381 | 11,612,160 | 71 |
| lcm‑optimal Pareto variant (`pareto.pkl` B) | 7 | 965,201 | **7,741,440** | 71 |
| {0,1,2,3,6} (5 extras) | 5 | 706,213 | 11,612,160 | 69 |
| {0,2,3,5,6} | 4 | 913,711 | 11,612,160 | 69 |
| M4e R5_nice_reduced | 5 | 3,061,997 | 7,741,440 | 68 |
| best 3‑extra form {0,2,3} | 3 | 31,765,985 | 11,612,160 | 66 |

Note the pattern: more extras ⇒ more affine freedom ⇒ smaller coefficients; the winner uses all 7 pool extras with tiny coefficients (at the printed 133,381 stage: u₄=97/6−(7/6)√2, u₅=256/3−(23/6)√2, u₆=−3+√2; ζ(5) = 109987/1536 − (133381/28672)√2; L(χ8,5) = −57754/63 + (37464/5)√2). 11,612,160 = 2¹²·3⁴·5·7.

## 4. The exact relations are the real discovery

u(1−b⁴) — the pool element M4e never tested — is what unlocks the big win: its relation supplies the √2·L‑entangled direction that kills the large coefficients. 1−b⁴ participates with nonzero coefficient in the optimum.

## 5. Verification status (printed evidence)

| check | result | status |
|---|---|---|
| LLL planted relation [1,√2,1+√2] | recovered [1,1,−1] | **PASS (printed)** |
| 14 dep residuals @1200 digits | 2e‑1310 … 7.2e‑1306 | **printed** |
| target relation residual | 6.076e‑1289 | **printed** |
| M4e's own lchi_final_1200 point ∈ our affine space | exact membership; stats (3,061,997 / 7,741,440 / 68) = M4e's reduced form | **PASS (printed)** — anchors the whole exact setup |
| value of the 443,867 form − W5 (5400 bits), unreduced and reduced | 1.04398e‑1295 | **printed** |
| exact symbol of the winner | holds **by construction** (target‑in‑span printed True; kv rows certified zero‑symbol; each u_z exactly symbol‑zero by exact rref) | construction + asserts printed; verify_nice‑style rerun pending |
| exact (3,2) cobracket | Li₅‑atoms reduce to 0 uniformly; M4e printed it for all our z except 1−b⁴ | rerun pending for the final form |
| winner (80,354) value vs W5 at ≥1590 digits | **NOT run** (inherits ~1e‑1300 iff deps exact) | **open** — rerun below |
| deps re‑confirmed at 1500 digits | **NOT run** | **open** — `overcomplete.py 1500` |
| drop‑L test | **NOT run**; note u₅,u₆ relations carry ∓(416/5)√2·L, so L is load‑bearing | **open** |
| `w5_closed_form.tex` + pdflatex; M4 vs M4_800.txt | **NOT produced** | **open** |
| progress.log per‑step lines | only the copy line was appended; steps tracked via `transcript.json`/`step_NN.py` | noted |

## 6. Winner's data location

`state_opt2.pkl` = `{score:(80354,11612160,71), v, tt}`; slot layout: indices 0–64/65–129 = p/q parts of the 65 S0 atoms (order of `nice_search.atoms(S0)`), 130/131 = L(χ8,5), 132–138/139–145 = p/q of Li₅ at the 7 extras in pool order {−b/2, 2b, (2−√2)/4, (2+√2)/4, b⁴, (4−3√2)/8, 1−b⁴}. `tt` = coordinates in the reduced nullspace basis of `subset_exact.pkl` subset (0,1,2,3,4,5,6).

## 7. Rerun commands (cwd = `new_results/M4f_w5_extras`, PY = `<project>/.venv/bin/python`)

```
nice $PY pool.py                 # step 1: pool.json (16 units, all pass symbol test)
nice $PY overcomplete.py 1200    # step 2: overcomplete_1200.json (14 deps, c0=-5806080)
nice $PY subset_exact.py         # step 3: subset_exact.pkl (76 closing subsets, exact)
nice $PY optimize2.py            # cascaded optimizer -> 443,867
nice $PY optimize3.py            # reweighted+restarts -> 164,273
nice $PY optimize4.py 1300       # -> 133,381 (state_opt.pkl)
nice $PY optimize5.py 1000       # Pareto A/B/C -> B: 965,201/7,741,440; M4e-membership cross-check
nice $PY optimize6.py 1000       # polish -> 80,354 (state_opt2.pkl); optimize7.py continues
nice $PY inspect_winner.py       # coefficient listing + 1.04e-1295 value check
```
Open‑item recipes: re‑confirm deps `nice $PY overcomplete.py 1500`; verify winner by adapting `verify_nice.py` to load `state_opt2.pkl` (exact symbol, exact (3,2) cobracket incl. Li₅(1−b⁴), fresh 1650‑digit values vs `W5_5400.json`); drop‑L test = `subset_exact.py` with columns 70–71 removed; LaTeX/M4 assembly per `display_nice.py` + M4c `planner_assemble_M4.py` against `M4a_decomposition/M4_800.txt`.

Scripts produced this session (work dir): `pool.py, m4f_common.py, overcomplete.py, subset_exact.py, height_lll.py, optimize.py, optimize2.py, optimize3.py, optimize4.py, optimize5.py, optimize6.py, optimize7.py, inspect_winner.py` (+ harness copies `step_01.py…`).
---

## Planner addendum (2026-10-22 22:40 CDT): part 2 of the run
Part 2 ended without a FINAL message, so the section "Session ended early" above is out of date. The following was
run in part 2 and printed, as recorded in transcript.json, turns 50–72.
- **Exports:** R5_best.json, R5_pareto.json and R5_min3.json, each with the reductions b², −1 and −1/√2 applied.

  | file | classical terms | max num | lcm denominator |
  |---|---|---|---|
  | R5_best.json | 70 | 80,354 | 11,612,160 |
  | R5_pareto.json | 70 | 965,201 | 7,741,440 |
  | R5_min3.json | 66 | 31,765,985 | 11,612,160 |

- **Verification of all three forms,** with the M4e verify_nice checks adapted to load a JSON:
  - exact full point symbol: PASS;
  - every Li5, including b⁴ and 1−b⁴, reduces to 0 in the (3,2) cobracket: PASS;
  - reduce(N5 − target) = 0: PASS;
  - fresh 1620-digit value of R5 + N5 − W5_5400: −2.1e-1596 for each form.
- **Dependencies:** `overcomplete.py` at 1200 and 1500 digits gives 14 dependencies at each precision. The vectors
  differ but span the same 14-dimensional space (combined rank 14). The target relation is identical, with
  c0 = −5806080 and residual 2.2e-1589 at 1500 digits.
- **M4 assembly** with R5_best: M4 − M4_800 = 1.2e-803.
- **Not run:** the drop-L test. w5_closed_form.tex compiles to 7 pages with one overfull vbox; it is superseded by
  M4g's supplement.
- **Planner, independent check** (`../M4e_w5_lchi/planner/planner_eval_terms.py`, 1600 digits): R5_best gives
  −2.1e-1596 against M4c's W5 and 2.8e-601 against M4a's W5.
- **Status.** R5_best is the form used in the final supplement:
  `../M4e_w5_lchi/supplement_W5.tex` and `paper_fragment.tex`. The top session also checked it independently, and it
  is logged in the LEDGER.
