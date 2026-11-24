# M4e — W5 via the missing weight‑5 constant L(χ8,5)

## Verdict

**Confirmed (weaker form) / strongest form refuted.** L(χ8,5) is provably the missing symbol‑invisible weight‑5 constant: it lies in the Li5 value‑span with the predicted 13‑denominator, and adding it as an atom removes M4d's 13 and shrinks the W5 relation ~50‑fold. **But** the planner's strong expectation fails: S0 Li5's + ζ(5) + L(χ8,5) do **not** close the value equation (checked with the *saturated* kernel at 1200 **and** 1500 digits); the target genuinely needs extra Li5 arguments. The final **R5_nice is verified** (exact symbol, exact (3,2) cobracket, 1596 digits vs certified W5, 374 digits vs M4c's R5), includes L(χ8,5) explicitly, and is ~100–7000× lower‑height than M4c's R5. The M4 reassembly checks to 1.2e‑803.

## 1. L(χ8,5) (done, 3‑way cross‑validated)

`8^{-5}[ζ(5,1/8)+ζ(5,7/8)−ζ(5,3/8)−ζ(5,5/8)]`:

```
L(χ8,5) = 0.995633856312967456342217717127517682183198240596321485345135329265557739460546…
```
| check | digits | result |
|---|---|---|
| Hurwitz vs direct series Σχ8(n)/n⁵ (10⁶ terms) | 80 | diff 5.0e‑36 (default `nsum` is *unreliable* here: off 7.4e‑6) |
| Hurwitz vs [Li5(ξ8)−Li5(ξ8³)−Li5(ξ8⁵)+Li5(ξ8⁷)]/(2√2) | 80 | full agreement (G(χ8)=2√2 checked) |
| mpmath Hurwitz vs flint `acb.zeta(5,a)` ball | 1650 | diff 4.24e‑1656, radius 2.08e‑1690 |

Saved `Lchi8_5.json` (1650 digits).

## 2. The decisive check: S0 + L does NOT close

- **Saturated the S0 symbol kernel** (`nice_sat_s0.py`, HNF of transposed integer symbol matrix, 65 atoms, rank 30, kernel 35, SNF invariants all 1 ⇒ saturated). **Index of M4d's RREF kernel in it = 967458816000 = 2^17·3^10·5^3.** This huge index invalidates the part‑1 S0+L negative (which used the non‑saturated kernel).
- **Reran with the saturated kernel at 1200 and 1500 digits** (`nice_lchi2.py`): round 0 finds only the expected S0‑internal dependencies (duplications {±b}→b², {±1/√2}→1/2, Li5(−1)=−15ζ(5)/16; the saturated basis exposes them as exact zero‑value kernel vectors, e.g. `16Li5(1/√2)+16Li5(−1/√2)−Li5(1/2)=0`). After dropping them: **no genuine relation** (height < ~1e22 at 1500 digits / 66 columns). 
- **Each extra Li5's symbol‑invisible part u_z is NOT in span(V0∪{L,√2L})** at 800 digits (5 of 6 extras; only b⁴ reduces, with height 571) (`lchi_second.py`). So the extras carry independent directions: the target needs them, and S0+L cannot close at findable height. (Motivically these reductions should exist but with height > 1e12 — that is exactly M4c's 1.4e15 regime.)

## 3. Test B (extras + L): closes, 52× nicer

`lchi_diag.py` / `lchi_testB.py` at 800/1200 digits, extras setup (71 atoms, saturated kernel 41):
| test | result |
|---|---|
| L(χ8,5) ∈ extras‑kernel value span? | **yes, coefficient 638976 = 2^14·3·13**, residual 7.5e‑905 |
| target ∈ span(extras kernel) without L (M4d) | W5‑coeff 100638720 = 2^13·3³·5·7·13, maxcoef 7.67e10 |
| target ∈ span(extras kernel ∪ {L,√2L}) | **W5‑coeff 1935360 = 2^11·3³·5·7**, maxcoef 1475349120, residual −2.0e‑1289 @1200 |

The 13 is gone (that was L's encoding); the **7 persists** (see §6).

## 4. R5_nice (the final closed form) — `closed_form_nice.md`

From the S0‑particular formulation (`lchi_final.py`; c0 = −967680 = 2^12·3³·5·7), expanded to atoms:
**R5_nice = N5 + Σ(68 classical atoms) + L(χ8,5)** — **73 terms**: 14 Li5 + 4 Li_{3,2} + 1 L(χ8,5) + 14 Li4·w1 + 28 w3·w2 + 12 pure products.
- **Li5 arguments (14):** S0 {b, −b, −b², 1/√2, 1/2, 1−1/√2, 2−√2, −1, b²} (note Li5(−1/√2) and ζ(5) get coefficient 0) **+ 5 extras** {−2+2√2, 1/2−√2/4, 1/2+√2/4, 1/2−√2/2, 1/2−3√2/8}. b⁴ is redundant (in span). **5 extras is the proven minimum** (the target has full support on the 5 independent u_z directions).
- **L(χ8,5) coefficient: −8/63 − (82052/315)√2** (in the LLL‑shortest representative; representative‑dependent, see §5).
- Height: max numerator **3.06e6**; lcm denominator **D = 1935360 = 2^11·3³·5·7**; primitive integer relation max coefficient **1.24e10**, W5 coefficient **1935360**. (M4c R5: 1.39e15 and 1.33e10 — i.e. **112× and 6900× smaller**.)

**Verification (`verify_nice.py`, at 1200 and 1500):** all PASS at both precisions —
1. exact full point symbol of R5_nice−N5 = symbol(W5)−symbol(N5);
2. every Li5 atom reduces to 0 in the (3,2) cobracket, and reduce(N5 − target(W5)) = 0 (M4c's identity re‑run);
3. R5_nice − W5(5400 bits) = **3.84e‑1596** (the certified‑W5 precision limit);
4. R5_nice − M4c R5 = **3.94e‑374**.

## 5. Uniqueness structure (`lchi_valrel.py`)

Exactly **2 Q(√2) value‑relations** entangle the solution: one among S0‑kernel values + u(−2+2√2) + u(1/2−√2/2), and one expressing **L(χ8,5) (and √2·L) via the Li5 atoms with denominators 2^11·3·13 and 2^12·3·13**. So L is expressible through Li5's (the predicted 13), the representative is unique only modulo these 2 relations, and both the 1200‑ and 1500‑digit representatives verify independently to 1596 digits. This is the honest sense of "stability": the relation space is fixed; LLL picks different (both exact) representatives.

## 6. Factor 7 (addendum point 4) — `n5_factor7.py`

All 4‑subsets of the 8 Li32 candidates were solved for the cobracket correction. **33/63 are 7‑free; the smallest denominator is 23** (e.g. {(1/√2,1),(1−b,1),(b,−1),(2−√2,1)}), vs the current 56 = 8·7. **Not adopted:** it merely swaps the prime 7 for the larger prime 23 in D, and would require a full re‑run with no clear height gain; the verified 56‑form is kept.

## 7. M4 assembly — `assemble_m4_nice.py`

M4 = R5_nice + R''_nice + W3⁽⁴⁾ + T (R''_nice from M4a w4pp_closed_form.json; W3⁽⁴⁾ = W3′ − ζ(3)/8 − π²λ/24 − λ³/12 from thm:M3; T = (π²+6λ²)/24):
```
M4 = 0.006027993112609122504516108250877605571477…
M4 − M4_800 (M4a full-integrand quadrature) = 1.22e-803   (the M4a precision limit)
```

## Rerun commands (work dir `new_results/M4e_w5_lchi/`, PY = venv python)
```
PY nice_sat_s0.py            # saturated S0 kernel + index  -> nice_ksat_s0.json
PY nice_lchi2.py 1200/1500   # S0 + L saturated kernel LLL (negative) -> lchi_sat_search_*.json
PY lchi_diag.py 800          # tests A/B/C -> lchi_diag_800.json
PY lchi_second.py 800        # u_z membership (2nd-constant check)
PY lchi_testB.py 1200        # test B + atom expansion -> lchi_testB_1200.json
PY lchi_min.py 1200          # minimal-extras form
PY lchi_final.py 1200/1500   # FINAL candidate + atom expansion -> lchi_final_*.json
PY verify_nice.py 1200/1500  # symbol + cobracket + 1596-digit + M4c-374 checks
PY lchi_valrel.py 1200       # value-relation / uniqueness structure
PY n5_factor7.py             # 7-free N5 search
PY display_nice.py           # -> closed_form_nice.md
PY assemble_m4_nice.py       # M4 vs M4_800
```

## Files produced
`Lchi8_5.json` (L(χ8,5), 1650 digits + flint ball), `nice_ksat_s0.json` (saturated kernel), `closed_form_nice.md` (**the R5_nice display**), `lchi_final_1200.json`/`lchi_final_1500.json` (atom coefficients), `lchi_diag_800.json`, `lchi_testB_1200.json`, `lchi_min_1200.json`, `lchi_sat_search_{1200,1500}.json`, and scripts `nice_sat_s0.py, nice_lchi.py, nice_lchi2.py, lchi_diag.py, lchi_second.py, lchi_testB.py, lchi_min.py, lchi_final.py, verify_nice.py, lchi_valrel.py, n5_factor7.py, display_nice.py, assemble_m4_nice.py`, `progress.log`; copies from M4d: `w5_*.py, nice_search.py, W5_5400.json, cobracket.json, nice_symbol.json, nice_atoms.pkl, nice_ksat.json`. Producing steps (auto‑saved): `step_28.py`–`step_58.py` (L(χ8,5) checks 28–32; kernel prep 33; S0+L LLL 34–35,41–42; diagnostics 36–39,43; saturation 40; test B 44; minimal/final forms 45–47; verification 48–49; factor‑7 50; display 51–53; M4 assembly 54–55; value‑relations 56; stats 57).

**Bottom line.** L(χ8,5) is the missing weight‑5 constant (verified: it is in the Li5 value‑span with the 13). The strong "S0+L closes" hope fails (saturated kernel, 1200+1500 digits) — the target intrinsically needs 5 extra Li5 arguments whose reduction to {ζ(5),L(χ8,5)} has astronomical height. With those extras + L(χ8,5) explicit, R5_nice is a verified, readable (73‑term, L(χ8,5)‑carrying) closed form ~100–7000× lower‑height than M4c's, and M4 reassembles to 1.2e‑803.
---

## Planner verification (2026-10-22 21:40 CDT)

**Verdict.** R5_nice is accepted as a numerical identity: it agrees with W5 to 1596 digits, and its full point symbol and (3,2) cobracket agree with W5 exactly. It is not an exact certificate.
- L(χ8,5) is confirmed as the missing weight-5 constant.
- The brief's stronger claim, that the nice set S0 plus L(χ8,5) must close, is withdrawn. Kimi's negative on the saturated kernel is consistent with the corrected count below.

**Why S0 + L(χ8,5) is not enough (correction to the brief).** The Lie coalgebra of this category is cogenerated in weights 1, 3, 5, two generators each:
- weight 1: λ and L;
- weight 3: ζ(3) and L(χ8,3);
- weight 5: ζ(5) and L(χ8,5).

So L₅ is 16-dimensional: 6 brackets of weight-1 generators, 8 of the type [[f₃,f₁],f₁] or [f₃,[f₁,f₁]], and the 2 weight-5 generators.

Two things set the brief's claim apart from the truth:
- **Symbol zero is a weaker condition than the brief assumed.** It removes only the 6 pure weight-1 directions.
- **The (3,2) cobracket fixes only 2 of the 8 f₃-directions.** These are f₃ ∧ [f₁,f₁].

That leaves about 6 f₃-directions besides ζ(5) and L(χ8,5). Only the 4⊗1 cobracket sees them, through the ζ(3)/L(χ8,3) part of the Li₄'s. The S0 Li₅'s do not span them, which is why Kimi needed 5 extra arguments. The brief's "S0 + L must close" overlooked these directions.

**Checks run by the planner:**
1. Fresh rerun of `verify_nice.py 1500` (`planner/rerun_verify_nice_1500.out`). Every check passes:
   - exact symbol;
   - exact (3,2) cobracket;
   - R5_nice − W5 = 3.8e-1596;
   - R5_nice − M4c R5 = 3.9e-374.
2. Independent mpmath evaluation (`planner/planner_eval_R5nice.py 1600`; output `planner/planner_eval_R5nice_1600.out`).
   - Method: direct series for every Li_k and Li_{3,2}; N5 hard-coded from the brief; L(χ8,5) from mpmath Hurwitz ζ. No imports from the M4c/M4d/M4e code.
   - L(χ8,5) − Lchi8_5.json = 2.4e-1641.
   - R5_nice − W5 (M4c, 5400 bits) = −2.1e-1596, for both the 1200 and the 1500 representative.
   - R5_nice − W5 (M4a W5_600, an independent quadrature) = 2.8e-601.
   - rep1200 − rep1500 = −1.4e-1638, i.e. they differ by an exact value relation.
   - Also R5_nice = M4c's R5 to 1516 digits.
3. Reduced representative (`planner/planner_reduce.py`, giving `planner/R5_nice_reduced.json`). It applies two exact substitutions: Li5(b²) → 16Li5(b) + 16Li5(−b) (duplication), and Li5(−1) → −(15/16)ζ(5).
   - **72 terms:**
     - 12 Li5's: W4's seven arguments b, −b, −b², 1/√2, 1/2, 1−1/√2, 2−√2, plus five extras 2b, (2−√2)/4, (2+√2)/4, −b/2, (4−3√2)/8;
     - ζ(5) and L(χ8,5);
     - 4 Li_{3,2};
     - 54 products.
   - The extras all satisfy "x and 1−x are both S-units".
   - The 1200 representative is nicer than the 1500 one: max numerator 3.06e6 versus 2.76e7.
   - **Correction to §4 above.** The displayed coefficients have lcm denominator 7741440 = 2^13·3^3·5·7. The figure 1935360 there is the multiplier of W5 in the integer relation.
   - Coefficients after reduction (1200 representative):
     - Li5(b): 391577/252 + (1001239/2016)√2
     - Li5(−b): −683777/252 + (1272659/2016)√2
     - ζ(5): −1099999/21504 − (426661/32256)√2
     - L(χ8,5): −8/63 − (82052/315)√2
     - The other Li5 coefficients are as in closed_form_nice.md.
4. M4 assembly. R5_nice equals M4c's R5 to more than 1500 digits, so the earlier planner check (M4 − M4_800 = 1.2e-803) carries over. Kimi's own assembly gives 1.22e-803.

**Caveats:**
- closed_form_nice.md shows the unreduced 14-Li5 form. Its LaTeX does not compile as written: the first row of each aligned block lacks `&`, and it contains `\text{\operatorname{Li}_5}`.
- The coefficient of L(χ8,5) depends on the representative, because L(χ8,5) lies in the Li5 span through a relation with denominator 2^11·3·13.

**Open:**
- Is there a choice of the 5 extra arguments with much smaller coefficients? Kimi's "5 is the minimum" concerns the number of extras, not their choice.
- The factor 7 from N5's 56. Kimi's 7-free N5 alternatives have denominator 23, so they are not better.
- An exact certificate.

## Planner update (2026-10-22 22:25 CDT): final supplement and fragment
- The table now uses M4f's form `../M4f_w5_extras/R5_best.json` (copy: `planner/R5_best_M4f.json`).
  - It has 14 Li5's: W4's 7 arguments plus 2b, (2−√2)/4, (2+√2)/4, −b/2, (4−3√2)/8, b⁴, 1−b⁴.
  - It also has ζ(5), L(χ8,5), 54 products and N5: 74 terms in total.
  - Max numerator 80,354, against 3,061,997 for the 12-Li5 form. lcm denominator 11,612,160 = 2^12·3^4·5·7.
- Verification:
  - Kimi M4f: exact full point symbol PASS, exact (3,2) cobracket PASS (including Li5(b⁴) and Li5(1−b⁴)), and a fresh 1620-digit value giving R5 + N5 − W5 = −2.1e-1596. The 14 dependencies span the same space at 1200 and 1500 digits.
  - Planner, independent (`planner/planner_eval_terms.py ../../M4f_w5_extras/R5_best.json 1600`): −2.1e-1596 against M4c's W5 and 2.8e-601 against M4a's W5.
  - M4g's `check_supplement.py`: exact round-trip of all 70 + 4 terms, and a 430-digit value of 3.0e-427.
- Files:
  - `supplement_W5.tex` / `.pdf` (3 pages; tectonic, 0 errors, 0 overfull boxes);
  - `paper_fragment.tex` (the top session's three wording fixes, plus the 14-argument sentence).
  - The 12-Li5 versions are kept as `*_12li5*`.
