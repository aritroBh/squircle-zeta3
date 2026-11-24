# OCT3: does the regular octagon carry ζ(3) beyond the local edge law?

**Short answer.** The prescribed hunt (PSLQ/LLL on P‴(1) against the reduced weight-3 S-unit basis) finds **no relation** — over {1,√2} up to coefficient size ≈2·10⁷, over the natural cyclotomic extension {1,√2,√(2+√2),√2·√(2+√2)} up to ≈9·10³, and for the (structurally motivated) quotient P‴(1)/√(2+√2) over {1,√2} up to ≈8.7·10¹² (520-digit data). Meanwhile the lower derivatives **do** close exactly over this basis: P′(1) (residual 8.8·10⁻²⁸⁰) and P″(1) (residual 1.2·10⁻²⁷⁴). A Galois argument shows that *if* P‴(1) has any closed form in this span, its total ζ(3)-coefficient must be **rational**, hence **≠ C4** (since C4=(54−39√2)√(2+√2) ∉ ℚ(√2)). So the octagon question stays **open**, with a precise basis, precise bounds, and a concrete suspect for the missing constant (a weight-3 Clausen/16th-root-of-unity value: the weight-3 Borel rank of ℚ(ζ₁₆)⁺ is 4, our S-unit Li₃ span supplies only 3).

All numbers below were printed by the code runs (scripts auto-saved in the work directory as `step_01.py … step_29.py`; artifacts `consts31.json, consts62.json, basis34.json, pder.json, hidata.json`). This message is the content of `result.md`.

---

## 1. Data identification and verification (PROVED/NUMERICAL)

The file `vals_m2_8_200.json` keys are the number m of normal directions (m=2…8). **The octagon is key `"4"`** (4 directions ⇒ 8 edges); the brief's "key 8" refers to the 16-gon and contains none of the stated values. Key `"4"` = [P(1), P′(1), P″(1), P‴(1)]:

| k | value (first 52 digits) | check |
|---|---|---|
| 0 | 2.53620268449771540274021154137584885279564001 | = 8(√2−1)√(2−√2) exactly (residual 3.8·10⁻²¹¹) |
| 1 | 2.75023713004966371870863005296385194252854498 | recomp. matches file to 3.9·10⁻²⁰⁰ |
| 2 | −1.45620511038120464465373069667959364567771524 | recomp. matches to 2.3·10⁻¹⁹⁹ |
| 3 | −4.39449630527605215997911011380609991382237745 | = stated P‴(1); recomp. matches to 1.6·10⁻¹⁹⁶ |

- P(1) identification (PROVED + NUMERICAL): on [0,π/4], Σ|⟨nᵢ,x⟩| = (1+√2)max(|x|,|y|)+min(|x|,|y|), so K₁ is the regular octagon of circumradius √2−1; perimeter 8(√2−1)√(2−√2) = 8(2−√2)/√(2+√2) = 2.53620268449771540274…
- Independent recomputation (NUMERICAL): with r(θ)=u(θ)⁻¹/ᵖ, u=Σ|cos(θ−iπ/4)|ᵖ, P(p)=8∫₀^{π/4}√(r²+r′²)dθ, Taylor-expanded to order 3 in t=p−1 and integrated (mpmath quad). At 210 digits all four derivatives match the file to 196–201 digits; at 520 digits they match the 280-digit run to 2.2·10⁻²⁷³ (P‴).
- C4 (NUMERICAL): 192 sin⁴(π/8)cos(π/8)(3sin²(π/8)−1) = −2.1329217491383995109743830724674630027123771… = (54−39√2)·√(2+√2) (difference 1.8·10⁻²⁷⁹ at 280 dps). **PROVED** C4 ∉ ℚ(√2): (a+b√2)²=2+√2 ⇒ a²+2b²=2, 2ab=1 ⇒ b²=(2±√2)/4 ∉ ℚ.

## 2. Task 1 — candidate set (built at 280 and 520 digits)

17 S-unit points z (all verified in [−1,1), exact in ℚ(√2)): −1, −√2/2, −b, (1−√2)/2, −b², (4−3√2)/8, b⁴, (2−√2)/4, b², 1−√2/2, b, 1/2, 2−√2, √2/2, 2b, (2+√2)/4, 12√2−16 (b=√2−1); plus ζ₃, π²log2, π²log(1+√2), the 4 cubic log monomials, π², the 3 quadratic log monomials, the 2 logs, 1 — each ×{1,√2}: 62 constants. Sanity identities printed at 280 digits: Li₃(½) formula (residual −1.4·10⁻²⁸¹), Landen at b (−6.9·10⁻²⁸¹), duplication Li₃(b²)=4(Li₃(b)+Li₃(−b)) (−1.4·10⁻²⁸¹).

Cross-ratio audit (PROVED, exact sympy, homogeneous coordinates): the cross-ratios of {0,1,∞,−(1+√2),2+√2} take 15 values, mapping to 8 representatives in [−1,1) = {−1, −√2/2, 1−√2, 1−√2/2, √2−1, 1/2, 2−√2, √2/2} ⊂ our 17 (orbits O1,O2,O3 only). The basis is a superset of what the letters require.

## 3. Task 2 — internal relations and the independent basis (NUMERICAL)

Exact LLL (python-flint `fmpz_mat.lll`) on the 62 constants, scale B=10²⁵⁵: **exactly 28 relations** (residuals ≤ 6·10⁻²⁷⁸, coefficients ≤ 8), cleanly separated from spurious vectors (coefficients ≈ 5·10⁶, residuals ≈ 10⁻²⁵⁰). They are 14 ℚ-relations ×{1,√2}: duplication (Li₃(b²)=4Li₃(b)+4Li₃(−b), Li₃(½)=4Li₃(√2/2)+4Li₃(−√2/2), Li₃(b⁴)=4Li₃(b²)+4Li₃(−b²)), Landen-orbit relations, and Li₃(−1)=−¾ζ₃, Li₃(½)=⁷⁄₈ζ₃−¹⁄₁₂π²log2+⅙log³2 — as predicted (orbit count 11 + 3 squaring = 14).

**Independent basis (34 = 17 × {1,√2})**: {Li₃(2b), Li₃((2+√2)/4), Li₃(12√2−16), ζ₃, π²log2, π²log(1+√2), log³2, log²2·log(1+√2), log2·log²(1+√2), log³(1+√2), π², log²2, log2·log(1+√2), log²(1+√2), log2, log(1+√2), 1} × {1,√2}. Exactly 3 genuine weight-3 constants — matching Borel rank (K₅ of ℚ(√2) with S={2}: 2+1). All 28 eliminated constants are explicitly expressed in the basis (max residual 5.8·10⁻²⁷⁸; e.g. Li₃(b) = ¹³⁄₁₆Li₃(2b)+¹¹⁄₈Li₃((2+√2)/4)−⁵⁄₁₆Li₃(12√2−16)−¹³¹⁄₁₂₈ζ₃+…). **Completeness**: no integer relation among the 34 with residual < 10⁻²⁵⁰ and coefficients ≤ 3.2·10⁷ (LLL, B=10²⁶⁰) — improving the old 10⁶ bound.

Weight-2 layer (for Task 5): 42 constants → 32 relations → 5 independent = {π², log²2, log2·log(1+√2), log²(1+√2), Li₂(12√2−16)} ×{1,√2} (rank 1, as Borel predicts); e.g. Li₂(b) = −¹⁄₂₀Li₂(12√2−16) + ¹⁹⁄₂₄₀π² + ½log2·log(1+√2) − ⁷⁄₁₀log²(1+√2).

## 4. Task 3 — the P‴(1) hunt: NEGATIVE (NUMERICAL)

| hunt | basis | digits/B | result |
|---|---|---|---|
| P‴ vs 34 | {1,√2} | 280 / 10²⁶⁰ | smallest |s| = 1.55·10⁻²⁵⁴ at maxcoeff 2.56·10⁷ (spurious; a true relation there would show |s| ≲ 10⁻²⁷⁰) ⇒ **no relation, maxcoeff ≲ 2·10⁷** |
| P‴ vs 68 | {1,√2,u,u√2}, u=√(2+√2) | 280 / 10²⁶⁰ | nothing below maxcoeff ≈ 9.3·10³ |
| P‴−C4·ζ₃ vs 34 and vs 68 | both | 280 / 10²⁶⁰ | nothing (2.1·10⁻²⁵⁸, maxcoeff 7856) |
| P‴/u vs 36 | {1,√2}, weight-3 triple + products + weight-2 + logs + 1 | 520 / 10⁴⁸⁰ | smallest |s| = 1.14·10⁻⁴⁶⁹ at maxcoeff 8.68·10¹² = the Diophantine floor (10^(480/37)) ⇒ **no relation, maxcoeff ≲ 8.7·10¹²** |

The divided hunt is structurally justified: P(1)/u = 8b² ∈ ℚ(√2), and P′, P″ are exactly u×(ℚ(√2)-span) (below) — yet P‴/u is not in the span to 8.7·10¹².

## 5. Task 4 — interpretation

- **(PROVED)** No {1,√2}-relation can ever reproduce the local-law ζ(3) term: C4=(54−39√2)√(2+√2) ∉ ℚ(√2). Any exact closed form for P‴(1) needs the field ℚ(ζ₁₆)⁺=ℚ(√2,u).
- **(CONJECTURE, but a rigorous implication given its hypotheses)** The span of our constants over F=ℚ(ζ₁₆)⁺ is Gal(F/ℚ)-stable (inversion/Landen close it), ζ₃ is Galois-fixed and is numerically linearly disjoint from the other basis constants over F. Hence if P‴(1) ∈ span_F, its ζ₃-coefficient q is unique and Galois-invariant, so **q ∈ ℚ** — necessarily **q ≠ C4**. I.e., the total ζ(3) carried by the octagon cannot be the local-law value alone: the non-local terms must leak additional ζ(3) (through reductions like Li₃(½)=⁷⁄₈ζ₃+…), making the total rational.
- **(NUMERICAL)** Despite that, no closed form was found in the tested spans to the bounds above. The likeliest explanation **(CONJECTURE)**: the weight-3 part of P‴(1) contains a constant outside the S-unit Li₃ span. Field count: Borel rank at weight 3 for ℚ(ζ₁₆)⁺ is r₁+r₂ = 4, our span supplies 3; the π/8 local geometry (sin/cos(π/8) factors, already visible in P′, P″) points to 16th-root-of-unity letters, whose weight-3 values are Clausen-type (e.g. Im Li₃(ζ₁₆ᵏ)) — the natural next basis element.
- **Answer to the title question: OPEN.** Precise basis and bounds are reported above; the ζ(3) coordinate q is undetermined, but it cannot equal C4 in any Galois-closed closed form — so either the octagon carries ζ(3) beyond the local edge law (with rational total q), or its weight-3 content is not S-unit Li₃ at all.

## 6. Task 5 — sanity checks: P′(1) and P″(1) DO close (NUMERICAL)

u := √(2+√2). Both relations are exact to the working precision and were re-verified by direct evaluation.

**P′(1)** (integer relation: a₀=1, maxcoeff 40; residual 8.81·10⁻²⁸⁰ at 280 dps):
P′(1) = u·[(40−28√2)·log2 + (24−16√2)·log(1+√2)]  — pure weight 1.

**P″(1)** (integer relation: a₀=−60, maxcoeff 7380; residual 1.17·10⁻²⁷⁴; recheck 1.8·10⁻²⁷⁶):
P″(1) = u·[(−7/30 + 13√2/60)π² + (123−87√2)log²2 + (88−62√2)log2·log(1+√2) + (104/5−68√2/5)log²(1+√2) + (−4/5+3√2/5)·Li₂(12√2−16) + (−80+56√2)log2 + (−48+32√2)log(1+√2)].

Note: contrary to the brief's guess, P″(1) is **not** π²-and-log only — it needs the genuine weight-2 S-unit constant Li₂(12√2−16) and has weight-1 terms. This also shows why the monolithic 69-dim hunts fail: P″'s relation (maxcoeff 7380) sat below the big-basis detection floor and only the weight-targeted basis exposed it — the same strategy that then *failed* for P‴ even at 520 digits/10⁴⁸⁰.

## 7. Evidence table

| claim | status | evidence (printed) |
|---|---|---|
| file key "4" = octagon; entries [P,P′,P″,P‴] | PROVED/NUM | P(1) vs 8(√2−1)√(2−√2): 3.8·10⁻²¹¹ |
| P^{(k)} recomputed (independent series method) | NUM | vs file: 1.35·10⁻²⁰¹, 3.9·10⁻²⁰⁰, 2.3·10⁻¹⁹⁹, 1.6·10⁻¹⁹⁶; 520- vs 280-digit runs: ≤2.2·10⁻²⁷³ |
| C4 = (54−39√2)√(2+√2) | NUM | difference 1.8·10⁻²⁷⁹ |
| C4 ∉ ℚ(√2) | PROVED | b²=(2±√2)/4 ∉ ℚ (exact) |
| 28 internal relations, basis of 34 independent | NUM | residuals ≤ 5.8·10⁻²⁷⁸; completeness: no relation < 10⁻²⁵⁰ with coeffs ≤ 3.2·10⁷ |
| cross-ratio reps = 8 ⊂ 17 | PROVED | exact sympy |
| P′(1) formula | NUM | residual 8.81·10⁻²⁸⁰ |
| P″(1) formula | NUM | residual 1.17·10⁻²⁷⁴, recheck 1.8·10⁻²⁷⁶ |
| no {1,√2} relation for P‴ | NUM | spurious floor maxcoeff 2.56·10⁷ at |s|=1.55·10⁻²⁵⁴ |
| no {1,√2,u,u√2} relation for P‴ | NUM | floor maxcoeff 9.3·10³ |
| no relation for P‴/u over {1,√2} | NUM | floor maxcoeff 8.68·10¹² at |s|=1.14·10⁻⁴⁶⁹ (520 dps, B=10⁴⁸⁰) |
| q ≠ C4 if a closed form exists in the span | CONJECTURE | Galois argument (§5) |
| missing constant = Clausen/16th-root weight-3 value | CONJECTURE | rank count 4 vs 3; π/8 letters |

## 8. Step files (auto-saved in the work directory; numbering as saved by the harness)

- `step_01.py` — load JSON, identify C4 and target.
- `step_02.py` — dump all keys; direct P(1) check (identifies octagon = key "4").
- `step_03.py` — series method; P…P‴ at 210 dps; comparison with file.
- `step_04.py` — 62 constants at 280 dps; sanity identities; saves `consts31/62.json`.
- `step_05…08.py` — LLL setup; diagnosis of the `round()`-on-mpf precision bug (residuals 10⁻¹⁷ → fixed with `mp.nint`).
- `step_09.py` — exact LLL (B=10²⁵⁵): 28 relations isolated.
- `step_10.py` — relations printed; rref.
- `step_11.py` — 28 explicit expressions; independent 34-basis; saves `basis34.json`.
- `step_12.py` — completeness (LLL + PSLQ; no relation ≤3.2·10⁷).
- `step_13.py` — P^{(k)} at 280 dps; saves `pder.json`.
- `step_14…16.py` — P‴ hunts over {1,√2} and ℚ(ζ₁₆)⁺; P‴−C4ζ₃ hunts (all negative).
- `step_17…18.py` — P′(1) relation found and printed.
- `step_19…20.py` — cross-ratio set of the 5 letters (exact) = 8 reps ⊂ 17.
- `step_21.py` — lean-basis P‴ hunt; weight-2 internal relations (5 independent).
- `step_22…24.py` — P″(1) hunt in weight-targeted basis; relation extracted (a₀=−60, maxcoeff 7380).
- `step_25.py` — P‴/u over {1,√2}×36 at 280 dps (negative ≤9·10⁶).
- `step_26…29.py` — 520-digit recomputation of all constants and P^{(k)} (`hidata.json`); final divided hunt at B=10⁴⁸⁰ (negative ≤8.7·10¹²).

**Bottom line.** The reduced independent weight-3 basis is exactly as Goncharov's argument predicts (3 genuine S-unit Li₃'s + ζ(3) + products), its internal relations were found and printed at 10⁻²⁷⁸ level, and independence holds past coefficient size 3·10⁷. Within that basis — and its ℚ(ζ₁₆)⁺ extension — P‴(1) admits no relation with coefficients up to 2·10⁷ ({1,√2}), 9·10³ (full cyclotomic), or 8.7·10¹² (P‴/u over {1,√2}), while P′(1) and P″(1) close exactly with tiny coefficients. Since C4∉ℚ(√2) and any closed form in this span forces a rational total ζ(3)-coefficient q≠C4, the data say: either the octagon's total ζ(3) is a rational multiple of ζ(3) exceeding the local edge law, or its weight-3 content involves a non-S-unit (Clausen/16th-root) constant. The question stays open; the basis, the relations, and the exclusion bounds above are the precise state of the art.