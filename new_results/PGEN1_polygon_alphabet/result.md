# PGEN1: the trilogarithm alphabet of a general polygon, and the cobracket of the far kink lines

**Deliverables:** `result.md` (full report), `verify_pgen1.py` (runs in **8.6 s**, returncode 0, output in `verify_output.txt`), `CR_presentations.pkl`, `delta_oct9d.pkl`, `CR_regular.pkl`, `field_tables.pkl`, `task3_L3.pkl`, `symbolproof_engine.py`, `progress.log`, and the session `step_01.py … step_127.py`.

Notation: on an edge, far kink line j has B_j(1+κ_j t) and crosses the edge's line at **p_j = −1/κ_j**. δ(N₃) = Σ c_i (z_i∧(1−z_i))⊗z_i ∈ ∧²(F^×)⊗F^× over ℚ (the (2,1)-part of the symbol, δ{z}₃ = {z}₂⊗z followed by {z}₂ ↦ z∧(1−z)). Over ℚ, ℚ(√2), ℚ(√3), ℚ(√5) the pre-Bloch group injects into ∧²F^× (B₂ = 0 rationally since r₂ = 0), so the wedge decides δ = 0 faithfully and is presentation-independent.

## Task 1 — the alphabet rule: CONFIRMED on every known closed form, PROVED at the symbol level for s ≤ 4

**(a) Hexagon (s = 1), PROVED.** CR of {0,1,∞,p} = one S₃-orbit {1−r, r, 1/r, 1/(1−r), (r−1)/r, r/(r−1)}; both of thm:hexagon's arguments x = 1−r and x/(x−1) lie in it (exact, 6 sample r). At r = 2, 1/2 the points become −1, 1/2, and δ = 0 for any coefficients — exactly why prop:nonlocal's non-local ζ(3) is canonical there.

**(b) Generic rational octagon (s = 2), PROVED (symbol level).** Per edge: p₁,p₂ = (−2/3,12/7), (−3,−9), (−3/5,9/5), (−4,−24); |CR_e| = 30,30,30,18. Every edge's symbol TW3_e lies in span{S(Li₃(z)) : z ∈ CR_e} + D (exact ranks; CR-span = 10,8,8,5); explicit **cross-ratio-only Li₃ presentations of all four edges** extracted (25,18,18,14 terms; TW3 − sum ∈ D exactly). **Dimensions explained:** dim Q₃ = (n³−n)/3 (Witt, free Lie algebra degree 3): 8 for the 3-letter edges {2,3,5}, 20 for the 4-letter edge {2,3,5,7} — the prime 7 enters through p₂ = 12/7. The OCT9d basis arguments need not be in CR_e (chosen for small invisible constants); some follow by duplication (−9/16 = −(3/4)², 4/9 = (2/3)², 7−4√3 = (2−√3)²).

**(c) Regular polygons, PROVED.** Octagon: CR_e = 15 points = orbit(√2) ∪ orbit(−√2/2) ∪ orbit(1/2); all 6 Li₃ points in it. Decagon: CR_e = single orbit {±φ, ±1/φ, φ², 1/φ²}, contains both points ±1/φ. 12-gon: |CR_e| = 63 over the 3-letter space {1+√3, √3, 2+√3} (with (1+√3)² = 2(2+√3)); all 8 DOD5 points in CR_e or duplication (7 directly, z₄ = z₁²); span{s(CR)} mod D = full 8-dim; all 8 symbols in it (exact).

**(d) Symbol-level proof for s ≤ 4, PROVED (exact flint integer ranks).** With recursion S(I(0;a₁..aₙ;1)) = Σᵢ S(âᵢ) ⊗ (a_{i+1}−a_i)/(a_i−a_{i−1}) (verified numerically at weights 2,3), formal letters {pᵢ, 1−pᵢ, pᵢ−pⱼ}: every convergent distinct-letter word symbol lies in span{s(cross-ratios)} + D:

| s | n | dim Q₃ | #CR | rank s(CR) mod D | #words | rank words mod D | words ⊂ span s(CR)+D |
|---|---|---|---|---|---|---|---|
| 1 | 2 | 2 | 6 | 2 | 3 | 2 | YES |
| 2 | 5 | 40 | 30 | 10 | 14 | 7 | YES |
| 3 | 9 | 240 | 90 | 30 | 39 | 17 | YES |
| 4 | 14 | 910 | 210 | 70 | 84 | 34 | YES |

and every weight-3 word reduces via shuffles to those + products (PROVED, s ≤ 4). With Goncharov (weight-3 MPLs reduce to Li₃) and Rudenko (they reduce to Li₃ at single cross-ratios), this proves the rule at symbol level for ≤ 4 far kink lines — covering all known examples. **Theorem candidate:** P'''(1) of any centrally symmetric convex polygon is a combination of Li₃ at the cross-ratios of the points where the kink lines cut the line of each edge, ζ(3), and products.

## Task 2 — δ(N₃) of the far kink lines

**Controls (PROVED):** δ{1/2}₃ = 0 (consistent with Li₃(1/2) = 7/8ζ(3)+products); Landen triple δ = 0 (7 rationals); hexagon arguments in CR; δ = 0 at r = 2, 1/2.

- **Generic octagon (PROVED, δ ≠ 0 on every edge; identical from the OCT9d basis and from the CR presentations):** edge 1: (2∧3)⊗(5976·2−6516·3+2448·5)/3125 + (2∧5)⊗(−6192·2−5328·3+6012·5)/15625 + (3∧5)⊗(−17568·2+14220·3−7111·5)/15625·(2/2, see full list in verify output) — edge 1 has 9 terms, edge 2: 8, edge 3: 6, edge 0: 24 over {2,3,5,7} (denominators 1681, 68921). Total δ ≠ 0 (edge scalars √41, 1, √13, √29 are ℚ-independent). **Geometric formula:** δ(N₃)|_e = Σ_{z∈CR_e} c_z (z∧(1−z))⊗z with the cross-ratio presentation coefficients; each lies in ker(∧³) = Lie³.
- **Regular octagon (PROVED, δ ≠ 0):** δ(N₃) = (r2∧u)⊗r2·(−5/2 + 15√2/8) + (r2∧u)⊗u·(−11/4 + 2√2); invariant under the whole relation lattice (Landen triples T1, T2 and duplication all have δ = 0, rechecked exactly).
- **Decagon (PROVED, δ = 0):** letters collapse to span{φ}, every wedge z∧(1−z) = 0. **Canonical pair (NUMERICAL):** N₃ = **(1/64 + 31√5/320)·ζ(3) + (875/768 − 125√5/192)·L(χ₅,3)** (val(N₃) − aζ(3) − bL₅ = −8.2e-41 at 120 dps). Key input, the golden analog of eq:chi3: **χ₃(1/φ) = (25√5/48)·L(χ₅,3) − (π²/12)log φ + (1/4)log³φ** — residuals **1.5e-151 (150 dps), 8.7e-252 (250 dps)**; eq:chi3 itself re-verified to 4.9e-150 (L(χ₈,3) = 0.958380454563094562…, and ℒ₃(√2−1)−ℒ₃(1−√2) = (4√2/3)L(χ₈,3) exactly, two precisions).
- **12-gon (PROVED, δ ≠ 0):** 9+9 terms over {1+√3, √3, 2+√3} (full list in verify output and result.md §2(d)); invariant under duplication z₄ = z₁² (exact).
- **Hypothesis test:** "regular 2m-gon's non-local part is δ-closed after summing over edges" — **FALSE for m = 4 (octagon) and m = 6 (12-gon); holds for m = 5 (decagon, trivially: golden letters collapse); trivially true m = 2, 3.**

## Task 3 — single-valued L₃

Convention: L₃(x) = Re Li₃(x) − log|x| Re Li₂(x) − (1/3)log²|x|·log|1−x| (Zagier; verified: Landen, inversion, duplication, L₃(1/2) = 7/8ζ(3), series cross-check).

- **Edges 1, 2, 3: OCT9c and OCT9d give identical L₃(ξ)** (diffs 5.7e-100, 2.6e-100, 1.7e-99): 0.02367462581917732…, 0.48743643627030…, −0.02102185222422….
- **Edge 0: the bases differ by exactly (80/206763)·ζ(3)** (pslq at 200 dps; 206763 = 3·41³). All three presentations (OCT9c, OCT9d, cross-ratio) have **identical symbols** (exact, 0 difference) and **rational ζ(3) pairwise L₃-differences** (pslq at 200 dps): Sc−Sd = −(31119817701814913/654492976038000)ζ(3), S_CR−Sd = −(6228044/5169075)ζ(3) — so they share one B₃-class; the shift is a rational ζ(3)-trade between the Li₃ basis and the invisible constants, the paper's "the rest of the ζ(3) depends on the basis", consistent with δ(N₃) ≠ 0 for edge 0. **Why the invariant exists:** L₃ descends to B₃ (satisfies all trilog functional equations), so the single-valued evaluation of the class is basis-independent up to the {1}-share; for edges 1-3 the two sessions' shares agree exactly.

## Evidence table

| check | label | residual/exactness |
|---|---|---|
| All controls | PROVED | exact |
| TW3_e ∈ span s(CR_e)+D, 4 edges | PROVED | exact ranks (10,8,8,5) |
| CR-only presentations (25,18,18,14 terms) | PROVED | TW3−sum ∈ D exactly |
| δ per edge: OCT9d basis ≡ CR basis | PROVED | exact, identical |
| regular octagon / decagon / 12-gon δ | PROVED | exact; δ = 0 only for decagon |
| δ invariant under duplication (12-gon) & Landen lattice (octagon) | PROVED | exact |
| alphabet proof s ≤ 4 + shuffle reduction | PROVED | exact flint ranks |
| decagon pair (a, b) | NUMERICAL | 8.2e-41 at 120 dps |
| golden χ₃ identity (25√5/48) | NUMERICAL | 1.5e-151 / 8.7e-252 |
| eq:chi3 control; D₈/L(χ₈,3) = 4√2/3 | NUMERICAL | 4.9e-150; exact 150 & 250 dps |
| L₃(ξ_e) OCT9c ≡ OCT9d, edges 1–3 | NUMERICAL | 5.7e-100 … 1.7e-99 |
| edge-0 trade (80/206763)ζ(3); rational q's | NUMERICAL | pslq 200 dps |

**Files that produced these:** `verify_pgen1.py` (all checks, 8.6 s), `symbolproof_engine.py` (s ≤ 4 proof), `CR_presentations.pkl` (exact CR presentations), `delta_oct9d.pkl`, and session steps: step_10–14 (machinery, controls, CR presentations, δ two ways), step_16–22 (regular-polygon geometry, δ, ∧³ checks), step_29–32 (12-gon/decagon CR sets and spans), step_34–44 (alphabet proof s ≤ 4, shuffle reduction), step_45–53 (decagon pair, golden identity, eq:chi3 control), step_54–86 (task-3 L₃ chain, edge-0 resolution), step_87–97 (result.md, verify script, final robustness).