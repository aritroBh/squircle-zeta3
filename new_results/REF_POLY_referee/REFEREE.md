# Referee report: thm:poly, thm:Cm, prop:edgelaw (main.tex draft 36, lines 1582-1626)

All code in this folder is the referee's own; nothing from the project is imported.
Python: `.venv/bin/python` (mpmath, sympy). The paper was only read.

## Verdicts

- **thm:poly (Theorem 23): PROVED AFTER REPAIRS.** Every number is right. Repairs P1, P2 are wording and one missing sentence.
- **thm:Cm (Theorem 24): PROVED AFTER REPAIRS.** C_m is right (exact). Repair C1: the mixed terms contribute C_m ζ(3) *plus* other constants; only their ζ(3) part is C_m ζ(3). Repair C2 (presentation).
- **prop:edgelaw (Proposition 7): PROVED AFTER REPAIRS.** Z_e is re-derived exactly, in full generality. The parallelogram and hexagon expansions are right (re-derived; see §2). Repairs E1-E4 are wording. The printed statement is the corrected one: it claims "all of its ζ(3)" only for fully local polygons, and the next paragraph gives the counterexample. It claims nothing that the E1 verifier refuted.

No printed coefficient of these three results was obtained by fitting. Integer relations appear in the notes only as checks or as first discovery; every printed formula now has an exact derivation (notes, and re-done here).

## 1. Statement checks (numerical, independent)

Method (`jet.py`). Polar perimeter integral P(p) = ∫ r √(1+(r_θ/r)²) dθ, r = F^{-1/p}, F = Σ_j w_j|⟨u_j,x⟩|^p (power convention: w = 1, u = v). At each node the integrand is a truncated power series in d = p−1 (generic series arithmetic: exp, log, sqrt, division), i.e. differentiated under the integral sign. Tanh-sinh on each sector between kink directions; the vanishing kink term at each end of a sector is computed from the offset to that end (no cancellation). Two runs: dps 75 / h = 1/64 and dps 90 / h = 1/128 (`log_jet_*.txt`, `jet_*.json`); `check.py` compares (`log_check.txt`).

Polygons (`polys.py`): diamond; parallelogram v = (1,1/2), (−1/3,1); regular hexagon, octagon, decagon (n_i = (cos iπ/m, sin iπ/m)); affine-regular hexagon v = M(1,0), M(0,1), M(1,−1), M = [[2,1/3],[−1/2,1]]; generic hexagon v = (1,0), (1/5,1), (2/3,−1/2); generic octagon v = (1,0), (1/5,1), (−2/3,1/2), (1/2,1/3).

| test | agreement |
|---|---|
| two jet runs, P, P′, P″, P‴, all 8 polygons | ≤ 1.8e-69 (≈ 69 digits) |
| P(1) vs polygon perimeter, all 8 | ≤ 4.5e-85 |
| thm:poly P(1) = 4m(sec−cos) and P′(1) = P(1)R₁, m = 2,3,4,5 | ≤ 4.9e-85 |
| hexagon P′(1) = 2√3 log 2 | 3.8e-85 |
| referee's general first-order formula (below), all 8 polygons | ≤ 4.9e-85 |
| parallelogram law P″ = π²𝓔/3, P‴ = 6𝓔(ζ(3)−π²/3): diamond, generic parallelogram | ≤ 2.5e-84 |
| hexagon law (P, P′, P″, P‴ from 2^{d/(1+d)}[…]): regular and generic affine-regular hexagon | ≤ 3.8e-85 |
| printed regular-hexagon P‴(1) | 2.6e-85 |
| Σ_e Z_e − C_m, m = 2,3,4,5 (Z_e from the printed formula and the polygon's own v_j) | ≤ 2.9e-101 (exact by hand too, below) |
| finite differences of the *direct* P(p) (19-point fit, δ = 1/200, `fd.py`) vs jet: P′, P″, P‴, all 8 polygons | P′ ≤ 2e-17, P″ ≤ 1.4e-14, P‴ ≤ 8e-12 (13-15 digits typical) |

The FD row is the check of differentiation under the integral sign; it uses no expansion.

**Integer relations (fully local polygons, `pslq_local.py`, jet at dps 140, h = 1/256; relations identical at 90 and 130 digits).**

| polygon | relation | ζ(3) coefficient − Σ_e Z_e |
|---|---|---|
| regular hexagon | 8P‴/√3 + 9ζ(3) − 3π²log2 − 16log³2 + 6π² + 96log²2 − 96log2 = 0 (residual 1.6e-136): exactly the printed closed form | 2.6e-141 |
| affine-regular hexagon \|x\|+\|y\|+\|x−y\| (edges 1/2, 1/2, √2/2) | on {1,√2}×{ζ(3), π²log2, log³2, π², log²2, log2, 1}, residual 5.7e-136 | 5.2e-141 |
| parallelogram (1,1/2), (−1/3,1) (L² = 73/49, 97/49) | on {1/L₁,1/L₂}×{ζ(3), π², 1}, residual 6.5e-133; the constant term is 0 | 2.1e-140 |

So for these fully local polygons the ζ(3) coefficient found by integer relation (against the constants the paper names) is Σ_e Z_e. For the hexagon this is −9√3/8 = C₃.

**Convention (`weighted.py`).** The same parallelogram K_1 with the weights 2, 3 held fixed (Σ w_j|⟨v_j/w_j,x⟩|^p) gives P′(1) = 4.771…, P″(1) = −0.3817…, P‴(1) = −0.2941…, against 0, 4.394…, −16.731… in the unweighted family. The parallelogram and hexagon expansions are statements about the unweighted family Σ|⟨v_j,x⟩|^p only. Z_e itself is convention-free (exact argument in §2).

## 2. Derivation checks

### thm:poly
- Notes (D15 §2): integration by parts in the polar angle; the remaining integrand is S₀(S₀+S₀″)/(…) ≡ 0 because S₀″ = −S₀ (straight edge). **By hand.** I checked the algebra and the boundary values; correct.
- Referee's own independent proof (general polygon, by hand): on an edge, x_p = y·F^{-1/p} with F = 1 + dσ(t) + O(d²), σ = Σ_j a_j log a_j. Since y″ = 0, the first variation of the length of y(1−dσ) is −[σ⟨τ,y⟩]₀¹. Hence for every centrally symmetric polygon in the power convention
  P′(1) = Σ_e (σ(V₀)⟨τ_e,V₀⟩ − σ(V₁)⟨τ_e,V₁⟩).
  Confirmed numerically on all 8 polygons (≤ 4.9e-85). For the regular 2m-gon it gives the printed R₁. "P′(1) = 0 only for m = 2": V = cot(π/2m) > 1 ≥ sin(iπ/m) for m ≥ 3, so every term of R₁ is ≥ 0 and the i = 1 term is > 0. Correct.
- Differentiation under the integral sign: argued in the notes (domination by t^{−δ}(1−t)^{−δ}·logs, E1 Lemma 1), **not stated in the paper** for this theorem (P is not analytic at p = 1). Confirmed by FD (above).

### thm:Cm
- Notes (D15c): exact sympy for the polar mixed coefficients G₂₁, G₁₂; collapse to six sec² / sec⁴ terms; exact integral. **Exact computer algebra.** Also checked numerically (10^-120).
- Referee: Σ_e Z_e = C_m **by hand**: for the regular 2m-gon with s = sin(π/2m), c = cos(π/2m): h = s, L = 2s²/c, L₀ = L₁ = L/2, ξ₀ = ξ₁ = R sin(π/m) = 2s². Then Z_e = 6·4s⁴·s²·(c/2s²)(6s²−2) = 24s⁴c(3s²−1), and 2m edges give C_m. With the exact general Z_e (next item) this re-proves thm:Cm independently of the polar computation.
- The ζ(3) part of a mixed term in t is 2ζ(3)∫coefficient (Beta lemma, re-checked in `exact_laws2.py`). The proof idea integrates the *polar* coefficients; this is legitimate because G₂₁dθ = g₂₁dt pointwise (the cubic coefficient does not change when the logarithms are shifted by functions free of them) — the paper should say the mixed terms are taken in the edge parameter t, as in the paragraph before the theorem (C2).
- **Issue C1.** The mixed terms do *not* contribute exactly C_m ζ(3). Exact computation (`misc.py`, regular hexagon): 6·Σ_e∫(g₂₁ℓ₀²ℓ₁ + g₁₂ℓ₀ℓ₁²)dt = −(9√3/8)ζ(3) + √3(131 − 11π²)/12. Their ζ(3) part is C₃ζ(3); they also carry π² and rationals.

### prop:edgelaw
- **Z_e: re-derived exactly, in full generality (`edge_sym.py`, `log_edge_sym.txt`).** Symbols for the two corner terms αt, β(1−t), two far terms B₁+C₁t and the remainder (Σa_j = 1), arbitrary corners (x₀,y₀), (x₁,y₁), the logs as free symbols, exact truncated series in d. Results: the coefficients of X₀²X₁ and X₀X₁² do not involve the far data B₁, C₁; and 12∫₀¹(g₂₁+g₁₂)dt − Z_e(printed) = 0 identically (after L² = |V₁−V₀|²). Exact computer algebra.
- The symbolic E_k was validated against the polar jet arc by arc, k = 0..3, on non-fully-local polygons (generic hexagon, generic octagon, regular octagon): ≤ 2.6e-45 at 45 digits (`log_e3_vs_polar.txt`). So the "local terms" really are the terms of the length element in the radial-projection parameter.
- Beta lemma (ζ(3) in ∫tⁿ(1−t)^k logᵃt logᵇ(1−t) only for (a,b) = (2,1), (1,2), coefficient 2B(n+1,k+1)): exact sympy, n ≤ 3, k ≤ 2, all (a,b) with a+b ≤ 3 (`log_exact_laws.txt`).
- Convention: weights held fixed shift every log by a constant; the cubic part of a polynomial does not change under constant shifts, so Z_e depends only on K_1 (by hand). The v_j are determined by K_1 (the jump of ∇ of the gauge across each kink line), so "Z_e depends only on K_1" is well posed.
- Fully local: "every parallelogram; a hexagon iff affine-regular; none for m ≥ 4": by hand (for m ≥ 4 an edge has ≥ 2 far kink lines, which are not parallel to each other; for a hexagon, V₁−V₀ ∥ V₂ and V₂−V₁ ∥ V₀ force V₁ = V₀+V₂). Correct.
- Parallelogram and hexagon expansions: in the notes by exact symbolic integration (E1 s04, s08); here checked to 84 digits on two parallelograms and two hexagons, and by integer relation. Exact re-derivation on four rational instances: total − law = 0 exactly at every order (§5).
- "Z_e vanishes as soon as either corner is flattened": ξ₀ = |⟨v_{j0},V₁⟩| → 0 as the jump v_{j0} → 0. Correct.

## 3. Hypotheses and wording

1. **"Fully local"** is defined precisely enough ("every edge is parallel to the kink lines through the corners it does not touch"), given the representation Σ|⟨v_j,x⟩| with pairwise non-parallel v_j. But "so that every term is local" is loose: on a fully local polygon the far logarithms are *constants* (log 2 for the regular hexagon), and terms with them remain; by the paper's own definition ("terms without these last logarithms") they are not local. What is true: no non-constant far logarithm occurs, so the expansion is a finite combination of Beta moments.
2. **"This is all of its ζ(3)"** / "C₂ and C₃ are the whole ζ(3)": true in the closed form (presentation) produced by the Beta lemma. As a statement about the number it presupposes that ζ(3) is linearly independent of the other constants (π², π²log2, log³2, …, over the field of the edge data), which is open. Line 1469 and Prop. order5 already carry this caveat; Prop. 7 and Thm 24 do not.
3. **Non-degeneracy:** the v_j must be pairwise non-parallel (then K_1 is a strictly convex 2m-gon). Stated in the sentence before Prop. 7; fine.
4. **Convention:** Prop. 7 is stated for Σ|⟨v_j,x⟩|^p (power convention). Z_e is convention-free (proved); the two expansions are not (weighted.py). The text after the proposition (line 1628) says this for P‴(1) in general. Fine, but the statement could say "in the family K_p = {Σ|⟨v_j,x⟩|^p ≤ 1}".
5. **The refuted form (E1 verify.md).** Not reprinted. The general claim is now restricted to fully local polygons, and line 1628 states the counterexample with the correct numbers (231/125 √5 − 15/64 √2 vs 72/125 √5 − 15/64 √2). One sentence before the proposition still reads like the refuted general claim (E4 below).

## 4. Repairs and sentences to change

- **P1 (thm:poly, title/statement).** The title says "at the polygon" but the theorem is for the regular 2m-gon only. Suggest: "First derivative at the regular polygon". Optional: add the general formula P′(1) = Σ_e(σ(V₀)⟨τ_e,V₀⟩ − σ(V₁)⟨τ_e,V₁⟩), σ(V) = Σ_j|⟨v_j,V⟩|log|⟨v_j,V⟩| (checked here to 85 digits on 8 polygons).
- **P2 (thm:poly, proof paragraph).** Add: "P is C^∞ but not analytic at p = 1; its derivatives may be taken under the integral sign, because for |p−1| ≤ δ < 1 the p-derivatives of the integrand are dominated by C t^{−δ}(1−t)^{−δ}(1+|log t|+|log(1−t)|)^k." The same sentence covers thm:Cm and prop:edgelaw.
- **C1 (thm:Cm, line 1593).** "the mixed terms contribute exactly C_m ζ(3) to P‴(1)" → "the ζ(3) part of the mixed terms of P‴(1) is exactly C_m ζ(3) (the mixed terms also contribute multiples of π² and rationals)".
- **C2 (thm:Cm, proof idea).** Add "(the mixed coefficients in the polar angle and in the edge parameter t agree after multiplying by dθ/dt, so their integrals agree)".
- **E1 (prop:edgelaw, line 1608).** "the local terms of e contribute exactly Z_e ζ(3) to P‴(1)" → "the ζ(3) part of the local terms of e is exactly Z_e ζ(3)".
- **E2 (prop:edgelaw, line 1611).** "If the polygon is fully local, this is all of its ζ(3)." → "If the polygon is fully local, the expansion is a finite combination of Beta moments, and in the resulting closed form Σ_eZ_e is the coefficient of ζ(3) (as a coefficient of the number, this assumes that ζ(3) is linearly independent of the other constants, which is open)." Same caveat for line 1626 "So C₂ and C₃ are the whole ζ(3) of P‴(1)…" → "…the whole ζ(3) of these closed forms."
- **E3 (definition, line 1604).** "so that every term is local" → "so that on every edge the far logarithms are constants and every non-constant term is local".
- **E4 (line 1604).** "So where does ζ(3) live? Not at a single corner …, but where two corners meet across the side that joins them." This is the refuted general form in words: for polygons that are not fully local the far kink lines add ζ(3) (line 1628). Suggest: "… but, in the local terms, where two corners meet across the side that joins them."
- **E5 (optional, statement).** Say "for K_p = {Σ_j|⟨v_j,x⟩|^p ≤ 1}" in the proposition, since the two expansions fail for fixed weights (here P‴(1) = −0.294 vs −16.731 for the same parallelogram).

## 5. Status of the exact re-derivation of the two expansions

`exact_laws2.py` integrates the exact E_k of `edge_sym.py` with an exact Beta table and compares with the printed expansions, exactly (sympy, radicals simplified). Fully symbolic corners were too slow (stopped after 25 min), so it was run on four rational instances: parallelograms with corners (1,0),(1/3,2) and (2,−1/2),(−1/5,3/2); affine-regular hexagons with V₀,V₁ = (1/2,0),(1/2,1/2) and (1,1/3),(2/5,3/2), V₂ = V₁−V₀. For every instance and every order d⁰…d³: total − law = 0 exactly (`log_exact_laws.txt`). Together with the notes' general symbolic derivation (E1 s04, s08) and the 84-digit checks, the two expansions are established. The parallelogram and hexagon P′(1) values (0 and P(1)log2) also follow from the referee's general first-order formula.

## Files

- `jet.py` (polar jet and direct perimeter), `polys.py`, `geom.py` (edge data, Z_e, first-order formula)
- `run_jet.py` → `log_jet_75_64.txt`, `log_jet_90_128.txt`, `jet_*.json`; `check.py` → `log_check.txt`
- `fd.py` → `log_fd.txt`
- `pslq_local.py` → `log_pslq_local.txt`, `log_pslq_parallelogram.txt`
- `edge_sym.py` → `log_edge_sym.txt`, `edge_E_general.pkl`
- `e3_vs_polar.py` → `log_e3_vs_polar.txt`
- `exact_laws2.py`, `exact_laws_lib.py` → `log_exact_laws.txt`
- `misc.py` → `log_misc.txt` (mixed terms of the regular hexagon)
- `weighted.py` → `log_weighted.txt`
- `t_diamond.py` (sanity)
