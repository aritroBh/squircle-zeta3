result = r'''# OCT7 — is "no double polylog survives" special to the REGULAR octagon?

**Question.** For a centrally symmetric octagon whose two non-local ratios per edge are NOT Landen-paired,
does P'''(1) need a genuine depth-2 constant? (Regular octagon: OCT6 showed Λ, swap reduce to single Li₃'s.)

**Answer (with the evidence below): YES, depth-2 survives generically.** On the rational octagon
K = { |x| + |y| + |2x+y| + |x+3y| ≤ 1 } (normals m = (1,0),(0,1),(2,1),(1,3), weights 1), the edge
V0=(1/4,0) → V1=(1/3,−1/9) has non-local ratios r1 = 4/3, r2 = 10/9 (κ1 = 1/3, κ2 = 1/9 > 0,
r1·r2 = 40/27 ≠ 1, x2 ≠ x1/(x1−1), x2 ≠ 1/x1). Its edge integral E_arc = ∫₀¹E₃ dt is exactly

  E_arc = (single Li₂/Li₃, logs, ζ2, ζ3 part) + (33·Λ1 + 24·Λ2 + 99·Λ3 + 60·Λ4)/3125

with Λ1 = Li_{2,1}(−1/3, 1/3), Λ2 = Li_{2,1}(−1/9, 3), Λ3 = Li_{2,1}(1/4, 2/5), Λ4 = Li_{2,1}(1/10, 5/2),
and all four Λ's survive LLL reduction tests against the natural single-polylog basis (60-digit smoke test
below; 400+ digit definitive test in step_06). The regular-octagon reduction is therefore a
**symmetry accident tied to Landen pairing r1·r2 = 1**.

## Task 1 — the octagon and its per-edge data (step_01, exact)

K = {Σ_j|⟨m_j,x⟩| ≤ 1}, normals (1,0),(0,1),(2,1),(1,3) — a rational centrally symmetric octagon
(polar of the zonotope Z = Σ[±m_j]; convex as a norm ball). Half-boundary edges (CCW), term values
a_j(t) = |⟨m_j, V0 + t(V1−V0)⟩| (two vanish at the endpoints = local; two non-local with ratios r_j):

| edge | V0 → V1 | A0 | A1 | B1 | κ1 | r1 | B2 | κ2 | r2 | \|y'\|² | r1r2 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 0 | (0,1/5) → (1/4,0) | 1/4 | 1/5 | 1/5 | 3/2 | 5/2 | 3/5 | −7/12 | 5/12 | 41/400 | 25/24 |
| **1** | **(1/4,0) → (1/3,−1/9)** | **1/9** | **1/4** | **1/4** | **1/3** | **4/3** | **1/2** | **1/9** | **10/9** | **25/1296** | **40/27** |
| 2 | (1/8,−1/4) → (1/3,−1/9) | 5/9 | 5/8 | 1/8 | 5/3 | 8/3 | 1/4 | −5/9 | 4/9 | 325/5184 | 32/27 |
| 3 | (0,−1/5) → (1/8,−1/4) | 1/8 | 1/5 | 1/5 | 1/4 | 5/4 | 3/5 | 1/24 | 25/24 | 29/1600 | 125/96 |

No edge is Landen-paired (r1r2 ≠ 1 everywhere), no x2 = x1/(x1−1), no x2 = 1/x1, no r = 1.
Checks: Σa_j(t) ≡ 1 on each edge (F[0] = 1, exact); working edge 1 has **rational |y'| = 5/36**
(9-12-15 triangle) — its E₃ coefficients are purely rational.

## Task 2 — exact E₃ and reduction (steps 03–05)

- step_03: exact E₃ for edge 1 (OCT5 builder, rational data; |y'| = 5/36). 164 log-monomials
  (83 local, 36+36 single-ratio, 9 mixed); point-verified vs numeric engine at t = 0.2,0.3,0.55,0.8,
  worst 3.2e-52.
- step_04: generalized mixed engine (OCT6 recursion structure; Jloc + K1d Jtab at x1 = −1/3, x2 = −1/9;
  Kölbig reduction of the A-atom, 6A(x) = 12Li₃(x)+6Li₃(x/(x−1))−6Li₂(x)ln(1−x)−ln³(1−x), residual 0/1.5e-67);
  **new exact D-constants D1–D10** (derived by the series identity F(a,b) = Li_{2,1}(a,b/a)/a, the reflection
  t↦1−t, and IBP relations κ1D5+κ2D6 = ℓ1ℓ2, κ2D7+2κ1D9 = ℓ1²ℓ2, κ1D8+2κ2D10 = ℓ1ℓ2²; each verified by
  quadrature at step_02, ≤ 7e-71). Verification (step_04): 30 K1-integrals worst 1.9e-63,
  20 mixed integrals worst 7.7e-62, JLOC 1.2e-66.
- step_05: exact assembly E_arc, 69 terms; **vs independent quadrature of ∫E₃: 8.4e-65**.

**Depth-2 content of E_arc (exact rational coefficients):**
(33·Λ1 + 24·Λ2 + 99·Λ3 + 60·Λ4)/3125, where
Λ1 = Li_{2,1}(x1, x2/x1), Λ2 = Li_{2,1}(x2, x1/x2) [the D1/D2 constants],
Λ3 = Li_{2,1}(β1, β2/β1), Λ4 = Li_{2,1}(β2, β1/β2) [the D3/D4 constants],
x_j = −κ_j = (−1/3, −1/9), β_j = κ_j/(1+κ_j) = (1/4, 1/10).

## Task 3 — reduction tests (smoke at 60 dps; definitive 400+ in step_06)

Iterative-PSLQ reducer (internal classical relations quotiented out; only target-involving relations accepted):
all four Λ's SURVIVE at 60 digits over a 67-element basis (logs 2,3,5; ζ2, ζ3; Li₂/Li₃ at
−1/3,−1/9,1/4,1/10,−1/2,−2/3,2/3,3/5; all weight-≤3 products), maxcoeff 1e9.
Definitive 400-digit test: see step_06 (below when run).

## Task 4 — control (regular octagon): see step_08 (below when run).

## Task 5 — structural criterion

Landen pairing ⟺ **r1·r2 = 1** ⟺ x2 = ι(x1) with ι(x) = −x/(1−x) ⟺ β1 = x2 (⟺ β2 = x1)
⟺ the two points share one anharmonic S₃-orbit. Only then do Λ3 = Λ2, Λ4 = Λ1 collapse,
leaving 2 depth-2 constants which (for the √2-orbit of the regular octagon) reduce to singles.
Generic edge: 4 distinct orbits-worth of constants; they survive.

## Files
step_01 (geometry), step_02 (Λ's + D-formulas), step_03 (exact E₃), step_04 (engine verification),
mixed_engine_gen.py (the generalized engine), step_05 (exact assembly), geometry.pkl, E3_edge1.pkl,
Earc_exact.pkl. Control + 400-digit tests: step_06+.
'''
with open("result.md","w") as f: f.write(result)
with open("progress.log","a") as f: f.write("RESULT.md written (draft, tasks 1-3 + engine)\n")
print("result.md written", len(result), "chars")
