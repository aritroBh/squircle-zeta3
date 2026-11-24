# DEC1 — the regular decagon's P'''(1) in closed form: the octagon pattern HOLDS (no depth-2 constants)

## Headline result

**P'''(1) = 60·L·S**, L = 2 sin(π/5)/(√5 φ²) (edge length; L² = (25−11√5)/10), φ = (1+√5)/2,
where S is the following **15-term** ℚ(√5)-combination (p+q√5 shown as p, q):

| constant | p | q |
|---|---|---|
| Li₃(1/φ)        | −75/32 | 81/64 |
| Li₃(−1/φ)       | −35/32 | 53/64 |
| ζ(3)            | 33/64 | −19/80 |
| π²/6            | −1/16 | −9/80 |
| log φ           | 5/2 | −1/2 |
| log 2           | 1 | 0 |
| log φ·π²/6      | 3/16 | −13/320 |
| log²φ           | −73/8 | 29/8 |
| log 2·π²/6      | 1/32 | 9/160 |
| log 2·log φ     | −5 | 1 |
| log²2           | −1 | 0 |
| log³φ           | 1099/96 | −1991/384 |
| log 2·log²φ     | 81/16 | −33/16 |
| log²2·log φ     | 5/4 | −1/4 |
| log³2           | 1/6 | 0 |

(rational part, coefficient of "1": exactly 0.) Verified against an independent tanh-sinh quadrature of the
whole exact integrand: **|60 L S − P'''(1)| = 1.1e-676** (P'''(1) computed to ~680 digits, error est. 1e-827;
also matches D15b's 200-digit value to 2e-101·10^0, all its digits). P'''(1) = −4.776200180882238416632…

**The octagon phenomenon is confirmed for m = 5: NO depth-2 constant survives.** The two double polylogs
Λ = Li₂,₁(−1/φ, −1/φ) and swap = Li₂,₁(1/φ², −φ) both reduce to single Li₃'s (verified 830 digits,
cross-validated independently to 700 digits), and the A-atoms A(±1/φ…), A(1/φ²) reduce PROVABLY
(orbit system of three proved function identities).

## 1. Edge data (Task 1) — and a structural surprise

K_p = {Σ_{i=1}^5 |⟨n_i,x⟩|^p ≤ 1}, n_i = (cos(iπ/5), sin(iπ/5)). K_1 is a regular decagon
(vertices concyclic, radius R = 1/(2φ² sin(π/5)); rotation by π/5 permutes the normals).
One edge: from V0 = (−h, R sin(π/10)) to V1 = (−h, −R sin(π/10)) (vertical edge x = −h), h = 1/(2φ) the
apothem, L = R/φ = 2R sin(π/10) its length. The five term values a_j(t) = w_j|⟨n_j, y(t)⟩| are affine:
- local: a₂ = A₀ t, a₃ = A₁(1−t), A₀ = A₁ = 1/(2φ²);
- non-local: a₁ = B₁(1+κ₁t), B₁ = 1/(2φ²), **κ₁ = φ−1 = 1/φ** (ratio r₁ = φ);
  a₄ = B₄(1+κ₄t), B₄ = 1/(2φ), **κ₄ = 1/φ−1 = −1/φ²** (ratio r₄ = 1/φ);
  a₅ = B₅ = 1/(2φ), **κ₅ = 0 — CONSTANT**.

**Why:** the edge's reflection symmetry maps terms to terms with reciprocal ratios; with three non-local
terms one must be fixed, and the fixed one is the normal perpendicular to the edge, whose term value is the
constant apothem h. So the decagon edge has genuinely only TWO non-local logs — Landen-paired
(x₁ = −κ₁ = −1/φ ↔ x₄ = −κ₄ = 1/φ², x₄ = x₁/(x₁−1)) — plus one constant letter log(2φ). The integral-type
inventory is therefore EXACTLY the octagon's (five mixed types), not a new triple-mixed family.

**Field:** all coefficients in ℚ(√5); the only radical is the overall factor L (the analog of the octagon's u).

**E3 exactly:** built by a 100-dps series engine with exact ℚ(√5) reconstruction (max denominator 384,
residual 1.6e-99): 55 log-monomials in the 5 term logs; in the reduced basis
{ℓ₀=log t, ℓ₁=log(1−t), w₁=log(1+t/φ), w₄=log(1−t/φ²), l₂=log 2, l_φ=log φ}: **261 monomials**
= 95 local + 70+70 single-ratio + 26 mixed (20 distinct mixed integrals: the five types
(a,b,c₁,c₄) ∈ {(0,0,1,1),(1,0,1,1),(0,1,1,1),(0,0,2,1),(0,0,1,2)}, n = 0..3 — the octagon's exact set).

P'''(1) = 3! · (10 edges) · ∫₀¹ E₃ = **60·L·∫₀¹ E₃n**. Checks: E₃(0.3) vs direct d-expansion 6.3e-82;
60 L∫E₃n vs 200-digit P'''(1) 9.3e-59 (60-dps quad); then P'''(1) to ~680 digits by tanh-sinh of the exact
integrand (err est. 1e-827; matches D15b's 200 digits, saved P3_660.pkl/txt).

## 2. Mixed integrals exactly (Task 2)

The decagon mixed engine (dec_engine.py, modelled on OCT6's): Jloc (40 keys) and K1d's symbolic J-table
(40 keys) instantiated at x₁ and x₄ (L₁ = −log(1−x₁) = −l_φ, L₁' = +l_φ), plus ten master constants
D1–D10 with closed forms in the named atoms, each derived by partial fractions / IBP and verified by
quadrature (≤1.8e-101):
- D1 = Λ/x₁ = ∫ log t·w₄/u₁, D2 = swap/x₄ = ∫ log t·w₁/u₄ (series identities, 6.3e-102);
- D3 = φ(l_φ Li₂(x₄)+swap), D4 = φ²(l_φ Li₂(x₁)−Λ) (t↦1−t reflection; 1.4e-101);
- D5 = φ(l_φ² − Li₂(1/φ) + Li₂(x₄)), D6 = φ²(Li₂(x₄) − Li₂(1/φ) + 2l_φ²) (log-splitting; 5.4e-102, 1.4e-101);
- D7 = 2φ²(l_φ³ − l_φ Li₂(1/φ) + Li₃(1/φ) − Li₃(x₄)), D8 = φ(2Li₃(1/φ) − l_φ³ − 2l_φ Li₂(x₄) − 2Li₃(x₄))
  (exact Nielsen antiderivatives; 1.8e-101, 1.3e-101);
- D9 = −φ l_φ³/2 + D7/(2φ), D10 = −φ² l_φ³/2 + φ D8/2 (IBP; 1.8e-102).
All **20 mixed integrals** (5 types × n=0..3) reduce by the IBP recursions to the named basis and match
quadrature: worst **8.4e-111**. The only depth-2 objects appearing anywhere: Λ and swap.

## 3. Reductions (Task 3)

**Golden Li₃ lattice** (verified 830 digits): the Landen 3-term identity
**Li₃(1/φ²) = 4Li₃(−1/φ) + 4Li₃(1/φ)** (R_a, 1.1e-830); the sum relation
Li₃(−1/φ)+Li₃(1/φ) = (ζ(3) − l_φ·π²/6)/5 + l_φ³/6 (2.4e-831); Lewin's classical
Li₃(1/φ²) = 4(ζ(3) − l_φ π²/6)/5 + 2l_φ³/3 (1.3e-831). Individually Li₃(−1/φ), Li₃(1/φ) do NOT reduce
(no relation vs {ζ(3), l_φ π²/6, l_φ³} at 830 digits) — two independent constants, like the octagon's six reps.
Golden Li₂'s all reduce: Li₂(1/φ²) = 2z₂/5 − l_φ², Li₂(1/φ) = 3z₂/5 − l_φ², Li₂(−1/φ) = −2z₂/5 + l_φ²/2.

**A-atoms — PROVED** (not just LLL). With OCT6's three function identities (Landen-sum L2, inversion L1,
the difference identity DIFF — each proved there by differentiation + pinning), the 4×4 orbit system
{L2@x₄, DIFF@x₄, L2@(1/φ), L1@φ} (with the Li₃ inversion Li₃(−x) = Li₃(−1/x) − z₂ ln x − ln³x/6,
itself proved by differentiating and pinning at x=1, verified at x = φ, 2, 3 to 1e-100) solves:
  A(−1/φ) = 4ζ(3)/5 − 2 l_φ z₂/5 + 2 Li₃(−1/φ)
  A(1/φ²) = ζ(3) − 3 l_φ z₂/5 + 4 Li₃(−1/φ) + 3 Li₃(1/φ)
The orbit solutions equal these exactly modulo the golden Li₃ lattice (symbolic difference 0), and match
the direct series (residuals 6.5e-103, 8.1e-103 at 100-digit eval; the LLL forms match series to 5.2e-830 / 2.9e-831).

**Depth-2 constants — reduce, NUMERICAL (LLL) + independent cross-validation:**
  Λ = Li₂,₁(−1/φ,−1/φ) = 2ζ(3)/5 + 4 l_φ z₂/5 − l_φ³ + 3Li₃(−1/φ) + Li₃(1/φ)
  swap = Li₂,₁(1/φ²,−φ) = 3ζ(3)/5 − 7 l_φ z₂/5 + 2 l_φ³ + Li₃(−1/φ) + Li₃(1/φ)
LLL at 830 digits (residuals 3.6e-832 / 1.8e-830). Cross-validated: Λ via the stuffle
Λ + Li₁₂(−1/φ,−1/φ) = Li₂(−1/φ)Li₁(−1/φ) − Li₃(1/φ²) with Li₁₂ from an independent double series
(3.1e-831); both Λ and swap via independent 700-dps quadratures of D1/D2 (1.3e-701 / 3.0e-702);
and finally by the whole-assembly check (1.1e-676).

## 4. Assembly and the ζ(3) report (Task 4)

Full assembly S = Σ (p+q√5)·letters·INT = 580 symbolic terms → after the reductions, the **15-term table**
above (basis l₂, l_φ, π²/6, ζ(3), Li₃(±1/φ)). Field: **ℚ(√5)**, one overall radical L.
- **Local ζ(3) coefficient = C₅, exactly.** The Jloc-only part carries ζ(3)-coefficient
  60 L (5/32 − √5/8) = √(10−2√5)·(105√5−255)/32 = **C₅** (symbolic difference 0; equals E1's Σ_e Z_e
  with Z_e = 6αβ(h²/L)(3(αu+βv)−2), α=β=1/(2φ²), u=v=1/2; matches D15b's 120-digit PSLQ value in all digits).
- **Total ζ(3) coefficient (this named basis):** 60 L (33/64 − 19√5/80) = L(495/16 − 57√5/4).
- **Non-local ζ(3) per edge:** N_e = (T − C₅)/10 = L·(69/32 − 27√5/40) = 3L(115−36√5)/160 ≈ 0.1299057313.
  (The local/non-local split is canonical — the local part carries no Li₃'s; the absolute total is
  **basis-dependent**: e.g. eliminating Li₃(1/φ), Li₃(−1/φ) via the sum-relation would change it.)

## 5. Conjecture for the paper (Task 5)

Data points for regular 2m-gons K_p = {Σ_{i=1}^m |⟨n_i,x⟩|^p ≤ 1}:
- **m = 2** (square): P'''(1) = 6√2(ζ(3) − π²/3) — elementary (E1, PROVED; checked here vs D15b, 2.3e-100).
- **m = 3** (hexagon): P'''(1) = (√3/8)(16log³2 + 3π²log2 − 96log²2 + 96log2 − 9ζ(3) − 6π²) — elementary (E1, PROVED).
- **m = 4** (octagon): P'''(1) = 48u·S, 35-term ℚ(√2)-combination of log2, log(1+√2), π², ζ(3), Li₂/Li₃ at
  six Landen-closed points; NO depth-2 (OCT6, PROVED assembly; Λ/swap reduce).
- **m = 5** (decagon, THIS WORK): P'''(1) = 60L·S, 15-term ℚ(√5)-combination of log2, log φ, π², ζ(3),
  Li₃(±1/φ); NO depth-2.

**What the pattern suggests (honest):**
(C1) For every regular 2m-gon, P'''(1) = 3!(2m)·L·S_m with S_m a combination over the maximal real
cyclotomic field ℚ(ζ_{2m})⁺ of: logs of algebraic numbers in the field, π², ζ(3), and SINGLE (depth-1)
polylogs Li₂, Li₃ at algebraic points in the field. The local ζ(3) coefficient is exactly the edge-law C_m.
(C2) **No multiple polylog of depth ≥ 2 survives in any regular 2m-gon.** Evidence: m = 4, 5 with two
Landen-paired non-local logs per edge — all Li₂,₁/Li₁₂/Li₁₁ reduce. CAVEAT (important): m = 5 did NOT test
the genuinely new triple-mixed integrals ∫tⁿ w_i w_j w_k — the decagon's third non-local term is constant,
so its inventory equals the octagon's. The first case with ≥ 3 NONCONSTANT non-local logs is **m = 6**
(four per edge, two reciprocal pairs), where Li₁,₁,₁-type constants can appear; the pattern predicts they
still reduce, but that is untested. (For m even, no term is fixed by the edge reflection, cf. §1.)
(C3) The number of independent Li₃ constants equals the number of Landen/orbit classes of the edge points
(1 for m=4's six reps after relations; 2 for m=5's {±1/φ}); the logs needed are logs of the field's
term-value units (log 2, log(1+√2) for m=4; log 2, log φ for m=5).

## 6. Evidence table

| check | residual |
|---|---|
| geometry: ratios {φ, 1/φ, 1}, a₅ constant (50 dps) | ≤ 6.7e-52 |
| exact E3n coefficients ℚ(√5) reconstruction | 1.6e-99 (max denom 384) |
| E3(0.3) vs direct d-series (80 dps) | 6.3e-82 |
| 60 L ∫E3n vs 200-digit P'''(1) (60-dps quad) | 9.3e-59 |
| P'''(1) tanh-sinh gold (P3_660), vs D15b 200 d | 2.0e-101 (all of its digits) |
| Jloc (40), Jtab at x₁,x₄ (40 each) vs quad (110 dps) | 2.9e-107, 3.2e-110, 9.1e-110 |
| D1–D10 closed forms vs quad (100–110 dps) | ≤ 1.8e-101 |
| 20 reduced mixed integrals vs quad (110 dps) | 8.4e-111 |
| A(−1/φ), A(1/φ²): LLL vs series | 5.2e-830, 2.9e-831 |
| A-atoms: orbit system (PROVED) vs LLL forms | exact (symbolic 0 mod lattice) |
| golden Li₃ lattice R_a, sum-rel, Lewin | 1.1e-830, 2.4e-831, 1.3e-831 |
| Λ, swap LLL reductions vs 830-d series | 4.9e-832, 3.6e-831 |
| Λ stuffle cross-check (Li₁₂ indep. series) | 3.1e-831 |
| Λ, swap via independent 700-d quad (D1, D2) | 1.3e-701, 3.0e-702 |
| local ζ(3) = C₅ | exact 0 |
| **final 15-term P'''(1) vs 680-digit gold** | **1.1e-676** |
| m=2 closed form vs D15b | 2.3e-100 |

## 7. Files (work dir DEC1_decagon)

step_23 (geometry numeric, ratios/constant term), step_25/26 (exact geometry, C₅ vs D15b), step_30/31
(100-dps series engine E3n_num100.pkl), step_33/34 (E3(0.3) check), step_35 (60 L ∫ vs 200 d), step_36
(ℚ(√5) reconstruction E3n_exact_pq.pkl), step_37 (basis change, mons.pkl: 261 monomials), step_39 (early
result.md), step_40 (680-digit tanh-sinh gold P3_660.pkl/txt), step_42 (D-constants + closed-form checks),
step_43 (Li₂,₁ series; D1 series check; D4 reflection), step_45 (830-d atoms atoms830.pkl), step_46/47/48
(LLL reductions rels_lll*.pkl; Li₃ lattice), step_50 (Li₃ inversion), step_52 (A-atom orbit system PROOF),
step_54 (stuffle cross-checks, li12_830.pkl), step_55 (D1/D2 700-d cross-checks), step_56/57 (mixed engine,
20 integrals 8.4e-111; dec_engine.py), step_58 (Jloc/Jtab verification), step_59 (assembly S_assembly.pkl,
7.3e-108), step_60 (reductions S_reduced.pkl), step_62 (canonical 15-term table final_table_pq.pkl),
step_64 (D3/D4 explicit), step_65/66 (ζ(3) analysis), step_67 (final 1.1e-676), step_68 (m=2 check,
gold vs D15b). progress.log, result.md.
