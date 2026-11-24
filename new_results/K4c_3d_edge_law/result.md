# K4c: a 3D local law for the ζ(3) term of smoothed polyhedra — final report

**Folder:** `research/08_theory/K4c_3d_edge_law/`. Engine: copied from K4b (`k4_engine.py`, `k4_quad.py`) plus a new generic driver `k4c_driver.py` and `pslq_lib.py`. dps 70–80 production, tanh–sinh h up to 1/96. Labels: **PROVED** / **NUMERICAL** / **CONJECTURE** / **FAILED-IN-BUDGET**.

---

## 0. Required first step — engine check (PROVED machinery, NUMERICAL checks)

`step_01.py` (h=1/32, dps 70): octahedron [d⁰]S = 4√3, residual **3.9e-60**; cuboctahedron [d⁰]S = (9+3√3)/4, residual **1.4e-60** (region A vs 3√3/4: 4.3e-61; region B vs 9/4: 9.5e-61). Engine works.

**New driver development (`step_04*.py`, `k4c_driver.py`):** I integrate F(u) = Φ^{−2/(1+d)}√(1+W) over the O_h chamber {x≥y≥z≥0} (φ∈[0,π/4], θ∈[atan(1/sinφ),π/2], mult 48), with generic sign-pattern breakpoints per normal. Two discoveries: (i) the outer φ-integral must be **split at every latitude where a breakpoint θ_C(φ) enters the chamber** (φ₀ = atan(1/2) for cubo/model B) — otherwise an interior log-cusp caps accuracy at ~1e-5; (ii) after splitting, the accuracy bottleneck is the **outer** grid, not the inner one (inner refinement 1/64→1/96 changes nothing). Final recipe: outer 1/48–1/96, inner 1/64, WMIN = 1e-66, LMIN = 1e-58. Validation (`step_04b/d`): octa k=0..3 residuals 4.9e-62/9.1e-61/1.2e-58/8.0e-57 vs the exact H1 expansion; cubo k=0..4 residuals 4.9e-60/3.0e-60/1.4e-58/1.0e-56/4.7e-55 vs K4b's certified `tetra_Sdk_final.json`. **Driver PROVED in the K4b sense (a priori bounds) and NUMERICALLY validated to 55–60 digits.**

---

## 1. Task 1 — the two new configurations and their [dᵏ]S (k=0 exact, k≤3 ≥40 digits)

### Model A — 6 face-diagonal normals (±1,±1,0)/√2 etc. → **K₁ = tetrakis hexahedron** (PROVED)

Z = Σ[−n_j,n_j] is the A₃≅D₃ zonotope = **truncated octahedron with edge 2** (24 V, 14 F: 6 squares + 8 regular hexagons). K₁ = Z° is the Catalan **tetrakis hexahedron**: 14 vertices (6 valence-4 at (±1/(2√2),0,0)-type; 8 valence-6 at (±1,±1,±1)/(3√2)), 24 congruent isosceles-triangle faces (area √5/36 each), 36 edges: **24 of length 1/(2√2)** (apex–base) + **12 of length √2/3** (base–base), and — Catalan property, verified two independent ways — **all dihedral angles = arccos(−4/5) = 2.498091544796509**. S(1) = 24·√5/36 = **2√5/3** (sympy exact; engine residual 4.0e-62).

Coefficients (h=1/96, errors from h32→h48→h64→h96 trend ratios ~5e5–4e6 per halving):

| k | [dᵏ]S (model A) | err est | digits |
|---|---|---|---|
| 0 | 1.490711984999859797606115779154184156960412239741017149514 | 4.1e-62 (= 2√5/3) | 61 |
| 1 | 4.310068730207151345980382494723116602408017857946782637185 | 1.9e-49 | 48 |
| 2 | 2.175841230811692460507624858609765605772129518494031707693 | 9.7e-48 | 47 |
| 3 | **−2.330156241024368053549557284389126408582825209655730929953** | 2.0e-46 | 45 |
| 4 | 0.004482674056558090243652254900043497274132120333173616794591 | 1.6e-45 | 44 |

NUMERICAL (45+ digits), files `modelA_h{32,48,64,96}.json`.

### Model B — 3 coordinate + 4 body-diagonal normals → **K₁ = polar of the cube⊕rhombic-dodecahedron zonotope** (PROVED)

Z = cube(edge 2) ⊕ RD (7 zones): 32 vertices, 18 facets (6 squares + 12 hexagons) — combinatorially the **chamfered cube** (24 vertices (4.6.6), 8 (6.6.6); facet directions from the 9 arrangement lines: 3 axial → squares, 6 face-diagonal → hexagons). K₁ = Z°: **18 vertices** — 6 valence-4 at ±a e_i, a = (4√3−3)/13, and 12 valence-6 at permutations of (0,±c,±c), c = (2√3−3)/2 — **32 triangular faces** (24 type-1 with sides L₁,L₁,L₂; 8 type-2 equilateral with side L₂), **48 edges**: **24 axis–diag edges** of length L₁ = √((2805−1608√3)/338) = 0.24241329462710659 with dihedral **θ₁ = arccos(−(283+48√3)/433) = 2.5784569556559880**, and **24 diag–diag edges** of length L₂ = √(42−24√3)/2 = (2√6−3√2)/2 = 0.32816939922353554 with dihedral θ₂, cos θ₂ = (−1152+665√3)/(3√(6889627−3977728√3)), cos²θ₂ = (1497+768√3)/3897, **θ₂ = 2.5901203971413503**. S(1) = **−36 + 21√3 + (9/13)√(35515−20504√3)** = 1.07576493049084851987376558708427009279065569782681829776786 (sympy exact; engine matches to ~60 digits).

Coefficients (outer 1/48; **an outer-1/96 rerun gave identical values to all 58 printed digits** for k=1..4, so error ≲ 1e-55; k=0 exact):

| k | [dᵏ]S (model B) | err est | digits |
|---|---|---|---|
| 0 | 1.075764930490848519873765587084270092790655697826818297768 | < 1e-59 (exact above) | 59 |
| 1 | 3.516495437145298347373890733252633379646687669001453640008 | < 1e-55 | 55 |
| 2 | 2.354480678056587151344461624642939462487721480650698816991 | < 1e-55 | 55 |
| 3 | **−1.898327187347440896619891901482403412469397080520072254548** | < 1e-55 | 55 |
| 4 | −0.4126500481804702236949455301550953680647740676969757841329 | < 1e-55 | 55 |

NUMERICAL (~55 digits), files `modelB_h{32,48,96}.json`. (Tiny imaginary artifacts ~1e-112 at breakpoint nodes — K4b-type rounding, 45 orders below the error floor; removed before saving.)

---

## 2. Task 2 — normalization and coefficient identification (partial)

- **Model A (NUMERICAL, strong):** [d]S/[d⁰]S = q_A satisfies **15 q_A = 13 − 18 ln 2 + 39 ln 3**, i.e. **q_A = 13/15 + (1/5) ln(3¹³/2⁶)** — PSLQ (`pslq_lib.py`, own Bailey–Broadhurst implementation, unit-tested) residual **4.7e-53** at dps 70, relation height 39 in a 4-dim lattice (spurious level ~1e17, so significant). With **c_A = 2^{−6/5}·3^{13/5}**, S(1+d) = c_A^{d/(1+d)}·Σt_k dᵏ gives **t₁ = (13/15)S₀** (verified: t₁/S₀ = 0.8666… = 13/15 to 48 digits). This mirrors cubo's c=4, t₁ = S₀/2.
- **Model B:** q_B = [d]S/[d⁰]S = 3.2688325650667910058383550480616858794415364290159 was **not** identified: no small-height relation on {1, ln2, ln3, ln5, ln7, ln13, ln(2±√3), ln(4±√3), ln(1+√2), ln(√2+√3), Catalan G, Cl₂-type} (all found relations have height ≥ 1e8, i.e. spurious at 65 digits). Model B's [d]S apparently already leaves the small-log world.
- **[d²], [d³] full identification: FAILED-IN-BUDGET.** On the hinted bases ({1,√2,√3,√6}⊗{1,π²,ζ(3)}, also with √5, √10, and ln-power extensions, with and without the c_A normalization, and on S_k/√5 rational-log spans), every candidate relation at 45–55 digits has residual **above** the data error (spurious height range for 12–15-dim lattices) or is not reproducible across the independent h48/h64/h96 grids. Concretely: model A would need ≳60 digits at 12-dim; model B most likely needs its (unidentified) q_B log-constants in the basis first. **The ζ(3) coefficients of [d³]S for models A and B were therefore not isolated.** (For the record: octa −64√3/9 and cubo −√3/6 are planner-verified.)

---

## 3. Task 3 — the edge law (framework + two exact anchors; test incomplete)

- Counts-only law 24Z_e+12Z_v was already **REFUTED** (K4b): octa/cubo have identical E, V ratios but c_ζ = −64√3/9 ≈ −12.317 vs −√3/6 ≈ −0.289.
- Assembled edge data (all PROVED): octa 12×(△△, θ=arccos(−1/3)=1.9106332362, ℓ=√2); cubo 24×(△□, θ=arccos(1/√3)=0.9553166181, ℓ=√(3/8)); model A 24×(ℓ=1/(2√2)) + 12×(ℓ=√2/3), all θ=arccos(−4/5); model B 24×(ℓ₁, θ₁) + 24×(ℓ₂, θ₂) with θ₁, θ₂ above.
- Under the ansatz [d³]S|ζ(3) = Σ ℓ_e F(θ_e, face/vertex data): the two exact anchors give F(θ_oct, △△) = (−64√3/9)/(12√2) = **−8√6/27 ≈ −0.7257747** and F(θ_cubo, △□) = (−√3/6)/(24√(3/8)) = **−√2/72 ≈ −0.0196419**, ratio 64√3/3 ≈ 36.95. Since θ_oct = 2θ_cubo (tan(θ_oct/2) = tan θ_cubo = √2) yet F changes 37-fold while face types change △△→△□, the data are consistent with the 2D E1 structure (Z_e depends on the two end-corner sharpnesses α, β multiplicatively): **any 3D local law must depend on face/vertex data beyond the dihedral angle** (CONJECTURE, anchored on 2 exact values).
- Models A (two edge orbits at the same θ) and B (two orbits at distinct θ) were designed to separate these dependencies — but their c_ζ values could not be extracted within budget (§2), so the 4-point law fit and its residuals remain open.

---

## 4. Evidence table

| claim | residual / digits | status | evidence |
|---|---|---|---|
| octa [d⁰]S = 4√3 | 3.9e-60 | NUMERICAL (engine PROVED) | step_01.py |
| cubo [d⁰]S = (9+3√3)/4 | 1.4e-60 | NUMERICAL | step_01.py |
| driver vs octa exact k=1..3 | ≤ 8.0e-57 | NUMERICAL | step_04b.py |
| driver vs K4b certified cubo k=0..4 | ≤ 4.7e-55 | NUMERICAL | step_04b/04d.py |
| A: K₁ = tetrakis hexahedron; S(1) = 2√5/3; edges 24×1/(2√2)+12×√2/3; θ = arccos(−4/5) | exact (sympy), engine k=0 4.0e-62 | PROVED | step_02/03*.py, modelA_h96.json |
| B: K₁ = polar of cube⊕RD (chamfered-cube dual); S(1) = −36+21√3+(9/13)√(35515−20504√3); L₁, L₂, θ₁, θ₂ as above | exact (sympy), engine k=0 ~60 digits | PROVED | step_02/03*.py, modelB_h48.json |
| A: [dᵏ]S k=1..3 to 45–48 digits | errs 1.9e-49/9.7e-48/2.0e-46 | NUMERICAL | step_05/06/07.py |
| B: [dᵏ]S k=1..3 to ~55 digits | outer 1/96 ≡ 1/48 to 58 digits | NUMERICAL | step_08/09.py |
| A: q_A = 13/15 + (1/5)ln(3¹³/2⁶); c_A = 2^{−6/5}3^{13/5} ⇒ t₁ = 13S₀/15 | PSLQ resid 4.7e-53 | NUMERICAL (strong) | step_10.py, pslq_lib.py |
| A,B: [d²],[d³] on hinted bases; ζ(3) coefficients | all candidates residual > data error | FAILED-IN-BUDGET | step_10.py runs |
| 3D law must use face/vertex data, not θ alone | anchored on 2 exact c_ζ | CONJECTURE | this report §3 |

## 5. Files

`step_01.py` (engine checks), `step_02.py` (hull/polar analysis), `step_03.py/03b/03c` (exact sympy geometry A & B), `step_04b/04c/04d.py` (driver v2 + validation), `k4c_driver.py` (production driver with outer-split logic), `step_05/06/07.py` (model A grids 1/32→1/96), `step_08/09.py` (model B grids), `pslq_lib.py` (tested PSLQ), `step_10.py` + inline runs (q-analysis, normalization, coefficient PSLQ), data `modelA_h{32,48,64,96}.json`, `modelB_h{32,48,96}.json`, `K1V_{A,B}.npy`, `progress.log` (full audit trail).

**Bottom line:** Task 1 complete (both polyhedra identified exactly; [dᵏ]S to 45–55 digits with exact k=0). Task 2 half-complete (model A normalized via c_A = 2^{−6/5}3^{13/5} with t₁ = 13S₀/15; model B's q unidentified; [d²]/[d³] and the two ζ(3) coefficients need ≳60 digits for A and a wider constant basis for B). Task 3: machinery and exact edge data for all four polyhedra in place; the two exact anchors already force face/vertex-dependence in any local law; final 4-point fit awaits the two missing c_ζ values — a straightforward continuation: model A at h=1/128 (~55 digits) plus PSLQ with the q_A-logs, and a q_B hunt over weight-2 constants before attacking its [d³].