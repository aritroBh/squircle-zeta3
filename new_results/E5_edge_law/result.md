# E5: the ζ(5) edge law for smoothed polygons

Work directory: `research/08_theory/E5_edge_law/` (steps in `steps/step_NN.py`, logs `steps/*_out.txt`, `progress.log`). Builds on E1 (`10_frontier/E1_edge_law`): K_p = {Σ_j w_j|⟨n_j,x⟩|^p ≤ 1}, d = p−1, edge data (α, β, L, h, s0, u = −s0/L, v = 1−u), |x_p′| = F^{−1/p}√((L−σΦ)²+h²Φ²), P^{(k)}(1) = k!Σ_e∫₀¹E_k dt. **PROVED** = exact symbolic identity; **NUMERICAL** = high-precision agreement (residual printed).

## 1. The ζ(5) edge law (PROVED)

The local (corner–corner) terms of an edge e contribute Z_e⁽⁵⁾ζ(5) + Z_e^{(2,3)}ζ(2)ζ(3) to P⁽⁵⁾(1). With w = αu+βv, e = αβuv, q = α²u²+β²v², Q = α²+β², p = αβ, s = α+β, r = h²/L²:

$$Z_e^{(5)} = 360\,\alpha\beta\,\tfrac{h^2}{L}\Bigl[20w^3-20ew-32q-48e+15w-2 + r\bigl(-15Qw+15pw+5ps+8Q-12p\bigr)\Bigr]$$
$$Z_e^{(2,3)} = 120\,\alpha\beta\,\tfrac{h^2}{L}\Bigl[-60ew+12q+72e-12w+2 + r\bigl(-45pw+15ps-3Q+18p\bigr)\Bigr]$$

**Derivation chain (each step verified).** (i) `step_01.py`: E1's s01_edge_sym.py extended to N=6; E5 = 475 terms, log-degree 5; cross coefficients g_ab = [l0^a l1^b]E5 (a+b=5) are degree-5 polynomials in t. (ii) `step_02.py`: full E5 with lumped τ_k, γ_k (k=1..5): 1910 terms, counts by log-degree {1:8,2:37,3:70,4:74,5:56}, 35 (τ,γ) patterns; g_ab are τ/γ-free (any τ/γ factor consumes a d-slot) and equal the reduced ones; full E5 validated vs Cauchy d-series on a decagon edge (≤3.1e-17, `step_09.py`). (iii) `step_03.py`: Beta table to weight 5 (ζ4=2ζ2²/5, γ cancels; max numeric error 4.8e-35 over 210 entries) with rules [ζ5]J = 24B(n+1,k+1) iff a+b=5,a,b≥1 else 0; [ζ2ζ3]J = −12B iff (a,b)=(3,2)/(2,3) else 0 — hence Z_e⁽⁵⁾ = 5!·24∫(g41+g32+g23+g14) = 2880∫G and Z_e^{(2,3)} = −5!·12∫(g32+g23) = −1440∫(g32+g23), plain polynomial integrals; μ-shifts and lower-degree terms give only weight ≤4 constants × log-monomials, so the laws are convention-independent. (iv) `step_04/05.py`: symmetrized closed forms above, law−raw difference exactly 0, reversal-invariant (α,u)↔(β,v) exactly.

**Parallelograms (PROVED; fully local ⇒ local = total, α=β=1):**
P⁽⁵⁾(1)|_ζ(5) = Σ_e 360(h²/L)(1−4uv−r),  P⁽⁵⁾(1)|_ζ(2)ζ(3) = Σ_e 120(h²/L)(2−12uv−3r).
Complete closed forms per edge (power convention, `step_12.py`): F_2 = (h²/L)ζ2; F_3 = (h²/L)(ζ3−2ζ2) (E1 ✓); F_4 = ζ2²[(9/5)(h²/L)+(42/5)(h²s0/L²)−(21/10)(h⁴/L³)+(42/5)(h²s0²/L³)] + 3(h²/L)(ζ2−ζ3); F_5 = ζ2²[−(36/5)(h²/L)−(168/5)(h²s0/L²)+(42/5)(h⁴/L³)−(168/5)(h²s0²/L³)] + ζ5[3(h²/L)+12(h²s0/L²)−3(h⁴/L³)+12(h²s0²/L³)] + ζ2ζ3[2(h²/L)+12(h²s0/L²)−3(h⁴/L³)+12(h²s0²/L³)] + (h²/L)(6ζ3−4ζ2); P^{(k)}(1) = k!Σ_4 edges F_k.

**Regular 2m-gons (PROVED symbolic in h=π/(2m)):** Σ_eZ_e⁽⁵⁾ = C_m⁽⁵⁾ = 2880m s⁴c(65s⁶−63s⁴+17s²−1), Σ_eZ_e^{(2,3)} = C_m^{(2,3)} = −960m c s⁴(45s⁶−27s⁴+6s²−1) (difference exactly 0 after c²=1−s²); per edge Z_e⁽⁵⁾ = 1440s⁴c(65s⁶−63s⁴+17s²−1), Z_e^{(2,3)} = −480cs⁴(45s⁶−27s⁴+6s²−1). Diamond m=2: −90√2 ζ(5) − 210√2 ζ(2)ζ(3) exactly.

**Non-local caveat (same as ζ(3)):** on a general polygon the τ/γ and nonconstant-log terms integrate to convergent Beta series of weight ≤5 that may add more ζ(5)/ζ(2)ζ(3); the defect is exhibited numerically (`step_09.py`, decagon: local-reduced vs full E5 differ by 0.04–0.15). The series stop exactly for parallelograms and affine-regular hexagons.

## 2. Evidence table

| Check | Object | Result | Residual | Digits |
|---|---|---|---|---|
| E5 reduced build (step_01) | 475 terms, deg 5; g_ab deg-5 polys | exact | 0 | — |
| Full E5 τ/γ (step_02) | 1910 terms; g_ab τ/γ-free, == reduced | exact | 0 | — |
| Full E5 numeric (step_09) | decagon edge, 5 t-values | matches Cauchy d-series | ≤ 3.1e-17 | 60 dps |
| Beta table (step_03) | 210 moments vs mp.quad; ζ5/ζ2ζ3 rules | max err 4.8e-35; rules hold everywhere | 4.8e-35 | 40 |
| Law vs raw integrals (step_05) | Z_e⁽⁵⁾, Z_e^{(2,3)} | exact | 0 | — |
| Regular 2m-gon (step_05) | per-edge & Σ vs D17 C_m | exact (symbolic in h); m=2..8 numeric | ≤ 6.8e-49 | 50 |
| Diamond m=2 (step_05) | −90√2, −210√2 | exact | 0 | — |
| Sheared parallelogram (step_06) | exact P^{(2..5)}; ζ5,ζ2ζ3 coords == ΣZ_e; E1 P‴ law | exact | 0 | — |
| Independent numeric (step_07) | per-t Cauchy + tanh-sinh, no series used | P⁽⁵⁾ diff | 7.8e-57 | 100 dps |
| Rhombus 1,2 exact (step_08) | ζ5 = 1152√5/25, ζ2ζ3 = −768√5/25 == 4Z_e | exact | 0 | — |
| Rhombus numeric (step_10) | same machine | P⁽⁵⁾ diff | 6.9e-57 | 100 dps |
| Rhombus PSLQ (step_11) | P⁽⁵⁾/√5 on {1,ζ2,ζ3,ζ2²,ζ2ζ3,ζ5} | relation (−125,0,−38400,57600,−4608,−3840,5760); ζ5 coordinate 1152√5/25 = Σ_eZ_e⁽⁵⁾ | 3.1e-57 | 100 dps |
| General parallelogram F_k (step_12) | F_2, F_3 == E1; F_4, F_5 vs sheared exact | exact | 0 | — |

Key exact values — sheared parallelogram ±(2,0),±(1,2): P⁽⁵⁾(1) = ζ5(−6912√5/25 − 172800√13/28561) + ζ2ζ3(−10752√5/25 − 821760√13/28561) + ζ2²(119808√5/125 + 1262592√13/28561) + ζ3(4608√5/5 + 23040√13/169) + ζ2(−3072√5/5 − 15360√13/169). Rhombus |x|^p+|y/2|^p≤1: P⁽⁵⁾(1) = √5(1152ζ5/25 − 768ζ2ζ3/25 − 4608ζ2²/125 + 2304ζ3/5 − 1536ζ2/5).

## 3. Files that produced these

`steps/step_01.py` (E5 reduced N=6, cross coefficients, 2.9s) · `step_02.py` (full E5 with τ₁..₅,γ₁..₅, 9.4s) · `step_03.py` (weight-5 Beta table + ζ5/ζ2ζ3 rules, 159s) · `step_04.py` (raw edge-law integrals ∫G, ∫(g32+g23)) · `step_05.py` (closed forms; regular 2m-gon symbolic-in-h proof; m=2..8; diamond) · `step_06.py` (exact sheared parallelogram, 44s) · `step_07.py` (independent 100-dps numeric, adaptive-radius branch-tracked per-t Cauchy + tanh-sinh with node-reuse halving, 925s) · `step_08.py` (exact rhombus, 45s) · `step_09.py` (full-E5 decagon validation + non-local defect) · `step_10.py` (rhombus numeric, 310s) · `step_11.py` (PSLQ) · `step_12.py` (general parallelogram F_k, 48s). Outputs: `steps/*_out.txt`, `*.pkl`, `*_num.json`; `result.md`; `progress.log`.

## 4. Statement (task item 4)

For **every parallelogram**, the ζ(5) coefficient of P⁽⁵⁾(1) is Σ_e 360(h²/L)(1−4uv−h²/L²) and the ζ(2)ζ(3) coefficient is Σ_e 120(h²/L)(2−12uv−3h²/L²) — **PROVED** (items 1–3 hold; parallelograms are fully local so the local integrand is the whole integrand, and the exact integration of the symbolic E5 is a proof, independently confirmed to 55 digits and by PSLQ). For **general polygons**, Z_e⁽⁵⁾ and Z_e^{(2,3)} as boxed in §1 are the canonical **local** ζ(5)/ζ(2)ζ(3) edge laws — with E1's ζ(3) caveat: non-local terms may add more; whenever all edges are parallel to the kinks they do not touch (parallelograms, affine-regular hexagons), the local laws are the whole story.