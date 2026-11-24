result = r'''# E5: the ζ(5) edge law for smoothed polygons

Folder: `research/08_theory/E5_edge_law/`. Read after E1 (`research/10_frontier/E1_edge_law/result.md`),
whose notation and lemmas are used: K_p = {Σ_j w_j|<n_j,x>|^p ≤ 1}, d = p−1, edge data
(α, β, L, h, s0, u = −s0/L, v = 1−u), radial-projection integrand |x_p′| = F^{−1/p}√((L−σΦ)²+h²Φ²),
P^{(k)}(1) = k! Σ_e ∫₀¹ E_k(t) dt, "fully local" = every edge parallel to the kinks it does not touch.
Labels: **PROVED** = exact symbolic identity (residual exactly 0, printed); **NUMERICAL** = high-precision
agreement (residual printed).

## 0. Bottom line

1. **The ζ(5) edge law (PROVED).** The local (corner–corner) terms of an edge e contribute
   Z_e^(5)·ζ(5) + Z_e^(2,3)·ζ(2)ζ(3) to P^{(5)}(1). With w = αu+βv, e = αβuv, q = α²u²+β²v²,
   Q = α²+β², p = αβ, s = α+β, r = h²/L²:
   $$Z_e^{(5)} = 360\,\alpha\beta\,\tfrac{h^2}{L}\Bigl[20w^3 - 20ew - 32q - 48e + 15w - 2
      + r\bigl(-15Qw + 15pw + 5ps + 8Q - 12p\bigr)\Bigr],$$
   $$Z_e^{(2,3)} = 120\,\alpha\beta\,\tfrac{h^2}{L}\Bigl[-60ew + 12q + 72e - 12w + 2
      + r\bigl(-45pw + 15ps - 3Q + 18p\bigr)\Bigr].$$
   Derivation and proofs in §1–§2 (`steps/step_01..05.py`). Both are invariant under edge reversal
   (α,u)↔(β,v) (exact). Same caveat as E1's ζ(3) law: on a general polygon the non-local terms may add
   more (§4); the law is the canonical local ζ(5)/ζ(2)ζ(3) content of an edge.
2. **Parallelograms (PROVED).** Fully local: local = total, α = β = 1, so
   $$P^{(5)}(1)\big|_{\zeta(5)} = \sum_e 360\,\tfrac{h^2}{L}\bigl(1 - 4uv - r\bigr),\qquad
     P^{(5)}(1)\big|_{\zeta(2)\zeta(3)} = \sum_e 120\,\tfrac{h^2}{L}\bigl(2 - 12uv - 3r\bigr).$$
   Complete P^{(5)}(1) (and P^{(4)}(1)) for every parallelogram in §2.2 (`steps/step_12.py`).
3. **Regular 2m-gons (PROVED, symbolic in h = π/(2m)).** Σ_e Z_e^(5) and Σ_e Z_e^(2,3) equal D17's
   C_m^(5) = 2880 m s⁴c(65s⁶−63s⁴+17s²−1) and C_m^(2,3) = −960 m c s⁴(45s⁶−27s⁴+6s²−1)
   (s = sin h, c = cos h): difference exactly 0 after c² = 1−s² reduction; numeric for m = 2..8
   ≤ 6.8e-49. Diamond m = 2: −90√2 ζ(5) − 210√2 ζ(2)ζ(3) exactly (`steps/step_05.py`).
   Per edge: Z_e^(5) = 1440 s⁴c(65s⁶−63s⁴+17s²−1), Z_e^(2,3) = −480 c s⁴(45s⁶−27s⁴+6s²−1).
4. **Non-square parallelograms, independent checks (NUMERICAL 55 digits + PSLQ).**
   - Sheared parallelogram ±(2,0), ±(1,2): exact P^{(5)}(1) in §3; ζ5- and ζ2ζ3-coordinates equal the
     law sums with difference exactly 0. Independent numeric (per-t Cauchy d-derivative + tanh-sinh,
     no symbolic series used): |num − exact| = 7.8e-57 (P^{(5)}), 2.7e-59 (P^{(4)}), 4.7e-61 (P^{(3)}),
     2.5e-63 (P^{(2)}) (`steps/step_06.py`, `steps/step_07.py`, dps 100).
   - Rhombus |x|^p+|y/2|^p ≤ 1: exact P^{(5)}(1) = √5(1152 ζ5/25 − 768 ζ2ζ3/25 − 4608 ζ2²/125
     + 2304 ζ3/5 − 1536 ζ2/5); numeric |num − exact| = 6.9e-57; PSLQ of P^{(5)}(1)/√5 on
     {1,ζ2,ζ3,ζ2²,ζ2ζ3,ζ5} recovers the exact integer relation (−125,0,−38400,57600,−4608,−3840,5760),
     residual 3.1e-57; ζ(5) coordinate 1152√5/25 = Σ_e Z_e^(5), ζ(2)ζ(3) coordinate −768√5/25 =
     Σ_e Z_e^(2,3) (`steps/step_08.py`, `step_10.py`, `step_11.py`).

## 1. Derivation (extends E1 §1–§2 to weight 5)

- `steps/step_01.py`: E1's s01_edge_sym.py extended to N=6 (d⁰..d⁵), reduced (τ=γ=0). E5 has 475
  terms, log-degree 5; cross coefficients g_ab = [l0^a l1^b]E5, a+b=5, are polynomials in t of degree 5.
- `steps/step_02.py`: full N=6 run with the lumped non-local symbols τ_k = Σ' a_i ℓ_i^k,
  γ_k = Σ' a_i′ ℓ_i^k (k = 1..5): E5 has 1910 terms, counts by log-degree {1:8, 2:37, 3:70, 4:74,
  5:56}, 35 distinct (τ,γ) patterns. g41, g32, g23, g14 are τ/γ-free and equal the reduced ones
  (printed): any τ/γ factor consumes one of the five d-slots, so it cannot multiply a degree-5 l0/l1
  monomial. Full E5 validated numerically against the Cauchy d-series on a regular-decagon edge
  (non-fully-local): ≤ 3.1e-17 at five t-values (`steps/step_09.py`).
- `steps/step_03.py`: Beta table J(n,k,a,b) = ∫₀¹ tⁿ(1−t)^k log^a t log^b(1−t) dt, a+b≤5, exact in
  {1, ζ2, ζ3, ζ2², ζ2ζ3, ζ5} (ζ(4) = 2ζ2²/5); Euler γ cancels identically; numeric max error 4.8e-35
  over all 210 entries. Rules (from the log-Beta expansion, verified on the whole table):
  [ζ5]J = 24 B(n+1,k+1) iff a+b=5, a,b≥1, else 0; [ζ2ζ3]J = −12 B(n+1,k+1) iff (a,b) ∈ {(3,2),(2,3)},
  else 0. Hence for polynomial c(t), [ζ5]∫₀¹c(t)X^aY^b = 24∫c (a+b=5, a,b≥1) and
  [ζ2ζ3]∫₀¹c(t)X^aY^b = −12∫c ((a,b) = (3,2)/(2,3)).
- Consequence: with G = g41+g32+g23+g14,
  Z_e^(5) = 5!·24∫G = 2880∫G,  Z_e^(2,3) = −5!·12∫(g32+g23) = −1440∫(g32+g23) — plain polynomial
  integrals. The shifts l0 = log t+μ0, l1 = log(1−t)+μ1 and all lower-degree terms contribute only
  weight ≤4 constants times log-monomials: never bare ζ5 or ζ2ζ3. So the two laws depend only on
  (α,β,u,v,h,L), not on the weight convention (fixed weights or power convention).

## 2. Closed forms

### 2.1 The edge law (§0.1)
∫G and ∫(g32+g23) are αβ(h²/L⁴)×(polynomials in L,α,β,s0,h²) (`steps/step_04_out.txt`); symmetrised
with u = −s0/L, v = 1+s0/L to the boxed formulas; law−raw difference exactly 0 (`steps/step_05.py`).

### 2.2 Whole P^{(k)}(1) for every parallelogram (power convention Σ|<v_j,x>|^p, α=β=1, μ=0)
Per edge, F_k = ∫₀¹E_k dt (`steps/step_12.py`; F_2, F_3 reproduce E1):
F_2 = (h²/L)ζ2;  F_3 = (h²/L)(ζ3−2ζ2);
F_4 = ζ2²[(9/5)(h²/L) + (42/5)(h²s0/L²) − (21/10)(h⁴/L³) + (42/5)(h²s0²/L³)] + 3(h²/L)(ζ2−ζ3);
F_5 = ζ2²[−(36/5)(h²/L) − (168/5)(h²s0/L²) + (42/5)(h⁴/L³) − (168/5)(h²s0²/L³)]
      + ζ5[3(h²/L) + 12(h²s0/L²) − 3(h⁴/L³) + 12(h²s0²/L³)]
      + ζ2ζ3[2(h²/L) + 12(h²s0/L²) − 3(h⁴/L³) + 12(h²s0²/L³)] + (h²/L)(6ζ3−4ζ2);
and P^{(k)}(1) = k! Σ_{4 edges} F_k. The ζ5/ζ2ζ3 brackets are 5!^{-1}×(the §0.2 law), as they must be.
No rational (ζ-free) term appears in any F_k.

## 3. The sheared parallelogram ±(2,0),±(1,2) (exact, `steps/step_06.py`)

P^{(0)} = 2√5+2√13; P^{(1)} = 0; P^{(2)} = ζ2(64√5/25 + 64√13/169);
P^{(3)} = ζ2(−384√5/25 − 384√13/169) + ζ3(192√5/25 + 192√13/169)  [= 6S(ζ3−π²/3), E1's law ✓];
P^{(4)} = ζ2²(−29952√5/625 − 315648√13/142805) + ζ2(2304√5/25 + 2304√13/169) + ζ3(−2304√5/25 − 2304√13/169);
P^{(5)} = ζ5(−6912√5/25 − 172800√13/28561) + ζ2ζ3(−10752√5/25 − 821760√13/28561)
          + ζ2²(119808√5/125 + 1262592√13/28561) + ζ3(4608√5/5 + 23040√13/169) + ζ2(−3072√5/5 − 15360√13/169).
ζ5 and ζ2ζ3 coordinates = Σ_e Z_e^(5), Σ_e Z_e^(2,3): difference exactly 0.

## 4. The non-local caveat

On a non-fully-local edge the τ_k, γ_k and the nonconstant parts log|1+(r_i−1)t| of the non-local logs
integrate to convergent series of Beta moments of weight ≤5, which can add further ζ(5)/ζ(2)ζ(3) beyond
Σ_e Z_e — exactly E1's ζ(3) situation. Whenever every edge is parallel to the kink lines it does not
touch (parallelograms, affine-regular hexagons), the series stop and §2.2/E1-style closed forms are the
whole P^{(5)}(1). For regular 2m-gons the local law reproduces D17's perturbative ζ(5) coefficients
exactly for all m (§0.3); for m ≥ 4 the phrase "ζ(5) coefficient of P^{(5)}(1)" carries E1's
basis-choice caveat for the cyclotomic remainder.

## 5. Files

steps/step_01.py (E5 reduced, cross coefficients; 2.9 s), step_02.py (full E5 with τ,γ; 9.4 s),
step_03.py (Beta table to weight 5, rules; 159 s), step_04.py (raw edge-law integrals),
step_05.py (closed forms + regular 2m-gon + diamond proofs), step_06.py (exact sheared parallelogram;
44 s), step_07.py (independent 55-digit numeric sheared parallelogram; 925 s), step_08.py (exact
rhombus), step_09.py (full-E5 numeric validation on decagon edge), step_10.py (numeric rhombus; 310 s),
step_11.py (PSLQ), step_12.py (general parallelogram F_k; 48 s). Outputs: steps/*_out.txt, *.pkl/json.
progress.log: run history.
'''
open('result.md','w').write(result)
open('progress.log','a').write("step_12: general parallelogram F_2..F_5 (power convention); F_2=z2 h^2/L, F_3=(z3-2z2)h^2/L (E1), F_4/F_5 verified vs sheared exact. result.md finalized.\n")
print("result.md finalized")
