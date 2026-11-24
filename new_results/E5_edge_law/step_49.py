result = r'''# E5: the ζ(5) edge law for smoothed polygons

Folder: `research/08_theory/E5_edge_law/`. Read after E1 (`research/10_frontier/E1_edge_law/result.md`)
whose notation (α, β, u, v, h, L, s0, radial projection, fully local) is used throughout.
Labels: **PROVED** = exact symbolic identity (or printed exact check); **NUMERICAL** = high-precision
agreement, residual printed.

## 0. Bottom line

1. **The ζ(5) edge law (PROVED).** For K_p = {Σ_j w_j|<n_j,x>|^p ≤ 1}, d = p−1, the local
   (corner–corner) terms of an edge e contribute Z_e^(5) ζ(5) + Z_e^(2,3) ζ(2)ζ(3) to P^{(5)}(1), with
   w = αu+βv, e = αβuv, q = α²u²+β²v², Q = α²+β², p = αβ, s = α+β, r = h²/L²:
   $$Z_e^{(5)} = 360\,\alpha\beta\,\frac{h^2}{L}\Bigl[\,20w^3 - 20ew - 32q - 48e + 15w - 2
      \;+\; r\,\bigl(-15Qw + 15pw + 5ps + 8Q - 12p\bigr)\Bigr],$$
   $$Z_e^{(2,3)} = 120\,\alpha\beta\,\frac{h^2}{L}\Bigl[\,-60ew + 12q + 72e - 12w + 2
      \;+\; r\,\bigl(-45pw + 15ps - 3Q + 18p\bigr)\Bigr].$$
   (`steps/step_04.py`, `steps/step_05.py`: identical to the raw integrals of the cross coefficients,
   difference 0; reversal-invariant under (α,u)↔(β,v).) Same caveat as ζ(3): for a general polygon the
   non-local terms may add more ζ(5) (see §4).
2. **Parallelograms (PROVED).** Parallelograms are fully local (α=β=1, no non-local normals), so
   local = total: for every parallelogram, with u+v=1, uv and r=h²/L² per edge,
   $$P^{(5)}(1)\big|_{\zeta(5)} = \sum_e 360\,\frac{h^2}{L}\bigl(1 - 4uv - r\bigr),\qquad
     P^{(5)}(1)\big|_{\zeta(2)\zeta(3)} = \sum_e 120\,\frac{h^2}{L}\bigl(2 - 12uv - 3r\bigr).$$
3. **Regular 2m-gons (PROVED, symbolic in h=π/(2m)).** Σ_e Z_e^(5) = C_m^(5) and Σ_e Z_e^(2,3) =
   C_m^(2,3) of D17, i.e. 2880 m sin⁴h cos h(65 sin⁶h−63 sin⁴h+17 sin²h−1) and
   −960 m cos h sin⁴h(45 sin⁶h−27 sin⁴h+6 sin²h−1), difference 0 after c²=1−s² reduction; m=2..8 numeric
   ≤ 6.8e-49; diamond m=2: −90√2 ζ(5) − 210√2 ζ(2)ζ(3) exactly (`steps/step_05.py`).
4. **Sheared parallelogram ±(2,0),±(1,2) (PROVED + NUMERICAL 55 digits).** Exact:
   P^{(5)}(1) = ζ5·(−6912√5/25 − 172800√13/28561) + ζ2ζ3·(−10752√5/25 − 821760√13/28561)
   + ζ2²·(119808√5/125 + 1262592√13/28561) + ζ3·(4608√5/5 + 23040√13/169) + ζ2·(−3072√5/5 − 15360√13/169);
   ζ5 and ζ2ζ3 coordinates equal Σ_e Z_e^(5), Σ_e Z_e^(2,3) with difference 0 (`steps/step_06.py`).
   Independent numeric (per-t Cauchy d-derivative + tanh-sinh, no use of the symbolic series):
   |num − exact| = 7.8e-57 for P^{(5)}, ≤ 2.7e-59 for P^{(4)} (`steps/step_07.py`, dps 100).

## 1. Method (extends E1 §1–§2 to weight 5)

- `steps/step_01.py`: E1's s01_edge_sym.py at N=6 (d⁰..d⁵), reduced (τ=γ=0). E5 has 475 terms,
  log-degree 5; cross coefficients g_ab = [l0^a l1^b]E5 (a+b=5) are polynomials in t of degree 5.
- `steps/step_02.py`: full N=6 run with lumped τ_k, γ_k (k=1..5). E5: 1910 terms, term counts by
  log-degree {1:8, 2:37, 3:70, 4:74, 5:56}, 35 distinct (τ,γ) patterns. The four cross coefficients
  g41, g32, g23, g14 are τ/γ-free and equal the reduced ones (structure: any τ/γ factor consumes one
  of the 5 d-slots, so it cannot multiply a degree-5 l0/l1 monomial).
- `steps/step_03.py`: Beta table J(n,k,a,b) = ∫₀¹ tⁿ(1−t)^k log^a t log^b(1−t) dt, a+b≤5, exact in
  {1, ζ2, ζ3, ζ2², ζ2ζ3, ζ5} (ζ4 = 2ζ2²/5); Euler γ cancels identically; numeric max error 4.8e-35
  over 210 entries. Rules proved by the log-Beta expansion and verified on the whole table:
  [ζ5]J = 24 B(n+1,k+1) iff a+b=5, a,b≥1, else 0; [ζ2ζ3]J = −12 B(n+1,k+1) iff (a,b)∈{(3,2),(2,3)},
  else 0. Hence for polynomial c(t): [ζ5]∫c(t)X^aY^b = 24∫c, [ζ2ζ3]∫c(t)(X³Y² or X²Y³) = −12∫c.
- Consequence: Z_e^(5) = 5!·24·∫(g41+g32+g23+g14) = 2880∫G, Z_e^(2,3) = −5!·12·∫(g32+g23) =
  −1440∫(g32+g23); both are plain polynomial integrals; shifts l0=log t+μ0, l1=log(1−t)+μ1 and all
  lower-weight terms contribute only weight ≤4 constants times log-monomials — never bare ζ5 or ζ2ζ3.

## 2. The edge law and its forms

Raw integrals (`steps/step_04_out.txt`): ∫G and ∫(g32+g23) are αβh²/L⁴ × (polynomials in L,α,β,s0,h²).
Symmetrised with u=−s0/L, v=1+s0/L into the boxed formulas of §0.1. Checks there:
law−raw difference 0 (exact); reversal invariance exact; regular-polygon per-edge values
1440 s⁴c(65s⁶−63s⁴+17s²−1) and −480 c s⁴(45s⁶−27s⁴+6s²−1) exactly.

Parallelogram specialization (α=β=1): w=1, e=uv, q=1−2uv, Q=2, p=1, s=2 gives
A5 = 1−4uv, B5 = −1, A23 = 2−12uv, B23 = −3 (boxed in §0.2).

## 3. Checks

(a) regular 2m-gon: §0.3. (b) diamond: §0.3. (c) non-square parallelograms: §0.4 and
`steps/step_06_out.txt` (also reproduces E1's P‴(1) = 6S(ζ3−π²/3) exactly, S = Σ h²/L).
Independent numeric for the sheared parallelogram: §0.4 (per-t Cauchy values E_k(t) match the
symbolic E_k to ≤ 2e-49 on a grid including t=0.01..0.99; branch tracking + adaptive radius; tanh-sinh
tail analysis in `steps/step_07.py`).

## 4. The non-local caveat

For a general polygon, E5's terms involving τ_k, γ_k and the nonconstant parts of the non-local logs
l_i(t) = μ_i + log|1+(r_i−1)t| integrate to convergent series of Beta moments of weight ≤5; their
weight-5 constants can add further ζ(5)/ζ(2)ζ(3) beyond Σ_e Z_e — exactly E1's ζ(3) situation.
Whenever every edge is parallel to the kink lines it does not touch (parallelograms, affine-regular
hexagons), the series stop and §0.2/E1-style closed forms are the whole P^{(5)}(1).

## Files

steps/step_01.py (E5 reduced, cross coefficients), step_02.py (full E5 with τ,γ), step_03.py (Beta
table to weight 5 + rules), step_04.py (raw edge-law integrals), step_05.py (closed forms + regular
2m-gon + diamond proofs), step_06.py (exact sheared parallelogram), step_07.py (independent 55-digit
numeric), step_08.py (rhombus 1,2 exact + PSLQ), step_09.py (full-E5 numeric validation on a
non-fully-local edge). progress.log has the run history.
'''
open('result.md','w').write(result)
open('progress.log','a').write("step_06: exact sheared parallelogram P^(5); z5/z2z3 coords == law sums (diff 0); E1 P3 law reproduced. step_07: independent numeric |diff| 7.8e-57 (P5). result.md written.\n")
print("result.md written, progress.log appended")
