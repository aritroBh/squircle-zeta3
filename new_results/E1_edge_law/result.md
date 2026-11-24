# E1: the general edge law for ζ(3) at the diamond end (p → 1)

Folder: `research/10_frontier/E1_edge_law/`. Read first: LEDGER (D15, D15b, D15c), `D15c_edge_zeta3/result.md`
and `sym_G.py`, `D15b_polygon_law/poly_coeffs.py`, `R2_rhombus_ab/result.md`.
Every number below is printed by a script in this folder (log file named in brackets).
Labels: **PROVED** = exact symbolic computation or written proof; **NUMERICAL** = high-precision agreement
(residual printed); **CONJECTURE** = pattern without proof.

## 0. Bottom line

1. **The edge law (PROVED, every centrally symmetric convex polygon, any number of sides, any normals and
   weights).** For K_p = {Σ_j w_j |⟨n_j,x⟩|^p ≤ 1},
   $$P'''(1)\big|_{\zeta(3)}=\sum_{\text{edges } e} Z_e,\qquad
   Z_e = 6\,\alpha_e\beta_e\,\frac{h_e^2}{L_e}\Bigl(3(\alpha_e u_e+\beta_e v_e)-2\Bigr).$$
   L_e is the edge length, h_e the distance from the centre to the edge line, u_e and v_e = 1 − u_e the
   fractions into which the foot of the perpendicular from the centre cuts the edge, and α_e, β_e ∈ (0,1]
   are the "sharpness" of the two end corners (definitions in §1). Equivalent forms in terms of the two
   exterior angles, or the two base angles and the two local weights, are in §2.3. It reduces exactly to
   C_m (regular 2m-gon) and to 24a²b²/r³ (rhombus) (`s14_log.txt`).
2. **Where ζ(3) can enter (PROVED, all m).** Parametrise each arc of ∂K_p by the flat edge of K_1 it
   projects to (radial projection; §1). Then the d³ coefficient of the perimeter element is a polynomial
   in t times logarithms of affine functions of t. It contains no square roots and no log of the distance
   to the centre. So
   P'''(1) = ζ(3)·Σ_e Z_e + A + ζ(2)·B,
   where A and B are explicit, geometrically convergent series. Each term of these series is a polynomial
   with rational coefficients in the edge data and in logarithms of geometric lengths. ζ(3) occurs in
   none of them (§3). The single-corner log³ and log² terms give only rationals and ζ(2). The terms
   with logs that do not vanish on the edge give only rationals and ζ(2), term by term.
3. **Closed forms (PROVED).** If every edge is parallel to the kink lines it does not touch, the series
   stop and P'''(1) is elementary plus ζ(3). That happens exactly for parallelograms and for
   affine-regular hexagons. With S = Σ_e h_e²/L_e and the power convention Σ|⟨v_j,x⟩|^p:
   - every parallelogram: P(p) = P(1) + (π²S/6)d² + S(ζ(3) − π²/3)d³ + O(d⁴). This generalises R2.
   - every affine-regular hexagon: P(p) = 2^{1−1/p}[P(1) + (π²S/24)d² − (S/8)(ζ(3) + 2π²/3)d³] + O(d⁴).
   - D15's hexagon formula P‴(1) = (√3/8)(16log³2 + 3π²log2 − 96log²2 + 96log2 − 9ζ(3) − 6π²) is now a
     **theorem** (exact difference 0, `s04_log.txt`). Before this it was a 170-digit PSLQ result.
4. **Tests (NUMERICAL).** These cover the two requested hexagons (three weights 1, 2, 3; three unequally
   spaced normals), an affine-regular hexagon, octagons, and a degeneration. They are checked against
   independent totals: D15b's polar coefficient integrals agree to ≤ 4.5e-58 at 60 digits, and a direct
   perimeter with finite differences agrees to ≤ 1.5e-16 (hexagons A and B) and ≤ 3e-15 (H_p). For the affine-regular hexagon, PSLQ on the
   independent polar value isolates the ζ(3) coefficient −3/2 − 3√2/16 at 80 and at 120 digits. That is
   exactly the edge-law value (§4).
5. **What "ζ(3)-free" means (important caveat; NUMERICAL examples, §5).** For a general polygon the
   remainder A + ζ(2)B has no ζ(3) *term by term* in the canonical expansion, which is centred on the
   configuration where the other kink lines are parallel to the edge. But a series of rational numbers can
   still add up to a multiple of ζ(3) plus elementary terms. Example (`s12_log.txt`): the series for
   ∫₀¹log²t·log(1+t)dt, with every term rational apart from log 2, equals 4log2 − 6 + π²/6 + (3/2)ζ(3).
   Concretely, when a hexagon is squeezed into a rhombus, P‴(1) is continuous but the edge-law
   coefficient jumps by 264√5/125. The jump is carried by the hexagon's "ζ(3)-free" remainder.
   So the edge law is the canonical ζ(3) coefficient on each combinatorial type of polygon. It is an
   unambiguous number-theoretic coefficient whenever the remainder is elementary: parallelograms,
   affine-regular hexagons, and in particular the regular hexagon. For regular 2m-gons with m ≥ 4,
   whether C_m equals the ζ(3) coefficient in a chosen basis of cyclotomic constants depends on that basis.
6. **Corrections to earlier records (§7).** The ledger's "|C_m| → 9π⁵/(16m⁴)" is wrong: C_m = −3π⁴/m³ +
   25π⁶/(8m⁵) + O(m⁻⁷). D15c's "single-corner log³ carries ζ(3)" is an artefact of the polar-angle
   bookkeeping.

---

## 1. Set-up and the radial-projection formula (PROVED)

**The polygon.** N(x) = Σ_{j=1}^m w_j|⟨n_j,x⟩|, with pairwise non-parallel unit normals n_j and w_j > 0.
K_1 = {N ≤ 1} is a centrally symmetric convex 2m-gon. Its vertices lie on the m kink lines ⟨n_j,x⟩ = 0,
two vertices on each line. Every centrally symmetric convex 2m-gon arises this way, uniquely: it is
the polar of a zonotope. The p-family is K_p = {Σ_j w_j e^{λ_j d}|⟨n_j,x⟩|^p ≤ 1} with d = p − 1:
- λ = 0 is the task's convention (fixed weights);
- λ_j = log w_j is the "power" convention Σ_j|⟨v_j,x⟩|^p with v_j = w_j n_j (R2, D15c's rhombus).

In code, kink j is given by its direction b_j and n_j = (−sin b_j, cos b_j), as in poly_coeffs.py.

**Edge data.** Take an edge e from V0 to V1 (counter-clockwise), with V0 on kink line j0 and V1 on kink
line j1.
- L = |V1 − V0|; τ = (V1 − V0)/L; h = distance from 0 to the edge line.
- σ(t) = ⟨y(t),τ⟩ = s0 + Lt along y(t) = (1−t)V0 + tV1. The foot F of the perpendicular is at σ = 0, so
  u = −s0/L and v = 1 − u.
- a_j(t) = w_j|⟨n_j,y(t)⟩|. These are affine in t, they satisfy Σ_j a_j(t) = N(y) = 1, and a_{j0}(0) =
  a_{j1}(1) = 0.
- **α = a_{j0}(1) = w_{j0}|⟨n_{j0},V1⟩|** is the share of corner V0's term in the norm at the far vertex.
  Likewise **β = a_{j1}(0)**. So 0 < α, β ≤ 1, with equality only when m = 2.
- ℓ_j(t) = log|⟨n_j,y(t)⟩| + λ_j.

**Lemma 1 (exact arclength).** x_p(t) = y(t)·F(t)^{−1/p} with F = Σ_j a_j e^{dℓ_j} is the radial
projection of y(t) onto ∂K_p. It maps the edge onto the arc of ∂K_p between the rays through V0 and V1,
and
$$|x_p'(t)| = F^{-1/p}\sqrt{(L-\sigma\Phi)^2+h^2\Phi^2},\qquad \Phi=\frac{\sum_j a_j'\,e^{d\ell_j}}{\sum_j a_j\,e^{d\ell_j}}.$$
*Proof.* F_p(y) = Σ w_j e^{λ_j d}|⟨n_j,y⟩|^p is homogeneous of degree p, so y/F_p^{1/p} lies on
∂K_p. Next, dF/dt = pΣ_j a_j' e^{dℓ_j}, because a_jℓ_j' = a_j'. So x_p' = F^{−1/p}(y′ − Φy), with
|y′|² = L², ⟨y,y′⟩ = σL and |y|² = h² + σ². ∎

The perimeter is P(p) = Σ_e ∫₀¹|x_p′(t)|dt for every p > 0. The form (L − σΦ)² + h²Φ² is positive
definite in (1, Φ): its determinant is L²h² > 0. Every d-derivative of the integrand is therefore
bounded by C·t^{−δ}(1−t)^{−δ}(1 + |log t| + |log(1−t)|)^k for |d| ≤ δ < 1. Hence P is C^∞ at p = 1 and
P^{(k)}(1) = k!·Σ_e ∫₀¹E_k(t)dt, where E_k is the d^k coefficient of |x_p′|. (Smooth, not analytic: R2
proved that the Taylor series has radius 0.)

**Lemma 2 (structure; PROVED, `s01_edge_sym.py`).** E_k(t) is a polynomial of total degree ≤ k in the
logarithms ℓ_j(t). Its coefficients are polynomials in t, whose own coefficients are polynomials in α, β,
L, 1/L, h, s0 and the affine coefficients of the other a_j. At p = 1, E_0 = L.
Contrast with the polar angle θ used in D15–D15c. There log|sin(θ − b_j)| = ℓ_j − λ_j − ν(t), with
ν = log|y(t)| = ½log(h² + σ²). This log of the distance to the centre has complex branch points
t_± = (−s0 ± ih)/L (the isotropic directions). The reparametrisation θ(t) does not depend on p, so
D_k(θ)dθ = E_k(t)dt pointwise.
(NUMERICAL: ≤ 4.6e-40 at 40 digits on hexagons A and B and the regular decagon; `s09_log.txt`.
Integrated edge by edge at 30 digits: ≤ 8.2e-27 for k = 3; `s02_log.txt`.)

## 2. The cross coefficients and the edge law (PROVED)

### 2.1 g₂₁, g₁₂
Write ℓ_{j0} = log t + μ0 and ℓ_{j1} = log(1−t) + μ1. The cubic part of a polynomial does not change
when its variables are shifted by constants. So the coefficient of log²t·log(1−t) in E_3 equals the
coefficient g₂₁ of ℓ_{j0}²ℓ_{j1}, and similarly for g₁₂.

The top-degree part of E_3 is the degree-3 Taylor part, in X_j = dℓ_j, of
Ψ = F(X)^{−2}·√((L·F(X) − σN(X))² + h²N(X)²), where F = Σ a_j e^{X_j} and N = Σ a_j′e^{X_j}. (The
explicit 1/p and the factor F^{d/p} only lower the log-degree.) Its X_{j0}²X_{j1} coefficient is taken
at X_j = 0 for every other j. There only Σ_{j∉{j0,j1}} a_j = 1 − αt − β(1−t) and
Σ_{j∉{j0,j1}} a_j′ = β − α enter. **So g₂₁ and g₁₂ do not depend on the other kink terms at all.**
(This is the structural reason for D15c's observation.) Exactly, with σ = s0 + Lt (`s01_log.txt`):

g₂₁ = αβ[24L³αt³ − 18L³αt² − 6L³t² + 4L³t + 18L²αs0t² − 12L²αs0t − 4L²s0t + 2L²s0 + 6Lαh²t − 3Lαh² − Lh² − 3αh²s0]/(2L²),

g₁₂ = −αβ[24L³βt³ − 36L³βt² + 12L³βt + 6L³t² − 4L³t + 18L²βs0t² − 24L²βs0t + 6L²βs0 + 4L²s0t − 2L²s0 + 6Lβh²t − 6Lβh² + Lh² − 3βh²s0]/(2L²).

These are cubic polynomials in t. Integrating them:
$$\int_0^1 g_{21}=-\frac{\alpha\beta h^2(L+3\alpha s_0)}{2L^2},\qquad \int_0^1 g_{12}=\frac{\alpha\beta h^2(3L\beta-L+3\beta s_0)}{2L^2},\qquad
\int_0^1(g_{21}+g_{12})=\frac{\alpha\beta h^2}{2L^2}\bigl(3\beta\sigma_1-3\alpha\sigma_0-2L\bigr),$$
with σ0 = s0 and σ1 = s0 + L. The Beta lemma (§3) gives the ζ(3) part 2ζ(3)∫(g₂₁+g₁₂) per edge, and
P‴ = 6Σ∫E_3, so Z_e = 12∫(g₂₁+g₁₂). This is the boxed formula of §0.

**Link to D15c (NUMERICAL, `s09_log.txt`, `s10_log.txt`).** I made D15c's symbolic G₂₁, G₁₂
(`sym_G.py`) weighted by s_i → w_is_i and c_i → w_ic_i. Then G₂₁(θ(t))·dθ/dt = g₂₁(t) and
G₁₂·dθ/dt = g₁₂ pointwise, to ≤ 6.7e-41. And 12∫_arc(G₂₁ + G₁₂)dθ = Z_e to ≤ 1.9e-40 on every edge of
hexagons A and B, the regular decagon, a tangential hexagon and a tangential octagon. This is forced by
uniqueness: the coefficient of log²t·log(1−t) beside functions analytic on [0,1] is unique.

### 2.2 Independence of the p-dependence of the weights
Z_e depends only on K_1, not on λ. Any smooth positive w_j(p) with w_j(1) = w_j shifts the ℓ_j by
constants plus O(d), and that leaves the cubic part unchanged. Examples: hexagon A has Z identical in
the fixed and power conventions, while P‴ differs (−1.46553 vs −0.40717; `s06_log_60.txt`). The same
holds for H_p (`s04_log.txt`).

### 2.3 Equivalent geometric forms (identities PROVED; checked NUMERICALLY ≤ 1.9e-40, `s10_log.txt`, `s16_log.txt`)
- **Pure geometry of K_1.** α_e = L_e sin ε_0/(2h_{e⁻}) and β_e = L_e sin ε_1/(2h_{e⁺}). Here ε_0, ε_1 are
  the exterior angles at V0, V1, and h_{e∓} are the distances from the centre to the neighbouring edge
  lines. This follows from w_{j0}n_{j0} = ±(∇N|_e − ∇N|_{e⁻})/2. So
  Z_e = (3/2)·L_e h_e² sin ε_0 sin ε_1/(h_{e⁻}h_{e⁺})·(3(α_eu_e + β_ev_e) − 2).
  **Z_e is proportional to sin ε_0·sin ε_1, so it vanishes as soon as either end stops being a corner.**
  This is the "interference of the two corners" made quantitative.
- **Length, two angles, local weights.** Let ψ_0, ψ_1 be the angles of the triangle O V0 V1 at V0, V1.
  Then α = w_{j0}L sin ψ_0 and β = w_{j1}L sin ψ_1, and
  Z_e = 6w_{j0}w_{j1}L³ sin³ψ_0 sin³ψ_1/sin²(ψ_0+ψ_1)·(3L sinψ_0 sinψ_1(w_{j0}cosψ_0 + w_{j1}cosψ_1)/sin(ψ_0+ψ_1) − 2).
- **Tangential polygons (all h_e equal; the centre is the incentre).**
  Z_e = (3/2)L_e sin ε_0 sin ε_1(3sin²(ε_0/2) + 3sin²(ε_1/2) − 2) (exact, `s14_log.txt`).

### 2.4 Reductions (exact, `s14_log.txt`)
| case | edge data | result |
|---|---|---|
| regular 2m-gon, D15 normalisation | α=β=R sin(π/m), R=tan(π/2m), u=v=½ | 2m·Z_e − C_m = 0, C_m = 48m sin⁴(π/2m)cos(π/2m)(3sin²(π/2m)−1) |
| rhombus \|x/a\|+\|y/b\| ≤ 1 | α=β=1, L=r, h=ab/r, u=a²/r² | 4Z_e − 24a²b²/r³ = 0 |
| any parallelogram | α=β=1 | Z_e = 6h²/L |
| any affine-regular hexagon | α=β=½ (affine invariant) | Z_e = −(3/4)h²/L |

The sign of Z_e is the sign of 3(αu + βv) − 2. Hexagon A has one pair of edges with Z_e > 0 and two
pairs with Z_e < 0 (`s17_log.txt`).

## 3. The non-cross cancellation, for all m (PROVED)

**Theorem.** On each edge put N_i(t) = log(1 + (r_i − 1)t) for every kink i ∉ {j0, j1}, with
r_i = ⟨n_i,V1⟩/⟨n_i,V0⟩ > 0. Then
E_3 = Σ c_{n,a,b,J}·t^n·log^a t·log^b(1−t)·Π_i N_i^{J_i}, with a + b + |J| ≤ 3.
Each coefficient c lies in R_e = ℚ[α, β, L, 1/L, h, s0, (affine data of the other a_i), μ0, μ1, μ_i],
where μ_i = log|⟨n_i,V0⟩| + λ_i. Moreover
$$6\int_0^1E_3\,dt = Z_e\,\zeta(3) + A_e + \zeta(2)B_e .$$
A_e and B_e are absolutely convergent series. Each term is an element of R_e[log r_i] times powers of
(r_i − 1) or (1/r_i − 1), and ζ(3) occurs in no term. The series terminate when every r_i = 1.

*Proof.*
(i) **Beta lemma** (`s03_beta_table.py`; verified for all n ≤ 4, k ≤ 2, a+b ≤ 3, to 4e-39). We have
∫₀¹t^n(1−t)^k log^a t log^b(1−t)dt = ∂_x^a∂_y^b B(x,y) at (n+1, k+1), and
log B(x+u, y+v) − log B(x,y) = Σ_j [ψ^{(j−1)}(x)u^j + ψ^{(j−1)}(y)v^j − ψ^{(j−1)}(x+y)(u+v)^j]/j!.
At integers, ψ = −γ + H, ψ′ = ζ(2) − H^{(2)} and ψ″ = −2ζ(3) + 2H^{(3)}. Euler's γ cancels
identically, because u + v − (u+v) = 0. The ζ(3) part of the cubic terms is
−(2ζ(3)/6)[u³ + v³ − (u+v)³] = ζ(3)(u²v + uv²). The pure u³ and v³ parts telescope: ψ″(x) − ψ″(x+y)
is rational. After exponentiation, ζ(3) therefore appears, through total degree 3, only in the u²v
and uv² coefficients. It enters ∂_x^a∂_y^b B only for (a,b) ∈ {(2,1),(1,2)}, with coefficient
2B(n+1,k+1). For a + b ≤ 2 the value is in ℚ + ℚζ(2).
For (a,b) = (3,0) or (0,3) it is rational.
(ii) **Local terms (J = 0).** The coefficients are polynomials in t, so the sums are finite. The ζ(3)
part is 2∫(c_{2,1} + c_{1,2}) = 2∫(g₂₁ + g₁₂) = Z_e/6.
(iii) **Non-local terms (|J| ≥ 1).** Here a + b ≤ 2. Since r > 0, we have min(|r−1|, |1/r−1|) < 1.
Expand N_i either as Σ_q (−1)^{q+1}((r_i−1)t)^q/q or as
log r_i + Σ_q (−1)^{q+1}((1/r_i−1)(1−t))^q/q. Both converge uniformly on [0,1]. Termwise integration
is justified by domination, and every term is a Beta moment with a + b ≤ 2, which lies in ℚ + ℚζ(2). ∎

**Why the polar-angle bookkeeping showed ζ(3) in the non-cross terms (D15c §3).** In θ-coordinates every
log carries −ν(t) = −½log(h² + σ²). For example, a single-corner (log sin)³ term contains
−3·log²t·ν, a "cross term" between the corner and the isotropic points t_±. Such integrals give
polylogarithms at complex arguments and contain ζ(3); D15c printed ∫₀^{π/2}log³sin θ dθ ∋ −(3π/4)ζ(3).
All of these cancel identically, because D_3dθ = E_3dt contains no ν at all. In words: **the distance
to the centre drops out of the perimeter element once the arc is parametrised by its flat edge.** The
resummation caveat (∫log²t/(1+t) = (3/2)ζ(3)) is answered as follows. In the radial parametrisation the
local coefficients are polynomials, so the Beta rule is exact. The only non-polynomial factors are logs
of affine functions whose zeros lie on the real line outside [0,1], and the theorem expands them about
the nearer endpoint.

## 4. Numerical tests

**4.1 Totals by three independent routes (`s05_log.txt`, `s06_log_30.txt`, `s06_log_60.txt`).**
- Hexagon A: kink directions 0, −π/3, −2π/3 (regular normals), weights 1, 2, 3, fixed-weight convention.
  Exact data (`s17_log.txt`): (α, β) = (1/3, 3/5), (3/4, 2/3), (2/5, 1/4) and r = 5/3, 3/4, 4/5.
  **Z = 90√39/2197 − 90√57/6859 − 26√7/343 = −0.0437912799017308550242690047576.**
- Hexagon B: kink directions 0, 1, 2.2 (unequally spaced), weights 1, fixed-weight convention.
  **Z = −1.99721421038240647936943175483.**

| polygon | P‴(1), decomposed Zζ(3)+Q0+Q2ζ(2)+NL | vs D15b polar, 30 / 60 digits | vs direct perimeter + FD (h=0.003, K=8) |
|---|---|---|---|
| hexagon A | −1.4655259166759924661521072656959599745 | 7.6e-27 / 3.2e-58 | 1.5e-16 |
| hexagon B | −5.285825448806008251345266145529503989022 | 1.4e-28 / 3.9e-58 | 8.5e-17 |
| hexagon A, power convention | −0.4071685118319488838268519800913010133889 | 6.9e-27 / 4.5e-58 | — |

The remainders are P‴ − Zζ(3) = Q0 + Q2ζ(2) + NL:
- Hexagon A: Q0 = −0.32125207948487103447, Q2 = −0.64428749704302872237, NL = −0.03182377415659795681.
- Hexagon B: Q0 = 5.61501303528409662877, Q2 = −5.17270404346017770122, NL = 0.00868374339972020572.

The NL series is built only from Beta moments with a + b ≤ 2, and it agrees with direct quadrature of the
non-local terms to ≤ 7e-31 per edge at 30 digits. The ζ(3) coefficient of the exact local part equals Z
to the working precision: ≤ 1.3e-60 at 60 digits. So the remainder is ζ(3)-free in the canonical sense,
and this was checked at two precisions against an independent total.

**4.2 Affine-regular hexagon H_p = {|x|^p + |y|^p + |x−y|^p ≤ 1}: ζ(3) isolated by PSLQ on independent
data (`s07_log.txt`).**
P‴(1) was computed with D15b's polar code alone, which does not use this note's derivation. PSLQ was run
on {P‴, ζ(3), √2ζ(3), π², √2π², π²log2, √2π²log2, log³2, √2log³2, log²2, √2log²2, log2, √2log2, 1, √2}.
It found the same integer relation at 80 and at 120 digits (residuals 3.4e-80 and 1.6e-119):
16P‴ = −(24+3√2)ζ(3) − (16+2√2)π² + (8+√2)π²log2 + (32+16√2)log³2 − (192+96√2)log²2 + (192+96√2)log2.
So the ζ(3) coefficient is −3/2 − 3√2/16, exactly the edge law, and the remainder P‴ − Zζ(3) has zero
ζ(3) in this basis. The fixed-weight version, with the √2-weighted normal, gives the same ζ(3)
coefficient (64P‴ relation, residual 8e-119). This matches the exact closed forms in `s04_log.txt`
(difference 0).

**4.3 General m (`s13_log_40.txt`).** These polygons have two non-local kinks per edge, so products
N_iN_j occur.
- Regular octagon: Z = C_4 = −2.132921749138399510974383; P‴(1) decomposed = −4.39449630527605215997911011381.
  This agrees with D15b's 200-digit value to 1.1e-25 at 25 digits and with the polar integral to
  1.1e-37 at 40 digits.
- Generic octagon (kinks 0, 0.7, 1.6, 2.5; weights 1, 1.5, 0.8, 1.2): Z = −1.809866442920175809070811;
  decomposed vs polar agree to 1.6e-38.

**4.4 Family laws vs polar integrals (`s15_log.txt`, 45 digits).** A sheared parallelogram, a second
parallelogram and T(regular hexagon) with T = [[1, 0.4], [0.2, 1.3]] all agree with the laws of §6 to
≤ 2.5e-41.

## 5. What "ζ(3)-free" can and cannot mean (NUMERICAL examples)

- **A series of rationals can sum to ζ(3).** `s12_log.txt`: I = ∫₀¹log²t·log(1+t)dt is a non-local term
  of the kind that occurs on an edge with r = 2. Its canonical expansion about t = 1,
  log 2·J₀ − Σ_k J_k/(k2^k), where every J_k = ∫log²t(1−t)^k dt is rational, gives
  0.2206081438273991022409508947459064476684. That equals 4log2 − 6 + π²/6 + (3/2)ζ(3) to all 40
  digits shown.
- **Degeneration hexagon → rhombus (`s11_log.txt`, `s11b_log.txt`).** Take kinks 0, π/2, π/2+ε with
  weights 1. As ε → 0, K_p tends to the rhombus family {|y|^p + 2|x|^p ≤ 1}.
  - P‴(1) is continuous: P‴ − P‴_rhombus = 7.5e-2, 2.0e-3, 4.0e-5, 9.3e-8 at ε = 1e-5, 1e-7, 1e-9,
    1e-12. The ratio to ε·log³ε is −4.9, −4.7, −4.5, −4.4.
  - The edge law jumps: Z_hex → −24√5/125 (error 2.46ε), whereas Z_rhombus = 48√5/25. On the two long
    edges, β changes from ½ to 1 when the two kink terms merge.
  - The jump is 264√5/125 = 4.72257556847955583882. The hexagon's canonical remainder therefore tends
    to R_rhombus + (264√5/125)ζ(3): the hidden ζ(3) emerges as r → 0, where the expansion variable
    r − 1 reaches the edge of its disc of convergence.

So the edge law is the canonical ζ(3) coefficient on a fixed combinatorial type of polygon, expanded
about the "parallel" configuration. It is a basis-free number-theoretic coefficient exactly when the
remainder is elementary (§6). For regular 2m-gons with m ≥ 4, each edge is parallel to at most one other
kink line, so non-local logs are unavoidable. C_m is then the canonical value, and its agreement with a
PSLQ coefficient depends on the chosen basis of cyclotomic constants.

## 6. Closed forms for fully local polygons (PROVED; `s04_log.txt`, `s08_log.txt`)

Fully local means r_i = 1 on every edge: each edge is parallel to all kink lines it does not touch. For
m = 2 this holds always. For m = 3 it holds exactly for affine-regular hexagons (Q − P ∥ R and its
cyclic versions force Q = P + R). For m ≥ 4 it is impossible. In the power convention, put
S = Σ_{all edges} h_e²/L_e and d = p − 1. Per edge (exact):
- parallelogram: ∫E₂ = ζ(2)h²/L and ∫E₃ = (h²/L)(ζ(3) − 2ζ(2));
- affine-regular hexagon: ∫E₃ = L(λ³/6 − λ² + λ) + (h²/L)(ζ(2)(λ/4 − ½) − ζ(3)/8), with λ = log 2.

Hence:
- **Parallelogram law.** P′(1) = 0, P″(1) = (π²/3)S, P‴(1) = 6(ζ(3) − π²/3)S, with
  S = (A²/2)(L₁⁻³ + L₂⁻³) for area A and side lengths L₁, L₂. For the rhombus this is exactly R2:
  4π²a²b²/(3r³) and (24a²b²/r³)(ζ(3) − π²/3), difference 0.
- **Affine-regular hexagon law.** P(p) = 2^{1−1/p}·Q(p) + O(d⁴), with Q(1) = P(1), Q′(1) = 0,
  Q″(1) = (π²/12)S and Q‴(1) = −(3/4)S(ζ(3) + 2π²/3).
  - The regular hexagon (S = 3√3/2) gives D15's P′, P″, P‴ exactly.
  - H_p = {|x|^p+|y|^p+|x−y|^p ≤ 1} (S = 2 + √2/4, P(1) = 2 + √2) gives
    P‴(1) = (2+√2)(log³2 − 6log²2 + 6log2) − ((8+√2)/16)(3ζ(3) + 2π² − π²log2).
- **Fixed-weight convention.** There are extra log-of-weight terms but the same ζ(3). For example, the
  rhombus {|y|^p/2 + |x|^p ≤ 1} has P‴(1) = (96√5/25)ζ(3) + ζ(2)(−192√5/25 + 192√5log2/125) −
  96√5log2/5 − 576√5log²2/25 − 448√5log³2/125.

## 7. Corrections to earlier records

1. LEDGER/D15c: "|C_m| → 9π⁵/(16m⁴)" is wrong. Exact series: C_m = −3π⁴/m³ + 25π⁶/(8m⁵) + O(m⁻⁷)
   (`s14_log.txt`: m³C_m = −263.0729, −291.9269, −292.2243 at m = 10, 100, 1000, against
   −3π⁴ = −292.2273).
2. D15c's "cross = total, numerical to ≤ 8e-95 for m = 2..8" compared 12m∫D₃ with P‴(1). That checks the
   integral, not the ζ(3) content. The cancellation is now PROVED (§3) in the canonical sense for all m,
   and as an absolute statement for m = 2, 3 regular and all fully local polygons.
3. D15c's rule "[ζ(3)]∫g·ℓ₀²ℓ₁ = 2ζ(3)∫g for analytic g" is exact for polynomial g, which is what the
   radial parametrisation produces. For general analytic g, "the ζ(3) part" is not well defined (§5).
4. D15/D15c: "the single-corner log³ term carries ζ(3)" (∫log³sin θ) is a polar-angle artefact (§3). In
   the edge parametrisation, single-corner terms give rationals only for log³ and ℚ + ℚζ(2) for log².
5. D15's hexagon P″(1) and P‴(1), and R2's rhombus P″ and P‴, are re-derived exactly here (difference 0).
   D15's hexagon values were PSLQ results before; they are now PROVED.

## 8. For the paper (plain words, Aritro's voice)

Where does ζ(3) come from when we round off a polygon? Let's look at one side at a time. Each side runs
from one corner to the next, and as p moves away from 1 the side bends a little near both of its
corners. If we follow the bent side using the straight side it came from, everything the perimeter
needs turns out to be a polynomial times logarithms: the log of the distance to the left corner, the
log of the distance to the right corner, and logs of distances to the polygon's other fold lines, which
never vanish on this side. A single corner's logarithm, even cubed, only ever gives rational numbers
and π². ζ(3) shows up only in the mixed terms, the square of one corner's log times the other corner's
log, and those integrate in closed form. A side of length L at distance h from the centre contributes
Z = 6αβ(h²/L)(3(αu + βv) − 2), where α and β measure how sharp its two corners are and u, v say where
the perpendicular from the centre lands on the side. Z is proportional to the product of the sines of the
two corners' turning angles, so ζ(3) really needs both corners: flatten either one and that side's ζ(3)
disappears. Adding one Z per side gives all the ζ(3) in P‴(1); for the regular 2m-gon this is exactly
C_m, and for the rhombus it is 24a²b²/r³. For every parallelogram and every affinely regular hexagon
the other logarithms are constants, and we get exact formulas, for example P‴(1) = 6(ζ(3) − π²/3)Σh²/L
for any parallelogram. For other polygons, everything else is a convergent series of rational numbers,
logarithms and π², which contains no ζ(3) term by term. We found, though, that such a series can still
add up to a multiple of ζ(3) (this happens when a hexagon is squeezed into a rhombus), so for those
polygons the edge law tells us exactly where ζ(3) enters, not that it can never hide elsewhere.

## 9. Open
- Closed form of the non-local part in polylogarithms at the letters r_i. That would give a natural basis
  in which to compare C_m (m ≥ 4) with a PSLQ coefficient.
- The same mechanism at higher order: ζ(5) at d⁵ should come only from the quintic cross terms
  (ψ⁗ telescopes the same way), and ζ(3)·log terms appear at d⁴. This is a CONJECTURE with an obvious
  proof route: the Lemma 2 and Beta-lemma argument.
- Smooth (non-polygonal) convex bodies: is there a continuum limit of Σ_e Z_e?

## Files
- `s01_edge_sym.py`: exact E_0..E_3 (`s01_E.pkl`), g₂₁, g₁₂, their integrals (`s01_out.txt`).
- `edge_num.py`: numeric edge data, E_k(t), the edge law, direct polar perimeter and finite differences.
- `s02_validate.py`: E_k vs D15b's polar D_k per edge; Z vs C_m and the rhombus.
- `s03_beta_table.py`: exact Beta table and the ζ(3) rule.
- `s04_exact_local.py`: exact closed forms (`s04_closed_forms.txt`).
- `s05_fd.py`: direct perimeter and finite differences.
- `s06_generic_hex.py`: decomposition of hexagons A and B (`s06_out_30.json`, `s06_out_60.json`).
- `s07_pslq_H.py`: PSLQ on independent polar data (`s07_pslq_H.json`).
- `s08_families.py`: parallelogram and affine-regular hexagon laws.
- `s09_theta_check.py`: pointwise link to D15c's G₂₁, G₁₂.
- `s10_geometric_forms.py`, `s16_angle_form.py`: equivalent forms.
- `s11_degeneration.py`, `s11b_degeneration_small.py`: the jump 264√5/125.
- `s12_basis_example.py`: rational series = ζ(3).
- `s13_general_m.py`: octagons.
- `s14_reductions.py`: exact reductions and the C_m asymptotics.
- `s15_family_check.py`: family laws vs polar integrals.
- `s17_exact_Z_hexA.py`: exact Z for hexagon A.
- `run_all.sh`: reruns everything into `*_log.txt`.
