# H1: the octahedron end of the superellipsoid, |x|^p + |y|^p + |z|^p = 1 as p → 1

Folder: `research/10_frontier/H1_octahedron_end/`. Every number below was printed by a script in this
folder (Python = `.venv/bin/python`, mpmath 1.3.0, sympy 1.14, python-flint 0.9). Notation:

* d = p − 1, ε = 1/p, η = 1 − 1/p = d/(1+d). Note that η and d agree to first order.
* S(p) = surface area. S(1) = 4√3 is the regular octahedron.
* S_n(p) = surface ((n−1)-volume) of {Σ|x_i|^p = 1} in ℝⁿ. S_n(1) = 2ⁿ√n/(n−1)! is the cross-polytope.
* ζ(2) = π²/6. "Weight" counts ζ(j) as j.

## 0. Answers in one place

| # | statement | status |
|---|---|---|
| T1 | exact representation S_n(p) = 2ⁿ ε^{n−1} Γ(ε)ⁿ/Γ(1+(n−1)ε) · E[(Σ_i G_i^{2η})^{1/2}], with G_i iid Gamma(ε,1) | **PROVED** (this is equivalent to the Naor–Romik surface/cone density formula, §8) |
| T2 | S_n(1+d) has a complete asymptotic expansion in powers of d (both signs of d, no log d terms); its Taylor coefficients come from a finite Gamma-function formula | **PROVED** |
| T3 | structure: the d^k coefficient of S_n/S_n(1) is a polynomial with rational coefficients in ζ(2), ζ(3), …, ζ(k), of weight ≤ k. There are no MZVs, no Cl₂(π/3), no L(χ₋₃, s), no log 2 or log 3, no Catalan constant, no dihedral angle, and no √2 or √3 beyond the prefactor S_n(1) | **PROVED** |
| T4 | **octahedron:** S(1+d) = 4√3 [1 + d + (1 − π²/9)d² + (1 + π²/9 − (16/9)ζ(3))d³ + (1 − 2π²/9 + (32/9)ζ(3) − 16π⁴/405)d⁴ + O(d⁵)]. Exact terms are computed through d¹² | **PROVED**; checked against direct quadrature to 33 digits and by an independent 70-digit PSLQ |
| T5 | **log form (n = 3):** log[(2−p)S(p)/(4√3)] = −(π²/9)η² − (16/9)ζ(3)η³ − (37π⁴/810)η⁴ − (160/27)ζ(5)η⁵ + (16ζ(3)²/81 − 572π⁶/76545)η⁶ + … | **PROVED** (exact through η¹²) |
| T6 | **all n:** log[S_n/S_n(1)] = Σ_k σ_k(n)d^k/k − ((n−1)(n²−2n−1)/(2n))ζ(2)η² − ((n²−1)((n−1)³−2)/(3n²))ζ(3)η³ + W₄(n)η⁴ + O(η⁵), with σ_k and W₄ given in §3 | **PROVED** (derivation in §3; exact checks for n = 2, …, 12) |
| T7 | n = 2 reproduces the diamond law P/√2 = 4 + (π²/6)d² + (ζ(3) − π²/3)d³ + … exactly, through d¹⁰ (D11's a₀…a₁₀) | **PROVED** (exact comparison) |
| T8 | the series diverges. The nearest Borel singularities are at t = (−log(n−1) ± iπ)/2, and they come from the **edge** (ridge) configurations, where one barycentric coordinate goes to 0. For n = 2 this gives ±iπ/2 (known). For n = 3 it gives (−log 2 ± iπ)/2, with \|t\| = √(π² + log²2)/2 = 1.608575… | **NUMERICAL** (Padé–Borel, n = 2…5, 80 terms) + heuristic derivation; zero radius is a **CONJECTURE** for n ≥ 3 (proved for n = 2 in R2) |

**The questions in the brief:**

1. **Does ζ(3) appear, at which order, and with what coefficient?** Yes. It first appears at **d³**, in
   S₃ = 4√3(1 + π²/9 − 16ζ(3)/9). The coefficient of ζ(3) in S is **−64√3/9** (−16/9 relative to S(1)).
   * S₁ = 4√3 = 6.928203230275509174109785366023489467771221015241522512…
   * S₂ = 4√3(1 − π²/9) = −0.669421780076565987998202567942442659217994592036021046824…
   * S₃ = 4√3(1 + π²/9 − 16ζ(3)/9) = −0.279673127275128620615069083260231597445440926451507814613…
2. **Do edges and vertices contribute different constants?** No, not at any finite order (T3). At the octahedron end
   there is no boundary layer. The deformation is a global bending of each flat face, and every coefficient is a
   log-moment of the face, which produces only single zeta values. The constants that would signal edges or vertices,
   such as Cl₂(π/3), L(χ₋₃, 2), log 3, π/√3 or arccos(1/3), never occur (T3). A 68-digit PSLQ that allows log 2, log 3,
   √3, π√3, √3·L(χ₋₃,2), Catalan, π log 2 and π²√3 returns the pure {1, π², ζ(3)} relations. Two finer statements
   hold: (a) is proved (§5); (b) is numerical (§6).
   * (a) In 3D ζ(3) arrives through three channels. One of them is the skewness term E[m₁m₂m₃], a genuinely
     three-coordinate quantity. It vanishes in 2D, and in 3D it is a signed balance of edge regions (negative)
     against vertex regions (positive).
   * (b) The edges control the **large-order** behaviour: the dominant Borel singularity comes from configurations
     with one coordinate → 0, while the vertex direction (two → 0) gives a farther singularity, log 2 ± iπ.
   Note on wording: on the octahedron, an edge is where *one* coordinate vanishes (x₃ = 0, i.e. t₃ = 0) and a vertex
   is where *two* vanish. "Three vanish" never happens on the surface.
3. **How does this compare with the 2D law?** ζ(2) appears at d², and ζ(3) first appears at d³, in both 2D and 3D.
   The differences:
   * (i) In 3D S′(1) = S(1) ≠ 0. In 2D the diamond is a critical point, P′(1) = 0.
   * (ii) The ζ(2) and ζ(3) terms change sign: P/P(1) = 1 + (π²/24)d² + (ζ(3)/4 − π²/12)d³ + … against the octahedron's −π²/9 and −16ζ(3)/9.
   * (iii) The new skewness channel above.
   * (iv) The rational parts of all 3D coefficients equal 1, because of the exact factor 1/(2−p).
4. **The n-dimensional cross-polytope:** T6. The ζ(3) coefficient −(n²−1)((n−1)³−2)/(3n²) never vanishes for an
   integer n ≥ 2. It is positive only for n = 2.

## 1. Exact representation (T1, PROVED)

**Octant formula.** Take the octant x, y, z ≥ 0 and put t₁ = x^p, t₂ = y^p, t₃ = z^p = 1 − t₁ − t₂ on the simplex Δ.
Then dA = |∇F|/|F_z| dx dy, with F = x^p + y^p + z^p, and x = t₁^ε, dx = ε t₁^{ε−1} dt₁. This gives

  S(p) = 8 ε² ∫_Δ (t₁t₂t₃)^{−η} (t₁^{2η} + t₂^{2η} + t₃^{2η})^{1/2} dt₁ dt₂.  (1)

In ℝⁿ: S_n(p) = 2ⁿ ε^{n−1} ∫_{Δ_{n−1}} (Π t_i)^{−η} (Σ t_i^{2η})^{1/2} dσ. At p = 1 this gives 2ⁿ√n/(n−1)!.

**Gamma variables.** The integrand of (1) equals (Π t_i)^{ε−1} times h(t) = (Σ t_i^{2η})^{1/2}, and h is
homogeneous of degree η. For G_i iid Gamma(ε) we have t = G/ΣG ~ Dirichlet(ε, …, ε), and ΣG is independent of t
(Lukacs). Hence

  ∫_Δ (Πt)^{ε−1} h dσ = Γ(ε)ⁿ/Γ(nε+η) · E[h(G)],  with nε + η = 1 + (n−1)ε.

That is T1. Checks: at p = 2, G_i = Z_i²/2 gives the sphere 2π^{n/2}/Γ(n/2) exactly; at p = 1 it gives 2ⁿ√n/(n−1)!.

**Centred form (used in Route B).** Write m_i = log t_i − mean_j log t_j; these log-ratios are invariant under scaling t ↦ λt.
Then (Πt)^{−η}(Σt^{2η})^{1/2} = √n (Πt)^{−(n−1)η/n} f(η, m), where

  f(η, m) = ( (1/n) Σ_i e^{2η m_i} )^{1/2} = exp( ½ K(2η) ),

and K is the cumulant generating function of the n numbers m_i. The same Dirichlet → Gamma step with weight (Πt)^{a},
a = −(n−1)η/n, gives

  S_n/S_n(1) = (n−1)! ε^{n−1} Γ(1+a)ⁿ/Γ(1+(n−1)ε) · B_n(η),  B_n(η) = E[f(η, m)],  m_i = log G_i − mean log G,  G_i iid Gamma(1+a).  (2)

The first factor is the exact "flat-face" (Dirichlet-normaliser) part. B_n is the "bending factor".

## 2. Expansion and structure (T2, T3, PROVED)

**Route A (uncentred).**
* Write (ΣX_i)^{1/2} = √n (1+x)^{1/2}, with x = (1/n)Σ(G_i^{2η} − 1) ≥ −1.
* For every K, |(1+x)^{1/2} − Σ_{k≤K} C(½,k)x^k| ≤ C_K|x|^{K+1} on x ≥ −1.
* |G^{2η} − 1| ≤ 2|η||log G| max(1, G^{2η}), and E[|log G|^m G^c] < ∞ for c > −ε. So E|x|^{K+1} = O(η^{K+1})
  uniformly near ε = 1, on both sides.
* Each E[x^k] is a finite combination of the analytic functions Γ(ε + 2jη)/Γ(ε).

Hence

  S_n/S_n(1) = (n−1)! ε^{n−1}/Γ(n − (n−1)η) · Σ_{k≤K} C(½,k) n^{−k} T_k + O(η^{K+1}),
  T_k = Σ_{k₁+…+k_n=k} (k!/Πk_i!) Π Φ_{k_i},  Φ_k = Σ_{j=0}^{k} C(k,j)(−1)^{k−j} Γ(1 + (2j−1)η).  (3)

Φ_k = O(η^k) is a k-th finite difference. The remainder bound gives T2: an asymptotic power series in η, hence in d,
with no logarithms. The n = 2 case of (3) is the D11/R2 Gamma-operator formula for the diamond.

**Route B (centred, independent code).** In (2), expand f in η as a polynomial in the m_i. Use the cumulants of log G:
κ₁ irrelevant (shift invariance), and κ_j = ψ^{(j−1)}(1+a) = Σ_r ψ^{(j−1+r)}(1)a^r/r! with ψ^{(m)}(1) = (−1)^{m+1}m!ζ(m+1).

**Proof of T3.**
* The η-expansion of f is Σ_k η^k (homogeneous polynomial of degree k in m).
* The expectation of a degree-k homogeneous polynomial in centred independent variables is a sum over set partitions
  (blocks of size ≥ 2) of products of cumulants κ_b. The total block size is k.
* κ_b(1+a) contributes ζ(b+r)η^r, of weight b + r. So [η^k] log B_n is a homogeneous weight-k polynomial in single zeta values.
* The flat-face factor equals Rat_n(η) · exp Λ_n(η), with:
  * Rat_n = (1−η)^{n−1}/Π_{i=1}^{n−1}(1 − (n−1)η/i), which is rational;
  * Λ_n = −Σ_{j≥2} ζ(j)(n−1)^j(1 − n^{1−j})η^j/j, which has weight j at η^j.
* Euler's γ cancels identically: n·logΓ(1+a) and logΓ(1−(n−1)η) have the same linear term.
* Substituting η = d/(1+d) mixes orders downward. So the d^k coefficient has weight ≤ k. ∎

**Agreement of the routes.**
* Route A (`01_series_routeA.py`) and Route B (`02_series_routeB.py`) agree exactly, coefficient by coefficient: through
  order 10 for n = 3 and order 6 for n = 2 (`06_crosscheck_routes.py`, all differences 0).
* The numerical high-order version of Route A (`05_high_order_numeric.py`, 80 orders, 150 digits) matches the exact
  values to the 1e-42 to 1e-47 comparison precision.

## 3. General n (T6, PROVED)

Write
* A₂(n) = −(n−1)(n²−2n−1)/(2n),
* A₃(n) = −(n²−1)((n−1)³−2)/(3n²),
* σ_k(n) = Σ_{i=1}^{n−2}((n−1−i)/i)^k, so that log Rat_n = Σ_k σ_k d^k/k, i.e. Rat_n = Π_{i=1}^{n−2} i/(i − (n−1−i)d).

Then:

  log[S_n(p)/S_n(1)] = Σ_k σ_k d^k/k + A₂ζ(2)η² + A₃ζ(3)η³ + W₄η⁴ + O(η⁵),
  W₄ = −((n−1)/n²)ζ(2)² + ((n−1)(n²−3n+4)/n³)ζ(4) − ((n−1)⁴(n³−1)/(4n³))ζ(4).

*Derivation.* Λ_n is exact, and log B_n = η²E[μ₂] + (2/3)η³E[μ₃] + η⁴(E[κ̂₄]/3 + ½Var μ₂) + O(η⁵), where:
* μ_r are the empirical moments of the m_i, and κ̂₄ = μ₄ − 3μ₂²;
* E[μ₂] = ((n−1)/n)ψ′(1+a), with ψ′(1+a) = ζ(2) − 2ζ(3)a + 3ζ(4)a²;
* E[μ₃] = ((n−1)(n−2)/n²)ψ″(1+a);
* at a = 0: Var μ₂ = 2ζ(2)²(n−1)/n² + 6ζ(4)(n−1)²/n³, and E[κ̂₄] = 6ζ(4)(n−1)(n²−6n+6)/n³ − 6ζ(2)²(n−1)/n².

Collecting terms gives log B_n = ((n−1)/n)ζ(2)η² + (2(n²−1)/(3n²))ζ(3)η³ + [−((n−1)/n²)ζ(2)² + ((n−1)(n²−3n+4)/n³)ζ(4)]η⁴ + O(η⁵).
`08_logform_general_n.py` and `13_final_numbers.py` confirm each of these coefficients exactly for n = 2, …, 12
(exact Route A data, every difference 0).

Explicit d-coefficients (c_k := [d^k] S_n/S_n(1)):
* **c₁ = σ₁ = (n−1)(H_{n−1} − 1)**: 0, 1, 5/2, 13/3, 77/12, … for n = 2, 3, 4, 5, 6.
* c₂ = (σ₁² + σ₂)/2 + A₂ζ(2).
* c₃ = σ₁³/6 + σ₁(σ₂/2 + A₂ζ(2)) + σ₃/3 − 2A₂ζ(2) + A₃ζ(3).

| n | c₂ | c₃ |
|---|---|---|
| 2 | π²/24 | ζ(3)/4 − π²/12 |
| 3 | 1 − π²/9 | 1 + π²/9 − 16ζ(3)/9 |
| 4 | 21/4 − 7π²/16 | 85/8 − 7π²/32 − 125ζ(3)/16 |
| 5 | 130/9 − 14π²/15 | 1210/27 − 98π²/45 − 496ζ(3)/25 |

(n = 6, …, 12 are in `final_numbers.out`.)

For n = 2 the formula gives P/P(1) = 1 + (π²/24)d² + (ζ(3)/4 − π²/12)d³ + …, i.e. P/√2 = 4 + (π²/6)d² + (ζ(3) − π²/3)d³.
That is the verified diamond law. `07_check_n2_vs_D11.py` checks a₀…a₁₀ against D11's exact list: all differences 0.

## 4. The octahedron, n = 3 (T4, T5)

* Rat₃ = (1−η)/(1−2η) = 1/(1−d) = 1/(2−p). So (2−p)S(p)/(4√3) = exp(Σ L_k η^k), with these exact coefficients (`15_logform_n3.py`):
  * L₂ = −π²/9 and L₃ = −16ζ(3)/9;
  * L₄ = −37π⁴/810 and L₅ = −160ζ(5)/27;
  * L₆ = 16ζ(3)²/81 − 572π⁶/76545;
  * L₇ = −(136/3645)π⁴ζ(3) − (64/243)π²ζ(5) − (512/27)ζ(7).
* Up to η⁵ every L_k is a single ζ(k). For n = 3 the ζ(2)ζ(3) term at η⁵ cancels exactly: in log B_n its coefficient is
  4(n−1)(n−3)/n³ for n ≤ 12 (observation, from the printed table).
* d-coefficients of S/(4√3) (`13_final_numbers.py`):
  * c₄ = 1 − 2π²/9 + 32ζ(3)/9 − 16π⁴/405;
  * c₅ = 1 + 2π²/9 + 16π⁴/135 − 64ζ(3)/9 + 16π²ζ(3)/81 − 160ζ(5)/27 = 2.388781136862079097857575551381978370430703788…;
  * c₆ = −1.15634388622839768456700575000381699543…, c₇ = −18.17446862343036166675596234699…, c₈ = 53.73447568952023481505698111…
* Numerically, c₀…c₄ = 1, 1, −0.0966227112321509576482767777640167928126332674711989584903722,
  −0.0403673388293499941734799538118943018806756966896879424736662,
  −0.767525708869791359265497680078141717215540880109430897925472.

## 5. Where ζ(3) comes from (the honest edge/vertex picture)

In log S at η³ (n = 3), the coefficient −16/9 = −48/27 splits into three channels:

| channel | n = 3 | n = 2 | origin |
|---|---|---|---|
| flat-face normaliser Λ_n | −64/27 | −1/4 | Γ(1+a)ⁿ/Γ(1+(n−1)ε) |
| variance shift E[μ₂] = ((n−1)/n)ψ′(1+a) | +24/27 | +1/2 | ψ′(1+a) = ζ(2) − 2ζ(3)a |
| skewness (2/3)E[μ₃] | −8/27 | 0 | ψ″(1) = −2ζ(3), weight (n−1)(n−2)/n² |

* For n = 3, μ₃ = (1/3)Σm_i³ = m₁m₂m₃ (because Σm_i = 0). This is the one term that needs three coordinates.
* Near an edge (t₃ → 0) m₃ → −∞ and m₁, m₂ → +∞, so m₁m₂m₃ < 0. Near a vertex (t₂, t₃ → 0) it is > 0.
* Its mean is (2/9)ψ″(1) = −(4/9)ζ(3), so the edges win the balance.
* The ζ(2) at d² comes from Var(log G) = ψ′(1) (a cross-covariance of log t_i, log t_j on the face) and from Λ_n.
* None of these can be attributed to a region of the octahedron in a canonical way. Below any finite order there is
  no small length scale: the d^k term is ∫_Δ (polynomial in log t_i) dσ, which converges absolutely.
* This is the opposite of the cube end p → ∞ (07_push/S3). There, faces, edges and vertices separate into boundary
  layers and the vertex brings log 3, arcsin(1/√3) and Li₂(−(√2−1)⁴).

## 6. Large order: the edges dominate (T8, NUMERICAL / CONJECTURE)

**Heuristic.**
* The only factorial growth is in B_n; Rat_n and Λ_n have finite radius.
* a_k = E[c_k(m)], where c_k is the degree-k Taylor coefficient of f(η, m) = F(η·m).
* Along a ray l = L·v (L → ∞), the density of the log-Gamma variables decays like e^{−ρ(v)L}. Here ρ is the number of
  coordinates sent to −∞, since the right tail is super-exponentially thin.
* So the Borel transform Σ a_k t^k/k! is singular at t = ρ(v)·(ηL at the first complex zero of Σ_i e^{2ηm_i}).
* Edge ray (one coordinate → −∞, ρ = 1): (n−1)e^{2s/n} + e^{−2(n−1)s/n} = 0 with s = ηL, i.e. e^{2s} = −1/(n−1). That gives
  **t_edge = (−log(n−1) ± iπ)/2**.
* For n = 3 the family l = (0, −βL, −L), β ∈ [0, 1], traces a curve of singularities τ(β) = (1+β)·½log x(β), with
  x + x^{1−β} + 1 = 0. It runs from the edge point τ(0) = −0.34657 + 1.57080i to the vertex point τ(1) = log 2 + iπ
  = 0.69315 + 3.14159i. |τ| rises from 1.6086 at β = 0 to a maximum of 3.249 near β = 0.65, then falls slightly to
  3.2172 at β = 1 (`12_singular_curve.py`). Its minimum is at the edge endpoint, so that endpoint is the nearest
  singularity.

**Numerics** (`10_borel_growth.py`, Borel transform of 80 coefficients, [40/40] Padé; output `borel_growth.out`):

| n | predicted t_edge | nearest Padé–Borel pole ([30/30] → [38/38] → [40/40]) |
|---|---|---|
| 2 | ±1.570796i | ±1.578842i → ±1.575978i → ±1.575555i (Re ≈ 2·10⁻⁴) |
| 3 | −0.346574 ± 1.570796i | −0.34811 ± 1.57971i → −0.34732 ± 1.57672i → −0.34737 ± 1.57606i |
| 4 | −0.549306 ± 1.570796i | −0.55324 ± 1.58046i → −0.55158 ± 1.57680i → −0.55111 ± 1.57610i |
| 5 | −0.693147 ± 1.570796i | −0.69797 ± 1.57995i → −0.69647 ± 1.57648i → −0.69599 ± 1.57594i |

* The next poles line up along the predicted cut. For n = 3: −0.3498 + 1.5920i, −0.3538 + 1.6194i, moving up and to
  the left like τ(β).
* The d-series gives the same points.
* With t_edge fixed, the amplitude fit (`11_growth_exponent.py`) gives a_k/k! ≈ Re[Z k^{−3/2} t_edge^{−k}]:
  * n = 2: the exponent goes to −1.495 at k = 78 and |Z| goes to 1.408 ≈ √2, matching D11's √2·n!(2/π)ⁿn^{−3/2};
  * n ≥ 3: |Z_k|k^{3/2} is constant to ±3% (n = 3, η-series: 1.33–1.37), but the local exponent oscillates
    between about −2.0 and −0.9 (for n = 3, between −1.8 and −1.2), probably because of the nearby cut. The exponent
    for n ≥ 3 is therefore not determined.
* Consequences:
  * The series in d has radius 0 (CONJECTURE for n ≥ 3; for n = 2 it is R2's theorem).
  * The Borel radius is |t_edge| = √(π² + log²(n−1))/2, which is 1.608575 for the octahedron.
  * There is no singularity on ℝ, which suggests Borel summability for both p > 1 and p < 1 (CONJECTURE).

## 7. Evidence table

| quantity | method | agreement | script / output |
|---|---|---|---|
| S(p) at d ∈ {0, ±0.1, ±0.03, ±0.01, ±0.003, ±0.001} | N1 (graph over the "z largest" sixth of an octant, nested tanh–sinh) vs N2 (simplex form (1)), 40 digits | ≤ 1.0e−39 | `03_direct_quadrature.py`, `direct_dps40.out` |
| S(2) = 4π, S(1) = 4√3 | N1 | 40 digits | same |
| series (exact, K = 12) vs N1 | 4√3 Σ_{k≤12} c_k d^k | 1.1e−33 (d = ±0.001), 1.6–1.9e−27 (±0.003), 0.8–1.3e−20 (±0.01); 1.4–2.0e−7 at ±0.1 (asymptotic) | `04_compare_series_vs_direct.py`, `compare_n3_K12_dps40.out` |
| consistency of c₁…c₄ with the data | symmetric peeled remainders at d = ±10⁻³ (lower c_k exact) | avg − c_{j+1} = −4.0e−8, −7.7e−7, +2.4e−6, −1.2e−6 ≈ c_{j+3}·10⁻⁶, as expected | same |
| Route A ≡ Route B | exact sympy | identical through order 10 (n = 3) and 6 (n = 2) | `06_crosscheck_routes.py` |
| n = 2 vs D11 diamond a₀…a₁₀ | exact | all differences 0 | `07_check_n2_vs_D11.py` |
| c₁…c₄ with no Gamma algebra | pointwise η-expansion of (1), 2-D tanh–sinh at 70 digits, then PSLQ | c₂: 9c₂ − 9 + π² = 0; c₃: −9c₃ + 9 + π² − 16ζ(3) = 0; c₄: max coefficient 1440; at 50 and 60 digits, numeric − exact = 3.6e−71, 2.5e−70, 9.0e−69 | `09_pslq_independent.py`, `pslq_independent.out` |
| absence of edge/vertex constants | PSLQ on {c, 1, π², ζ(3), log 2, log 3, √3, π√3, √3·L(χ₋₃,2), G, π log 2, π²√3} at 68 digits | same relations; every extra constant has coefficient 0 | same |
| general-n η², η³, η⁴ (log B_n), c₁ = (n−1)(H_{n−1}−1), c₂, c₃ | exact | all differences 0, n = 2…12 | `08_logform_general_n.py`, `13_final_numbers.py` |
| n = 3, 4, 5 series vs a non-series evaluation of T1 (2-D nested quadrature of the Gaussian-linearised form) | 32-digit quadrature | n = 3: 1.6–7.8e−15, n = 4: 1.5–7.0e−15, n = 5: 0.9–4.5e−15 at d = ±0.01, ±0.001. This is the quadrature's own accuracy: for n = 3 it equals the gap to the 40-digit direct value | `14_ndim_check.py`, `ndim_check.out` |
| Borel singularities | Padé–Borel, 80 terms | table in §6 | `10_borel_growth.py`, `11_growth_exponent.py`, `12_singular_curve.py` |

## 8. Prior literature (checked 2026-10-21, web search; full texts not read)

* T1 is classical in substance. Naor–Romik, *Projecting the surface measure of the sphere of ℓ_pⁿ* (Ann. IHP 2003),
  give dσ/dμ = n Vol(B_pⁿ)|∇‖x‖_p| relative to the cone measure. Combined with the Schechtman–Zinn / Rachev–Rüschendorf
  Gamma representation of the cone measure, this is exactly T1. I checked the equivalence by hand:
  n Vol(B_pⁿ) = 2ⁿε^{n−1}Γ(ε)ⁿ/Γ(nε), and E[(ΣG)^η] = Γ(1+(n−1)ε)/Γ(nε) (Lukacs) converts one form into the other.
* 2D: P′(1) = 0 and P″(1) = (√2/3)π² are known (Montilla–Rodríguez–Linares 2607.04592; Rodríguez–Montilla 2607.09048,
  whose abstract covers only 2D perimeters); see the LEDGER.
* I found no expansion of the superellipsoid or ℓ_p-sphere surface area about p = 1 in any dimension ≥ 3. That covers
  S′(1) = S(1), the ζ(2) and ζ(3) terms, T3, T6, and the edge-controlled Borel singularity.
* Searches: MathWorld "Supersphere" (only the multiple-integral formula), and general ℓ_p surface-area and affine
  surface-area literature (nothing at p → 1).
* Status: "to our knowledge new", with the usual residual risk that the full texts were not read.

## 9. Scripts and outputs (all in this folder)

| script | what it does |
|---|---|
| `01_series_routeA.py` | exact Taylor coefficients via (3) (args: K, n…) → `routeA_n*_K*.json` |
| `02_series_routeB.py` | exact coefficients via the cumulant/centred form (2) → `routeB_n*_K*.json` |
| `03_direct_quadrature.py` | direct surface area, methods N1 and N2 → `direct_S_dps40.json` |
| `04_compare_series_vs_direct.py` | series vs direct, peeled remainders → `compare_n3_K12_dps40.out` |
| `05_high_order_numeric.py` | 80-order numerical coefficients → `hi_coeffs_n{2,3,4,5}_K80.json` |
| `06_crosscheck_routes.py` | Route A ≡ Route B ≡ numeric |
| `07_check_n2_vs_D11.py` | n = 2 against D11's exact diamond coefficients |
| `08_logform_general_n.py` | general-n η², η³ theorem; log B_n table → `logform_general_n.out`, `logB_coeffs.json` |
| `09_pslq_independent.py` | Gamma-free 70-digit log-moment integrals + PSLQ → `pslq_independent.out` |
| `10_borel_growth.py`, `11_growth_exponent.py`, `12_singular_curve.py` | large-order analysis → `borel_growth.out`, `growth_exponent.out`, `singular_curve.out` |
| `13_final_numbers.py` | W₄ check, general-n c₁…c₃, 60-digit S_k → `final_numbers.out` |
| `14_ndim_check.py` | n = 3, 4, 5 non-series check → `ndim_check.out` |
| `15_logform_n3.py` | exact log form for n = 3 → `logform_n3.out` |

## 10. Open questions

1. Prove zero radius and Borel summability for n ≥ 3. The R2/D11 argument for n = 2 should extend: the Borel
   transform is an integral over the ray family, singular at the edge endpoint. Also get the exact large-order constant
   (amplitude and exponent) at the endpoint of the singular curve τ(β).
2. Close log B_n in general n at η⁵ and beyond. The ζ(2)ζ(3) coefficient 4(n−1)(n−3)/n³ is so far only an
   observation, made for n ≤ 12.
3. Explain the single-zeta pattern L_k ∝ ζ(k) (k ≤ 5) of the n = 3 log form. Is there a hidden Γ-product formula
   up to that order?
4. Find a geometric reason for c₁ = (n−1)(H_{n−1} − 1). This is the first variation of the cross-polytope's area in p:
   it vanishes only for n = 2.
5. Do general semi-axes |x/a|^p + |y/b|^p + |z/c|^p = 1 work as in R2? Expect a weighted Dirichlet with
   ψ-values at shifted arguments, still giving only single zetas.
