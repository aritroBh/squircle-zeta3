import Mathlib
import Zeta3.Wave1.T09_OctSeries
import Zeta3.Wave4.X05_DiamondA4A5
/-!
# Y05 — the octahedron end of the superellipsoid (paper label `thm:oct`)

Target id: Y05.
Paper labels: `thm:oct` (main.tex lines 1975–1990, setup line 1973, remarks lines 1992–1993),
the superellipsoid `|x|^p+|y|^p+|z|^p=1` of Section `sec:surfaces` (main.tex line 1909),
`thm:pw` for `N = 3` (main.tex lines 1183–1195, the factor `Y_3 = 1/(2-p)`), and the detailed
proofs: paper/ieee/proofs/sec_g7_octa.tex (Lemmas `g7:lem:both`, `g7:lem:X2`, `g7:lem:w23`,
`g7:lem:dcoef`, the machine step `g7:sec:machine`) and paper/ieee/proofs/sec_pureweight.tex,
Section `sec:pw:n` (equation `eq:pw:Sn`, Lemma `lem:pw:dirichlet`, the centring and
equation `eq:pw:factor`, the flat-face factor, Lemmas `lem:pw:X`, `lem:pw:moments`).

Quoted verbatim.

main.tex 1909:

  The three-dimensional squircle is the superellipsoid $|x|^p+|y|^p+|z|^p=1$
  (Fig.~\ref{fig:superellipsoid}). Its surface area $S(p)$ goes from $4\sqrt3$ (the octahedron,
  $p=1$) through $4\pi$ (the sphere) to $24$ (the cube).

main.tex 1973 (setup):

  Let's write $d=p-1$, as at the diamond, and $\tau=1-\frac1p=\frac d{1+d}$, and $[d^k]S$ for
  the coefficient of $d^k$ in $S(1+d)$.

main.tex 1975–1990 (`thm:oct`):

  As $d\to0$ from either side,
  \frac{S(1+d)}{4\sqrt3}=1+d+\Bigl(1-\frac{\pi^2}{9}\Bigr)d^2+\Bigl(1+\frac{\pi^2}{9}
  -\frac{16\zeta(3)}{9}\Bigr)d^3
  +\Bigl(1-\frac{2\pi^2}{9}+\frac{32\zeta(3)}{9}-\frac{16\pi^4}{405}\Bigr)d^4+O(d^5),
  and for every $k$, $[d^k]S$ is $4\sqrt3$ times a polynomial with rational coefficients in
  $\zeta(2),\dots,\zeta(k)$ of weight at most $k$. Moreover
  \log\frac{(2-p)\,S(p)}{4\sqrt3}=-\frac{\pi^2}{9}\tau^2-\frac{16\zeta(3)}{9}\tau^3
  -\frac{37\pi^4}{810}\tau^4-\frac{160\zeta(5)}{27}\tau^5
  +\Bigl(\frac{16}{81}\zeta(3)^2-\frac{572}{76545}\pi^6\Bigr)\tau^6+O(\tau^7),
  where the coefficient of $\tau^k$ is a polynomial in $\zeta(2),\dots,\zeta(k)$, homogeneous of
  weight $k$.

main.tex 1990 (text after the theorem, the face integral):

  The area of the first-octant face is, exactly, $(1-\tau)^2=p^{-2}$ times the integral over the
  triangle $t_1,t_2\ge0$, $t_1+t_2\le1$ of
  $(t_1t_2t_3)^{-\tau}(t_1^{2\tau}+t_2^{2\tau}+t_3^{2\tau})^{1/2}$, with $t_i=x_i^p$,
  $t_3=1-t_1-t_2$ and $\tau=1-\frac1p$.

sec_pureweight.tex, `eq:pw:Sn` (with `n = 3`):

  S_n(p)=2^n\varepsilon^{n-1}\int_\Delta\Bigl(\prod_{i=1}^nt_i\Bigr)^{-\tau}
  \Bigl(\sum_{i=1}^nt_i^{2\tau}\Bigr)^{1/2}dt_1\cdots dt_{n-1},
  $\Delta=\{t_i\ge0,\ \sum_{i<n}t_i\le1\}$. The integral converges since $\tau<1$.

sec_pureweight.tex, graph-area derivation of `eq:pw:Sn`:

  That part is the graph $x_n=(1-\sum_{i<n}x_i^p)^{1/p}$ over
  $\{x_i\ge0,\ \sum_{i<n}x_i^p\le1\}$. ... the area element is
  $\sqrt{1+\sum_{i<n}\bigl(\partial_ix_n\bigr)^2}\,dx$ ...

sec_pureweight.tex, centring and Lemma `lem:pw:dirichlet` (`n = 3`, `a = -2τ/3`):

  S_n(p)=2^n\sqrt n\,\varepsilon^{n-1}\frac{\Gamma(1+a)^n}{\Gamma(n-(n-1)\tau)}\,B_n(\tau),
  \qquad B_n(\tau)=\E\bigl[f(\tau,m(G))\bigr],
  where $G_1,\dots,G_n$ are independent with density $g^a\eu^{-g}/\Gamma(1+a)$ on $(0,\infty)$
  and $m_i(G)=\log G_i-\frac1n\sum_j\log G_j$;
  $f(\tau,m)=\bigl(\frac1n\sum_i\eu^{2\tau m_i}\bigr)^{1/2}$.

sec_pureweight.tex, flat-face factor:

  $\Gamma(1+a)^n/\Gamma(1-(n-1)\tau)=\exp\bigl(\sum_{k\ge2}\lambda_k\zeta(k)\tau^k\bigr)$ with
  $\lambda_k=-\frac{(n-1)^k(1-n^{1-k})}k\in\QQ$.

sec_g7_octa.tex, proof of Lemma `g7:lem:w23`:

  \E[Y]=\frac{\Gamma\bigl(1+a+2\tau(1-\frac1n)\bigr)\,\Gamma\bigl(1+a-\frac{2\tau}n\bigr)^{n-1}}
  {\Gamma(1+a)^n}.

sec_g7_octa.tex, Lemma `lem:pw:moments` / machine step (`|ν| = s`, `c_v = ν_v - s/3`):

  \E\bigl[\eu^{2\tau\sum_{l\in S}m_{i_l}}\bigr]=\prod_{v=1}^n
  \frac{\Gamma\bigl(1+a+2\tau c_v\bigr)}{\Gamma(1+a)}

sec_g7_octa.tex, Lemma `g7:lem:X2`:

  \Bigl|B_n(\tau)-\sum_{j\le J}\binom{1/2}j\E[X^j]\Bigr|\le C_J'\,|\tau|^{2J+2}
  \qquad(|\tau|\le\tau_{J+1}).

sec_g7_octa.tex, Lemma `g7:lem:w23` (`n = 3`): $w_2(3)=-\frac23Z_2$ and
$w_3(3)=-\frac{16}9Z_3$, that is $-\frac{\pi^2}9$ and $-\frac{16\zeta(3)}9$ at $\tau^2$ and
$\tau^3$.

Formalization choices.
* `ζ(n)` is `X05.z n = (riemannZeta n).re` (imported from `Zeta3.Wave4.X05_DiamondA4A5`), the
  real part of Mathlib's complex `riemannZeta`; for `n ≥ 2` it is real and equals
  `∑ 1/k^n` (`zeta_nat_eq_tsum_of_gt_one`). `ζ(2)`, `ζ(4)`, `ζ(6)` are written as powers of `π`,
  as the paper does.
* The truncated series of `thm:oct` through order 4 are `T09.dSeries (X05.z 3) d` and
  `T09.logSeries (X05.z 3) τ` (imported from `Zeta3.Wave1.T09_OctSeries`, which also proved
  that the `d`-series follows from the `τ`-series); they are exactly the displays quoted above,
  truncated after `d^4` resp. `τ^4`.
* DEFINITION OF `S`. `S p` is DEFINED as the paper's exact face integral `eq:pw:Sn` for `n = 3`:
  `S(p) = 8 (1/p)^2 ∫_Δ (t₁t₂t₃)^{-τ} (t₁^{2τ}+t₂^{2τ}+t₃^{2τ})^{1/2} dt₁dt₂`, `t₃ = 1-t₁-t₂`,
  `τ = 1 - 1/p`, `Δ = {t₁,t₂ ≥ 0, t₁+t₂ ≤ 1}` (Lebesgue measure on `ℝ × ℝ`, Bochner integral).
  stmt_1 says the integrand is integrable for every `p > 0`, so the integral is not a junk
  value. stmt_2 says `S` equals `8` times the classical graph-area integral
  `∬ √(1 + (∂₁x₃)² + (∂₂x₃)²)` of the first-octant piece `x₃ = (1-x₁^p-x₂^p)^{1/p}`
  (the paper's derivation of `eq:pw:Sn`). Boundary points of `Δ` (where `0^x`-conventions of
  `rpow` and `deriv`-junk could matter) form a null set.
* `τ = 1 - 1/p`, `a = -2τ/3` (`n = 3`). The Gamma-variable expectation
  `E_τ[F] = Γ(1+a)^{-3} ∫_{(0,∞)^3} (g₁g₂g₃)^a e^{-(g₁+g₂+g₃)} F(g) dg` is `gammaE τ F`.
  `Ymean τ g = (g₁^{2τ}+g₂^{2τ}+g₃^{2τ}) / (3 (g₁g₂g₃)^{2τ/3})` is the paper's
  `Y = (1/3) Σ e^{2τ m_i}` (since `e^{2τ m_i} = g_i^{2τ}/(g₁g₂g₃)^{2τ/3}`), and the bending
  factor is `B τ = E_τ[√Y]`. `X = Y - 1`.
* `O(·)` claims are explicit: `∃ C δ, 0 < δ ∧ ∀ … , |var| ≤ δ → |error| ≤ C |var|^k`, for both
  signs of the variable ("as `d → 0` from either side"). In the `τ`-statements `p > 0` ranges
  over a two-sided neighbourhood of `1` via `|1 - 1/p| ≤ δ`. The `log` is `Real.log`; near
  `p = 1` its argument is positive (it is `≈ 1`), so the statements are not true by the
  `log 0 = 0` convention (the `τ²` term is nonzero).
* `Y_3 = 1/(2-p)`: `(2-p) S(p)/(4√3) = Γ(1+a)^3/Γ(1-2τ) · B(τ)` follows from stmt_3 by
  `Γ(3-2τ) = (2-2τ)(1-2τ)Γ(1-2τ)` and `1-2τ = (2-p)/p`.

Statements.
* stmt_1: integrability of the face integrand on `Δ` for every `p > 0` ("The integral converges
  since τ<1"; both signs of `τ`, sec_g7_octa.tex proof of `g7:lem:both`).
* stmt_2: `S(p)` equals `8 ×` the graph-area integral of the first-octant piece, `p > 0`.
* stmt_3: Dirichlet–Gamma factorization with centring:
  `S(p) = 8√3 (1/p)² Γ(1+a)³/Γ(3-2τ) · B(τ)`, `p > 0`.
* stmt_4: flat-face factor: `log(Γ(1+a)³/Γ(1-2τ)) = Σ_{k≥2} λ_k ζ(k) τ^k`,
  `λ_k = -2^k(1-3^{1-k})/k`, for `|τ| < 1/2` (HasSum).
* stmt_5: `E[Y] = Γ(1+2τ/3) Γ(1-4τ/3)² / Γ(1-2τ/3)³` for `|τ| < 3/4` (proof of `g7:lem:w23`).
* stmt_6: `E[Y²] = (3 Γ(1+2τ) Γ(1-2τ)² + 6 Γ(1-2τ)) / (9 Γ(1-2τ/3)³)` for `|τ| < 1/2`
  (Lemma `lem:pw:moments` with `j = 2`: multiplicity vectors `2e_i` and `e_i+e_j`).
* stmt_7: Lemma `g7:lem:X2` with `J = 2`: `|B - (1 + E[X]/2 - E[X²]/8)| ≤ C τ⁶`.
* stmt_8: Lemma `g7:lem:w23` (`n = 3`), the hand-proved part:
  `log((2-p)S/(4√3)) = -π²/9 τ² - 16ζ(3)/9 τ³ + O(τ⁴)`.
* stmt_9: the hand-proved `d`-coefficients `d⁰,…,d³` (sec_g7_octa.tex "What is proved here"):
  `S(1+d)/(4√3) = 1 + d + (1-π²/9)d² + (1+π²/9-16ζ(3)/9)d³ + O(d⁴)`.
* stmt_10: the target claim: the `d`-series through `d⁴` with `O(d⁵)` AND the `log` series
  through `τ⁴` with `O(τ⁵)`.
* stmt_11: `thm:oct` at full strength for its analytic part: the `d`-series through `d⁴` with
  `O(d⁵)` AND the `log` series through `τ⁶` with `O(τ⁷)` (`ζ(5)`, `ζ(3)²`, `π⁶` terms).

NOT COVERED:
* That `S(p)` as defined is the surface area (2-dimensional Hausdorff measure) of
  `{|x|^p+|y|^p+|z|^p = 1}`: the 8-fold symmetry reduction, the null coordinate planes, and the
  identification of the graph-area integral with Hausdorff measure. Only the graph-area formula
  for the first-octant piece is stated (stmt_2).
* The weight claims: "for every `k`, `[d^k]S` is `4√3` times a polynomial with rational
  coefficients in `ζ(2),…,ζ(k)` of weight at most `k`", and "the coefficient of `τ^k` is a
  polynomial in `ζ(2),…,ζ(k)`, homogeneous of weight `k`" (formal-polynomial statements for all
  `k`; also the all-orders existence of the expansion, Lemma `g7:lem:both` for general `K`).
* The general-`J` form of Lemma `g7:lem:X2`, the general Lemma `lem:pw:moments`, and the
  formal Lemma `g7:lem:Eformal`; only `J = 2`, `j ≤ 2` are stated.
* Other dimensions `N ≠ 3` (cross-polytope, `w_2(n)`, `w_3(n)`, the `ζ(3)`-coefficient
  `-(N²-1)((N-1)³-2)/(3N²)`), the values `S(2) = 4π`, `S(∞) = 24`, the remark
  `[d^1]S = S(1)`, the `τ⁷` and higher coefficients, and the divergence / Borel statements of
  the text after the theorem.
-/
-- BEGIN STATEMENTS
namespace Y05

/-- The triangle `Δ = {t₁, t₂ ≥ 0, t₁ + t₂ ≤ 1}` in `ℝ × ℝ`. -/
def Δ : Set (ℝ × ℝ) := {t | 0 ≤ t.1 ∧ 0 ≤ t.2 ∧ t.1 + t.2 ≤ 1}

/-- The integrand of `eq:pw:Sn` (`n = 3`):
`(t₁t₂t₃)^{-τ} (t₁^{2τ}+t₂^{2τ}+t₃^{2τ})^{1/2}`, with `t₃ = 1 - t₁ - t₂`. -/
noncomputable def faceIntegrand (τ : ℝ) (t : ℝ × ℝ) : ℝ :=
  (t.1 * t.2 * (1 - t.1 - t.2)) ^ (-τ) *
    Real.sqrt (t.1 ^ (2 * τ) + t.2 ^ (2 * τ) + (1 - t.1 - t.2) ^ (2 * τ))

/-- The surface area `S(p)` of `|x|^p+|y|^p+|z|^p = 1`, DEFINED by the paper's exact formula
`eq:pw:Sn` (`n = 3`): `S(p) = 8 (1/p)² ∫_Δ faceIntegrand (1 - 1/p)`. -/
noncomputable def S (p : ℝ) : ℝ :=
  8 * (1 / p) ^ 2 * ∫ t in Δ, faceIntegrand (1 - 1 / p) t

/-- The first-octant piece of the surface as a graph: `x₃ = (1 - x₁^p - x₂^p)^{1/p}`. -/
noncomputable def zgraph (p x₁ x₂ : ℝ) : ℝ := (1 - x₁ ^ p - x₂ ^ p) ^ (1 / p)

/-- Its domain `{x₁, x₂ ≥ 0, x₁^p + x₂^p ≤ 1}`. -/
def D (p : ℝ) : Set (ℝ × ℝ) := {x | 0 ≤ x.1 ∧ 0 ≤ x.2 ∧ x.1 ^ p + x.2 ^ p ≤ 1}

/-- The graph-area element `√(1 + (∂₁x₃)² + (∂₂x₃)²)`. -/
noncomputable def areaElt (p : ℝ) (x : ℝ × ℝ) : ℝ :=
  Real.sqrt (1 + (deriv (fun u => zgraph p u x.2) x.1) ^ 2
    + (deriv (fun v => zgraph p x.1 v) x.2) ^ 2)

/-- `Y = (1/3) Σ_i e^{2τ m_i}` written in the Gamma variables `g = (g₁, g₂, g₃)`:
`(g₁^{2τ}+g₂^{2τ}+g₃^{2τ}) / (3 (g₁g₂g₃)^{2τ/3})`. -/
noncomputable def Ymean (τ : ℝ) (g : ℝ × ℝ × ℝ) : ℝ :=
  (g.1 ^ (2 * τ) + g.2.1 ^ (2 * τ) + g.2.2 ^ (2 * τ)) /
    (3 * (g.1 * g.2.1 * g.2.2) ^ (2 * τ / 3))

/-- The expectation over independent `G₁, G₂, G₃` with density `g^a e^{-g}/Γ(1+a)` on
`(0,∞)`, `a = -2τ/3`. -/
noncomputable def gammaE (τ : ℝ) (F : ℝ × ℝ × ℝ → ℝ) : ℝ :=
  (Real.Gamma (1 - 2 * τ / 3) ^ 3)⁻¹ *
    ∫ g in Set.Ioi (0 : ℝ) ×ˢ Set.Ioi (0 : ℝ) ×ˢ Set.Ioi (0 : ℝ),
      (g.1 * g.2.1 * g.2.2) ^ (-2 * τ / 3) * Real.exp (-(g.1 + g.2.1 + g.2.2)) * F g

/-- The bending factor `B₃(τ) = E[√Y]`. -/
noncomputable def B (τ : ℝ) : ℝ := gammaE τ (fun g => Real.sqrt (Ymean τ g))

/-- sec_pureweight.tex after `eq:pw:Sn` ("The integral converges since τ<1"; both signs of `τ`
in the proof of `g7:lem:both`): for every `p > 0` the face integrand is integrable on `Δ`. -/
def stmt_1 : Prop :=
  ∀ p : ℝ, 0 < p → MeasureTheory.IntegrableOn (faceIntegrand (1 - 1 / p)) Δ

/-- sec_pureweight.tex, derivation of `eq:pw:Sn`: for `p > 0`, `S(p)` is `8` times the
graph-area integral `∬_{D} √(1 + (∂₁x₃)² + (∂₂x₃)²)` of `x₃ = (1-x₁^p-x₂^p)^{1/p}`
(and that integrand is integrable). -/
def stmt_2 : Prop :=
  ∀ p : ℝ, 0 < p →
    MeasureTheory.IntegrableOn (areaElt p) (D p) ∧ S p = 8 * ∫ x in D p, areaElt p x

/-- sec_pureweight.tex, Lemma `lem:pw:dirichlet` and the centring (`n = 3`): for `p > 0`,
`τ = 1 - 1/p`, `a = -2τ/3`: `S(p) = 8√3 (1/p)² Γ(1+a)³ / Γ(3-2τ) · B₃(τ)`. -/
def stmt_3 : Prop :=
  ∀ p : ℝ, 0 < p →
    S p = 8 * Real.sqrt 3 * (1 / p) ^ 2 * Real.Gamma (1 - 2 * (1 - 1 / p) / 3) ^ 3
      / Real.Gamma (3 - 2 * (1 - 1 / p)) * B (1 - 1 / p)

/-- sec_pureweight.tex, the flat-face factor (`n = 3`): for `|τ| < 1/2`,
`log(Γ(1-2τ/3)³/Γ(1-2τ)) = Σ_{k≥2} λ_k ζ(k) τ^k` with `λ_k = -2^k (1 - 3^{1-k})/k`
(Euler's `γ` cancels). -/
def stmt_4 : Prop :=
  ∀ τ : ℝ, |τ| < 1 / 2 →
    HasSum (fun k : ℕ => -(2 ^ (k + 2) * (1 - 1 / (3 : ℝ) ^ (k + 1))) / ((k + 2 : ℕ) : ℝ)
        * X05.z (k + 2) * τ ^ (k + 2))
      (Real.log (Real.Gamma (1 - 2 * τ / 3) ^ 3 / Real.Gamma (1 - 2 * τ)))

/-- sec_g7_octa.tex, proof of Lemma `g7:lem:w23` (`n = 3`):
`E[Y] = Γ(1+a+4τ/3) Γ(1+a-2τ/3)² / Γ(1+a)³ = Γ(1+2τ/3) Γ(1-4τ/3)² / Γ(1-2τ/3)³`,
for `|τ| < 3/4`. -/
def stmt_5 : Prop :=
  ∀ τ : ℝ, |τ| < 3 / 4 →
    gammaE τ (Ymean τ) =
      Real.Gamma (1 + 2 * τ / 3) * Real.Gamma (1 - 4 * τ / 3) ^ 2 / Real.Gamma (1 - 2 * τ / 3) ^ 3

/-- sec_pureweight.tex, Lemma `lem:pw:moments` for `j = 2` (`n = 3`; multiplicity vectors
`2eᵢ`, `c = (4/3,-2/3,-2/3)`, and `eᵢ+eⱼ`, `c = (1/3,1/3,-2/3)`): for `|τ| < 1/2`,
`E[Y²] = (3 Γ(1+2τ) Γ(1-2τ)² + 6 Γ(1)² Γ(1-2τ)) / (9 Γ(1-2τ/3)³)`. -/
def stmt_6 : Prop :=
  ∀ τ : ℝ, |τ| < 1 / 2 →
    gammaE τ (fun g => Ymean τ g ^ 2) =
      (3 * Real.Gamma (1 + 2 * τ) * Real.Gamma (1 - 2 * τ) ^ 2 + 6 * Real.Gamma (1 - 2 * τ))
        / (9 * Real.Gamma (1 - 2 * τ / 3) ^ 3)

/-- sec_g7_octa.tex, Lemma `g7:lem:X2` with `J = 2` (`n = 3`), both signs of `τ`:
`|B₃(τ) - (1 + E[X]/2 - E[X²]/8)| ≤ C |τ|⁶` for `|τ| ≤ δ`, `X = Y - 1`. -/
def stmt_7 : Prop :=
  ∃ C δ : ℝ, 0 < δ ∧ ∀ τ : ℝ, |τ| ≤ δ →
    |B τ - (1 + gammaE τ (fun g => Ymean τ g - 1) / 2
        - gammaE τ (fun g => (Ymean τ g - 1) ^ 2) / 8)| ≤ C * |τ| ^ 6

/-- sec_g7_octa.tex, Lemma `g7:lem:w23` for `n = 3` with Lemma `g7:lem:both` (`K = 3`), the
hand-proved part: `log((2-p)S(p)/(4√3)) = -π²/9 τ² - 16ζ(3)/9 τ³ + O(τ⁴)`, both signs. -/
def stmt_8 : Prop :=
  ∃ C δ : ℝ, 0 < δ ∧ ∀ p : ℝ, 0 < p → |1 - 1 / p| ≤ δ →
    |Real.log ((2 - p) * S p / (4 * Real.sqrt 3))
        - (-(Real.pi ^ 2 / 9) * (1 - 1 / p) ^ 2 - 16 * X05.z 3 / 9 * (1 - 1 / p) ^ 3)|
      ≤ C * |1 - 1 / p| ^ 4

/-- sec_g7_octa.tex ("the coefficients of `d⁰,d¹,d²,d³` in dimension 3" by hand):
`S(1+d)/(4√3) = 1 + d + (1-π²/9)d² + (1+π²/9-16ζ(3)/9)d³ + O(d⁴)` as `d → 0`, both sides. -/
def stmt_9 : Prop :=
  ∃ C δ : ℝ, 0 < δ ∧ ∀ d : ℝ, |d| ≤ δ →
    |S (1 + d) / (4 * Real.sqrt 3)
        - (1 + d + (1 - Real.pi ^ 2 / 9) * d ^ 2
            + (1 + Real.pi ^ 2 / 9 - 16 * X05.z 3 / 9) * d ^ 3)| ≤ C * |d| ^ 4

/-- `thm:oct` through order 4 (the target claim): as `d → 0` from either side,
`S(1+d)/(4√3) = 1 + d + (1-π²/9)d² + (1+π²/9-16ζ(3)/9)d³
  + (1-2π²/9+32ζ(3)/9-16π⁴/405)d⁴ + O(d⁵)` (`= T09.dSeries`), and
`log((2-p)S(p)/(4√3)) = -π²/9 τ² - 16ζ(3)/9 τ³ - 37π⁴/810 τ⁴ + O(τ⁵)` (`= T09.logSeries`). -/
def stmt_10 : Prop :=
  (∃ C δ : ℝ, 0 < δ ∧ ∀ d : ℝ, |d| ≤ δ →
    |S (1 + d) / (4 * Real.sqrt 3) - T09.dSeries (X05.z 3) d| ≤ C * |d| ^ 5) ∧
  (∃ C δ : ℝ, 0 < δ ∧ ∀ p : ℝ, 0 < p → |1 - 1 / p| ≤ δ →
    |Real.log ((2 - p) * S p / (4 * Real.sqrt 3)) - T09.logSeries (X05.z 3) (1 - 1 / p)|
      ≤ C * |1 - 1 / p| ^ 5)

/-- `thm:oct` at full strength (its analytic part): the `d`-series through `d⁴` with `O(d⁵)`,
and `log((2-p)S(p)/(4√3)) = -π²/9 τ² - 16ζ(3)/9 τ³ - 37π⁴/810 τ⁴ - 160ζ(5)/27 τ⁵
  + (16/81 ζ(3)² - 572/76545 π⁶) τ⁶ + O(τ⁷)`, both signs of `τ`. -/
def stmt_11 : Prop :=
  (∃ C δ : ℝ, 0 < δ ∧ ∀ d : ℝ, |d| ≤ δ →
    |S (1 + d) / (4 * Real.sqrt 3) - T09.dSeries (X05.z 3) d| ≤ C * |d| ^ 5) ∧
  (∃ C δ : ℝ, 0 < δ ∧ ∀ p : ℝ, 0 < p → |1 - 1 / p| ≤ δ →
    |Real.log ((2 - p) * S p / (4 * Real.sqrt 3))
        - (T09.logSeries (X05.z 3) (1 - 1 / p) - 160 * X05.z 5 / 27 * (1 - 1 / p) ^ 5
          + (16 / 81 * X05.z 3 ^ 2 - 572 / 76545 * Real.pi ^ 6) * (1 - 1 / p) ^ 6)|
      ≤ C * |1 - 1 / p| ^ 7)

end Y05
-- END STATEMENTS

open MeasureTheory Real Set Filter Topology

/-! ## stmt_1: integrability of the face integrand -/

namespace Y05

lemma Δ_meas : MeasurableSet Δ := by
  unfold Δ
  simp only [Set.ofPred_and]
  exact (measurableSet_le measurable_const measurable_fst).inter
    ((measurableSet_le measurable_const measurable_snd).inter
      (measurableSet_le (measurable_fst.add measurable_snd) measurable_const))

lemma Δ_sub : Δ ⊆ Icc (0 : ℝ) 1 ×ˢ Icc (0 : ℝ) 1 := by
  rintro ⟨a, b⟩ ⟨ha, hb, hab⟩
  exact ⟨⟨ha, by linarith⟩, ⟨hb, by linarith⟩⟩

lemma Δ_fin : volume Δ ≠ ⊤ := by
  refine ne_top_of_le_ne_top ?_ (measure_mono Δ_sub)
  rw [Measure.volume_eq_prod, Measure.prod_prod]
  simp

lemma face_meas (τ : ℝ) : Measurable (faceIntegrand τ) := by
  unfold faceIntegrand
  fun_prop

/-- bound for `τ ≤ 0` -/
lemma face_bd_neg {τ : ℝ} (hτ : τ ≤ 0) {t : ℝ × ℝ} (ht : t ∈ Δ) :
    ‖faceIntegrand τ t‖ ≤ 2 := by
  obtain ⟨ha, hb, hab⟩ := ht
  set a := t.1; set b := t.2
  have hc : 0 ≤ 1 - a - b := by linarith
  set c := 1 - a - b with hcdef
  have hnn : 0 ≤ faceIntegrand τ t := by unfold faceIntegrand; positivity
  rw [Real.norm_eq_abs, abs_of_nonneg hnn]
  unfold faceIntegrand
  rw [← hcdef]
  rcases eq_or_lt_of_le hτ with h0 | hneg
  · rw [h0]; simp only [neg_zero, Real.rpow_zero, mul_zero, one_mul]
    rw [Real.sqrt_le_iff]; norm_num
  by_cases hz : a = 0 ∨ b = 0 ∨ c = 0
  · have : a * b * c = 0 := by
      rcases hz with h | h | h <;> simp [h]
    rw [this, Real.zero_rpow (by linarith)]; norm_num
  push Not at hz
  obtain ⟨ha0, hb0, hc0⟩ := hz
  have pa : 0 < a := lt_of_le_of_ne ha (Ne.symm ha0)
  have pb : 0 < b := lt_of_le_of_ne hb (Ne.symm hb0)
  have pc : 0 < c := lt_of_le_of_ne hc (Ne.symm hc0)
  have e2 : ∀ x : ℝ, 0 < x → x ^ (2 * τ) = ((x ^ (-τ)) ^ 2)⁻¹ := fun x hx => by
    rw [← Real.rpow_natCast, ← Real.rpow_mul hx.le, ← Real.rpow_neg hx.le]; congr 1; push_cast; ring
  have b1 : ∀ x : ℝ, 0 < x → x ≤ 1 → x ^ (-τ) ≤ 1 := fun x hx hx1 =>
    Real.rpow_le_one hx.le hx1 (by linarith)
  rw [Real.mul_rpow (by positivity) pc.le, Real.mul_rpow pa.le pb.le, e2 a pa, e2 b pb, e2 c pc]
  have qa := b1 a pa (by linarith); have qb := b1 b pb (by linarith)
  have qc := b1 c pc (by linarith)
  have ra := Real.rpow_pos_of_pos pa (-τ); have rb := Real.rpow_pos_of_pos pb (-τ)
  have rc := Real.rpow_pos_of_pos pc (-τ)
  generalize a ^ (-τ) = α at *
  generalize b ^ (-τ) = β at *
  generalize c ^ (-τ) = γ at *
  rw [show α * β * γ * √((α ^ 2)⁻¹ + (β ^ 2)⁻¹ + (γ ^ 2)⁻¹) =
      √((β * γ) ^ 2 + (α * γ) ^ 2 + (α * β) ^ 2) by
    rw [show (β * γ) ^ 2 + (α * γ) ^ 2 + (α * β) ^ 2 =
        (α * β * γ) ^ 2 * ((α ^ 2)⁻¹ + (β ^ 2)⁻¹ + (γ ^ 2)⁻¹) by field_simp,
      Real.sqrt_mul (by positivity), Real.sqrt_sq (by positivity)]]
  rw [Real.sqrt_le_iff]
  refine ⟨by norm_num, ?_⟩
  have h1 : (β * γ) ^ 2 ≤ 1 := by
    rw [mul_pow]; exact mul_le_one₀ (pow_le_one₀ rb.le qb) (by positivity) (pow_le_one₀ rc.le qc)
  have h2 : (α * γ) ^ 2 ≤ 1 := by
    rw [mul_pow]; exact mul_le_one₀ (pow_le_one₀ ra.le qa) (by positivity) (pow_le_one₀ rc.le qc)
  have h3 : (α * β) ^ 2 ≤ 1 := by
    rw [mul_pow]; exact mul_le_one₀ (pow_le_one₀ ra.le qa) (by positivity) (pow_le_one₀ rb.le qb)
  linarith

/-- the three dominating pieces for `0 < τ < 1` -/
noncomputable def domA (τ : ℝ) (t : ℝ × ℝ) : ℝ := t.1 ^ (-τ) * t.2 ^ (-τ)
noncomputable def domG (τ : ℝ) (t : ℝ × ℝ) : ℝ := t.1 ^ (-τ) * (1 - t.2) ^ (-τ)

lemma face_bd_pos {τ : ℝ} (hτ0 : 0 < τ) (hτ1 : τ < 1) {t : ℝ × ℝ} (ht : t ∈ Δ) :
    ‖faceIntegrand τ t‖ ≤ 6 * (domA τ t + domG τ (t.2, t.2 + t.1) + domG τ (t.1, t.1 + t.2)) := by
  obtain ⟨ha, hb, hab⟩ := ht
  unfold domA domG
  simp only
  set a := t.1; set b := t.2
  have hc : 0 ≤ 1 - a - b := by linarith
  have e1 : 1 - (b + a) = 1 - a - b := by ring
  have e2 : 1 - (a + b) = 1 - a - b := by ring
  rw [e1, e2]
  set c := 1 - a - b with hcdef
  have hnn : 0 ≤ faceIntegrand τ t := by unfold faceIntegrand; positivity
  rw [Real.norm_eq_abs, abs_of_nonneg hnn]
  unfold faceIntegrand
  rw [← hcdef]
  have hs : √(a ^ (2 * τ) + b ^ (2 * τ) + c ^ (2 * τ)) ≤ 2 := by
    have := Real.rpow_le_one ha (by linarith : a ≤ 1) (by linarith : 0 ≤ 2 * τ)
    have := Real.rpow_le_one hb (by linarith : b ≤ 1) (by linarith : 0 ≤ 2 * τ)
    have := Real.rpow_le_one hc (by linarith : c ≤ 1) (by linarith : 0 ≤ 2 * τ)
    rw [Real.sqrt_le_iff]; constructor <;> linarith
  have h3 : ∀ x : ℝ, 1 / 3 ≤ x → x ^ (-τ) ≤ 3 := fun x hx => by
    calc x ^ (-τ) ≤ (1 / 3 : ℝ) ^ (-τ) := Real.rpow_le_rpow_of_nonpos (by norm_num) hx (by linarith)
      _ = (3 : ℝ) ^ τ := by rw [one_div, Real.inv_rpow (by norm_num), Real.rpow_neg (by norm_num), inv_inv]
      _ ≤ (3 : ℝ) ^ (1 : ℝ) := Real.rpow_le_rpow_of_exponent_le (by norm_num) hτ1.le
      _ = 3 := Real.rpow_one 3
  have pa := Real.rpow_nonneg ha (-τ); have pb := Real.rpow_nonneg hb (-τ)
  have pc := Real.rpow_nonneg hc (-τ)
  have hP : (a * b * c) ^ (-τ) ≤ 3 * (a ^ (-τ) * b ^ (-τ) + b ^ (-τ) * c ^ (-τ) + a ^ (-τ) * c ^ (-τ)) := by
    rw [Real.mul_rpow (by positivity) hc, Real.mul_rpow ha hb]
    have pab := mul_nonneg pa pb; have pbc := mul_nonneg pb pc; have pac := mul_nonneg pa pc
    by_cases hc3 : 1 / 3 ≤ c
    · have := h3 c hc3; nlinarith
    by_cases ha3 : 1 / 3 ≤ a
    · have := h3 a ha3; nlinarith
    have hb3 : 1 / 3 ≤ b := by linarith
    have := h3 b hb3; nlinarith
  have hP0 : 0 ≤ (a * b * c) ^ (-τ) := Real.rpow_nonneg (by positivity) _
  calc (a * b * c) ^ (-τ) * √(a ^ (2 * τ) + b ^ (2 * τ) + c ^ (2 * τ))
      ≤ (3 * (a ^ (-τ) * b ^ (-τ) + b ^ (-τ) * c ^ (-τ) + a ^ (-τ) * c ^ (-τ))) * 2 :=
        mul_le_mul hP hs (Real.sqrt_nonneg _) (by positivity)
    _ = _ := by ring

lemma pow_int {τ : ℝ} (hτ : τ < 1) : IntegrableOn (fun x : ℝ => x ^ (-τ)) (Icc 0 1) := by
  rw [← intervalIntegrable_iff_integrableOn_Icc_of_le zero_le_one]
  exact intervalIntegral.intervalIntegrable_rpow' (by linarith)

lemma pow_int' {τ : ℝ} (hτ : τ < 1) : IntegrableOn (fun x : ℝ => (1 - x) ^ (-τ)) (Icc 0 1) := by
  rw [← intervalIntegrable_iff_integrableOn_Icc_of_le zero_le_one]
  have := (intervalIntegral.intervalIntegrable_rpow' (a := 1 - 0) (b := 1 - 1)
    (by linarith : -1 < -τ)).comp_sub_left 1
  simpa using this

lemma box_int {f g : ℝ → ℝ} (hf : IntegrableOn f (Icc 0 1)) (hg : IntegrableOn g (Icc 0 1)) :
    Integrable ((Icc (0 : ℝ) 1 ×ˢ Icc (0 : ℝ) 1).indicator (fun z : ℝ × ℝ => f z.1 * g z.2)) := by
  rw [integrable_indicator_iff (measurableSet_Icc.prod measurableSet_Icc)]
  show IntegrableOn _ _ (volume.prod volume)
  rw [IntegrableOn, ← Measure.prod_restrict]
  exact Integrable.mul_prod hf hg

end Y05

theorem Y05.proof_1 : Y05.stmt_1 := by
  intro p hp
  set τ := 1 - 1 / p with hτdef
  have hτ1 : τ < 1 := by have : 0 < 1 / p := by positivity
                         linarith
  rcases le_or_gt τ 0 with hτ0 | hτ0
  · refine Integrable.mono' (integrableOn_const Δ_fin (C := (2 : ℝ)))
      (face_meas τ).aestronglyMeasurable ?_
    exact (ae_restrict_iff' Δ_meas).mpr (Eventually.of_forall fun t ht => face_bd_neg hτ0 ht)
  · have iA : Integrable ((Icc (0 : ℝ) 1 ×ˢ Icc (0 : ℝ) 1).indicator (domA τ)) :=
      box_int (pow_int hτ1) (pow_int hτ1)
    have iG : Integrable ((Icc (0 : ℝ) 1 ×ˢ Icc (0 : ℝ) 1).indicator (domG τ)) :=
      box_int (pow_int hτ1) (pow_int' hτ1)
    have iG1 := (measurePreserving_prod_add (volume : Measure ℝ) (volume : Measure ℝ)
      ).integrable_comp_of_integrable iG
    have iG2 := (measurePreserving_prod_add_swap (volume : Measure ℝ) (volume : Measure ℝ)
      ).integrable_comp_of_integrable iG
    have hsum : IntegrableOn (fun t => 6 * (domA τ t + domG τ (t.2, t.2 + t.1)
        + domG τ (t.1, t.1 + t.2))) Δ := by
      have hA : IntegrableOn (domA τ) Δ := by
        have := iA.integrableOn (s := Δ)
        refine this.congr_fun (fun t ht => ?_) Δ_meas
        rw [indicator_of_mem (Δ_sub ht)]
      have h1 : IntegrableOn (fun t : ℝ × ℝ => domG τ (t.1, t.1 + t.2)) Δ := by
        refine (iG1.integrableOn (s := Δ)).congr_fun (fun t ht => ?_) Δ_meas
        obtain ⟨ha, hb, hab⟩ := ht
        simp only [Function.comp]
        rw [indicator_of_mem]
        exact ⟨⟨ha, by linarith⟩, ⟨by linarith, hab⟩⟩
      have h2 : IntegrableOn (fun t : ℝ × ℝ => domG τ (t.2, t.2 + t.1)) Δ := by
        refine (iG2.integrableOn (s := Δ)).congr_fun (fun t ht => ?_) Δ_meas
        obtain ⟨ha, hb, hab⟩ := ht
        simp only [Function.comp]
        rw [indicator_of_mem]
        exact ⟨⟨hb, by linarith⟩, ⟨by linarith, by linarith⟩⟩
      exact ((hA.add h2).add h1).const_mul 6
    refine Integrable.mono' hsum (face_meas τ).aestronglyMeasurable ?_
    exact (ae_restrict_iff' Δ_meas).mpr
      (Eventually.of_forall fun t ht => face_bd_pos hτ0 hτ1 ht)


/-! ## stmt_4: the flat-face factor -/

namespace Y05

/-- two-factor log-Gamma series (the argument of X05.proof_3 with `|u|,|v|,|u+v| < 1`). -/
lemma logGamma2 {u v r : ℝ} (hr : r < 1) (hu : |u| ≤ r) (hv : |v| ≤ r) (hη : |u + v| ≤ r) :
    HasSum (fun k : ℕ => (-1 : ℝ) ^ (k + 2) * X05.z (k + 2) / ((k + 2 : ℕ) : ℝ) *
        (u ^ (k + 2) + v ^ (k + 2) - (u + v) ^ (k + 2)))
      (Real.log (Real.Gamma (1 + u) * Real.Gamma (1 + v) / Real.Gamma (1 + (u + v)))) := by
  set η := u + v with hηdef
  have huv : u + v = η := rfl
  have hr0 : 0 ≤ r := (abs_nonneg u).trans hu
  have hu' : |u| < 1 := by linarith
  have hv' : |v| < 1 := by linarith
  have hη' : |η| < 1 := by linarith
  have hApos : 0 < Real.Gamma (1 + u) * Real.Gamma (1 + v) / Real.Gamma (1 + η) := by
    have := Real.Gamma_pos_of_pos (show 0 < 1 + u by cases abs_lt.mp hu'; linarith)
    have := Real.Gamma_pos_of_pos (show 0 < 1 + v by cases abs_lt.mp hv'; linarith)
    have := Real.Gamma_pos_of_pos (show 0 < 1 + η by cases abs_lt.mp hη'; linarith)
    positivity
  set f : ℕ × ℕ → ℝ := fun p => (-1 : ℝ) ^ (p.2 + 2) * (1 / ((p.1 : ℝ) + 1) ^ (p.2 + 2)) /
      ((p.2 + 2 : ℕ) : ℝ) * (u ^ (p.2 + 2) + v ^ (p.2 + 2) - η ^ (p.2 + 2)) with hfdef
  have hfs : Summable f := by
    have hg1 : Summable (fun k : ℕ => 1 / ((k : ℝ) + 1) ^ 2) :=
      (X05.z_hasSum (by norm_num)).summable
    have hg2 : Summable (fun j : ℕ => 3 * r ^ j) :=
      (summable_geometric_of_lt_one hr0 hr).mul_left 3
    refine Summable.of_norm_bounded (hg1.mul_of_nonneg hg2 (fun k => by positivity)
      (fun j => by positivity)) ?_
    rintro ⟨k, j⟩
    simp only [hfdef, Real.norm_eq_abs]
    have hk1 : (1 : ℝ) ≤ (k : ℝ) + 1 := by have := k.cast_nonneg (α := ℝ); linarith
    rw [abs_mul, abs_div, abs_mul, abs_pow, abs_neg, abs_one, one_pow, one_mul,
      abs_of_pos (by positivity : (0 : ℝ) < 1 / ((k : ℝ) + 1) ^ (j + 2)),
      abs_of_pos (by positivity : (0 : ℝ) < ((j + 2 : ℕ) : ℝ))]
    have b1 : 1 / ((k : ℝ) + 1) ^ (j + 2) ≤ 1 / ((k : ℝ) + 1) ^ 2 := by
      apply one_div_le_one_div_of_le (by positivity)
      exact pow_le_pow_right₀ hk1 (by omega)
    have pw : ∀ x : ℝ, |x| ≤ r → |x| ^ (j + 2) ≤ r ^ j := fun x hx =>
      (pow_le_pow_left₀ (abs_nonneg x) hx _).trans
        (by rw [pow_add]; exact mul_le_of_le_one_right (by positivity)
              (pow_le_one₀ hr0 hr.le))
    have b2 : |u ^ (j + 2) + v ^ (j + 2) - η ^ (j + 2)| ≤ 3 * r ^ j := by
      have pu := pw u hu; have pv := pw v hv; have pη := pw η hη
      calc |u ^ (j + 2) + v ^ (j + 2) - η ^ (j + 2)|
          ≤ |u ^ (j + 2)| + |v ^ (j + 2)| + |η ^ (j + 2)| := by
            have := abs_sub (u ^ (j + 2) + v ^ (j + 2)) (η ^ (j + 2))
            have := abs_add_le (u ^ (j + 2)) (v ^ (j + 2))
            linarith
        _ ≤ 3 * r ^ j := by rw [abs_pow, abs_pow, abs_pow]; linarith
    have hj : (1 : ℝ) ≤ ((j + 2 : ℕ) : ℝ) := by
      push_cast; have := j.cast_nonneg (α := ℝ); linarith
    calc 1 / ((k : ℝ) + 1) ^ (j + 2) / ((j + 2 : ℕ) : ℝ) *
          |u ^ (j + 2) + v ^ (j + 2) - η ^ (j + 2)|
        ≤ 1 / ((k : ℝ) + 1) ^ 2 / 1 * (3 * r ^ j) := by gcongr
      _ = _ := by simp
  have hf := hfs.hasSum
  have hT := hf.prod_fiberwise (fun k => X05.log_series_k huv hu' hv' hη' k)
  have hf' : HasSum (f ∘ Equiv.prodComm ℕ ℕ) (∑' p, f p) :=
    (Equiv.prodComm ℕ ℕ).hasSum_iff.mpr hf
  have hZ := hf'.prod_fiberwise (fun j => by
    have := ((X05.z_hasSum (show 1 < j + 2 by omega)).mul_left ((-1 : ℝ) ^ (j + 2))).div_const
      ((j + 2 : ℕ) : ℝ) |>.mul_right (u ^ (j + 2) + v ^ (j + 2) - η ^ (j + 2))
    exact this)
  have hlim1 := (tendsto_add_atTop_iff_nat 1).mpr hT.tendsto_sum_nat
  have hlim2 : Tendsto (fun n : ℕ => ∑ k ∈ Finset.range (n + 1), (Real.log (1 + η / ((k : ℝ) + 1))
      - Real.log (1 + u / ((k : ℝ) + 1)) - Real.log (1 + v / ((k : ℝ) + 1)))) atTop
      (𝓝 (Real.log (Real.Gamma (1 + u) * Real.Gamma (1 + v) / Real.Gamma (1 + η)))) := by
    have := ((Real.continuousAt_log hApos.ne').tendsto).comp
      (X05.R_tendsto huv hu' hv' hη' (by cases abs_lt.mp hη'; linarith))
    refine this.congr fun n => ?_
    simp only [Function.comp]
    exact X05.log_R hu' hv' hη' (n + 1)
  have heq := tendsto_nhds_unique hlim1 hlim2
  rw [heq] at hZ
  convert hZ using 1

end Y05

theorem Y05.proof_4 : Y05.stmt_4 := by
  intro τ hτ
  set u := -2 * τ / 3 with hu
  have hr : 2 * |τ| < 1 := by linarith
  have h1 : |u| ≤ 2 * |τ| := by
    rw [hu, abs_div, abs_mul]; norm_num; linarith [abs_nonneg τ]
  have h2 : |u + u| ≤ 2 * |τ| := by
    rw [show u + u = -4 * τ / 3 by rw [hu]; ring, abs_div, abs_mul]; norm_num
    linarith [abs_nonneg τ]
  have h3 : |u + u + u| ≤ 2 * |τ| := by
    rw [show u + u + u = -2 * τ by rw [hu]; ring, abs_mul]; norm_num
  have A := Y05.logGamma2 hr h1 h1 h2
  have B := Y05.logGamma2 hr h2 h1 h3
  have hG : ∀ x : ℝ, |x| ≤ 2 * |τ| → 0 < Real.Gamma (1 + x) := fun x hx =>
    Real.Gamma_pos_of_pos (by cases abs_lt.mp (lt_of_le_of_lt hx hr); linarith)
  have gu := hG u h1; have g2 := hG _ h2; have g3 := hG _ h3
  have hlog : Real.log (Real.Gamma (1 - 2 * τ / 3) ^ 3 / Real.Gamma (1 - 2 * τ)) =
      Real.log (Real.Gamma (1 + u) * Real.Gamma (1 + u) / Real.Gamma (1 + (u + u))) +
      Real.log (Real.Gamma (1 + (u + u)) * Real.Gamma (1 + u) / Real.Gamma (1 + (u + u + u))) := by
    rw [← Real.log_mul (by positivity) (by positivity)]
    congr 1
    rw [show 1 - 2 * τ / 3 = 1 + u by rw [hu]; ring,
      show 1 - 2 * τ = 1 + (u + u + u) by rw [hu]; ring]
    generalize Real.Gamma (1 + u) = a at *
    generalize Real.Gamma (1 + (u + u)) = b at *
    generalize Real.Gamma (1 + (u + u + u)) = c at *
    field_simp
  rw [hlog]
  convert A.add B using 1
  ext k
  set n := k + 2
  have key : ∀ x : ℝ, (-1 : ℝ) ^ n * (-x) ^ n = x ^ n := fun x => by
    rw [← mul_pow, neg_one_mul, neg_neg]
  have e1 : u + u = -(2 * (2 * τ / 3)) := by rw [hu]; ring
  have e2 : u + u + u = -(2 * τ) := by rw [hu]; ring
  have e0 : u = -(2 * τ / 3) := by rw [hu]; ring
  have hn : ((n : ℕ) : ℝ) ≠ 0 := by positivity
  calc -(2 ^ n * (1 - 1 / (3 : ℝ) ^ (k + 1))) / ((n : ℕ) : ℝ) * X05.z n * τ ^ n
      = X05.z n / ((n : ℕ) : ℝ) * (3 * (2 * τ / 3) ^ n - (2 * τ) ^ n) := by
        have h3 : (3 : ℝ) ^ n = 3 * 3 ^ (k + 1) := by ring
        have h3n : (3 : ℝ) ^ (k + 1) ≠ 0 := by positivity
        rw [div_pow, mul_pow, h3]
        field_simp
        ring
    _ = X05.z n / ((n : ℕ) : ℝ) * (3 * ((-1 : ℝ) ^ n * u ^ n) - (-1 : ℝ) ^ n * (u + u + u) ^ n) := by
        rw [e2, e0, key, key]
    _ = _ := by ring

/-! ## stmt_5, stmt_6: Gamma moments -/

namespace Y05

/-- the positive octant -/
abbrev Q : Set (ℝ × ℝ × ℝ) := Ioi (0 : ℝ) ×ˢ Ioi (0 : ℝ) ×ˢ Ioi (0 : ℝ)

/-- a Gamma monomial `e^{-(g₁+g₂+g₃)} g₁^x g₂^y g₃^z` -/
noncomputable def mono (x y z : ℝ) (g : ℝ × ℝ × ℝ) : ℝ :=
  Real.exp (-(g.1 + g.2.1 + g.2.2)) * (g.1 ^ x * g.2.1 ^ y * g.2.2 ^ z)

lemma gamma1 {x : ℝ} (hx : -1 < x) :
    IntegrableOn (fun t : ℝ => Real.exp (-t) * t ^ x) (Ioi 0) ∧
      ∫ t in Ioi (0 : ℝ), Real.exp (-t) * t ^ x = Real.Gamma (1 + x) := by
  have h : 0 < 1 + x := by linarith
  have e : 1 + x - 1 = x := by ring
  have := Real.GammaIntegral_convergent h
  rw [e] at this
  refine ⟨this, ?_⟩
  rw [Real.Gamma_eq_integral h, e]

lemma mono_split (x y z : ℝ) :
    mono x y z = fun g => (Real.exp (-g.1) * g.1 ^ x) *
      ((fun h : ℝ × ℝ => (Real.exp (-h.1) * h.1 ^ y) * (Real.exp (-h.2) * h.2 ^ z)) g.2) := by
  funext g
  simp only [mono, neg_add, Real.exp_add]
  ring

lemma mono_int {x y z : ℝ} (hx : -1 < x) (hy : -1 < y) (hz : -1 < z) :
    IntegrableOn (mono x y z) Q ∧
      ∫ g in Q, mono x y z g = Real.Gamma (1 + x) * Real.Gamma (1 + y) * Real.Gamma (1 + z) := by
  obtain ⟨ix, vx⟩ := gamma1 hx
  obtain ⟨iy, vy⟩ := gamma1 hy
  obtain ⟨iz, vz⟩ := gamma1 hz
  have iyz : IntegrableOn (fun h : ℝ × ℝ => (Real.exp (-h.1) * h.1 ^ y) *
      (Real.exp (-h.2) * h.2 ^ z)) (Ioi (0 : ℝ) ×ˢ Ioi (0 : ℝ)) (volume.prod volume) := by
    rw [IntegrableOn, ← Measure.prod_restrict]
    exact Integrable.mul_prod iy iz
  rw [mono_split]
  constructor
  · show IntegrableOn _ _ (volume.prod volume)
    rw [IntegrableOn, ← Measure.prod_restrict]
    exact Integrable.mul_prod ix iyz
  · refine Eq.trans (setIntegral_prod_mul (μ := (volume : Measure ℝ))
      (ν := (volume : Measure (ℝ × ℝ)))
      (fun t : ℝ => Real.exp (-t) * t ^ x)
      (fun h : ℝ × ℝ => (Real.exp (-h.1) * h.1 ^ y) * (Real.exp (-h.2) * h.2 ^ z))
      (Ioi (0 : ℝ)) (Ioi (0 : ℝ) ×ˢ Ioi (0 : ℝ))) ?_
    rw [vx]
    refine Eq.trans (congrArg _ (setIntegral_prod_mul (μ := (volume : Measure ℝ))
      (ν := (volume : Measure ℝ))
      (fun t : ℝ => Real.exp (-t) * t ^ y) (fun t : ℝ => Real.exp (-t) * t ^ z)
      (Ioi (0 : ℝ)) (Ioi (0 : ℝ)))) ?_
    rw [vy, vz]
    ring


lemma Q_meas : MeasurableSet Q :=
  measurableSet_Ioi.prod (measurableSet_Ioi.prod measurableSet_Ioi)

lemma rp_pow {t : ℝ} (ht : 0 < t) (c : ℝ) (n : ℕ) : t ^ (c * n) = (t ^ c) ^ n := by
  rw [Real.rpow_mul ht.le, Real.rpow_natCast]

lemma rp_neg {t : ℝ} (ht : 0 < t) (c : ℝ) (n : ℕ) : t ^ (-(c * n)) = ((t ^ c) ^ n)⁻¹ := by
  rw [Real.rpow_neg ht.le, rp_pow ht]

lemma pt5 (τ : ℝ) {g : ℝ × ℝ × ℝ} (hg : g ∈ Q) :
    (g.1 * g.2.1 * g.2.2) ^ (-2 * τ / 3) * Real.exp (-(g.1 + g.2.1 + g.2.2)) * Ymean τ g =
      1 / 3 * (mono (2 * τ / 3) (-(4 * τ / 3)) (-(4 * τ / 3)) g
        + mono (-(4 * τ / 3)) (2 * τ / 3) (-(4 * τ / 3)) g
        + mono (-(4 * τ / 3)) (-(4 * τ / 3)) (2 * τ / 3) g) := by
  obtain ⟨ha, hb, hc⟩ := hg
  simp only [mem_Ioi] at ha hb hc
  unfold Ymean mono
  set c := 2 * τ / 3 with hcdef
  have e1 : -2 * τ / 3 = -(c * (1 : ℕ)) := by rw [hcdef]; push_cast; ring
  have e2 : 2 * τ = c * (3 : ℕ) := by rw [hcdef]; push_cast; ring
  have e4 : -(4 * τ / 3) = -(c * (2 : ℕ)) := by rw [hcdef]; push_cast; ring
  rw [Real.mul_rpow (mul_pos ha hb).le hc.le, Real.mul_rpow ha.le hb.le,
    Real.mul_rpow (mul_pos ha hb).le hc.le, Real.mul_rpow ha.le hb.le, e1, e2, e4]
  simp only [rp_pow ha, rp_pow hb, rp_pow hc, rp_neg ha, rp_neg hb, rp_neg hc]
  have qa := Real.rpow_pos_of_pos ha c
  have qb := Real.rpow_pos_of_pos hb c
  have qc := Real.rpow_pos_of_pos hc c
  generalize hA : g.1 ^ c = A at *
  generalize hB : g.2.1 ^ c = B at *
  generalize hC : g.2.2 ^ c = C at *
  field_simp

lemma pt6 (τ : ℝ) {g : ℝ × ℝ × ℝ} (hg : g ∈ Q) :
    (g.1 * g.2.1 * g.2.2) ^ (-2 * τ / 3) * Real.exp (-(g.1 + g.2.1 + g.2.2)) * Ymean τ g ^ 2 =
      1 / 9 * (mono (2 * τ) (-(2 * τ)) (-(2 * τ)) g
        + mono (-(2 * τ)) (2 * τ) (-(2 * τ)) g
        + mono (-(2 * τ)) (-(2 * τ)) (2 * τ) g
        + 2 * (mono 0 0 (-(2 * τ)) g + mono 0 (-(2 * τ)) 0 g + mono (-(2 * τ)) 0 0 g)) := by
  obtain ⟨ha, hb, hc⟩ := hg
  simp only [mem_Ioi] at ha hb hc
  unfold Ymean mono
  set c := 2 * τ / 3 with hcdef
  have e1 : -2 * τ / 3 = -(c * (1 : ℕ)) := by rw [hcdef]; push_cast; ring
  have e2 : 2 * τ = c * (3 : ℕ) := by rw [hcdef]; push_cast; ring
  rw [Real.mul_rpow (mul_pos ha hb).le hc.le, Real.mul_rpow ha.le hb.le,
    Real.mul_rpow (mul_pos ha hb).le hc.le, Real.mul_rpow ha.le hb.le, e1, e2]
  simp only [rp_pow ha, rp_pow hb, rp_pow hc, rp_neg ha, rp_neg hb, rp_neg hc, Real.rpow_zero]
  have qa := Real.rpow_pos_of_pos ha c
  have qb := Real.rpow_pos_of_pos hb c
  have qc := Real.rpow_pos_of_pos hc c
  generalize hA : g.1 ^ c = A at *
  generalize hB : g.2.1 ^ c = B at *
  generalize hC : g.2.2 ^ c = C at *
  field_simp
  ring

lemma gE_eq {τ : ℝ} {F : ℝ × ℝ × ℝ → ℝ} {G : ℝ × ℝ × ℝ → ℝ}
    (h : ∀ g ∈ Q, (g.1 * g.2.1 * g.2.2) ^ (-2 * τ / 3) * Real.exp (-(g.1 + g.2.1 + g.2.2)) * F g
      = G g) :
    gammaE τ F = (Real.Gamma (1 - 2 * τ / 3) ^ 3)⁻¹ * ∫ g in Q, G g := by
  unfold gammaE
  rw [setIntegral_congr_fun Q_meas (fun g hg => h g hg)]

end Y05

theorem Y05.proof_5 : Y05.stmt_5 := by
  intro τ hτ
  obtain ⟨h1, h2⟩ := abs_lt.mp hτ
  have ha : -1 < 2 * τ / 3 := by linarith
  have hb : -1 < -(4 * τ / 3) := by linarith
  obtain ⟨i1, v1⟩ := mono_int ha hb hb
  obtain ⟨i2, v2⟩ := mono_int hb ha hb
  obtain ⟨i3, v3⟩ := mono_int hb hb ha
  have j12 : Integrable (fun g => mono (2 * τ / 3) (-(4 * τ / 3)) (-(4 * τ / 3)) g
      + mono (-(4 * τ / 3)) (2 * τ / 3) (-(4 * τ / 3)) g) (volume.restrict Q) := i1.add i2
  rw [gE_eq (fun g hg => pt5 τ hg), integral_const_mul, integral_add j12 i3, integral_add i1 i2,
    v1, v2, v3, show 1 + -(4 * τ / 3) = 1 - 4 * τ / 3 by ring]
  field_simp
  ring

theorem Y05.proof_6 : Y05.stmt_6 := by
  intro τ hτ
  obtain ⟨h1, h2⟩ := abs_lt.mp hτ
  have ha : -1 < 2 * τ := by linarith
  have hb : -1 < -(2 * τ) := by linarith
  have h0 : (-1 : ℝ) < 0 := by norm_num
  obtain ⟨i1, v1⟩ := mono_int ha hb hb
  obtain ⟨i2, v2⟩ := mono_int hb ha hb
  obtain ⟨i3, v3⟩ := mono_int hb hb ha
  obtain ⟨i4, v4⟩ := mono_int h0 h0 hb
  obtain ⟨i5, v5⟩ := mono_int h0 hb h0
  obtain ⟨i6, v6⟩ := mono_int hb h0 h0
  set a := 2 * τ
  set b := -(2 * τ)
  have j12 : Integrable (fun g => mono a b b g + mono b a b g) (volume.restrict Q) := i1.add i2
  have j123 : Integrable (fun g => mono a b b g + mono b a b g + mono b b a g)
      (volume.restrict Q) := j12.add i3
  have j45 : Integrable (fun g => mono 0 0 b g + mono 0 b 0 g) (volume.restrict Q) := i4.add i5
  have j456 : Integrable (fun g => mono 0 0 b g + mono 0 b 0 g + mono b 0 0 g)
      (volume.restrict Q) := j45.add i6
  have j456' : Integrable (fun g => 2 * (mono 0 0 b g + mono 0 b 0 g + mono b 0 0 g))
      (volume.restrict Q) := j456.const_mul 2
  rw [gE_eq (fun g hg => pt6 τ hg), integral_const_mul, integral_add j123 j456',
    integral_const_mul, integral_add j12 i3, integral_add i1 i2,
    integral_add j45 i6, integral_add i4 i5, v1, v2, v3, v4, v5, v6,
    show 1 + b = 1 - 2 * τ by ring, add_zero, Real.Gamma_one]
  have : 0 < Real.Gamma (1 - 2 * τ / 3) := Real.Gamma_pos_of_pos (by linarith)
  field_simp
  ring

/-! ## stmt_3: the Dirichlet–Gamma factorization -/

namespace Y05

/-- the gammaE integrand of `B` -/
noncomputable def Fb (τ : ℝ) (g : ℝ × ℝ × ℝ) : ℝ :=
  (g.1 * g.2.1 * g.2.2) ^ (-2 * τ / 3) * Real.exp (-(g.1 + g.2.1 + g.2.2)) *
    Real.sqrt (Ymean τ g)

lemma Fb_form (τ : ℝ) {g : ℝ × ℝ × ℝ} (hg : g ∈ Q) :
    Fb τ g = Real.exp (-(g.1 + g.2.1 + g.2.2)) * (g.1 * g.2.1 * g.2.2) ^ (-τ) *
      Real.sqrt (g.1 ^ (2 * τ) + g.2.1 ^ (2 * τ) + g.2.2 ^ (2 * τ)) / Real.sqrt 3 := by
  obtain ⟨ha, hb, hc⟩ := hg
  simp only [mem_Ioi] at ha hb hc
  have hP : 0 < g.1 * g.2.1 * g.2.2 := by positivity
  unfold Fb Ymean
  generalize g.1 * g.2.1 * g.2.2 = P at *
  rw [Real.sqrt_div' _ (by positivity), Real.sqrt_mul (by norm_num),
    Real.sqrt_eq_rpow (P ^ (2 * τ / 3)), ← Real.rpow_mul hP.le,
    show -2 * τ / 3 = -τ + 2 * τ / 3 * (1 / 2) by ring, Real.rpow_add hP]
  have h1 := Real.rpow_pos_of_pos hP (2 * τ / 3 * (1 / 2))
  have h3 : (0 : ℝ) < Real.sqrt 3 := by positivity
  field_simp

lemma Fb_meas (τ : ℝ) : Measurable (Fb τ) := by
  unfold Fb Ymean; fun_prop

lemma Fb_nn (τ : ℝ) (g : ℝ × ℝ × ℝ) (hg : g ∈ Q) : 0 ≤ Fb τ g := by
  obtain ⟨ha, hb, hc⟩ := hg
  simp only [mem_Ioi] at ha hb hc
  unfold Fb; positivity

lemma Fb_bd (τ : ℝ) {g : ℝ × ℝ × ℝ} (hg : g ∈ Q) :
    Fb τ g ≤ mono 0 (-τ) (-τ) g + mono (-τ) 0 (-τ) g + mono (-τ) (-τ) 0 g := by
  have hg' := hg
  rw [Fb_form τ hg]
  obtain ⟨ha, hb, hc⟩ := hg
  simp only [mem_Ioi] at ha hb hc
  unfold mono
  rw [Real.mul_rpow (by positivity) hc.le, Real.mul_rpow ha.le hb.le]
  simp only [Real.rpow_zero]
  have e : ∀ x : ℝ, 0 < x → x ^ (2 * τ) = (x ^ τ) ^ 2 := fun x hx => by
    rw [← Real.rpow_natCast, ← Real.rpow_mul hx.le]; congr 1; push_cast; ring
  have ek : ∀ x : ℝ, 0 < x → x ^ (-τ) * x ^ τ = 1 := fun x hx => by
    rw [← Real.rpow_add hx]; simp
  rw [e _ ha, e _ hb, e _ hc]
  have pa := Real.rpow_pos_of_pos ha τ; have pb := Real.rpow_pos_of_pos hb τ
  have pc := Real.rpow_pos_of_pos hc τ
  have na := Real.rpow_pos_of_pos ha (-τ); have nb := Real.rpow_pos_of_pos hb (-τ)
  have nc := Real.rpow_pos_of_pos hc (-τ)
  have ka := ek _ ha; have kb := ek _ hb; have kc := ek _ hc
  have hs : Real.sqrt ((g.1 ^ τ) ^ 2 + (g.2.1 ^ τ) ^ 2 + (g.2.2 ^ τ) ^ 2) ≤
      g.1 ^ τ + g.2.1 ^ τ + g.2.2 ^ τ := by
    rw [Real.sqrt_le_iff]; constructor
    · positivity
    · nlinarith [mul_pos pa pb, mul_pos pb pc, mul_pos pa pc]
  have h3 : (1 : ℝ) ≤ Real.sqrt 3 := by
    rw [show (1 : ℝ) = Real.sqrt 1 by simp]; exact Real.sqrt_le_sqrt (by norm_num)
  have hE := Real.exp_pos (-(g.1 + g.2.1 + g.2.2))
  have hs0 := Real.sqrt_nonneg ((g.1 ^ τ) ^ 2 + (g.2.1 ^ τ) ^ 2 + (g.2.2 ^ τ) ^ 2)
  generalize Real.exp (-(g.1 + g.2.1 + g.2.2)) = E at *
  generalize g.1 ^ τ = A at *
  generalize g.2.1 ^ τ = B' at *
  generalize g.2.2 ^ τ = C at *
  generalize g.1 ^ (-τ) = A' at *
  generalize g.2.1 ^ (-τ) = Bn at *
  generalize g.2.2 ^ (-τ) = C' at *
  generalize Real.sqrt (A ^ 2 + B' ^ 2 + C ^ 2) = s at *
  rw [div_le_iff₀ (by linarith)]
  have : E * (A' * Bn * C') * s ≤ E * (A' * Bn * C') * (A + B' + C) :=
    mul_le_mul_of_nonneg_left hs (by positivity)
  have hR : E * (A' * Bn * C') * (A + B' + C) =
      E * (1 * Bn * C') + E * (A' * 1 * C') + E * (A' * Bn * 1) := by
    linear_combination (E * Bn * C') * ka + (E * A' * C') * kb + (E * A' * Bn) * kc
  have hR0 : 0 ≤ E * (1 * Bn * C') + E * (A' * 1 * C') + E * (A' * Bn * 1) := by positivity
  nlinarith

lemma Fb_int {τ : ℝ} (hτ : τ < 1) : IntegrableOn (Fb τ) Q := by
  have h0 : (-1 : ℝ) < 0 := by norm_num
  have h1 : (-1 : ℝ) < -τ := by linarith
  have i := ((mono_int h0 h1 h1).1.add (mono_int h1 h0 h1).1).add (mono_int h1 h1 h0).1
  refine Integrable.mono' i (Fb_meas τ).aestronglyMeasurable ?_
  refine (ae_restrict_iff' Q_meas).mpr (Eventually.of_forall fun g hg => ?_)
  rw [Real.norm_eq_abs, abs_of_nonneg (Fb_nn τ g hg)]
  exact Fb_bd τ hg

/-- the open triangle -/
def Δo : Set (ℝ × ℝ) := {t | 0 < t.1 ∧ 0 < t.2 ∧ 0 < 1 - t.1 - t.2}

lemma Δo_meas : MeasurableSet Δo := by
  unfold Δo
  simp only [Set.ofPred_and]
  exact (measurableSet_lt measurable_const measurable_fst).inter
    ((measurableSet_lt measurable_const measurable_snd).inter
      (measurableSet_lt measurable_const (by fun_prop)))

lemma Δo_ae : Δo =ᵐ[volume] Δ := by
  rw [ae_eq_set]
  constructor
  · have : Δo \ Δ = ∅ := by
      ext t; simp only [Δo, Δ, Set.mem_sdiff, mem_ofPred_eq, mem_empty_iff_false, iff_false, not_and,
        not_le]
      intro h1 h2; obtain ⟨a, b, c⟩ := h1; have := h2 a.le b.le; linarith
    rw [this]; simp
  · have hsub : Δ \ Δo ⊆ ({0} ×ˢ univ) ∪ (univ ×ˢ {0}) ∪
        ((fun z : ℝ × ℝ => (z.1, z.1 + z.2)) ⁻¹' (univ ×ˢ {1})) := by
      rintro ⟨a, b⟩ ⟨⟨ha, hb, hab⟩, hn⟩
      simp only [Δo, mem_ofPred_eq, not_and, not_lt] at hn
      simp only [mem_union, mem_prod, mem_singleton_iff, mem_univ, and_true, true_and,
        mem_preimage]
      by_cases h1 : a = 0
      · left; left; exact h1
      by_cases h2 : b = 0
      · left; right; exact h2
      right
      have := hn (lt_of_le_of_ne ha (Ne.symm h1)) (lt_of_le_of_ne hb (Ne.symm h2))
      linarith
    refine measure_mono_null hsub ?_
    refine measure_union_null (measure_union_null ?_ ?_) ?_
    · rw [Measure.volume_eq_prod, Measure.prod_prod]; simp
    · rw [Measure.volume_eq_prod, Measure.prod_prod]; simp
    · rw [Measure.volume_eq_prod, (measurePreserving_prod_add (volume : Measure ℝ)
        (volume : Measure ℝ)).measure_preimage
        (MeasurableSet.univ.prod (measurableSet_singleton 1)).nullMeasurableSet,
        Measure.prod_prod]
      simp

lemma Fb_scale (τ : ℝ) {r : ℝ} (hr : 0 < r) {t : ℝ × ℝ} (ht : t ∈ Δo) :
    Fb τ (r * t.1, r * t.2, r - r * t.1 - r * t.2) =
      r ^ (-2 * τ) * Real.exp (-r) * faceIntegrand τ t / Real.sqrt 3 := by
  obtain ⟨ha, hb, hc⟩ := ht
  have hq : (r * t.1, r * t.2, r - r * t.1 - r * t.2) ∈ Q := by
    refine ⟨?_, ?_, ?_⟩ <;> simp only [mem_Ioi]
    · positivity
    · positivity
    · nlinarith
  rw [Fb_form τ hq]
  simp only
  unfold faceIntegrand
  rw [show r - r * t.1 - r * t.2 = r * (1 - t.1 - t.2) by ring,
    show r * t.1 + r * t.2 + r * (1 - t.1 - t.2) = r by ring,
    show r * t.1 * (r * t.2) * (r * (1 - t.1 - t.2)) = r ^ 3 * (t.1 * t.2 * (1 - t.1 - t.2)) by
      ring,
    Real.mul_rpow (by positivity) (by positivity),
    Real.mul_rpow hr.le ha.le, Real.mul_rpow hr.le hb.le, Real.mul_rpow hr.le hc.le,
    show r ^ (2 * τ) * t.1 ^ (2 * τ) + r ^ (2 * τ) * t.2 ^ (2 * τ) +
      r ^ (2 * τ) * (1 - t.1 - t.2) ^ (2 * τ) = (r ^ τ) ^ 2 *
      (t.1 ^ (2 * τ) + t.2 ^ (2 * τ) + (1 - t.1 - t.2) ^ (2 * τ)) by
      rw [← Real.rpow_natCast, ← Real.rpow_mul hr.le]; push_cast; ring_nf,
    Real.sqrt_mul (by positivity), Real.sqrt_sq (by positivity),
    ← Real.rpow_natCast, ← Real.rpow_mul hr.le]
  have : r ^ ((3 : ℕ) * -τ : ℝ) * r ^ τ = r ^ (-2 * τ) := by
    rw [← Real.rpow_add hr]; congr 1; push_cast; ring
  rw [← this]
  ring

lemma Fb_zero {τ : ℝ} {r : ℝ} (hr : r ≤ 0) (u : ℝ × ℝ) :
    Q.indicator (Fb τ) (u.1, u.2, r - u.1 - u.2) = 0 := by
  rw [indicator_of_notMem]
  rintro ⟨ha, hb, hc⟩
  simp only [mem_Ioi] at ha hb hc
  linarith

lemma inner_r (τ : ℝ) {r : ℝ} (hr : 0 < r) :
    ∫ u : ℝ × ℝ, Q.indicator (Fb τ) (u.1, u.2, r - u.1 - u.2) =
      r ^ (2 - 2 * τ) * Real.exp (-r) * (∫ t in Δ, faceIntegrand τ t) / Real.sqrt 3 := by
  have hs := Measure.integral_comp_smul (volume : Measure (ℝ × ℝ))
    (fun u : ℝ × ℝ => Q.indicator (Fb τ) (u.1, u.2, r - u.1 - u.2)) r
  have hfr : Module.finrank ℝ (ℝ × ℝ) = 2 := by simp
  rw [hfr] at hs
  have hpt : (fun x : ℝ × ℝ => Q.indicator (Fb τ) ((r • x).1, (r • x).2, r - (r • x).1 - (r • x).2))
      = Δo.indicator (fun t => r ^ (-2 * τ) * Real.exp (-r) / Real.sqrt 3 * faceIntegrand τ t) := by
    funext t
    simp only [Prod.smul_fst, Prod.smul_snd, smul_eq_mul]
    by_cases ht : t ∈ Δo
    · rw [indicator_of_mem ht, indicator_of_mem, Fb_scale τ hr ht]
      · ring
      obtain ⟨ha, hb, hc⟩ := ht
      refine ⟨?_, ?_, ?_⟩ <;> simp only [mem_Ioi]
      · positivity
      · positivity
      · nlinarith
    · rw [indicator_of_notMem ht, indicator_of_notMem]
      rintro ⟨ha, hb, hc⟩
      simp only [mem_Ioi] at ha hb hc
      apply ht
      refine ⟨pos_of_mul_pos_right ha hr.le, pos_of_mul_pos_right hb hr.le, ?_⟩
      have : 0 < r * (1 - t.1 - t.2) := by linarith
      exact pos_of_mul_pos_right this hr.le
  rw [hpt, integral_indicator Δo_meas, integral_const_mul, setIntegral_congr_set Δo_ae] at hs
  have hr2 : (0 : ℝ) < r ^ 2 := by positivity
  rw [abs_of_pos (inv_pos.mpr hr2), smul_eq_mul] at hs
  have e : r ^ (2 - 2 * τ) = r ^ 2 * r ^ (-2 * τ) := by
    rw [show 2 - 2 * τ = (2 : ℕ) + -2 * τ by push_cast; ring, Real.rpow_add hr,
      Real.rpow_natCast]
  rw [e]
  field_simp at hs ⊢
  linarith

lemma dirichlet {τ : ℝ} (hτ : τ < 1) :
    ∫ g in Q, Fb τ g = Real.Gamma (3 - 2 * τ) * (∫ t in Δ, faceIntegrand τ t) / Real.sqrt 3 := by
  set F := Q.indicator (Fb τ) with hF
  have hFm : Measurable F := (Fb_meas τ).indicator Q_meas
  have hFi : Integrable F := (integrable_indicator_iff Q_meas).mpr (Fb_int hτ)
  rw [← integral_indicator Q_meas]
  -- reassociate
  have hA := (measurePreserving_prodAssoc (volume : Measure ℝ) (volume : Measure ℝ)
    (volume : Measure ℝ))
  have e1 : ∫ g, F g = ∫ v : (ℝ × ℝ) × ℝ, F (v.1.1, v.1.2, v.2) ∂(volume.prod volume) := by
    rw [show (volume : Measure (ℝ × ℝ × ℝ)) = volume.prod (volume.prod volume) from rfl,
      ← hA.integral_comp' F]
    rfl
  have hF1i : Integrable (fun v : (ℝ × ℝ) × ℝ => F (v.1.1, v.1.2, v.2)) (volume.prod volume) :=
    hA.integrable_comp_of_integrable (g := F) hFi
  -- shear
  have hS : MeasurePreserving (fun w : (ℝ × ℝ) × ℝ => (w.1, w.2 - w.1.1 - w.1.2))
      (volume.prod volume) (volume.prod volume) := by
    refine MeasurePreserving.skew_product (MeasurePreserving.id _)
      (g := fun (u : ℝ × ℝ) (r : ℝ) => r - u.1 - u.2) (by fun_prop) (Eventually.of_forall fun u => ?_)
    have := map_sub_right_eq_self (volume : Measure ℝ) (u.1 + u.2)
    simpa [sub_sub] using this
  set G : (ℝ × ℝ) × ℝ → ℝ := fun v => F (v.1.1, v.1.2, v.2) with hG
  have hGm : Measurable G := hFm.comp (by fun_prop)
  have e2 : ∫ v, G v ∂(volume.prod volume) =
      ∫ w, G (w.1, w.2 - w.1.1 - w.1.2) ∂(volume.prod volume) := by
    conv_lhs => rw [← hS.map_eq]
    rw [integral_map hS.measurable.aemeasurable hGm.aestronglyMeasurable]
  have hHi : Integrable (fun w : (ℝ × ℝ) × ℝ => G (w.1, w.2 - w.1.1 - w.1.2))
      (volume.prod volume) := hS.integrable_comp_of_integrable hF1i
  rw [e1, e2, integral_prod_symm _ hHi]
  simp only [hG]
  have e3 : ∀ r : ℝ, ∫ u : ℝ × ℝ, F (u.1, u.2, r - u.1 - u.2) =
      (Ioi (0 : ℝ)).indicator (fun r => Real.exp (-r) * r ^ (3 - 2 * τ - 1) *
        ((∫ t in Δ, faceIntegrand τ t) / Real.sqrt 3)) r := by
    intro r
    by_cases hr : 0 < r
    · rw [indicator_of_mem (show r ∈ Ioi 0 from hr), hF, inner_r τ hr,
        show 3 - 2 * τ - 1 = 2 - 2 * τ by ring]
      ring
    · rw [indicator_of_notMem (show r ∉ Ioi 0 from hr), hF]
      simp only [Fb_zero (not_lt.mp hr), integral_zero]
  simp_rw [e3]
  rw [integral_indicator measurableSet_Ioi, integral_mul_const,
    ← Real.Gamma_eq_integral (by linarith)]
  ring

end Y05

theorem Y05.proof_3 : Y05.stmt_3 := by
  intro p hp
  set τ := 1 - 1 / p with hτdef
  have hτ1 : τ < 1 := by have : 0 < 1 / p := by positivity
                         linarith
  have hB : B τ = (Real.Gamma (1 - 2 * τ / 3) ^ 3)⁻¹ * ∫ g in Q, Fb τ g := rfl
  rw [hB, dirichlet hτ1]
  have g1 : 0 < Real.Gamma (1 - 2 * τ / 3) := Real.Gamma_pos_of_pos (by linarith)
  have g2 : 0 < Real.Gamma (3 - 2 * τ) := Real.Gamma_pos_of_pos (by linarith)
  have s3 : 0 < Real.sqrt 3 := by positivity
  unfold S
  rw [← hτdef]
  generalize Real.Gamma (1 - 2 * τ / 3) = a at *
  generalize Real.Gamma (3 - 2 * τ) = b at *
  field_simp

/-! ## stmt_2: the graph-area formula -/

namespace Y05

lemma dz1 {p x₁ x₂ : ℝ} (hp : 0 < p) (h1 : 0 < x₁) (hW : 0 < 1 - x₁ ^ p - x₂ ^ p) :
    HasDerivAt (fun u => zgraph p u x₂)
      (-(x₁ ^ (p - 1)) * (1 - x₁ ^ p - x₂ ^ p) ^ (1 / p - 1)) x₁ := by
  have h := ((Real.hasDerivAt_rpow_const (p := p) (Or.inl h1.ne')).const_sub 1).sub_const
    (x₂ ^ p)
  have h2 := h.rpow_const (p := 1 / p) (Or.inl hW.ne')
  unfold zgraph
  convert h2 using 1
  field_simp

lemma dz2 {p x₁ x₂ : ℝ} (hp : 0 < p) (h2 : 0 < x₂) (hW : 0 < 1 - x₁ ^ p - x₂ ^ p) :
    HasDerivAt (fun v => zgraph p x₁ v)
      (-(x₂ ^ (p - 1)) * (1 - x₁ ^ p - x₂ ^ p) ^ (1 / p - 1)) x₂ := by
  have h := (Real.hasDerivAt_rpow_const (p := p) (Or.inl h2.ne')).const_sub (1 - x₁ ^ p)
  have h3 := h.rpow_const (p := 1 / p) (Or.inl (by rw [sub_sub] at hW ⊢; linarith))
  unfold zgraph
  simp only [sub_sub] at h3 ⊢
  convert h3 using 1
  field_simp

/-- the coordinate map `t ↦ (t₁^{1/p}, t₂^{1/p})` -/
noncomputable def φ (p : ℝ) : ℝ × ℝ → ℝ × ℝ :=
  Prod.map (fun x : ℝ => x ^ (1 / p)) (fun x : ℝ => x ^ (1 / p))

noncomputable def φ' (p : ℝ) (t : ℝ × ℝ) : ℝ × ℝ →L[ℝ] ℝ × ℝ :=
  (ContinuousLinearMap.smulRight (1 : ℝ →L[ℝ] ℝ) (1 / p * t.1 ^ (1 / p - 1))).prodMap
    (ContinuousLinearMap.smulRight (1 : ℝ →L[ℝ] ℝ) (1 / p * t.2 ^ (1 / p - 1)))

lemma φ'_det (p : ℝ) (t : ℝ × ℝ) :
    (φ' p t).det = (1 / p * t.1 ^ (1 / p - 1)) * (1 / p * t.2 ^ (1 / p - 1)) := by
  unfold φ'
  rw [ContinuousLinearMap.det, ContinuousLinearMap.coe_prodMap, LinearMap.det_prodMap]
  simp

/-- the half-open triangle `t₁, t₂ > 0`, `t₁ + t₂ ≤ 1` -/
def Δh : Set (ℝ × ℝ) := {t | 0 < t.1 ∧ 0 < t.2 ∧ t.1 + t.2 ≤ 1}

lemma Δh_meas : MeasurableSet Δh := by
  unfold Δh
  simp only [Set.ofPred_and]
  exact (measurableSet_lt measurable_const measurable_fst).inter
    ((measurableSet_lt measurable_const measurable_snd).inter
      (measurableSet_le (measurable_fst.add measurable_snd) measurable_const))

lemma Δh_ae : Δh =ᵐ[volume] Δ := by
  rw [ae_eq_set]
  have h0 := (ae_eq_set.mp Δo_ae).2
  constructor
  · have : Δh \ Δ = ∅ := by
      ext t
      simp only [Δh, Δ, Set.mem_sdiff, mem_ofPred_eq, mem_empty_iff_false, iff_false]
      rintro ⟨⟨a, b, c⟩, h⟩
      exact h ⟨a.le, b.le, c⟩
    rw [this]; simp
  · refine measure_mono_null (fun t ht => ?_) h0
    obtain ⟨h1, h2⟩ := ht
    refine ⟨h1, fun h3 => h2 ?_⟩
    obtain ⟨a, b, c⟩ := h3
    exact ⟨a, b, by linarith⟩

lemma Δo_Δh : Δo =ᵐ[volume] Δh := Δo_ae.trans Δh_ae.symm

lemma φ_deriv {p : ℝ} (t : ℝ × ℝ) (ht : t ∈ Δh) : HasFDerivWithinAt (φ p) (φ' p t) Δh t := by
  obtain ⟨a, b, _⟩ := ht
  have h1 := (Real.hasDerivAt_rpow_const (x := t.1) (p := 1 / p) (Or.inl a.ne')).hasFDerivAt
  have h2 := (Real.hasDerivAt_rpow_const (x := t.2) (p := 1 / p) (Or.inl b.ne')).hasFDerivAt
  exact (h1.prodMap t h2).hasFDerivWithinAt

lemma φ_inj {p : ℝ} (hp : 0 < p) : InjOn (φ p) Δh := by
  rintro ⟨a, b⟩ ⟨ha, hb, _⟩ ⟨c, d⟩ ⟨hc, hd, _⟩ h
  simp only [φ, Prod.map, Prod.mk.injEq] at h ha hb hc hd ⊢
  have hp' : (1 / p) ≠ 0 := by positivity
  constructor
  · rw [← Real.rpow_rpow_inv ha.le hp', ← Real.rpow_rpow_inv hc.le hp', h.1]
  · rw [← Real.rpow_rpow_inv hb.le hp', ← Real.rpow_rpow_inv hd.le hp', h.2]

/-- `D' = {x₁, x₂ > 0, x₁^p + x₂^p ≤ 1}` -/
def D' (p : ℝ) : Set (ℝ × ℝ) := {x | 0 < x.1 ∧ 0 < x.2 ∧ x.1 ^ p + x.2 ^ p ≤ 1}

lemma pw_inv {p x : ℝ} (hp : 0 < p) (hx : 0 ≤ x) : (x ^ (1 / p)) ^ p = x := by
  rw [← Real.rpow_mul hx, one_div, inv_mul_cancel₀ hp.ne', Real.rpow_one]

lemma φ_image {p : ℝ} (hp : 0 < p) : φ p '' Δh = D' p := by
  ext ⟨x, y⟩
  simp only [mem_image, φ, Prod.map, Prod.mk.injEq, D', Δh, mem_ofPred_eq, Prod.exists]
  constructor
  · rintro ⟨a, b, ⟨ha, hb, hab⟩, rfl, rfl⟩
    refine ⟨Real.rpow_pos_of_pos ha _, Real.rpow_pos_of_pos hb _, ?_⟩
    rw [pw_inv hp ha.le, pw_inv hp hb.le]; exact hab
  · rintro ⟨hx, hy, hxy⟩
    refine ⟨x ^ p, y ^ p, ⟨Real.rpow_pos_of_pos hx _, Real.rpow_pos_of_pos hy _, hxy⟩, ?_, ?_⟩
    · rw [← Real.rpow_mul hx.le, mul_one_div_cancel hp.ne', Real.rpow_one]
    · rw [← Real.rpow_mul hy.le, mul_one_div_cancel hp.ne', Real.rpow_one]

lemma D'_ae (p : ℝ) : D' p =ᵐ[volume] D p := by
  rw [ae_eq_set]
  constructor
  · have : D' p \ D p = ∅ := by
      ext t
      simp only [D', D, Set.mem_sdiff, mem_ofPred_eq, mem_empty_iff_false, iff_false]
      rintro ⟨⟨a, b, c⟩, h⟩
      exact h ⟨a.le, b.le, c⟩
    rw [this]; simp
  · have hsub : D p \ D' p ⊆ ({0} ×ˢ univ) ∪ (univ ×ˢ {0}) := by
      rintro ⟨a, b⟩ ⟨⟨ha, hb, hab⟩, hn⟩
      simp only [D', mem_ofPred_eq, not_and] at hn
      simp only [mem_union, mem_prod, mem_singleton_iff, mem_univ, and_true, true_and]
      by_contra hc
      push Not at hc
      exact hn (lt_of_le_of_ne ha (Ne.symm hc.1)) (lt_of_le_of_ne hb (Ne.symm hc.2)) hab
    refine measure_mono_null hsub (measure_union_null ?_ ?_)
    · rw [Measure.volume_eq_prod, Measure.prod_prod]; simp
    · rw [Measure.volume_eq_prod, Measure.prod_prod]; simp

lemma area_pt {p : ℝ} (hp : 0 < p) {t : ℝ × ℝ} (ht : t ∈ Δo) :
    |(φ' p t).det| • areaElt p (φ p t) = (1 / p) ^ 2 * faceIntegrand (1 - 1 / p) t := by
  obtain ⟨ha, hb, hc⟩ := ht
  rw [φ'_det, smul_eq_mul]
  set τ := 1 - 1 / p with hτ
  have hx1 : (t.1 ^ (1 / p)) ^ p = t.1 := pw_inv hp ha.le
  have hx2 : (t.2 ^ (1 / p)) ^ p = t.2 := pw_inv hp hb.le
  have hW : 0 < 1 - (t.1 ^ (1 / p)) ^ p - (t.2 ^ (1 / p)) ^ p := by rw [hx1, hx2]; exact hc
  have d1 := (dz1 hp (Real.rpow_pos_of_pos ha (1 / p)) hW).deriv
  have d2 := (dz2 hp (Real.rpow_pos_of_pos hb (1 / p)) hW).deriv
  unfold areaElt
  simp only [φ, Prod.map]
  rw [d1, d2, hx1, hx2]
  have e1 : ∀ x : ℝ, 0 < x → (x ^ (1 / p)) ^ (p - 1) = x ^ τ := fun x hx => by
    rw [← Real.rpow_mul hx.le, hτ]; congr 1; field_simp
  rw [e1 _ ha, e1 _ hb, show 1 / p - 1 = -τ by rw [hτ]; ring]
  set c := 1 - t.1 - t.2 with hcdef
  unfold faceIntegrand
  rw [← hcdef]
  have sq : ∀ x : ℝ, 0 < x → x ^ (2 * τ) = (x ^ τ) ^ 2 := fun x hx => by
    rw [← Real.rpow_natCast, ← Real.rpow_mul hx.le]; congr 1; push_cast; ring
  have ki : ∀ x : ℝ, 0 < x → x ^ (-τ) * x ^ τ = 1 := fun x hx => by
    rw [← Real.rpow_add hx]; simp
  rw [sq _ ha, sq _ hb, sq _ hc, Real.mul_rpow (by positivity) hc.le, Real.mul_rpow ha.le hb.le]
  have k1 := ki _ ha; have k2 := ki _ hb; have k3 := ki _ hc
  have na := Real.rpow_pos_of_pos ha (-τ); have nb := Real.rpow_pos_of_pos hb (-τ)
  have nc := Real.rpow_pos_of_pos hc (-τ)
  have pa := Real.rpow_pos_of_pos ha τ; have pb := Real.rpow_pos_of_pos hb τ
  have pc := Real.rpow_pos_of_pos hc τ
  have hj : 1 / p * t.1 ^ (-τ) * (1 / p * t.2 ^ (-τ)) =
      (1 / p) ^ 2 * (t.1 ^ (-τ) * t.2 ^ (-τ)) := by ring
  rw [hj, abs_of_pos (by positivity)]
  generalize t.1 ^ τ = A at *
  generalize t.2 ^ τ = B' at *
  generalize c ^ τ = C at *
  generalize t.1 ^ (-τ) = A' at *
  generalize t.2 ^ (-τ) = Bn at *
  generalize c ^ (-τ) = C' at *
  have hrad : 1 + (-A * C') ^ 2 + (-B' * C') ^ 2 = C' ^ 2 * (A ^ 2 + B' ^ 2 + C ^ 2) := by
    linear_combination (-1 : ℝ) * (C' * C + 1) * k3
  rw [hrad, Real.sqrt_mul (by positivity), Real.sqrt_sq nc.le]
  ring

end Y05

theorem Y05.proof_2 : Y05.stmt_2 := by
  intro p hp
  have hint : IntegrableOn (fun t => (1 / p) ^ 2 * faceIntegrand (1 - 1 / p) t) Δo :=
    ((Y05.proof_1 p hp).mono_set (show Δo ⊆ Δ from fun t ht =>
      ⟨ht.1.le, ht.2.1.le, by linarith [ht.2.2]⟩)).const_mul _
  have hcongr : ∫ t in Δo, |(φ' p t).det| • areaElt p (φ p t) =
      ∫ t in Δo, (1 / p) ^ 2 * faceIntegrand (1 - 1 / p) t :=
    setIntegral_congr_fun Δo_meas (fun t ht => area_pt hp ht)
  have hint2 : IntegrableOn (fun t => |(φ' p t).det| • areaElt p (φ p t)) Δh :=
    (hint.congr_fun (fun t ht => (area_pt hp ht).symm) Δo_meas).congr_set_ae Δo_Δh.symm
  have hI := integrableOn_image_iff_integrableOn_abs_det_fderiv_smul (μ := volume) Δh_meas
    (fun t ht => φ_deriv (p := p) t ht) (φ_inj hp) (areaElt p)
  have hV := integral_image_eq_integral_abs_det_fderiv_smul (μ := volume) Δh_meas
    (fun t ht => φ_deriv (p := p) t ht) (φ_inj hp) (areaElt p)
  rw [φ_image hp] at hI hV
  refine ⟨(hI.mpr hint2).congr_set_ae (D'_ae p).symm, ?_⟩
  rw [← setIntegral_congr_set (D'_ae p), hV, ← setIntegral_congr_set Δo_Δh, hcongr,
    integral_const_mul, setIntegral_congr_set Δo_ae]
  unfold S
  ring

/-! ## stmt_7: the `J = 2` Taylor bound -/

namespace Y05

lemma sqrt_tay {x : ℝ} (hx : 0 ≤ x) :
    |Real.sqrt (1 + x) - (1 + x / 2 - x ^ 2 / 8)| ≤ x ^ 3 := by
  have hs2 : Real.sqrt (1 + x) ^ 2 = 1 + x := Real.sq_sqrt (by linarith)
  have hs1 : 1 ≤ Real.sqrt (1 + x) := by
    rw [Real.le_sqrt (by norm_num) (by linarith)]; linarith
  generalize Real.sqrt (1 + x) = s at *
  have hx' : x = s ^ 2 - 1 := by linarith
  rw [hx']
  have h1 : 0 ≤ s - 1 := by linarith
  have e : s - (1 + (s ^ 2 - 1) / 2 - (s ^ 2 - 1) ^ 2 / 8) = (s - 1) ^ 3 * ((s + 3) / 8) := by
    ring
  rw [e, abs_of_nonneg (mul_nonneg (pow_nonneg h1 3) (by linarith)),
    show (s ^ 2 - 1) ^ 3 = (s - 1) ^ 3 * (s + 1) ^ 3 by ring]
  have h3 : (s + 3) / 8 ≤ (s + 1) ^ 3 := by nlinarith
  exact mul_le_mul_of_nonneg_left h3 (pow_nonneg h1 3)

lemma exp_sub_le' {a b : ℝ} (h : b ≤ a) : |Real.exp a - Real.exp b| ≤ |a - b| * (Real.exp a + Real.exp b) := by
  have k := Real.add_one_le_exp (b - a)
  rw [Real.exp_sub] at k
  have ea := Real.exp_pos a; have eb := Real.exp_pos b
  have k2 : (b - a + 1) * Real.exp a ≤ Real.exp b := by rwa [le_div_iff₀ ea] at k
  have hab : Real.exp b ≤ Real.exp a := Real.exp_le_exp.mpr h
  rw [abs_of_nonneg (by linarith), abs_of_nonneg (by linarith)]
  nlinarith

lemma exp_sub_le (a b : ℝ) : |Real.exp a - Real.exp b| ≤ |a - b| * (Real.exp a + Real.exp b) := by
  rcases le_total b a with h | h
  · exact exp_sub_le' h
  · rw [abs_sub_comm, abs_sub_comm a b, add_comm]; exact exp_sub_le' h

lemma alg7 {q1 q2 q3 E D e : ℝ} (h1 : 0 < q1) (h2 : 0 < q2) (h3 : 0 < q3) (he : 0 < e)
    (_hD : 0 ≤ D) (u1 : q1 ≤ E) (u2 : q2 ≤ E) (u3 : q3 ≤ E)
    (v1 : q1⁻¹ ≤ E) (v2 : q2⁻¹ ≤ E) (v3 : q3⁻¹ ≤ E)
    (d12 : |q1 - q2| ≤ D) (d23 : |q2 - q3| ≤ D) (d13 : |q1 - q3| ≤ D) :
    |(q1 * q2 * q3)⁻¹ * e * (Real.sqrt ((q1 ^ 3 + q2 ^ 3 + q3 ^ 3) / (3 * (q1 * q2 * q3)))
      - (1 + ((q1 ^ 3 + q2 ^ 3 + q3 ^ 3) / (3 * (q1 * q2 * q3)) - 1) / 2
        - ((q1 ^ 3 + q2 ^ 3 + q3 ^ 3) / (3 * (q1 * q2 * q3)) - 1) ^ 2 / 8))|
      ≤ E ^ 3 * e * ((3 / 2) * D ^ 2 * E ^ 4) ^ 3 := by
  have hE0 : 0 ≤ E := le_trans h1.le u1
  generalize hY : (q1 ^ 3 + q2 ^ 3 + q3 ^ 3) / (3 * (q1 * q2 * q3)) = Y
  have hX : Y - 1 = (q1 + q2 + q3) * ((q1 - q2) ^ 2 + (q2 - q3) ^ 2 + (q1 - q3) ^ 2) *
      (q1 * q2 * q3)⁻¹ / 6 := by
    rw [← hY]; field_simp; ring
  have hP : (q1 * q2 * q3)⁻¹ ≤ E ^ 3 := by
    rw [mul_inv, mul_inv, pow_three]
    have := mul_le_mul v1 v2 (by positivity) hE0
    calc q1⁻¹ * q2⁻¹ * q3⁻¹ ≤ (E * E) * E := mul_le_mul this v3 (by positivity) (by positivity)
      _ = E * (E * E) := by ring
  have sq : ∀ d : ℝ, |d| ≤ D → d ^ 2 ≤ D ^ 2 := fun d hd => by
    rw [← sq_abs]; exact pow_le_pow_left₀ (abs_nonneg d) hd 2
  have hS : (q1 - q2) ^ 2 + (q2 - q3) ^ 2 + (q1 - q3) ^ 2 ≤ 3 * D ^ 2 := by
    have := sq _ d12; have := sq _ d23; have := sq _ d13; linarith
  have hX0 : 0 ≤ Y - 1 := by rw [hX]; positivity
  have hXb : Y - 1 ≤ (3 / 2) * D ^ 2 * E ^ 4 := by
    rw [hX]
    have hq : q1 + q2 + q3 ≤ 3 * E := by linarith
    calc (q1 + q2 + q3) * ((q1 - q2) ^ 2 + (q2 - q3) ^ 2 + (q1 - q3) ^ 2) * (q1 * q2 * q3)⁻¹ / 6
        ≤ (3 * E) * (3 * D ^ 2) * E ^ 3 / 6 := by gcongr
      _ = _ := by ring
  have key : |Real.sqrt Y - (1 + (Y - 1) / 2 - (Y - 1) ^ 2 / 8)| ≤ (Y - 1) ^ 3 := by
    have := sqrt_tay hX0; rwa [show 1 + (Y - 1) = Y by ring] at this
  have hw : 0 < (q1 * q2 * q3)⁻¹ * e := by positivity
  rw [abs_mul, abs_of_pos hw]
  calc (q1 * q2 * q3)⁻¹ * e * |Real.sqrt Y - (1 + (Y - 1) / 2 - (Y - 1) ^ 2 / 8)|
      ≤ (q1 * q2 * q3)⁻¹ * e * (Y - 1) ^ 3 := mul_le_mul_of_nonneg_left key hw.le
    _ ≤ E ^ 3 * e * ((3 / 2) * D ^ 2 * E ^ 4) ^ 3 := by
      apply mul_le_mul (mul_le_mul_of_nonneg_right hP he.le)
        (pow_le_pow_left₀ hX0 hXb 3) (pow_nonneg hX0 3) (by positivity)

lemma M6 {M : ℝ} (hM : 0 ≤ M) : M ^ 6 ≤ 12 ^ 6 * Real.exp (M / 2) := by
  have := Real.add_one_le_exp (M / 12)
  have h : M / 12 ≤ Real.exp (M / 12) := by linarith
  have h6 := pow_le_pow_left₀ (by positivity) h 6
  rw [← Real.exp_nat_mul] at h6
  calc M ^ 6 = 12 ^ 6 * (M / 12) ^ 6 := by ring
    _ ≤ 12 ^ 6 * Real.exp (((6 : ℕ) : ℝ) * (M / 12)) := by gcongr
    _ = 12 ^ 6 * Real.exp (M / 2) := by congr 2; push_cast; ring

lemma ekap {x κ : ℝ} (hx : 0 < x) : Real.exp (κ * |Real.log x|) ≤ x ^ κ + x ^ (-κ) := by
  have a := Real.rpow_pos_of_pos hx κ
  have b := Real.rpow_pos_of_pos hx (-κ)
  rcases abs_cases (Real.log x) with ⟨h, _⟩ | ⟨h, _⟩
  · rw [h, Real.rpow_def_of_pos hx κ, mul_comm]; linarith
  · rw [h, Real.rpow_def_of_pos hx (-κ), show κ * -Real.log x = Real.log x * -κ by ring]
    linarith

/-- the one-variable dominating function -/
noncomputable def hh (x : ℝ) : ℝ := Real.exp (-x) * x ^ (7 / 10 : ℝ) + Real.exp (-x) * x ^ (-(7 / 10) : ℝ)

lemma hh_int : IntegrableOn hh (Ioi 0) := by
  have a := (gamma1 (show (-1 : ℝ) < 7 / 10 by norm_num)).1
  have b := (gamma1 (show (-1 : ℝ) < -(7 / 10) by norm_num)).1
  exact a.add b

lemma hh3_int : IntegrableOn (fun g : ℝ × ℝ × ℝ => hh g.1 * (hh g.2.1 * hh g.2.2)) Q := by
  have iyz : IntegrableOn (fun h : ℝ × ℝ => hh h.1 * hh h.2) (Ioi (0 : ℝ) ×ˢ Ioi (0 : ℝ))
      (volume.prod volume) := by
    rw [IntegrableOn, ← Measure.prod_restrict]
    exact Integrable.mul_prod hh_int hh_int
  show IntegrableOn _ _ (volume.prod volume)
  rw [IntegrableOn, ← Measure.prod_restrict]
  exact Integrable.mul_prod (f := hh) (g := fun h : ℝ × ℝ => hh h.1 * hh h.2) hh_int iyz

lemma rp_q {x : ℝ} (hx : 0 < x) (y : ℝ) : x ^ y = Real.exp (y * Real.log x) := by
  rw [Real.rpow_def_of_pos hx, mul_comm]

lemma pt7 {τ : ℝ} (hτ : |τ| ≤ 1 / 70) {g : ℝ × ℝ × ℝ} (hg : g ∈ Q) :
    |(g.1 * g.2.1 * g.2.2) ^ (-2 * τ / 3) * Real.exp (-(g.1 + g.2.1 + g.2.2)) *
      (Real.sqrt (Ymean τ g) - (1 + (Ymean τ g - 1) / 2 - (Ymean τ g - 1) ^ 2 / 8))|
      ≤ 216 * 12 ^ 6 * τ ^ 6 * (hh g.1 * (hh g.2.1 * hh g.2.2)) := by
  obtain ⟨ha, hb, hc⟩ := hg
  simp only [mem_Ioi] at ha hb hc
  set c := 2 * τ / 3 with hcdef
  have hY : Ymean τ g = ((Real.exp (c * Real.log g.1)) ^ 3 + (Real.exp (c * Real.log g.2.1)) ^ 3
      + (Real.exp (c * Real.log g.2.2)) ^ 3) / (3 * (Real.exp (c * Real.log g.1) *
        Real.exp (c * Real.log g.2.1) * Real.exp (c * Real.log g.2.2))) := by
    unfold Ymean
    rw [Real.mul_rpow (mul_pos ha hb).le hc.le, Real.mul_rpow ha.le hb.le, ← hcdef,
      rp_q ha, rp_q hb, rp_q hc, rp_q ha, rp_q hb, rp_q hc]
    have k : ∀ x : ℝ, Real.exp (2 * τ * x) = Real.exp (c * x) ^ 3 := fun x => by
      rw [← Real.exp_nat_mul]; congr 1; push_cast; rw [hcdef]; ring
    rw [k, k, k]
  have hw : (g.1 * g.2.1 * g.2.2) ^ (-2 * τ / 3) = (Real.exp (c * Real.log g.1) *
        Real.exp (c * Real.log g.2.1) * Real.exp (c * Real.log g.2.2))⁻¹ := by
    rw [Real.mul_rpow (mul_pos ha hb).le hc.le, Real.mul_rpow ha.le hb.le,
      rp_q ha, rp_q hb, rp_q hc, mul_inv, mul_inv, ← Real.exp_neg, ← Real.exp_neg,
      ← Real.exp_neg, hcdef]
    congr 3 <;> ring
  rw [hY, hw]
  set l1 := Real.log g.1
  set l2 := Real.log g.2.1
  set l3 := Real.log g.2.2
  set M := |l1| + |l2| + |l3| with hM
  have hM0 : 0 ≤ M := by positivity
  set E := Real.exp (|c| * M) with hE
  have hc0 : 0 ≤ |c| := abs_nonneg c
  have up : ∀ l : ℝ, |l| ≤ M → Real.exp (c * l) ≤ E := fun l hl => by
    apply Real.exp_le_exp.mpr
    calc c * l ≤ |c * l| := le_abs_self _
      _ = |c| * |l| := abs_mul c l
      _ ≤ |c| * M := mul_le_mul_of_nonneg_left hl hc0
  have upi : ∀ l : ℝ, |l| ≤ M → (Real.exp (c * l))⁻¹ ≤ E := fun l hl => by
    rw [← Real.exp_neg, ← mul_neg]; apply up; rwa [abs_neg]
  have b1 : |l1| ≤ M := by have := abs_nonneg l2; have := abs_nonneg l3; linarith
  have b2 : |l2| ≤ M := by have := abs_nonneg l1; have := abs_nonneg l3; linarith
  have b3 : |l3| ≤ M := by have := abs_nonneg l2; have := abs_nonneg l1; linarith
  have dd : ∀ l l' : ℝ, |l| ≤ M → |l'| ≤ M → |l| + |l'| ≤ M →
      |Real.exp (c * l) - Real.exp (c * l')| ≤ 2 * |c| * M * E := fun l l' hl hl' hll => by
    calc |Real.exp (c * l) - Real.exp (c * l')|
        ≤ |c * l - c * l'| * (Real.exp (c * l) + Real.exp (c * l')) := exp_sub_le _ _
      _ ≤ (|c| * M) * (E + E) := by
          apply mul_le_mul _ (add_le_add (up l hl) (up l' hl')) (by positivity) (by positivity)
          rw [← mul_sub, abs_mul]
          exact mul_le_mul_of_nonneg_left ((abs_sub _ _).trans hll) hc0
      _ = 2 * |c| * M * E := by ring
  have hall := alg7 (Real.exp_pos (c * l1)) (Real.exp_pos (c * l2)) (Real.exp_pos (c * l3))
    (Real.exp_pos (-(g.1 + g.2.1 + g.2.2))) (by positivity : (0 : ℝ) ≤ 2 * |c| * M * E)
    (up l1 b1) (up l2 b2) (up l3 b3) (upi l1 b1) (upi l2 b2) (upi l3 b3)
    (dd l1 l2 b1 b2 (by have := abs_nonneg l3; linarith))
    (dd l2 l3 b2 b3 (by have := abs_nonneg l1; linarith))
    (dd l1 l3 b1 b3 (by have := abs_nonneg l2; linarith))
  refine hall.trans ?_
  -- now pure estimates
  have hEpow : E ^ 3 * ((3 / 2) * (2 * |c| * M * E) ^ 2 * E ^ 4) ^ 3 =
      216 * |c| ^ 6 * M ^ 6 * Real.exp (21 * |c| * M) := by
    rw [hE, show 21 * |c| * M = ((21 : ℕ) : ℝ) * (|c| * M) by push_cast; ring, Real.exp_nat_mul]
    ring
  have hc6 : |c| ^ 6 ≤ τ ^ 6 := by
    rw [hcdef, abs_div, abs_mul, show |(2 : ℝ)| = 2 by norm_num, show |(3 : ℝ)| = 3 by norm_num,
      div_pow, mul_pow, show |τ| ^ 6 = τ ^ 6 from Even.pow_abs ⟨3, rfl⟩ τ]
    have := pow_nonneg (abs_nonneg τ) 6
    rw [show |τ| ^ 6 = τ ^ 6 from Even.pow_abs ⟨3, rfl⟩ τ] at this
    linarith
  have h21 : 21 * |c| * M ≤ M / 5 := by
    have : |c| = 2 * |τ| / 3 := by
      rw [hcdef, abs_div, abs_mul]; norm_num
    rw [this]; nlinarith [abs_nonneg τ]
  have hprod : Real.exp (-(g.1 + g.2.1 + g.2.2)) * Real.exp (7 / 10 * M) ≤
      hh g.1 * (hh g.2.1 * hh g.2.2) := by
    have e1 : Real.exp (-(g.1 + g.2.1 + g.2.2)) * Real.exp (7 / 10 * M) =
        (Real.exp (-g.1) * Real.exp (7 / 10 * |l1|)) * ((Real.exp (-g.2.1) *
          Real.exp (7 / 10 * |l2|)) * (Real.exp (-g.2.2) * Real.exp (7 / 10 * |l3|))) := by
      simp only [← Real.exp_add]; congr 1; rw [hM]; ring
    rw [e1]
    unfold hh
    simp only [← mul_add]
    have k1 := ekap (κ := 7 / 10) ha
    have k2 := ekap (κ := 7 / 10) hb
    have k3 := ekap (κ := 7 / 10) hc
    gcongr
  have hm6 := M6 hM0
  have hτ6 : 0 ≤ τ ^ 6 := by positivity
  set e := Real.exp (-(g.1 + g.2.1 + g.2.2))
  have he : 0 < e := Real.exp_pos _
  calc E ^ 3 * e * ((3 / 2) * (2 * |c| * M * E) ^ 2 * E ^ 4) ^ 3
      = e * (216 * |c| ^ 6 * M ^ 6 * Real.exp (21 * |c| * M)) := by rw [← hEpow]; ring
    _ ≤ e * (216 * τ ^ 6 * (12 ^ 6 * Real.exp (M / 2)) * Real.exp (M / 5)) := by
        gcongr
    _ = 216 * 12 ^ 6 * τ ^ 6 * (e * Real.exp (7 / 10 * M)) := by
        rw [show 7 / 10 * M = M / 2 + M / 5 by ring, Real.exp_add]; ring
    _ ≤ 216 * 12 ^ 6 * τ ^ 6 * (hh g.1 * (hh g.2.1 * hh g.2.2)) := by gcongr

end Y05

namespace Y05

/-- the (unnormalized) Gamma weight -/
noncomputable def wt (τ : ℝ) (g : ℝ × ℝ × ℝ) : ℝ :=
  (g.1 * g.2.1 * g.2.2) ^ (-2 * τ / 3) * Real.exp (-(g.1 + g.2.1 + g.2.2))

lemma gE_wt (τ : ℝ) (F : ℝ × ℝ × ℝ → ℝ) :
    gammaE τ F = (Real.Gamma (1 - 2 * τ / 3) ^ 3)⁻¹ * ∫ g in Q, wt τ g * F g := rfl

lemma wt_mono (τ : ℝ) {g : ℝ × ℝ × ℝ} (hg : g ∈ Q) :
    mono (-2 * τ / 3) (-2 * τ / 3) (-2 * τ / 3) g = wt τ g := by
  obtain ⟨ha, hb, hc⟩ := hg
  simp only [mem_Ioi] at ha hb hc
  unfold mono wt
  rw [Real.mul_rpow (mul_pos ha hb).le hc.le, Real.mul_rpow ha.le hb.le]
  ring

end Y05

theorem Y05.proof_7 : Y05.stmt_7 := by
  have hcont : ContinuousAt Real.Gamma 1 :=
    (Real.differentiableAt_Gamma (fun m h => by
      have := m.cast_nonneg (α := ℝ); linarith)).continuousAt
  obtain ⟨ε, hε, hεb⟩ := Metric.continuousAt_iff.mp hcont (1 / 2) (by norm_num)
  set K := ∫ g in Q, hh g.1 * (hh g.2.1 * hh g.2.2)
  refine ⟨8 * (216 * 12 ^ 6) * K, min (1 / 70) (ε / 2), lt_min (by norm_num) (by linarith),
    fun τ hτ => ?_⟩
  have hτ1 : |τ| ≤ 1 / 70 := hτ.trans (min_le_left _ _)
  have hτ2 : |τ| ≤ ε / 2 := hτ.trans (min_le_right _ _)
  have hG : 1 / 2 ≤ Real.Gamma (1 - 2 * τ / 3) := by
    have := hεb (x := 1 - 2 * τ / 3) (by
      rw [Real.dist_eq, show 1 - 2 * τ / 3 - 1 = -(2 * τ / 3) by ring, abs_neg, abs_div,
        abs_mul]
      norm_num; linarith)
    rw [Real.dist_eq, Real.Gamma_one] at this
    linarith [(abs_lt.mp this).1]
  obtain ⟨t1, t2⟩ := abs_le.mp hτ1
  have hτl : |τ| < 1 / 2 := by linarith
  have hτl' : |τ| < 3 / 4 := by linarith
  -- integrability
  have ha : (-1 : ℝ) < -2 * τ / 3 := by linarith
  obtain ⟨i0', v0'⟩ := mono_int ha ha ha
  have i0 : IntegrableOn (wt τ) Q := i0'.congr_fun (fun g hg => wt_mono τ hg) Q_meas
  have v0 : ∫ g in Q, wt τ g = Real.Gamma (1 - 2 * τ / 3) ^ 3 := by
    rw [← setIntegral_congr_fun Q_meas (fun g hg => wt_mono τ hg), v0',
      show 1 + -2 * τ / 3 = 1 - 2 * τ / 3 by ring]
    ring
  have i1 : IntegrableOn (fun g => wt τ g * Real.sqrt (Ymean τ g)) Q := Fb_int (by linarith)
  have iY : IntegrableOn (fun g => wt τ g * Ymean τ g) Q := by
    have h1 : -1 < 2 * τ / 3 := by linarith
    have h2 : -1 < -(4 * τ / 3) := by linarith
    have := (((mono_int h1 h2 h2).1.add (mono_int h2 h1 h2).1).add
      (mono_int h2 h2 h1).1).const_mul (1 / 3)
    exact IntegrableOn.congr_fun this (fun g hg => (pt5 τ hg).symm) Q_meas
  have iY2 : IntegrableOn (fun g => wt τ g * Ymean τ g ^ 2) Q := by
    have h1 : -1 < 2 * τ := by linarith
    have h2 : -1 < -(2 * τ) := by linarith
    have h0 : (-1 : ℝ) < 0 := by norm_num
    have := ((((mono_int h1 h2 h2).1.add (mono_int h2 h1 h2).1).add
      (mono_int h2 h2 h1).1).add ((((mono_int h0 h0 h2).1.add (mono_int h0 h2 h0).1).add
      (mono_int h2 h0 h0).1).const_mul 2)).const_mul (1 / 9)
    exact IntegrableOn.congr_fun this (fun g hg => (pt6 τ hg).symm) Q_meas
  have i2 : IntegrableOn (fun g => wt τ g * (Ymean τ g - 1)) Q :=
    (iY.sub i0).congr_fun (fun g _ => by simp only [Pi.sub_apply]; ring) Q_meas
  have i3 : IntegrableOn (fun g => wt τ g * (Ymean τ g - 1) ^ 2) Q :=
    ((iY2.sub (iY.const_mul 2)).add i0).congr_fun
      (fun g _ => by simp only [Pi.sub_apply, Pi.add_apply]; ring) Q_meas
  -- linearity
  set f : ℝ × ℝ × ℝ → ℝ := fun g => wt τ g * (Real.sqrt (Ymean τ g) -
    (1 + (Ymean τ g - 1) / 2 - (Ymean τ g - 1) ^ 2 / 8)) with hfdef
  have hf : f = fun g => wt τ g * Real.sqrt (Ymean τ g) - wt τ g
      - (wt τ g * (Ymean τ g - 1)) / 2 + (wt τ g * (Ymean τ g - 1) ^ 2) / 8 := by
    funext g; rw [hfdef]; ring
  have hint : ∫ g in Q, f g = (∫ g in Q, wt τ g * Real.sqrt (Ymean τ g)) - (∫ g in Q, wt τ g)
      - (∫ g in Q, wt τ g * (Ymean τ g - 1)) / 2 + (∫ g in Q, wt τ g * (Ymean τ g - 1) ^ 2) / 8 := by
    have e1 := integral_add (μ := volume.restrict Q)
      (f := fun g => wt τ g * Real.sqrt (Ymean τ g) - wt τ g - wt τ g * (Ymean τ g - 1) / 2)
      (g := fun g => wt τ g * (Ymean τ g - 1) ^ 2 / 8) ((i1.sub i0).sub (i2.div_const 2))
      (i3.div_const 8)
    have e2 := integral_sub (μ := volume.restrict Q)
      (f := fun g => wt τ g * Real.sqrt (Ymean τ g) - wt τ g)
      (g := fun g => wt τ g * (Ymean τ g - 1) / 2) (i1.sub i0) (i2.div_const 2)
    have e3 := integral_sub (μ := volume.restrict Q)
      (f := fun g => wt τ g * Real.sqrt (Ymean τ g)) (g := fun g => wt τ g) i1 i0
    rw [hf]
    simp only at e1 e2 e3 ⊢
    rw [e1, e2, e3, integral_div, integral_div]
  have hGp : 0 < Real.Gamma (1 - 2 * τ / 3) := by linarith
  have hrep : B τ - (1 + gammaE τ (fun g => Ymean τ g - 1) / 2
      - gammaE τ (fun g => (Ymean τ g - 1) ^ 2) / 8) =
      (Real.Gamma (1 - 2 * τ / 3) ^ 3)⁻¹ * ∫ g in Q, f g := by
    rw [hint, B, gE_wt, gE_wt, gE_wt, v0]
    have hGn : Real.Gamma (1 - 2 * τ / 3) ≠ 0 := hGp.ne'
    generalize Real.Gamma (1 - 2 * τ / 3) = G at *
    field_simp
    ring
  rw [hrep, abs_mul, abs_of_pos (by positivity)]
  have hbd : |∫ g in Q, f g| ≤ ∫ g in Q, 216 * 12 ^ 6 * τ ^ 6 * (hh g.1 * (hh g.2.1 * hh g.2.2)) := by
    rw [← Real.norm_eq_abs]
    refine norm_integral_le_of_norm_le (hh3_int.const_mul _) ?_
    refine (ae_restrict_iff' Q_meas).mpr (Eventually.of_forall fun g hg => ?_)
    rw [Real.norm_eq_abs]
    exact pt7 hτ1 hg
  rw [integral_const_mul] at hbd
  have hinv : (Real.Gamma (1 - 2 * τ / 3) ^ 3)⁻¹ ≤ 8 := by
    rw [inv_le_comm₀ (by positivity) (by norm_num)]
    calc 8⁻¹ = (1 / 2 : ℝ) ^ 3 := by norm_num
      _ ≤ _ := pow_le_pow_left₀ (by norm_num) hG 3
  calc (Real.Gamma (1 - 2 * τ / 3) ^ 3)⁻¹ * |∫ g in Q, f g|
      ≤ 8 * (216 * 12 ^ 6 * τ ^ 6 * K) := mul_le_mul hinv hbd (abs_nonneg _) (by norm_num)
    _ = 8 * (216 * 12 ^ 6) * K * |τ| ^ 6 := by
      rw [Even.pow_abs ⟨3, rfl⟩ τ]; ring

/-! ## stmt_8: the `τ²`, `τ³` coefficients -/

namespace Y05

lemma z2_eq : X05.z 2 = Real.pi ^ 2 / 6 := by
  unfold X05.z
  rw [show ((2 : ℕ) : ℂ) = 2 by norm_num, riemannZeta_two]
  rw [show (Real.pi : ℂ) ^ 2 / 6 = ((Real.pi ^ 2 / 6 : ℝ) : ℂ) by push_cast; ring, Complex.ofReal_re]

lemma z_le {n : ℕ} (hn : 2 ≤ n) : X05.z n ≤ 2 := by
  have h1 := X05.z_hasSum (show 1 < n by omega)
  have h2 := X05.z_hasSum (show 1 < 2 by norm_num)
  have hle : X05.z n ≤ X05.z 2 := by
    refine hasSum_le (fun k => ?_) h1 h2
    have hk : (1 : ℝ) ≤ (k : ℝ) + 1 := by have := k.cast_nonneg (α := ℝ); linarith
    exact one_div_le_one_div_of_le (by positivity) (pow_le_pow_right₀ hk hn)
  rw [z2_eq] at hle
  have := Real.pi_lt_d2
  have hp := Real.pi_pos
  nlinarith

lemma z_nn {n : ℕ} (hn : 2 ≤ n) : 0 ≤ X05.z n :=
  (X05.z_hasSum (show 1 < n by omega)).nonneg (fun k => by positivity)

lemma tail2 {f : ℕ → ℝ} {L A r : ℝ} (hf : HasSum f L) (hr0 : 0 ≤ r) (hr : r ≤ 1 / 2)
    (hb : ∀ k, |f (k + 2)| ≤ A * r ^ k) : |L - f 0 - f 1| ≤ 2 * A := by
  have hs : HasSum (fun k => f (k + 2)) (L - (f 0 + f 1)) := by
    have := (hasSum_nat_add_iff' 2).mpr hf
    simpa [Finset.sum_range_succ] using this
  have hA : 0 ≤ A := by
    have := (abs_nonneg _).trans (hb 0); simpa using this
  have hg : HasSum (fun k => A * r ^ k) (A * (1 - r)⁻¹) :=
    (hasSum_geometric_of_lt_one hr0 (by linarith)).mul_left A
  have hup := hasSum_le (fun k => (le_abs_self _).trans (hb k)) hs hg
  have hlo := hasSum_le (fun k => neg_le_of_abs_le (hb k)) hg.neg hs
  have hinv : (1 - r)⁻¹ ≤ 2 := by
    rw [inv_le_comm₀ (by linarith) (by norm_num)]; linarith
  have : A * (1 - r)⁻¹ ≤ 2 * A := by nlinarith
  rw [show L - f 0 - f 1 = L - (f 0 + f 1) by ring, abs_le]
  constructor <;> linarith

lemma ser4 {τ : ℝ} (hτ : |τ| ≤ 1 / 4) :
    |Real.log (Real.Gamma (1 - 2 * τ / 3) ^ 3 / Real.Gamma (1 - 2 * τ))
      - (-(2 * Real.pi ^ 2 / 9) * τ ^ 2 - 64 * X05.z 3 / 27 * τ ^ 3)| ≤ 64 * τ ^ 4 := by
  have hs := Y05.proof_4 τ (by linarith)
  have ht := tail2 hs (by positivity : (0 : ℝ) ≤ 2 * |τ|) (by linarith) (A := 32 * τ ^ 4)
    (fun k => by
    have hz := z_le (show 2 ≤ k + 2 + 2 by omega)
    have hz0 := z_nn (show 2 ≤ k + 2 + 2 by omega)
    have h3 : 0 ≤ 1 - 1 / (3 : ℝ) ^ (k + 2 + 1) := by
      have : 1 / (3 : ℝ) ^ (k + 2 + 1) ≤ 1 := by
        rw [div_le_one (by positivity)]; exact one_le_pow₀ (by norm_num)
      linarith
    have h3' : 1 - 1 / (3 : ℝ) ^ (k + 2 + 1) ≤ 1 := by
      have : 0 ≤ 1 / (3 : ℝ) ^ (k + 2 + 1) := by positivity
      linarith
    have hn : (1 : ℝ) ≤ ((k + 2 + 2 : ℕ) : ℝ) := by push_cast; linarith [k.cast_nonneg (α := ℝ)]
    rw [abs_mul, abs_mul, abs_div, abs_neg, abs_of_nonneg (mul_nonneg (by positivity) h3),
      abs_of_pos (by positivity : (0 : ℝ) < ((k + 2 + 2 : ℕ) : ℝ)), abs_of_nonneg hz0, abs_pow]
    calc 2 ^ (k + 2 + 2) * (1 - 1 / 3 ^ (k + 2 + 1)) / ((k + 2 + 2 : ℕ) : ℝ) * X05.z (k + 2 + 2)
          * |τ| ^ (k + 2 + 2)
        ≤ 2 ^ (k + 2 + 2) * 1 / 1 * 2 * |τ| ^ (k + 2 + 2) := by gcongr
      _ = 32 * τ ^ 4 * (2 * |τ|) ^ k := by
        rw [show τ ^ 4 = |τ| ^ 4 from (Even.pow_abs ⟨2, rfl⟩ τ).symm]; ring)
  rw [abs_le] at ht ⊢
  norm_num [z2_eq] at ht
  constructor <;> linarith [ht.1, ht.2]

lemma lg2 {u v R : ℝ} (hR : R ≤ 1 / 2) (hu : |u| ≤ R) (hv : |v| ≤ R) (huv : |u + v| ≤ R) :
    |Real.log (Real.Gamma (1 + u) * Real.Gamma (1 + v) / Real.Gamma (1 + (u + v)))
      - (-(X05.z 2) * (u * v) + X05.z 3 * (u * v * (u + v)))| ≤ 12 * R ^ 4 := by
  have hR0 : 0 ≤ R := (abs_nonneg u).trans hu
  have hs := logGamma2 (by linarith) hu hv huv
  have pw : ∀ x : ℝ, |x| ≤ R → ∀ n, |x ^ n| ≤ R ^ n := fun x hx n => by
    rw [abs_pow]; exact pow_le_pow_left₀ (abs_nonneg x) hx n
  have ht := tail2 hs hR0 hR (A := 6 * R ^ 4) (fun k => by
    have hz := z_le (show 2 ≤ k + 2 + 2 by omega)
    have hz0 := z_nn (show 2 ≤ k + 2 + 2 by omega)
    have hn : (1 : ℝ) ≤ ((k + 2 + 2 : ℕ) : ℝ) := by push_cast; linarith [k.cast_nonneg (α := ℝ)]
    have h3 : |u ^ (k + 2 + 2) + v ^ (k + 2 + 2) - (u + v) ^ (k + 2 + 2)| ≤ 3 * R ^ (k + 2 + 2) := by
      have a := pw u hu (k + 2 + 2); have b := pw v hv (k + 2 + 2); have c := pw _ huv (k + 2 + 2)
      have := abs_sub (u ^ (k + 2 + 2) + v ^ (k + 2 + 2)) ((u + v) ^ (k + 2 + 2))
      have := abs_add_le (u ^ (k + 2 + 2)) (v ^ (k + 2 + 2))
      linarith
    rw [abs_mul, abs_div, abs_mul, abs_pow, abs_neg, abs_one, one_pow, one_mul, abs_of_nonneg hz0,
      abs_of_pos (by positivity : (0 : ℝ) < ((k + 2 + 2 : ℕ) : ℝ))]
    calc X05.z (k + 2 + 2) / ((k + 2 + 2 : ℕ) : ℝ) *
          |u ^ (k + 2 + 2) + v ^ (k + 2 + 2) - (u + v) ^ (k + 2 + 2)|
        ≤ 2 / 1 * (3 * R ^ (k + 2 + 2)) := by gcongr
      _ = 6 * R ^ 4 * R ^ k := by ring)
  rw [abs_le] at ht ⊢
  norm_num at ht
  constructor <;> nlinarith [ht.1, ht.2]

end Y05

namespace Y05

lemma gam_lb : ∃ ε > 0, ∀ τ : ℝ, |τ| ≤ ε → 1 / 2 ≤ Real.Gamma (1 - 2 * τ / 3) := by
  have hcont : ContinuousAt Real.Gamma 1 :=
    (Real.differentiableAt_Gamma (fun m h => by
      have := m.cast_nonneg (α := ℝ); linarith)).continuousAt
  obtain ⟨ε, hε, hεb⟩ := Metric.continuousAt_iff.mp hcont (1 / 2) (by norm_num)
  refine ⟨ε / 2, by linarith, fun τ hτ => ?_⟩
  have := hεb (x := 1 - 2 * τ / 3) (by
    rw [Real.dist_eq, show 1 - 2 * τ / 3 - 1 = -(2 * τ / 3) by ring, abs_neg, abs_div,
      abs_mul]
    norm_num; linarith)
  rw [Real.dist_eq, Real.Gamma_one] at this
  linarith [(abs_lt.mp this).1]

lemma wt_pack {τ : ℝ} (hτ : |τ| < 1 / 2) :
    IntegrableOn (wt τ) Q ∧ (∫ g in Q, wt τ g) = Real.Gamma (1 - 2 * τ / 3) ^ 3 ∧
      IntegrableOn (fun g => wt τ g * Ymean τ g) Q := by
  obtain ⟨t1, t2⟩ := abs_lt.mp hτ
  have ha : (-1 : ℝ) < -2 * τ / 3 := by linarith
  obtain ⟨i0', v0'⟩ := mono_int ha ha ha
  refine ⟨i0'.congr_fun (fun g hg => wt_mono τ hg) Q_meas, ?_, ?_⟩
  · rw [← setIntegral_congr_fun Q_meas (fun g hg => wt_mono τ hg), v0',
      show 1 + -2 * τ / 3 = 1 - 2 * τ / 3 by ring]
    ring
  · have h1 : -1 < 2 * τ / 3 := by linarith
    have h2 : -1 < -(4 * τ / 3) := by linarith
    have := (((mono_int h1 h2 h2).1.add (mono_int h2 h1 h2).1).add
      (mono_int h2 h2 h1).1).const_mul (1 / 3)
    exact IntegrableOn.congr_fun this (fun g hg => (pt5 τ hg).symm) Q_meas

lemma EX_eq {τ : ℝ} (hτ : |τ| < 1 / 2) :
    gammaE τ (fun g => Ymean τ g - 1) = gammaE τ (Ymean τ) - 1 := by
  obtain ⟨i0, v0, iY⟩ := wt_pack hτ
  obtain ⟨t1, t2⟩ := abs_lt.mp hτ
  have hG : 0 < Real.Gamma (1 - 2 * τ / 3) := Real.Gamma_pos_of_pos (by linarith)
  rw [gE_wt, gE_wt]
  have e : (fun g => wt τ g * (Ymean τ g - 1)) = fun g => wt τ g * Ymean τ g - wt τ g := by
    funext g; ring
  rw [e, integral_sub iY i0, v0]
  generalize Real.Gamma (1 - 2 * τ / 3) = G at *
  field_simp

lemma EY_exp {τ : ℝ} (hτ : |τ| ≤ 1 / 8) :
    |Real.log (gammaE τ (Ymean τ)) - (4 / 3 * X05.z 2 * τ ^ 2 + 32 / 27 * X05.z 3 * τ ^ 3)|
      ≤ 384 * τ ^ 4 := by
  rw [Y05.proof_5 τ (by linarith)]
  set R := 2 * |τ| with hR
  have hR2 : R ≤ 1 / 2 := by linarith
  have ht0 := abs_nonneg τ
  have a1 : |2 * τ / 3| ≤ R := by rw [abs_div, abs_mul]; norm_num; linarith
  have a2 : |-(4 * τ / 3)| ≤ R := by rw [abs_neg, abs_div, abs_mul]; norm_num; linarith
  have a3 : |2 * τ / 3 + -(4 * τ / 3)| ≤ R := by
    rw [show 2 * τ / 3 + -(4 * τ / 3) = -(2 * τ / 3) by ring, abs_neg, abs_div, abs_mul]
    norm_num; linarith
  have a4 : |-(2 * τ / 3)| ≤ R := by rw [abs_neg]; exact a1
  have a5 : |-(2 * τ / 3) + -(2 * τ / 3)| ≤ R := by
    rw [show -(2 * τ / 3) + -(2 * τ / 3) = -(4 * τ / 3) by ring]; exact a2
  have h1 := lg2 hR2 a1 a2 a3
  have h2 := lg2 hR2 a4 a4 a5
  rw [show 1 + -(4 * τ / 3) = 1 - 4 * τ / 3 by ring,
    show 1 + (2 * τ / 3 + -(4 * τ / 3)) = 1 - 2 * τ / 3 by ring] at h1
  rw [show 1 + -(2 * τ / 3) = 1 - 2 * τ / 3 by ring,
    show 1 + (-(2 * τ / 3) + -(2 * τ / 3)) = 1 - 4 * τ / 3 by ring] at h2
  obtain ⟨t1, t2⟩ := abs_le.mp hτ
  have ga : 0 < Real.Gamma (1 - 2 * τ / 3) := Real.Gamma_pos_of_pos (by linarith)
  have gb : 0 < Real.Gamma (1 - 4 * τ / 3) := Real.Gamma_pos_of_pos (by linarith)
  have gc : 0 < Real.Gamma (1 + 2 * τ / 3) := Real.Gamma_pos_of_pos (by linarith)
  have hlog : Real.log (Real.Gamma (1 + 2 * τ / 3) * Real.Gamma (1 - 4 * τ / 3) ^ 2 /
      Real.Gamma (1 - 2 * τ / 3) ^ 3) =
      Real.log (Real.Gamma (1 + 2 * τ / 3) * Real.Gamma (1 - 4 * τ / 3) /
        Real.Gamma (1 - 2 * τ / 3)) -
      Real.log (Real.Gamma (1 - 2 * τ / 3) * Real.Gamma (1 - 2 * τ / 3) /
        Real.Gamma (1 - 4 * τ / 3)) := by
    rw [← Real.log_div (by positivity) (by positivity)]
    congr 1
    generalize Real.Gamma (1 - 2 * τ / 3) = a at *
    generalize Real.Gamma (1 - 4 * τ / 3) = b at *
    generalize Real.Gamma (1 + 2 * τ / 3) = c at *
    field_simp
  rw [hlog]
  have hR4 : R ^ 4 = 16 * τ ^ 4 := by
    rw [hR, mul_pow, show |τ| ^ 4 = τ ^ 4 from Even.pow_abs ⟨2, rfl⟩ τ]; norm_num
  rw [hR4] at h1 h2
  rw [abs_le] at h1 h2 ⊢
  constructor <;> nlinarith [h1.1, h1.2, h2.1, h2.2]

lemma EX_exp {τ : ℝ} (hτ : |τ| ≤ 1 / 8) :
    |gammaE τ (fun g => Ymean τ g - 1) - (4 / 3 * X05.z 2 * τ ^ 2 + 32 / 27 * X05.z 3 * τ ^ 3)|
      ≤ 484 * τ ^ 4 := by
  rw [EX_eq (by linarith)]
  have hV : 0 < gammaE τ (Ymean τ) := by
    rw [Y05.proof_5 τ (by linarith)]
    obtain ⟨t1, t2⟩ := abs_le.mp hτ
    have ga : 0 < Real.Gamma (1 - 2 * τ / 3) := Real.Gamma_pos_of_pos (by linarith)
    have gb : 0 < Real.Gamma (1 - 4 * τ / 3) := Real.Gamma_pos_of_pos (by linarith)
    have gc : 0 < Real.Gamma (1 + 2 * τ / 3) := Real.Gamma_pos_of_pos (by linarith)
    positivity
  have hL := EY_exp hτ
  rw [← Real.exp_log hV]
  set L := Real.log (gammaE τ (Ymean τ))
  set t := |τ| with ht
  have t0 : 0 ≤ t := abs_nonneg τ
  have e2 : τ ^ 2 = t ^ 2 := (sq_abs τ).symm
  have e4 : τ ^ 4 = t ^ 4 := (Even.pow_abs ⟨2, rfl⟩ τ).symm
  have e3 : |τ ^ 3| = t ^ 3 := by rw [abs_pow]
  have z2 := z_le (le_refl 2); have z3 := z_le (show 2 ≤ 3 by norm_num)
  have z2' := z_nn (le_refl 2); have z3' := z_nn (show 2 ≤ 3 by norm_num)
  set P := 4 / 3 * X05.z 2 * τ ^ 2 + 32 / 27 * X05.z 3 * τ ^ 3
  have hP : |P| ≤ 3 * t ^ 2 := by
    show |4 / 3 * X05.z 2 * τ ^ 2 + 32 / 27 * X05.z 3 * τ ^ 3| ≤ 3 * t ^ 2
    have := abs_add_le (4 / 3 * X05.z 2 * τ ^ 2) (32 / 27 * X05.z 3 * τ ^ 3)
    rw [abs_mul, abs_mul, abs_mul, abs_mul, e3, abs_of_nonneg z2', abs_of_nonneg z3',
      abs_of_pos (by norm_num : (0 : ℝ) < 4 / 3), abs_of_pos (by norm_num : (0 : ℝ) < 32 / 27),
      abs_of_nonneg (by positivity : 0 ≤ τ ^ 2)] at this
    have h1 : X05.z 2 * t ^ 2 ≤ 2 * t ^ 2 := by nlinarith
    have h2 : X05.z 3 * t ^ 3 ≤ 2 * t ^ 3 := by nlinarith [pow_nonneg t0 3]
    have h3 : t ^ 3 ≤ t ^ 2 / 8 := by nlinarith [pow_nonneg t0 2]
    have h1' : X05.z 2 * τ ^ 2 ≤ 2 * t ^ 2 := by rw [e2]; exact h1
    have := pow_nonneg t0 2
    linarith
  have hL1 : |L| ≤ 9 * t ^ 2 := by
    have : |L| ≤ |L - P| + |P| := by
      have := abs_add_le (L - P) P; rwa [sub_add_cancel] at this
    have h4 : 384 * τ ^ 4 ≤ 6 * t ^ 2 := by
      have h5 : t ^ 2 ≤ 1 / 64 := by nlinarith
      rw [e4, show t ^ 4 = t ^ 2 * t ^ 2 by ring]
      nlinarith [pow_nonneg t0 2]
    linarith
  have hL2 : |L| ≤ 1 := by nlinarith
  have hE := Real.abs_exp_sub_one_sub_id_le hL2
  have hsq : L ^ 2 ≤ 81 * τ ^ 4 := by
    rw [e4, ← sq_abs L]
    have := pow_le_pow_left₀ (abs_nonneg L) hL1 2
    nlinarith
  have : |Real.exp L - 1 - P| ≤ |Real.exp L - 1 - L| + |L - P| := by
    have := abs_add_le (Real.exp L - 1 - L) (L - P)
    rwa [show Real.exp L - 1 - L + (L - P) = Real.exp L - 1 - P by ring] at this
  have h4 : 0 ≤ τ ^ 4 := by positivity
  linarith

end Y05

namespace Y05

lemma algX {q1 q2 q3 E D : ℝ} (h1 : 0 < q1) (h2 : 0 < q2) (h3 : 0 < q3)
    (u1 : q1 ≤ E) (u2 : q2 ≤ E) (u3 : q3 ≤ E)
    (v1 : q1⁻¹ ≤ E) (v2 : q2⁻¹ ≤ E) (v3 : q3⁻¹ ≤ E)
    (d12 : |q1 - q2| ≤ D) (d23 : |q2 - q3| ≤ D) (d13 : |q1 - q3| ≤ D) :
    0 ≤ (q1 ^ 3 + q2 ^ 3 + q3 ^ 3) / (3 * (q1 * q2 * q3)) - 1 ∧
      (q1 ^ 3 + q2 ^ 3 + q3 ^ 3) / (3 * (q1 * q2 * q3)) - 1 ≤ (3 / 2) * D ^ 2 * E ^ 4 ∧
      (q1 * q2 * q3)⁻¹ ≤ E ^ 3 := by
  have hE0 : 0 ≤ E := le_trans h1.le u1
  generalize hY : (q1 ^ 3 + q2 ^ 3 + q3 ^ 3) / (3 * (q1 * q2 * q3)) = Y
  have hX : Y - 1 = (q1 + q2 + q3) * ((q1 - q2) ^ 2 + (q2 - q3) ^ 2 + (q1 - q3) ^ 2) *
      (q1 * q2 * q3)⁻¹ / 6 := by
    rw [← hY]; field_simp; ring
  have hP : (q1 * q2 * q3)⁻¹ ≤ E ^ 3 := by
    rw [mul_inv, mul_inv, pow_three]
    have := mul_le_mul v1 v2 (by positivity) hE0
    calc q1⁻¹ * q2⁻¹ * q3⁻¹ ≤ (E * E) * E := mul_le_mul this v3 (by positivity) (by positivity)
      _ = E * (E * E) := by ring
  have sq : ∀ d : ℝ, |d| ≤ D → d ^ 2 ≤ D ^ 2 := fun d hd => by
    rw [← sq_abs]; exact pow_le_pow_left₀ (abs_nonneg d) hd 2
  have hS : (q1 - q2) ^ 2 + (q2 - q3) ^ 2 + (q1 - q3) ^ 2 ≤ 3 * D ^ 2 := by
    have := sq _ d12; have := sq _ d23; have := sq _ d13; linarith
  refine ⟨by rw [hX]; positivity, ?_, hP⟩
  rw [hX]
  have hq : q1 + q2 + q3 ≤ 3 * E := by linarith
  calc (q1 + q2 + q3) * ((q1 - q2) ^ 2 + (q2 - q3) ^ 2 + (q1 - q3) ^ 2) * (q1 * q2 * q3)⁻¹ / 6
      ≤ (3 * E) * (3 * D ^ 2) * E ^ 3 / 6 := by gcongr
    _ = _ := by ring

lemma M4 {M : ℝ} (hM : 0 ≤ M) : M ^ 4 ≤ 8 ^ 4 * Real.exp (M / 2) := by
  have := Real.add_one_le_exp (M / 8)
  have h : M / 8 ≤ Real.exp (M / 8) := by linarith
  have h6 := pow_le_pow_left₀ (by positivity) h 4
  rw [← Real.exp_nat_mul] at h6
  calc M ^ 4 = 8 ^ 4 * (M / 8) ^ 4 := by ring
    _ ≤ 8 ^ 4 * Real.exp (((4 : ℕ) : ℝ) * (M / 8)) := by gcongr
    _ = 8 ^ 4 * Real.exp (M / 2) := by congr 2; push_cast; ring

lemma pt8 {τ : ℝ} (hτ : |τ| ≤ 1 / 70) {g : ℝ × ℝ × ℝ} (hg : g ∈ Q) :
    |wt τ g * (Ymean τ g - 1) ^ 2| ≤ 36 * 8 ^ 4 * τ ^ 4 * (hh g.1 * (hh g.2.1 * hh g.2.2)) := by
  unfold wt
  obtain ⟨ha, hb, hc⟩ := hg
  simp only [mem_Ioi] at ha hb hc
  set c := 2 * τ / 3 with hcdef
  have hY : Ymean τ g = ((Real.exp (c * Real.log g.1)) ^ 3 + (Real.exp (c * Real.log g.2.1)) ^ 3
      + (Real.exp (c * Real.log g.2.2)) ^ 3) / (3 * (Real.exp (c * Real.log g.1) *
        Real.exp (c * Real.log g.2.1) * Real.exp (c * Real.log g.2.2))) := by
    unfold Ymean
    rw [Real.mul_rpow (mul_pos ha hb).le hc.le, Real.mul_rpow ha.le hb.le, ← hcdef,
      rp_q ha, rp_q hb, rp_q hc, rp_q ha, rp_q hb, rp_q hc]
    have k : ∀ x : ℝ, Real.exp (2 * τ * x) = Real.exp (c * x) ^ 3 := fun x => by
      rw [← Real.exp_nat_mul]; congr 1; push_cast; rw [hcdef]; ring
    rw [k, k, k]
  have hw : (g.1 * g.2.1 * g.2.2) ^ (-2 * τ / 3) = (Real.exp (c * Real.log g.1) *
        Real.exp (c * Real.log g.2.1) * Real.exp (c * Real.log g.2.2))⁻¹ := by
    rw [Real.mul_rpow (mul_pos ha hb).le hc.le, Real.mul_rpow ha.le hb.le,
      rp_q ha, rp_q hb, rp_q hc, mul_inv, mul_inv, ← Real.exp_neg, ← Real.exp_neg,
      ← Real.exp_neg, hcdef]
    congr 3 <;> ring
  rw [hY, hw]
  set l1 := Real.log g.1
  set l2 := Real.log g.2.1
  set l3 := Real.log g.2.2
  set M := |l1| + |l2| + |l3| with hM
  have hM0 : 0 ≤ M := by positivity
  set E := Real.exp (|c| * M) with hE
  have hc0 : 0 ≤ |c| := abs_nonneg c
  have up : ∀ l : ℝ, |l| ≤ M → Real.exp (c * l) ≤ E := fun l hl => by
    apply Real.exp_le_exp.mpr
    calc c * l ≤ |c * l| := le_abs_self _
      _ = |c| * |l| := abs_mul c l
      _ ≤ |c| * M := mul_le_mul_of_nonneg_left hl hc0
  have upi : ∀ l : ℝ, |l| ≤ M → (Real.exp (c * l))⁻¹ ≤ E := fun l hl => by
    rw [← Real.exp_neg, ← mul_neg]; apply up; rwa [abs_neg]
  have b1 : |l1| ≤ M := by have := abs_nonneg l2; have := abs_nonneg l3; linarith
  have b2 : |l2| ≤ M := by have := abs_nonneg l1; have := abs_nonneg l3; linarith
  have b3 : |l3| ≤ M := by have := abs_nonneg l2; have := abs_nonneg l1; linarith
  have dd : ∀ l l' : ℝ, |l| ≤ M → |l'| ≤ M → |l| + |l'| ≤ M →
      |Real.exp (c * l) - Real.exp (c * l')| ≤ 2 * |c| * M * E := fun l l' hl hl' hll => by
    calc |Real.exp (c * l) - Real.exp (c * l')|
        ≤ |c * l - c * l'| * (Real.exp (c * l) + Real.exp (c * l')) := exp_sub_le _ _
      _ ≤ (|c| * M) * (E + E) := by
          apply mul_le_mul _ (add_le_add (up l hl) (up l' hl')) (by positivity) (by positivity)
          rw [← mul_sub, abs_mul]
          exact mul_le_mul_of_nonneg_left ((abs_sub _ _).trans hll) hc0
      _ = 2 * |c| * M * E := by ring
  obtain ⟨X0, Xb, hP⟩ := algX (Real.exp_pos (c * l1)) (Real.exp_pos (c * l2))
    (Real.exp_pos (c * l3))
    (up l1 b1) (up l2 b2) (up l3 b3) (upi l1 b1) (upi l2 b2) (upi l3 b3)
    (dd l1 l2 b1 b2 (by have := abs_nonneg l3; linarith))
    (dd l2 l3 b2 b3 (by have := abs_nonneg l1; linarith))
    (dd l1 l3 b1 b3 (by have := abs_nonneg l2; linarith))
  have hW0 : 0 ≤ (Real.exp (c * l1) * Real.exp (c * l2) * Real.exp (c * l3))⁻¹ := by positivity
  have hprod : Real.exp (-(g.1 + g.2.1 + g.2.2)) * Real.exp (7 / 10 * M) ≤
      hh g.1 * (hh g.2.1 * hh g.2.2) := by
    have e1 : Real.exp (-(g.1 + g.2.1 + g.2.2)) * Real.exp (7 / 10 * M) =
        (Real.exp (-g.1) * Real.exp (7 / 10 * |l1|)) * ((Real.exp (-g.2.1) *
          Real.exp (7 / 10 * |l2|)) * (Real.exp (-g.2.2) * Real.exp (7 / 10 * |l3|))) := by
      simp only [← Real.exp_add]; congr 1; rw [hM]; ring
    rw [e1]
    unfold hh
    simp only [← mul_add]
    have k1 := ekap (κ := 7 / 10) ha
    have k2 := ekap (κ := 7 / 10) hb
    have k3 := ekap (κ := 7 / 10) hc
    gcongr
  have hc4 : |c| ^ 4 ≤ τ ^ 4 := by
    rw [hcdef, abs_div, abs_mul, show |(2 : ℝ)| = 2 by norm_num, show |(3 : ℝ)| = 3 by norm_num,
      div_pow, mul_pow, show |τ| ^ 4 = τ ^ 4 from Even.pow_abs ⟨2, rfl⟩ τ]
    have := pow_nonneg (abs_nonneg τ) 4
    rw [show |τ| ^ 4 = τ ^ 4 from Even.pow_abs ⟨2, rfl⟩ τ] at this
    linarith
  have h15 : 15 * |c| * M ≤ M / 5 := by
    have : |c| = 2 * |τ| / 3 := by
      rw [hcdef, abs_div, abs_mul]; norm_num
    rw [this]; nlinarith [abs_nonneg τ]
  have hEpow : E ^ 3 * ((3 / 2) * (2 * |c| * M * E) ^ 2 * E ^ 4) ^ 2 =
      36 * |c| ^ 4 * M ^ 4 * Real.exp (15 * |c| * M) := by
    rw [hE, show 15 * |c| * M = ((15 : ℕ) : ℝ) * (|c| * M) by push_cast; ring, Real.exp_nat_mul]
    ring
  have hm4 := M4 hM0
  have hτ4 : 0 ≤ τ ^ 4 := by positivity
  set e := Real.exp (-(g.1 + g.2.1 + g.2.2))
  have he : 0 < e := Real.exp_pos _
  generalize (Real.exp (c * l1) ^ 3 + Real.exp (c * l2) ^ 3 + Real.exp (c * l3) ^ 3) /
    (3 * (Real.exp (c * l1) * Real.exp (c * l2) * Real.exp (c * l3))) - 1 = X at *
  generalize (Real.exp (c * l1) * Real.exp (c * l2) * Real.exp (c * l3))⁻¹ = W at *
  rw [abs_of_nonneg (by positivity)]
  calc W * e * X ^ 2 ≤ E ^ 3 * e * ((3 / 2) * (2 * |c| * M * E) ^ 2 * E ^ 4) ^ 2 :=
        mul_le_mul (mul_le_mul_of_nonneg_right hP he.le) (pow_le_pow_left₀ X0 Xb 2)
          (by positivity) (by positivity)
    _ = e * (36 * |c| ^ 4 * M ^ 4 * Real.exp (15 * |c| * M)) := by rw [← hEpow]; ring
    _ ≤ e * (36 * τ ^ 4 * (8 ^ 4 * Real.exp (M / 2)) * Real.exp (M / 5)) := by
        gcongr
    _ = 36 * 8 ^ 4 * τ ^ 4 * (e * Real.exp (7 / 10 * M)) := by
        rw [show 7 / 10 * M = M / 2 + M / 5 by ring, Real.exp_add]; ring
    _ ≤ 36 * 8 ^ 4 * τ ^ 4 * (hh g.1 * (hh g.2.1 * hh g.2.2)) := by gcongr

end Y05

namespace Y05

lemma EX2_bd : ∃ C δ : ℝ, 0 < δ ∧ ∀ τ : ℝ, |τ| ≤ δ →
    |gammaE τ (fun g => (Ymean τ g - 1) ^ 2)| ≤ C * τ ^ 4 := by
  obtain ⟨ε, hε, hG⟩ := gam_lb
  set K := ∫ g in Q, hh g.1 * (hh g.2.1 * hh g.2.2)
  refine ⟨8 * (36 * 8 ^ 4 * K), min (1 / 70) ε, lt_min (by norm_num) hε, fun τ hτ => ?_⟩
  have hτ1 : |τ| ≤ 1 / 70 := hτ.trans (min_le_left _ _)
  have hg := hG τ (hτ.trans (min_le_right _ _))
  rw [gE_wt, abs_mul, abs_of_pos (by positivity)]
  have hbd : |∫ g in Q, wt τ g * (Ymean τ g - 1) ^ 2| ≤
      ∫ g in Q, 36 * 8 ^ 4 * τ ^ 4 * (hh g.1 * (hh g.2.1 * hh g.2.2)) := by
    rw [← Real.norm_eq_abs]
    refine norm_integral_le_of_norm_le (hh3_int.const_mul _) ?_
    refine (ae_restrict_iff' Q_meas).mpr (Eventually.of_forall fun g hg => ?_)
    rw [Real.norm_eq_abs]
    exact pt8 hτ1 hg
  rw [integral_const_mul] at hbd
  have hinv : (Real.Gamma (1 - 2 * τ / 3) ^ 3)⁻¹ ≤ 8 := by
    rw [inv_le_comm₀ (by positivity) (by norm_num)]
    calc 8⁻¹ = (1 / 2 : ℝ) ^ 3 := by norm_num
      _ ≤ _ := pow_le_pow_left₀ (by norm_num) hg 3
  calc (Real.Gamma (1 - 2 * τ / 3) ^ 3)⁻¹ * |∫ g in Q, wt τ g * (Ymean τ g - 1) ^ 2|
      ≤ 8 * (36 * 8 ^ 4 * τ ^ 4 * K) := mul_le_mul hinv hbd (abs_nonneg _) (by norm_num)
    _ = 8 * (36 * 8 ^ 4 * K) * τ ^ 4 := by ring

lemma factor_id {p : ℝ} (hp : 0 < p) (hτ : |1 - 1 / p| ≤ 1 / 8) :
    (2 - p) * S p / (4 * Real.sqrt 3) =
      Real.Gamma (1 - 2 * (1 - 1 / p) / 3) ^ 3 / Real.Gamma (1 - 2 * (1 - 1 / p)) *
        B (1 - 1 / p) := by
  rw [Y05.proof_3 p hp]
  set τ := 1 - 1 / p with hτdef
  obtain ⟨t1, t2⟩ := abs_le.mp hτ
  have n1 : 1 - 2 * τ ≠ 0 := by intro h; linarith
  have n2 : 2 - 2 * τ ≠ 0 := by intro h; linarith
  have hG3 : Real.Gamma (3 - 2 * τ) = (2 - 2 * τ) * ((1 - 2 * τ) * Real.Gamma (1 - 2 * τ)) := by
    rw [show 3 - 2 * τ = (2 - 2 * τ) + 1 by ring, Real.Gamma_add_one n2,
      show 2 - 2 * τ = (1 - 2 * τ) + 1 by ring, Real.Gamma_add_one n1]
  have hkey : (2 - p) * (1 / p) ^ 2 = (2 - 2 * τ) * (1 - 2 * τ) / 2 := by
    rw [hτdef]; field_simp; ring
  have g1 : 0 < Real.Gamma (1 - 2 * τ) := Real.Gamma_pos_of_pos (by linarith)
  rw [hG3]
  have s3 : 0 < Real.sqrt 3 := by positivity
  generalize Real.Gamma (1 - 2 * τ / 3) = A at *
  generalize Real.Gamma (1 - 2 * τ) = G at *
  generalize B τ = b at *
  have : (2 - p) * (8 * Real.sqrt 3 * (1 / p) ^ 2 * A ^ 3 / ((2 - 2 * τ) * ((1 - 2 * τ) * G)) * b)
      / (4 * Real.sqrt 3) = ((2 - p) * (1 / p) ^ 2) *
        (8 * Real.sqrt 3 * A ^ 3 * b / ((2 - 2 * τ) * ((1 - 2 * τ) * G)) / (4 * Real.sqrt 3)) := by
    ring
  have n3 : 1 - τ ≠ 0 := by intro h; linarith
  rw [this, hkey]
  field_simp
  ring

end Y05

theorem Y05.proof_8 : Y05.stmt_8 := by
  obtain ⟨C7, δ7, hδ7, h7⟩ := Y05.proof_7
  obtain ⟨C2, δ2, hδ2, h2⟩ := EX2_bd
  set K1 := 242 + |C2| / 8 + |C7| with hK1
  have hK10 : 0 ≤ K1 := by positivity
  refine ⟨64 + K1 + 2 * (3 / 2 + K1) ^ 2, min (min δ7 δ2) (min (1 / 8) (1 / (4 * (K1 + 2)))),
    lt_min (lt_min hδ7 hδ2) (lt_min (by norm_num) (by positivity)), fun p hp hτ => ?_⟩
  rw [factor_id hp (hτ.trans ((min_le_right _ _).trans (min_le_left _ _)))]
  set τ := 1 - 1 / p with hτdef
  have hτ7 : |τ| ≤ δ7 := hτ.trans ((min_le_left _ _).trans (min_le_left _ _))
  have hτ2 : |τ| ≤ δ2 := hτ.trans ((min_le_left _ _).trans (min_le_right _ _))
  have hτ8 : |τ| ≤ 1 / 8 := hτ.trans ((min_le_right _ _).trans (min_le_left _ _))
  have hτK : |τ| ≤ 1 / (4 * (K1 + 2)) := hτ.trans ((min_le_right _ _).trans (min_le_right _ _))
  set t := |τ| with ht
  have t0 : 0 ≤ t := abs_nonneg τ
  have e2 : τ ^ 2 = t ^ 2 := (sq_abs τ).symm
  have e4 : τ ^ 4 = t ^ 4 := (Even.pow_abs ⟨2, rfl⟩ τ).symm
  have e6 : τ ^ 6 = t ^ 6 := (Even.pow_abs ⟨3, rfl⟩ τ).symm
  have e3 : |τ ^ 3| = t ^ 3 := by rw [abs_pow]
  have hA := h7 τ hτ7
  have hX2 := h2 τ hτ2
  have hX1 := EX_exp hτ8
  have hF := ser4 (show |τ| ≤ 1 / 4 by linarith)
  set E1 := gammaE τ (fun g => Ymean τ g - 1)
  set E2 := gammaE τ (fun g => (Ymean τ g - 1) ^ 2)
  set P := 4 / 3 * X05.z 2 * τ ^ 2 + 32 / 27 * X05.z 3 * τ ^ 3 with hPdef
  generalize B τ = b at hA ⊢
  clear_value P t E1 E2 τ
  have z2 := z_le (le_refl 2); have z3 := z_le (show 2 ≤ 3 by norm_num)
  have z2' := z_nn (le_refl 2); have z3' := z_nn (show 2 ≤ 3 by norm_num)
  have hP : |P| ≤ 3 * t ^ 2 := by
    rw [hPdef]
    have := abs_add_le (4 / 3 * X05.z 2 * τ ^ 2) (32 / 27 * X05.z 3 * τ ^ 3)
    rw [abs_mul, abs_mul, abs_mul, abs_mul, e3, abs_of_nonneg z2', abs_of_nonneg z3',
      abs_of_pos (by norm_num : (0 : ℝ) < 4 / 3), abs_of_pos (by norm_num : (0 : ℝ) < 32 / 27),
      abs_of_nonneg (by positivity : 0 ≤ τ ^ 2)] at this
    have h1 : X05.z 2 * t ^ 2 ≤ 2 * t ^ 2 := by nlinarith
    have h2 : X05.z 3 * t ^ 3 ≤ 2 * t ^ 3 := by nlinarith [pow_nonneg t0 3]
    have h3 : t ^ 3 ≤ t ^ 2 / 8 := by nlinarith [pow_nonneg t0 2]
    have h1' : X05.z 2 * τ ^ 2 ≤ 2 * t ^ 2 := by rw [e2]; exact h1
    have := pow_nonneg t0 2
    linarith
  -- y = B - 1
  have t1 : t ≤ 1 := by linarith
  have h64 : t ^ 6 ≤ t ^ 4 := pow_le_pow_of_le_one t0 t1 (by norm_num)
  have t4 : 0 ≤ t ^ 4 := by positivity
  have hy : |(b - 1) - P / 2| ≤ K1 * t ^ 4 := by
    have a1 : |C7 * |τ| ^ 6| ≤ |C7| * t ^ 4 := by
      rw [abs_mul, abs_pow, abs_abs, ← ht]
      exact mul_le_mul_of_nonneg_left h64 (abs_nonneg _)
    have a2 : |E2| ≤ |C2| * t ^ 4 := by
      rw [← e4]; exact hX2.trans (mul_le_mul_of_nonneg_right (le_abs_self _) (by positivity))
    have hA' : |b - (1 + E1 / 2 - E2 / 8)| ≤ |C7| * t ^ 4 := hA.trans ((le_abs_self _).trans a1)
    rw [abs_le] at hA' hX1 a2 ⊢
    rw [e4] at hX1
    rw [hK1]
    constructor <;> linarith [hA'.1, hA'.2, hX1.1, hX1.2, a2.1, a2.2]
  have hy2 : |b - 1| ≤ (3 / 2 + K1) * t ^ 2 := by
    have : |b - 1| ≤ |(b - 1) - P / 2| + |P / 2| := by
      have := abs_add_le ((b - 1) - P / 2) (P / 2); rwa [sub_add_cancel] at this
    rw [abs_div, abs_two] at this
    have h42 : t ^ 4 ≤ t ^ 2 := pow_le_pow_of_le_one t0 t1 (by norm_num)
    have := mul_le_mul_of_nonneg_left h42 hK10
    linarith
  have hsmall : (3 / 2 + K1) * t ^ 2 ≤ 1 / 4 := by
    have h1 : t * (4 * (K1 + 2)) ≤ 1 := by
      rw [le_div_iff₀ (by positivity)] at hτK; linarith
    have h2 : t ^ 2 ≤ t := by
      rw [pow_two]; exact (mul_le_mul_of_nonneg_right t1 t0).trans (by rw [one_mul])
    have h3 : (3 / 2 + K1) * t ^ 2 ≤ (K1 + 2) * t ^ 2 :=
      mul_le_mul_of_nonneg_right (by linarith) (by positivity)
    have h4 : (K1 + 2) * t ^ 2 ≤ (K1 + 2) * t := mul_le_mul_of_nonneg_left h2 (by linarith)
    linarith
  have hyB : |b - 1| ≤ 1 / 4 := hy2.trans hsmall
  have hBpos : 0 < b := by linarith only [(abs_le.mp hyB).1]
  have hlogB : |Real.log (b) - (b - 1)| ≤ 2 * (b - 1) ^ 2 := by
    have := Real.abs_log_sub_add_sum_range_le (x := -(b - 1))
      (by rw [abs_neg]; linarith only [hyB]) 1
    simp only [Finset.sum_range_one, zero_add, pow_one, Nat.cast_zero, div_one, abs_neg,
      sub_neg_eq_add, add_sub_cancel] at this
    rw [show -(b - 1) + Real.log (b) = Real.log (b) - (b - 1) by ring] at this
    refine this.trans ?_
    rw [div_le_iff₀ (by linarith only [hyB]), sq_abs]
    have := mul_le_mul_of_nonneg_left (show 1 ≤ 2 * (1 - |b - 1|) by linarith only [hyB])
      (sq_nonneg (b - 1))
    linarith only [this]
  have hsq : (b - 1) ^ 2 ≤ (3 / 2 + K1) ^ 2 * t ^ 4 := by
    rw [← sq_abs, show t ^ 4 = (t ^ 2) ^ 2 by ring, ← mul_pow]
    exact pow_le_pow_left₀ (abs_nonneg _) hy2 2
  obtain ⟨u1, u2⟩ := abs_le.mp (show |τ| ≤ 1 / 8 by rw [← ht]; exact hτ8)
  have g1 : 0 < Real.Gamma (1 - 2 * τ) := Real.Gamma_pos_of_pos (by linarith only [u2])
  have ga : 0 < Real.Gamma (1 - 2 * τ / 3) := Real.Gamma_pos_of_pos (by linarith only [u2])
  rw [Real.log_mul (by positivity) hBpos.ne']
  rw [e4] at hF
  rw [abs_le] at hF hlogB hy ⊢
  rw [hPdef, z2_eq] at hy
  constructor
  · linarith only [hsq, t4, hF.1, hlogB.1, hy.1]
  · linarith only [hsq, t4, hF.2, hlogB.2, hy.2]

/-! ## stmt_9: the `d`-coefficients through `d³` -/

namespace Y05

open Asymptotics in
lemma expstep4 (x G y : ℝ → ℝ) (hy : y =O[nhds 0] fun d : ℝ => d ^ 2)
    (hx : x =O[nhds 0] fun d : ℝ => d ^ 2) (hG : G =O[nhds 0] fun d : ℝ => d ^ 4)
    (hyx : ∀ d, y d = x d + G d) :
    (fun d => Real.exp (y d) - (1 + x d + x d ^ 2 / 2)) =O[nhds 0] fun d : ℝ => d ^ 4 := by
  have hy0 : Tendsto y (nhds 0) (nhds 0) :=
    hy.trans_tendsto (by simpa using ((continuous_pow 2).tendsto (0:ℝ)))
  have hy1 : ∀ᶠ d in nhds (0:ℝ), |y d| ≤ 1 := by
    filter_upwards [hy0.eventually (Metric.closedBall_mem_nhds (0:ℝ) one_pos)] with d hd
    simpa [Real.norm_eq_abs] using hd
  have h3 : (fun d => Real.exp (y d) - (1 + y d + y d ^ 2 / 2)) =O[nhds 0] fun d => y d ^ 3 := by
    refine IsBigO.of_bound 1 ?_
    filter_upwards [hy1] with d hd
    have := Real.exp_bound hd (n := 3) (by norm_num)
    simp [Finset.sum_range_succ, Nat.factorial] at this
    simp only [Real.norm_eq_abs, abs_pow]
    have h' : |Real.exp (y d) - (1 + y d + y d ^ 2 / 2)| ≤ |y d| ^ 3 * (4 / 18) := by
      convert this using 2 <;> norm_num <;> ring_nf
    nlinarith [pow_nonneg (abs_nonneg (y d)) 3]
  have h6 : (fun d => y d ^ 3) =O[nhds 0] fun d : ℝ => d ^ 6 := by
    have := hy.pow 3
    refine this.trans (isBigO_of_le _ (fun d => ?_))
    simp only [Real.norm_eq_abs, ← pow_mul]
    simp
  have h65 : (fun d : ℝ => d ^ 6) =O[nhds 0] fun d : ℝ => d ^ 4 := T09.hpow 4 2
  have hv : (fun d => y d + x d) =O[nhds 0] (1 : ℝ → ℝ) :=
    (hy.add hx).trans (T09.dpow_one 2 (by norm_num))
  have hGv : (fun d => G d * (1 + (y d + x d) / 2)) =O[nhds 0] fun d : ℝ => d ^ 4 := by
    have h1 : (fun d => 1 + (y d + x d) / 2) =O[nhds 0] (1 : ℝ → ℝ) := by
      refine (isBigO_const_const (1:ℝ) one_ne_zero _).add ?_
      exact (hv.const_mul_left (1/2)).congr_left (fun d => by ring)
    simpa using (hG.mul h1)
  have := (h3.trans (h6.trans h65)).add hGv
  refine this.congr_left (fun d => ?_)
  simp only [hyx d]
  ring

noncomputable def poly6 (c0 c1 c2 c3 c4 c5 c6 d : ℝ) : ℝ :=
  c0 + c1 * d + c2 * d ^ 2 + c3 * d ^ 3 + c4 * d ^ 4 + c5 * d ^ 5 + c6 * d ^ 6

lemma poly6_cont (c0 c1 c2 c3 c4 c5 c6 : ℝ) :
    Continuous (fun d : ℝ => poly6 c0 c1 c2 c3 c4 c5 c6 d) := by
  unfold poly6; fun_prop

noncomputable def R9 (P Z d : ℝ) : ℝ :=
  poly6
    ((1/1 : ℝ) * P ^ 0 * Z ^ 0 + (32/9 : ℝ) * P ^ 0 * Z ^ 1 + (-2/9 : ℝ) * P ^ 1 * Z ^ 0 + (1/162 : ℝ) * P ^ 2 * Z ^ 0)
    ((6/1 : ℝ) * P ^ 0 * Z ^ 0 + (32/3 : ℝ) * P ^ 0 * Z ^ 1 + (-8/9 : ℝ) * P ^ 1 * Z ^ 0 + (16/81 : ℝ) * P ^ 1 * Z ^ 1 + (1/81 : ℝ) * P ^ 2 * Z ^ 0)
    ((15/1 : ℝ) * P ^ 0 * Z ^ 0 + (64/9 : ℝ) * P ^ 0 * Z ^ 1 + (128/81 : ℝ) * P ^ 0 * Z ^ 2 + (-11/9 : ℝ) * P ^ 1 * Z ^ 0 + (16/81 : ℝ) * P ^ 1 * Z ^ 1 + (1/162 : ℝ) * P ^ 2 * Z ^ 0)
    ((20/1 : ℝ) * P ^ 0 * Z ^ 0 + (-80/9 : ℝ) * P ^ 0 * Z ^ 1 + (-4/9 : ℝ) * P ^ 1 * Z ^ 0)
    ((15/1 : ℝ) * P ^ 0 * Z ^ 0 + (-16/1 : ℝ) * P ^ 0 * Z ^ 1 + (4/9 : ℝ) * P ^ 1 * Z ^ 0)
    ((6/1 : ℝ) * P ^ 0 * Z ^ 0 + (-80/9 : ℝ) * P ^ 0 * Z ^ 1 + (4/9 : ℝ) * P ^ 1 * Z ^ 0)
    ((1/1 : ℝ) * P ^ 0 * Z ^ 0 + (-16/9 : ℝ) * P ^ 0 * Z ^ 1 + (1/9 : ℝ) * P ^ 1 * Z ^ 0)
    d / ((1 - d) * (1 + d) ^ 6)

lemma alg9 (Z d : ℝ) (h1 : 1 - d ≠ 0) (h2 : 1 + d ≠ 0) :
    (1 + (-(Real.pi ^ 2 / 9) * (1 - 1 / (1 + d)) ^ 2 - 16 * Z / 9 * (1 - 1 / (1 + d)) ^ 3)
      + (-(Real.pi ^ 2 / 9) * (1 - 1 / (1 + d)) ^ 2 - 16 * Z / 9 * (1 - 1 / (1 + d)) ^ 3) ^ 2 / 2)
      / (1 - d) - (1 + d + (1 - Real.pi ^ 2 / 9) * d ^ 2
        + (1 + Real.pi ^ 2 / 9 - 16 * Z / 9) * d ^ 3) = d ^ 4 * R9 (Real.pi ^ 2) Z d := by
  unfold R9 poly6
  field_simp
  ring

lemma R9_cont (P Z : ℝ) : ContinuousAt (fun d : ℝ => R9 P Z d) 0 := by
  unfold R9
  exact ContinuousAt.div (poly6_cont _ _ _ _ _ _ _).continuousAt
    (by fun_prop : Continuous fun d : ℝ => (1 - d) * (1 + d) ^ 6).continuousAt
    (by show (1 - (0:ℝ)) * (1 + 0) ^ 6 ≠ 0; norm_num)

lemma Fb_pos (τ : ℝ) {g : ℝ × ℝ × ℝ} (hg : g ∈ Q) : 0 < Fb τ g := by
  rw [Fb_form τ hg]
  obtain ⟨ha, hb, hc⟩ := hg
  simp only [mem_Ioi] at ha hb hc
  have h1 := Real.rpow_pos_of_pos ha (2 * τ)
  have h2 := Real.rpow_pos_of_pos hb (2 * τ)
  have h3 := Real.rpow_pos_of_pos hc (2 * τ)
  have hs : 0 < Real.sqrt (g.1 ^ (2 * τ) + g.2.1 ^ (2 * τ) + g.2.2 ^ (2 * τ)) :=
    Real.sqrt_pos.mpr (by positivity)
  have := Real.rpow_pos_of_pos (show 0 < g.1 * g.2.1 * g.2.2 by positivity) (-τ)
  positivity

lemma B_pos {τ : ℝ} (hτ : τ < 1) : 0 < B τ := by
  have hB : B τ = (Real.Gamma (1 - 2 * τ / 3) ^ 3)⁻¹ * ∫ g in Q, Fb τ g := rfl
  rw [hB]
  have hG : 0 < Real.Gamma (1 - 2 * τ / 3) := Real.Gamma_pos_of_pos (by linarith)
  apply mul_pos (by positivity)
  rw [setIntegral_pos_iff_support_of_nonneg_ae
    ((ae_restrict_iff' Q_meas).mpr (Eventually.of_forall fun g hg => Fb_nn τ g hg)) (Fb_int hτ)]
  have hs : Function.support (Fb τ) ∩ Q = Q := by
    ext g
    simp only [mem_inter_iff, Function.mem_support]
    exact ⟨fun h => h.2, fun hg => ⟨(Fb_pos τ hg).ne', hg⟩⟩
  rw [hs]
  have : volume Q = ⊤ := by
    rw [Measure.volume_eq_prod, Measure.prod_prod, Measure.volume_eq_prod, Measure.prod_prod]
    simp
  rw [this]; simp

end Y05

open Asymptotics in
theorem Y05.proof_9 : Y05.stmt_9 := by
  set z := X05.z 3
  set F : ℝ → ℝ := fun d => S (1 + d) / (4 * Real.sqrt 3) with hF
  set x : ℝ → ℝ := fun d => -(Real.pi ^ 2 / 9) * (1 - 1 / (1 + d)) ^ 2
    - 16 * z / 9 * (1 - 1 / (1 + d)) ^ 3 with hx
  set P3 : ℝ → ℝ := fun d => 1 + d + (1 - Real.pi ^ 2 / 9) * d ^ 2
    + (1 + Real.pi ^ 2 / 9 - 16 * z / 9) * d ^ 3 with hP3
  obtain ⟨C8, δ8, hδ8, h8⟩ := Y05.proof_8
  set η := min (1 / 16) (δ8 / 2) with hη
  have hη0 : 0 < η := lt_min (by norm_num) (by linarith)
  have hnear : ∀ᶠ d in nhds (0 : ℝ), |d| < η := by
    have := Metric.ball_mem_nhds (0 : ℝ) hη0
    filter_upwards [this] with d hd
    simpa [Real.dist_eq] using hd
  have hτb : ∀ d : ℝ, |d| < η → 0 < 1 + d ∧ |1 - 1 / (1 + d)| ≤ 2 * |d| := fun d hd => by
    have h16 : |d| < 1 / 16 := lt_of_lt_of_le hd (min_le_left _ _)
    obtain ⟨a, b⟩ := abs_lt.mp h16
    have hp : 0 < 1 + d := by linarith
    refine ⟨hp, ?_⟩
    rw [show 1 - 1 / (1 + d) = d / (1 + d) by field_simp; ring, abs_div, abs_of_pos hp,
      div_le_iff₀ hp]
    nlinarith [abs_nonneg d]
  have hpos : ∀ᶠ d in nhds (0 : ℝ), 0 < F d := by
    filter_upwards [hnear] with d hd
    obtain ⟨hp, hτ⟩ := hτb d hd
    have h16 : |d| < 1 / 16 := lt_of_lt_of_le hd (min_le_left _ _)
    have hτ8 : |1 - 1 / (1 + d)| ≤ 1 / 8 := by linarith
    have hfac := factor_id hp hτ8
    obtain ⟨u1, u2⟩ := abs_le.mp hτ8
    have g1 : 0 < Real.Gamma (1 - 2 * (1 - 1 / (1 + d))) := Real.Gamma_pos_of_pos (by linarith)
    have ga : 0 < Real.Gamma (1 - 2 * (1 - 1 / (1 + d)) / 3) :=
      Real.Gamma_pos_of_pos (by linarith)
    have hB := B_pos (show 1 - 1 / (1 + d) < 1 by linarith)
    have hpos' : 0 < (2 - (1 + d)) * S (1 + d) / (4 * Real.sqrt 3) := by
      rw [hfac]; positivity
    obtain ⟨a, b⟩ := abs_lt.mp h16
    have : (2 - (1 + d)) * S (1 + d) / (4 * Real.sqrt 3) = (1 - d) * F d := by
      rw [hF]; ring
    rw [this] at hpos'
    exact pos_of_mul_pos_right hpos' (by linarith)
  have hG : (fun d => Real.log ((1 - d) * F d) - x d) =O[nhds 0] fun d : ℝ => d ^ 4 := by
    refine IsBigO.of_bound (16 * |C8|) ?_
    filter_upwards [hnear] with d hd
    obtain ⟨hp, hτ⟩ := hτb d hd
    have hτδ : |1 - 1 / (1 + d)| ≤ δ8 := by
      have := lt_of_lt_of_le hd (min_le_right _ _); linarith
    have h := h8 (1 + d) hp hτδ
    have e : (2 - (1 + d)) * S (1 + d) / (4 * Real.sqrt 3) = (1 - d) * F d := by rw [hF]; ring
    rw [e] at h
    rw [Real.norm_eq_abs, Real.norm_eq_abs, abs_pow]
    calc |Real.log ((1 - d) * F d) - x d| ≤ C8 * |1 - 1 / (1 + d)| ^ 4 := h
      _ ≤ |C8| * (2 * |d|) ^ 4 := by
          apply mul_le_mul (le_abs_self _) (pow_le_pow_left₀ (abs_nonneg _) hτ 4)
            (by positivity) (abs_nonneg _)
      _ = 16 * |C8| * |d| ^ 4 := by ring
  have hxO : x =O[nhds 0] fun d : ℝ => d ^ 2 := by
    have hw : ContinuousAt (fun d : ℝ => -(Real.pi ^ 2 / 9) * ((1 + d)⁻¹) ^ 2
        - (16 * z / 9) * d * ((1 + d)⁻¹) ^ 3) 0 := by
      have := T09.cinv 1 one_ne_zero
      fun_prop
    have := T09.bdd hw (fun d => d ^ 2)
    refine this.congr' ?_ EventuallyEq.rfl
    filter_upwards [hnear] with d hd
    obtain ⟨hp, _⟩ := hτb d hd
    simp only [hx]
    rw [show 1 - 1 / (1 + d) = d * (1 + d)⁻¹ by field_simp; ring]
    ring
  have hy : (fun d : ℝ => Real.log ((1 - d) * F d)) =O[nhds 0] fun d : ℝ => d ^ 2 := by
    have := hxO.add (hG.trans (T09.hpow 2 2))
    refine this.congr_left (fun d => ?_)
    simp
  have hE := expstep4 _ _ _ hy hxO hG (fun d => by simp)
  have hT1 : (fun d : ℝ => (Real.exp (Real.log ((1 - d) * F d))
      - (1 + x d + x d ^ 2 / 2)) * (1 - d)⁻¹) =O[nhds 0] fun d : ℝ => d ^ 4 := by
    simpa using hE.mul (T09.cinv'.tendsto.isBigO_one ℝ)
  have hT2 : (fun d : ℝ => d ^ 4 * R9 (Real.pi ^ 2) z d) =O[nhds 0] fun d : ℝ => d ^ 4 :=
    T09.bdd (R9_cont _ _) _
  have hfin : (fun d => F d - P3 d) =O[nhds 0] fun d : ℝ => d ^ 4 := by
    refine (hT1.add hT2).congr' ?_ EventuallyEq.rfl
    filter_upwards [hnear, hpos] with d hd hFd
    have h16 : |d| < 1 / 16 := lt_of_lt_of_le hd (min_le_left _ _)
    obtain ⟨a, b⟩ := abs_lt.mp h16
    have h1 : 1 - d ≠ 0 := by intro h; linarith
    have h2 : 1 + d ≠ 0 := by intro h; linarith
    have hexp : Real.exp (Real.log ((1 - d) * F d)) = (1 - d) * F d :=
      Real.exp_log (mul_pos (by linarith) hFd)
    have hA := alg9 z d h1 h2
    show _ = F d - P3 d
    rw [hexp, ← hA]
    simp only [hx, hP3]
    field_simp
    ring
  obtain ⟨c, hc⟩ := hfin.bound
  obtain ⟨ε, hε, hεb⟩ := Metric.eventually_nhds_iff.mp hc
  refine ⟨c, ε / 2, by linarith, fun d hd => ?_⟩
  have := hεb (show dist d 0 < ε by rw [Real.dist_eq, sub_zero]; linarith)
  simpa [Real.norm_eq_abs, abs_pow, hF, hP3] using this

/-! ## stmt_10: the `τ⁴` coefficient and the `d`-series through `d⁴` -/

namespace Y05

lemma z4_eq : X05.z 4 = Real.pi ^ 4 / 90 := by
  unfold X05.z
  rw [show ((4 : ℕ) : ℂ) = 4 by norm_num, riemannZeta_four]
  rw [show (Real.pi : ℂ) ^ 4 / 90 = ((Real.pi ^ 4 / 90 : ℝ) : ℂ) by push_cast; ring,
    Complex.ofReal_re]

lemma tail3 {f : ℕ → ℝ} {L A r : ℝ} (hf : HasSum f L) (hr0 : 0 ≤ r) (hr : r ≤ 1 / 2)
    (hb : ∀ k, |f (k + 3)| ≤ A * r ^ k) : |L - f 0 - f 1 - f 2| ≤ 2 * A := by
  have hs : HasSum (fun k => f (k + 3)) (L - (f 0 + f 1 + f 2)) := by
    have := (hasSum_nat_add_iff' 3).mpr hf
    simpa [Finset.sum_range_succ] using this
  have hA : 0 ≤ A := by
    have := (abs_nonneg _).trans (hb 0); simpa using this
  have hg : HasSum (fun k => A * r ^ k) (A * (1 - r)⁻¹) :=
    (hasSum_geometric_of_lt_one hr0 (by linarith)).mul_left A
  have hup := hasSum_le (fun k => (le_abs_self _).trans (hb k)) hs hg
  have hlo := hasSum_le (fun k => neg_le_of_abs_le (hb k)) hg.neg hs
  have hinv : (1 - r)⁻¹ ≤ 2 := by
    rw [inv_le_comm₀ (by linarith) (by norm_num)]; linarith
  have : A * (1 - r)⁻¹ ≤ 2 * A := by nlinarith
  rw [show L - f 0 - f 1 - f 2 = L - (f 0 + f 1 + f 2) by ring, abs_le]
  constructor <;> linarith

/-- the truncated flat-face series -/
noncomputable def Φ4 (τ : ℝ) : ℝ :=
  -(4 / 3) * X05.z 2 * τ ^ 2 - 64 / 27 * X05.z 3 * τ ^ 3 - 104 / 27 * X05.z 4 * τ ^ 4

lemma ball_ev (ε : ℝ) (hε : 0 < ε) : ∀ᶠ τ in nhds (0 : ℝ), |τ| < ε := by
  filter_upwards [Metric.ball_mem_nhds (0 : ℝ) hε] with τ hτ
  simpa [Real.dist_eq] using hτ

open Asymptotics in
lemma flatO : (fun τ : ℝ => Real.log (Real.Gamma (1 - 2 * τ / 3) ^ 3 / Real.Gamma (1 - 2 * τ))
    - Φ4 τ) =O[nhds 0] fun τ => τ ^ 5 := by
  refine IsBigO.of_bound 128 ?_
  filter_upwards [ball_ev (1 / 4) (by norm_num)] with τ hτ
  have hs := Y05.proof_4 τ (by linarith)
  have ht := tail3 hs (by positivity : (0 : ℝ) ≤ 2 * |τ|) (by linarith) (A := 64 * |τ| ^ 5)
    (fun k => by
    have hz := z_le (show 2 ≤ k + 3 + 2 by omega)
    have hz0 := z_nn (show 2 ≤ k + 3 + 2 by omega)
    have h3 : 0 ≤ 1 - 1 / (3 : ℝ) ^ (k + 3 + 1) := by
      have : 1 / (3 : ℝ) ^ (k + 3 + 1) ≤ 1 := by
        rw [div_le_one (by positivity)]; exact one_le_pow₀ (by norm_num)
      linarith
    have h3' : 1 - 1 / (3 : ℝ) ^ (k + 3 + 1) ≤ 1 := by
      have : 0 ≤ 1 / (3 : ℝ) ^ (k + 3 + 1) := by positivity
      linarith
    have hn : (1 : ℝ) ≤ ((k + 3 + 2 : ℕ) : ℝ) := by push_cast; linarith [k.cast_nonneg (α := ℝ)]
    rw [abs_mul, abs_mul, abs_div, abs_neg, abs_of_nonneg (mul_nonneg (by positivity) h3),
      abs_of_pos (by positivity : (0 : ℝ) < ((k + 3 + 2 : ℕ) : ℝ)), abs_of_nonneg hz0, abs_pow]
    calc 2 ^ (k + 3 + 2) * (1 - 1 / 3 ^ (k + 3 + 1)) / ((k + 3 + 2 : ℕ) : ℝ) * X05.z (k + 3 + 2)
          * |τ| ^ (k + 3 + 2)
        ≤ 2 ^ (k + 3 + 2) * 1 / 1 * 2 * |τ| ^ (k + 3 + 2) := by gcongr
      _ = 64 * |τ| ^ 5 * (2 * |τ|) ^ k := by ring)
  rw [Real.norm_eq_abs, Real.norm_eq_abs, abs_pow]
  have e : Φ4 τ = -(2 ^ (0 + 2) * (1 - 1 / (3 : ℝ) ^ (0 + 1))) / ((0 + 2 : ℕ) : ℝ) * X05.z (0 + 2)
      * τ ^ (0 + 2) + -(2 ^ (1 + 2) * (1 - 1 / (3 : ℝ) ^ (1 + 1))) / ((1 + 2 : ℕ) : ℝ) *
      X05.z (1 + 2) * τ ^ (1 + 2) + -(2 ^ (2 + 2) * (1 - 1 / (3 : ℝ) ^ (2 + 1))) /
      ((2 + 2 : ℕ) : ℝ) * X05.z (2 + 2) * τ ^ (2 + 2) := by
    unfold Φ4; norm_num; ring
  rw [e]
  calc _ = |Real.log (Real.Gamma (1 - 2 * τ / 3) ^ 3 / Real.Gamma (1 - 2 * τ)) - _ - _ - _| := by
        congr 1; ring
    _ ≤ 2 * (64 * |τ| ^ 5) := ht
    _ = 128 * |τ| ^ 5 := by ring

/-- three-term `logGamma2` -/
lemma lg3 {u v R : ℝ} (hR : R ≤ 1 / 2) (hu : |u| ≤ R) (hv : |v| ≤ R) (huv : |u + v| ≤ R) :
    |Real.log (Real.Gamma (1 + u) * Real.Gamma (1 + v) / Real.Gamma (1 + (u + v)))
      - (-(X05.z 2) * (u * v) + X05.z 3 * (u * v * (u + v))
        + X05.z 4 / 4 * (u ^ 4 + v ^ 4 - (u + v) ^ 4))| ≤ 12 * R ^ 5 := by
  have hR0 : 0 ≤ R := (abs_nonneg u).trans hu
  have hs := logGamma2 (by linarith) hu hv huv
  have pw : ∀ x : ℝ, |x| ≤ R → ∀ n, |x ^ n| ≤ R ^ n := fun x hx n => by
    rw [abs_pow]; exact pow_le_pow_left₀ (abs_nonneg x) hx n
  have ht := tail3 hs hR0 hR (A := 6 * R ^ 5) (fun k => by
    have hz := z_le (show 2 ≤ k + 3 + 2 by omega)
    have hz0 := z_nn (show 2 ≤ k + 3 + 2 by omega)
    have hn : (1 : ℝ) ≤ ((k + 3 + 2 : ℕ) : ℝ) := by push_cast; linarith [k.cast_nonneg (α := ℝ)]
    have h3 : |u ^ (k + 3 + 2) + v ^ (k + 3 + 2) - (u + v) ^ (k + 3 + 2)| ≤
        3 * R ^ (k + 3 + 2) := by
      have a := pw u hu (k + 3 + 2); have b := pw v hv (k + 3 + 2)
      have c := pw _ huv (k + 3 + 2)
      have := abs_sub (u ^ (k + 3 + 2) + v ^ (k + 3 + 2)) ((u + v) ^ (k + 3 + 2))
      have := abs_add_le (u ^ (k + 3 + 2)) (v ^ (k + 3 + 2))
      linarith
    rw [abs_mul, abs_div, abs_mul, abs_pow, abs_neg, abs_one, one_pow, one_mul, abs_of_nonneg hz0,
      abs_of_pos (by positivity : (0 : ℝ) < ((k + 3 + 2 : ℕ) : ℝ))]
    calc X05.z (k + 3 + 2) / ((k + 3 + 2 : ℕ) : ℝ) *
          |u ^ (k + 3 + 2) + v ^ (k + 3 + 2) - (u + v) ^ (k + 3 + 2)|
        ≤ 2 / 1 * (3 * R ^ (k + 3 + 2)) := by gcongr
      _ = 6 * R ^ 5 * R ^ k := by ring)
  have e : -(X05.z 2) * (u * v) + X05.z 3 * (u * v * (u + v))
        + X05.z 4 / 4 * (u ^ 4 + v ^ 4 - (u + v) ^ 4) =
      (-1 : ℝ) ^ (0 + 2) * X05.z (0 + 2) / ((0 + 2 : ℕ) : ℝ) *
        (u ^ (0 + 2) + v ^ (0 + 2) - (u + v) ^ (0 + 2)) +
      (-1 : ℝ) ^ (1 + 2) * X05.z (1 + 2) / ((1 + 2 : ℕ) : ℝ) *
        (u ^ (1 + 2) + v ^ (1 + 2) - (u + v) ^ (1 + 2)) +
      (-1 : ℝ) ^ (2 + 2) * X05.z (2 + 2) / ((2 + 2 : ℕ) : ℝ) *
        (u ^ (2 + 2) + v ^ (2 + 2) - (u + v) ^ (2 + 2)) := by
    norm_num; ring
  rw [e]
  calc _ = |Real.log (Real.Gamma (1 + u) * Real.Gamma (1 + v) / Real.Gamma (1 + (u + v)))
        - _ - _ - _| := by congr 1; ring
    _ ≤ 2 * (6 * R ^ 5) := ht
    _ = 12 * R ^ 5 := by ring

open Asymptotics in
lemma lgO (a b : ℝ) : (fun τ : ℝ =>
    Real.log (Real.Gamma (1 + a * τ) * Real.Gamma (1 + b * τ) / Real.Gamma (1 + (a * τ + b * τ)))
      - (-(X05.z 2) * (a * τ * (b * τ)) + X05.z 3 * (a * τ * (b * τ) * (a * τ + b * τ))
        + X05.z 4 / 4 * ((a * τ) ^ 4 + (b * τ) ^ 4 - (a * τ + b * τ) ^ 4)))
      =O[nhds 0] fun τ => τ ^ 5 := by
  set M := |a| + |b| + 1 with hM
  have hM0 : 0 < M := by positivity
  refine IsBigO.of_bound (12 * M ^ 5) ?_
  filter_upwards [ball_ev (1 / (2 * M)) (by positivity)] with τ hτ
  have hR : M * |τ| ≤ 1 / 2 := by
    rw [lt_div_iff₀ (by positivity)] at hτ; nlinarith [abs_nonneg τ]
  have hu : |a * τ| ≤ M * |τ| := by
    rw [abs_mul]; exact mul_le_mul_of_nonneg_right (by rw [hM]; linarith [abs_nonneg b])
      (abs_nonneg τ)
  have hv : |b * τ| ≤ M * |τ| := by
    rw [abs_mul]; exact mul_le_mul_of_nonneg_right (by rw [hM]; linarith [abs_nonneg a])
      (abs_nonneg τ)
  have huv : |a * τ + b * τ| ≤ M * |τ| := by
    rw [← add_mul, abs_mul]
    exact mul_le_mul_of_nonneg_right (by rw [hM]; linarith [abs_add_le a b]) (abs_nonneg τ)
  have := lg3 hR hu hv huv
  rw [Real.norm_eq_abs, Real.norm_eq_abs, abs_pow]
  calc _ ≤ 12 * (M * |τ|) ^ 5 := this
    _ = 12 * M ^ 5 * |τ| ^ 5 := by ring

end Y05

namespace Y05

noncomputable def PY (τ : ℝ) : ℝ :=
  4 / 3 * X05.z 2 * τ ^ 2 + 32 / 27 * X05.z 3 * τ ^ 3 + 40 / 27 * X05.z 4 * τ ^ 4
noncomputable def PA (τ : ℝ) : ℝ := 4 * X05.z 2 * τ ^ 2 + 8 * X05.z 4 * τ ^ 4 - Φ4 τ
noncomputable def PB (τ : ℝ) : ℝ := -Φ4 τ
noncomputable def QQ (P : ℝ) : ℝ := P + P ^ 2 / 2
noncomputable def Yp (τ : ℝ) : ℝ :=
  QQ (PY τ) / 2 - (1 / 3 * QQ (PA τ) + 2 / 3 * QQ (PB τ) - 2 * QQ (PY τ)) / 8

noncomputable def W10 (Z τ : ℝ) : ℝ :=
  ((-4/1215 : ℝ) * Real.pi ^ 6 * Z ^ 0) * τ ^ 1 +
    ((-176/10935 : ℝ) * Real.pi ^ 4 * Z ^ 1) * τ ^ 2 +
    ((-128/2187 : ℝ) * Real.pi ^ 2 * Z ^ 2 + (-8/492075 : ℝ) * Real.pi ^ 8 * Z ^ 0) * τ ^ 3 +
    ((-2048/19683 : ℝ) * Real.pi ^ 0 * Z ^ 3 + (224/98415 : ℝ) * Real.pi ^ 6 * Z ^ 1) * τ ^ 4 +
    ((128/98415 : ℝ) * Real.pi ^ 4 * Z ^ 2 + (208/4428675 : ℝ) * Real.pi ^ 10 * Z ^ 0) * τ ^ 5 +
    ((-2048/177147 : ℝ) * Real.pi ^ 2 * Z ^ 3 + (5504/13286025 : ℝ) * Real.pi ^ 8 * Z ^ 1) * τ ^ 6 +
    ((-8192/531441 : ℝ) * Real.pi ^ 0 * Z ^ 4 + (512/531441 : ℝ) * Real.pi ^ 6 * Z ^ 2 + (-5696/1793613375 : ℝ) * Real.pi ^ 12 * Z ^ 0) * τ ^ 7 +
    ((28672/23914845 : ℝ) * Real.pi ^ 4 * Z ^ 3 + (256/119574225 : ℝ) * Real.pi ^ 10 * Z ^ 1) * τ ^ 8 +
    ((512/14348907 : ℝ) * Real.pi ^ 8 * Z ^ 2 + (-15872/16142520375 : ℝ) * Real.pi ^ 14 * Z ^ 0) * τ ^ 9 +
    ((-111104/48427561125 : ℝ) * Real.pi ^ 12 * Z ^ 1) * τ ^ 10 +
    ((-123008/2179240250625 : ℝ) * Real.pi ^ 16 * Z ^ 0) * τ ^ 11

lemma alg10 (τ : ℝ) :
    Φ4 τ + Yp τ - Yp τ ^ 2 / 2 - T09.logSeries (X05.z 3) τ = τ ^ 5 * W10 (X05.z 3) τ := by
  unfold Yp QQ PY PA PB Φ4 T09.logSeries W10
  rw [z2_eq, z4_eq]
  ring

lemma W10_cont (Z : ℝ) : Continuous (fun τ => W10 Z τ) := by
  unfold W10; fun_prop

open Asymptotics in
lemma O2_of {f w : ℝ → ℝ} (hw : Continuous w) (h : ∀ τ, f τ = τ ^ 2 * w τ) :
    f =O[nhds 0] fun τ : ℝ => τ ^ 2 :=
  (T09.bdd hw.continuousAt (fun τ => τ ^ 2)).congr_left (fun τ => (h τ).symm)

open Asymptotics in
lemma QO {P : ℝ → ℝ} (h : P =O[nhds 0] fun τ : ℝ => τ ^ 2) :
    (fun τ => QQ (P τ)) =O[nhds 0] fun τ : ℝ => τ ^ 2 := by
  have h2 : (fun τ => P τ ^ 2) =O[nhds 0] fun τ : ℝ => τ ^ 2 := by
    refine (h.pow 2).trans ?_
    exact (T09.hpow 2 2).congr_left (fun τ => by ring)
  unfold QQ
  exact h.add (h2.const_mul_left (1 / 2) |>.congr_left (fun τ => by ring))

open Asymptotics in
lemma PY_O : PY =O[nhds 0] fun τ : ℝ => τ ^ 2 :=
  O2_of (w := fun τ => 4 / 3 * X05.z 2 + 32 / 27 * X05.z 3 * τ + 40 / 27 * X05.z 4 * τ ^ 2)
    (by fun_prop) (fun τ => by unfold PY; ring)

open Asymptotics in
lemma PA_O : PA =O[nhds 0] fun τ : ℝ => τ ^ 2 :=
  O2_of (w := fun τ => 4 * X05.z 2 + 8 * X05.z 4 * τ ^ 2 +
      (4 / 3 * X05.z 2 + 64 / 27 * X05.z 3 * τ + 104 / 27 * X05.z 4 * τ ^ 2))
    (by fun_prop) (fun τ => by unfold PA Φ4; ring)

open Asymptotics in
lemma PB_O : PB =O[nhds 0] fun τ : ℝ => τ ^ 2 :=
  O2_of (w := fun τ => 4 / 3 * X05.z 2 + 64 / 27 * X05.z 3 * τ + 104 / 27 * X05.z 4 * τ ^ 2)
    (by fun_prop) (fun τ => by unfold PB Φ4; ring)

open Asymptotics in
lemma Yp_O : Yp =O[nhds 0] fun τ : ℝ => τ ^ 2 := by
  have a := QO PY_O; have b := QO PA_O; have c := QO PB_O
  have := ((a.const_mul_left (1 / 2)).sub ((((b.const_mul_left (1 / 3)).add
    (c.const_mul_left (2 / 3))).sub (a.const_mul_left 2)).const_mul_left (1 / 8)))
  refine this.congr_left (fun τ => ?_)
  unfold Yp; ring

open Asymptotics in
lemma expO {L P : ℝ → ℝ} (hP : P =O[nhds 0] fun τ : ℝ => τ ^ 2)
    (hLP : (fun τ => L τ - P τ) =O[nhds 0] fun τ : ℝ => τ ^ 5) :
    (fun τ => Real.exp (L τ) - (1 + QQ (P τ))) =O[nhds 0] fun τ : ℝ => τ ^ 5 := by
  have hy : L =O[nhds 0] fun τ : ℝ => τ ^ 2 :=
    (hP.add (hLP.trans (T09.hpow 2 3))).congr_left (fun τ => by ring)
  have := T09.expstep P (fun τ => L τ - P τ) L hy hP hLP (fun τ => by ring)
  refine this.congr_left (fun τ => ?_)
  unfold QQ; ring

open Asymptotics in
lemma LYO : (fun τ => Real.log (gammaE τ (Ymean τ)) - PY τ) =O[nhds 0] fun τ : ℝ => τ ^ 5 := by
  have h1 := lgO (2 / 3) (-4 / 3)
  have h2 := lgO (-2 / 3) (-2 / 3)
  refine (h1.sub h2).congr' ?_ EventuallyEq.rfl
  filter_upwards [ball_ev (1 / 8) (by norm_num)] with τ hτ
  obtain ⟨t1, t2⟩ := abs_lt.mp hτ
  show _ = Real.log (gammaE τ (Ymean τ)) - PY τ
  rw [Y05.proof_5 τ (by linarith)]
  rw [show 1 + 2 / 3 * τ = 1 + 2 * τ / 3 by ring, show 1 + -4 / 3 * τ = 1 - 4 * τ / 3 by ring,
    show 1 + (2 / 3 * τ + -4 / 3 * τ) = 1 - 2 * τ / 3 by ring,
    show 1 + -2 / 3 * τ = 1 - 2 * τ / 3 by ring,
    show 1 + (-2 / 3 * τ + -2 / 3 * τ) = 1 - 4 * τ / 3 by ring]
  have ga : 0 < Real.Gamma (1 - 2 * τ / 3) := Real.Gamma_pos_of_pos (by linarith)
  have gb : 0 < Real.Gamma (1 - 4 * τ / 3) := Real.Gamma_pos_of_pos (by linarith)
  have gc : 0 < Real.Gamma (1 + 2 * τ / 3) := Real.Gamma_pos_of_pos (by linarith)
  have hlog : Real.log (Real.Gamma (1 + 2 * τ / 3) * Real.Gamma (1 - 4 * τ / 3) ^ 2 /
      Real.Gamma (1 - 2 * τ / 3) ^ 3) =
      Real.log (Real.Gamma (1 + 2 * τ / 3) * Real.Gamma (1 - 4 * τ / 3) /
        Real.Gamma (1 - 2 * τ / 3)) -
      Real.log (Real.Gamma (1 - 2 * τ / 3) * Real.Gamma (1 - 2 * τ / 3) /
        Real.Gamma (1 - 4 * τ / 3)) := by
    rw [← Real.log_div (by positivity) (by positivity)]
    congr 1
    generalize Real.Gamma (1 - 2 * τ / 3) = a at *
    generalize Real.Gamma (1 - 4 * τ / 3) = b at *
    generalize Real.Gamma (1 + 2 * τ / 3) = c at *
    field_simp
  rw [hlog]
  unfold PY
  ring

noncomputable def Aval (τ : ℝ) : ℝ :=
  Real.Gamma (1 + 2 * τ) * Real.Gamma (1 - 2 * τ) ^ 2 / Real.Gamma (1 - 2 * τ / 3) ^ 3
noncomputable def Bval (τ : ℝ) : ℝ := Real.Gamma (1 - 2 * τ) / Real.Gamma (1 - 2 * τ / 3) ^ 3

open Asymptotics in
lemma LAO : (fun τ => Real.log (Aval τ) - PA τ) =O[nhds 0] fun τ : ℝ => τ ^ 5 := by
  have h1 := lgO 2 (-2)
  refine (h1.sub flatO).congr' ?_ EventuallyEq.rfl
  filter_upwards [ball_ev (1 / 8) (by norm_num)] with τ hτ
  obtain ⟨t1, t2⟩ := abs_lt.mp hτ
  show _ = Real.log (Aval τ) - PA τ
  rw [show 1 + -2 * τ = 1 - 2 * τ by ring, show 1 + (2 * τ + -2 * τ) = 1 by ring,
    Real.Gamma_one, div_one]
  have ga : 0 < Real.Gamma (1 - 2 * τ / 3) := Real.Gamma_pos_of_pos (by linarith)
  have gb : 0 < Real.Gamma (1 - 2 * τ) := Real.Gamma_pos_of_pos (by linarith)
  have gc : 0 < Real.Gamma (1 + 2 * τ) := Real.Gamma_pos_of_pos (by linarith)
  have hlog : Real.log (Aval τ) = Real.log (Real.Gamma (1 + 2 * τ) * Real.Gamma (1 - 2 * τ)) -
      Real.log (Real.Gamma (1 - 2 * τ / 3) ^ 3 / Real.Gamma (1 - 2 * τ)) := by
    rw [← Real.log_div (by positivity) (by positivity)]
    unfold Aval
    congr 1
    generalize Real.Gamma (1 - 2 * τ / 3) = a at *
    generalize Real.Gamma (1 - 2 * τ) = b at *
    field_simp
  rw [hlog]
  unfold PA
  ring

open Asymptotics in
lemma LBO : (fun τ => Real.log (Bval τ) - PB τ) =O[nhds 0] fun τ : ℝ => τ ^ 5 := by
  refine flatO.neg_left.congr_left (fun τ => ?_)
  unfold Bval PB
  rw [show Real.Gamma (1 - 2 * τ) / Real.Gamma (1 - 2 * τ / 3) ^ 3 =
    (Real.Gamma (1 - 2 * τ / 3) ^ 3 / Real.Gamma (1 - 2 * τ))⁻¹ by rw [inv_div], Real.log_inv]
  ring

open Asymptotics in
lemma EYO : (fun τ => gammaE τ (Ymean τ) - (1 + QQ (PY τ))) =O[nhds 0] fun τ : ℝ => τ ^ 5 := by
  refine (expO PY_O LYO).congr' ?_ EventuallyEq.rfl
  filter_upwards [ball_ev (1 / 8) (by norm_num)] with τ hτ
  obtain ⟨t1, t2⟩ := abs_lt.mp hτ
  have hV : 0 < gammaE τ (Ymean τ) := by
    rw [Y05.proof_5 τ (by linarith)]
    have ga : 0 < Real.Gamma (1 - 2 * τ / 3) := Real.Gamma_pos_of_pos (by linarith)
    have gb : 0 < Real.Gamma (1 - 4 * τ / 3) := Real.Gamma_pos_of_pos (by linarith)
    have gc : 0 < Real.Gamma (1 + 2 * τ / 3) := Real.Gamma_pos_of_pos (by linarith)
    positivity
  simp only [Real.exp_log hV]

open Asymptotics in
lemma AO : (fun τ => Aval τ - (1 + QQ (PA τ))) =O[nhds 0] fun τ : ℝ => τ ^ 5 := by
  refine (expO PA_O LAO).congr' ?_ EventuallyEq.rfl
  filter_upwards [ball_ev (1 / 8) (by norm_num)] with τ hτ
  obtain ⟨t1, t2⟩ := abs_lt.mp hτ
  have hV : 0 < Aval τ := by
    unfold Aval
    have ga : 0 < Real.Gamma (1 - 2 * τ / 3) := Real.Gamma_pos_of_pos (by linarith)
    have gb : 0 < Real.Gamma (1 - 2 * τ) := Real.Gamma_pos_of_pos (by linarith)
    have gc : 0 < Real.Gamma (1 + 2 * τ) := Real.Gamma_pos_of_pos (by linarith)
    positivity
  simp only [Real.exp_log hV]

open Asymptotics in
lemma BO : (fun τ => Bval τ - (1 + QQ (PB τ))) =O[nhds 0] fun τ : ℝ => τ ^ 5 := by
  refine (expO PB_O LBO).congr' ?_ EventuallyEq.rfl
  filter_upwards [ball_ev (1 / 8) (by norm_num)] with τ hτ
  obtain ⟨t1, t2⟩ := abs_lt.mp hτ
  have hV : 0 < Bval τ := by
    unfold Bval
    have ga : 0 < Real.Gamma (1 - 2 * τ / 3) := Real.Gamma_pos_of_pos (by linarith)
    have gb : 0 < Real.Gamma (1 - 2 * τ) := Real.Gamma_pos_of_pos (by linarith)
    positivity
  simp only [Real.exp_log hV]

lemma EX2_eq {τ : ℝ} (hτ : |τ| < 1 / 2) :
    gammaE τ (fun g => (Ymean τ g - 1) ^ 2) =
      1 / 3 * Aval τ + 2 / 3 * Bval τ - 2 * gammaE τ (Ymean τ) + 1 := by
  obtain ⟨i0, v0, iY⟩ := wt_pack hτ
  obtain ⟨t1, t2⟩ := abs_lt.mp hτ
  have iY2 : IntegrableOn (fun g => wt τ g * Ymean τ g ^ 2) Q := by
    have h1 : -1 < 2 * τ := by linarith
    have h2 : -1 < -(2 * τ) := by linarith
    have h0 : (-1 : ℝ) < 0 := by norm_num
    have := ((((mono_int h1 h2 h2).1.add (mono_int h2 h1 h2).1).add
      (mono_int h2 h2 h1).1).add ((((mono_int h0 h0 h2).1.add (mono_int h0 h2 h0).1).add
      (mono_int h2 h0 h0).1).const_mul 2)).const_mul (1 / 9)
    exact IntegrableOn.congr_fun this (fun g hg => (pt6 τ hg).symm) Q_meas
  have h6 := Y05.proof_6 τ hτ
  have hG : 0 < Real.Gamma (1 - 2 * τ / 3) := Real.Gamma_pos_of_pos (by linarith)
  have e2 : gammaE τ (fun g => Ymean τ g ^ 2) = 1 / 3 * Aval τ + 2 / 3 * Bval τ := by
    rw [h6]; unfold Aval Bval
    generalize Real.Gamma (1 - 2 * τ / 3) = G at *
    field_simp; ring
  rw [← e2]
  rw [gE_wt, gE_wt, gE_wt]
  have e : (fun g => wt τ g * (Ymean τ g - 1) ^ 2) =
      fun g => wt τ g * Ymean τ g ^ 2 - 2 * (wt τ g * Ymean τ g) + wt τ g := by
    funext g; ring
  have i2 := iY.const_mul 2
  have i3 : IntegrableOn (fun g => wt τ g * Ymean τ g ^ 2 - 2 * (wt τ g * Ymean τ g)) Q :=
    Integrable.sub iY2 i2
  have e1 := integral_add (μ := volume.restrict Q)
    (f := fun g => wt τ g * Ymean τ g ^ 2 - 2 * (wt τ g * Ymean τ g)) (g := fun g => wt τ g)
    i3 i0
  have e3 := integral_sub (μ := volume.restrict Q)
    (f := fun g => wt τ g * Ymean τ g ^ 2) (g := fun g => 2 * (wt τ g * Ymean τ g)) iY2 i2
  rw [e, e1, e3, integral_const_mul, v0]
  generalize Real.Gamma (1 - 2 * τ / 3) = G at *
  have hG' : G ^ 3 ≠ 0 := by positivity
  field_simp

end Y05

namespace Y05

open Asymptotics in
lemma h7O : (fun τ => B τ - (1 + gammaE τ (fun g => Ymean τ g - 1) / 2
    - gammaE τ (fun g => (Ymean τ g - 1) ^ 2) / 8)) =O[nhds 0] fun τ : ℝ => τ ^ 5 := by
  obtain ⟨C7, δ7, hδ7, h7⟩ := Y05.proof_7
  refine IsBigO.of_bound |C7| ?_
  filter_upwards [ball_ev (min δ7 1) (lt_min hδ7 one_pos)] with τ hτ
  have h1 : |τ| ≤ δ7 := (le_of_lt hτ).trans (min_le_left _ _)
  have h2 : |τ| ≤ 1 := (le_of_lt hτ).trans (min_le_right _ _)
  rw [Real.norm_eq_abs, Real.norm_eq_abs, abs_pow]
  calc _ ≤ C7 * |τ| ^ 6 := h7 τ h1
    _ ≤ |C7| * |τ| ^ 5 := by
      apply mul_le_mul (le_abs_self _) (pow_le_pow_of_le_one (abs_nonneg _) h2 (by norm_num))
        (by positivity) (abs_nonneg _)

open Asymptotics in
lemma yO : (fun τ => (B τ - 1) - Yp τ) =O[nhds 0] fun τ : ℝ => τ ^ 5 := by
  have := (h7O.add (EYO.const_mul_left (1 / 2))).sub
    ((((AO.const_mul_left (1 / 3)).add (BO.const_mul_left (2 / 3))).sub
      (EYO.const_mul_left 2)).const_mul_left (1 / 8))
  refine this.congr' ?_ EventuallyEq.rfl
  filter_upwards [ball_ev (1 / 8) (by norm_num)] with τ hτ
  have hτ' : |τ| < 1 / 2 := by linarith
  rw [EX_eq hτ', EX2_eq hτ']
  unfold Yp
  ring

open Asymptotics in
lemma yO2 : (fun τ => B τ - 1) =O[nhds 0] fun τ : ℝ => τ ^ 2 :=
  ((yO.trans (T09.hpow 2 3)).add Yp_O).congr_left (fun τ => by ring)

lemma log_quad {y : ℝ} (hy : |y| ≤ 1 / 2) :
    |Real.log (1 + y) - (y - y ^ 2 / 2)| ≤ 2 * |y| ^ 3 := by
  have := Real.abs_log_sub_add_sum_range_le (x := -y) (by rw [abs_neg]; linarith) 2
  have e : (∑ i ∈ Finset.range 2, (-y) ^ (i + 1) / ((i : ℝ) + 1)) + Real.log (1 - -y) =
      Real.log (1 + y) - (y - y ^ 2 / 2) := by
    simp [Finset.sum_range_succ]; ring
  rw [e, abs_neg] at this
  refine this.trans ?_
  rw [show |y| ^ (2 + 1) = |y| ^ 3 from rfl, div_le_iff₀ (by linarith)]
  have := pow_nonneg (abs_nonneg y) 3
  have := mul_le_mul_of_nonneg_left (show 1 ≤ 2 * (1 - |y|) by linarith) this
  linarith

open Asymptotics in
lemma logBO : (fun τ => Real.log (B τ) - ((B τ - 1) - (B τ - 1) ^ 2 / 2))
    =O[nhds 0] fun τ : ℝ => τ ^ 5 := by
  have ht : Tendsto (fun τ => B τ - 1) (nhds 0) (nhds 0) :=
    yO2.trans_tendsto (by simpa using ((continuous_pow 2).tendsto (0:ℝ)))
  have hev : ∀ᶠ τ in nhds (0 : ℝ), |B τ - 1| ≤ 1 / 2 := by
    filter_upwards [ht.eventually (Metric.closedBall_mem_nhds (0:ℝ) (by norm_num : (0:ℝ) < 1 / 2))]
      with τ hτ
    simpa [Real.norm_eq_abs] using hτ
  have h3 : (fun τ => Real.log (B τ) - ((B τ - 1) - (B τ - 1) ^ 2 / 2))
      =O[nhds 0] fun τ => (B τ - 1) ^ 3 := by
    refine IsBigO.of_bound 2 ?_
    filter_upwards [hev] with τ hτ
    have := log_quad hτ
    rw [show 1 + (B τ - 1) = B τ by ring] at this
    rw [Real.norm_eq_abs, Real.norm_eq_abs, abs_pow]
    exact this
  have h6 : (fun τ => (B τ - 1) ^ 3) =O[nhds 0] fun τ : ℝ => τ ^ 6 := by
    refine (yO2.pow 3).trans (isBigO_of_le _ (fun τ => ?_))
    simp only [Real.norm_eq_abs, ← pow_mul]
    simp
  exact h3.trans (h6.trans (T09.hpow 5 1))

open Asymptotics in
lemma quadO : (fun τ => ((B τ - 1) - (B τ - 1) ^ 2 / 2) - (Yp τ - Yp τ ^ 2 / 2))
    =O[nhds 0] fun τ : ℝ => τ ^ 5 := by
  have hb : (fun τ => 1 - ((B τ - 1) + Yp τ) / 2) =O[nhds 0] (1 : ℝ → ℝ) := by
    refine (isBigO_const_const (1:ℝ) one_ne_zero _).sub ?_
    exact (((yO2.add Yp_O).trans (T09.dpow_one 2 (by norm_num))).const_mul_left (1 / 2)).congr_left
      (fun τ => by ring)
  have := yO.mul hb
  simp only [Pi.one_apply, mul_one] at this
  refine this.congr_left (fun τ => ?_)
  ring

open Asymptotics in
lemma mainτ : (fun τ => Real.log (Real.Gamma (1 - 2 * τ / 3) ^ 3 / Real.Gamma (1 - 2 * τ) * B τ)
    - T09.logSeries (X05.z 3) τ) =O[nhds 0] fun τ : ℝ => τ ^ 5 := by
  have hW : (fun τ : ℝ => τ ^ 5 * W10 (X05.z 3) τ) =O[nhds 0] fun τ : ℝ => τ ^ 5 :=
    T09.bdd (W10_cont _).continuousAt _
  refine (((flatO.add logBO).add quadO).add hW).congr' ?_ EventuallyEq.rfl
  filter_upwards [ball_ev (1 / 8) (by norm_num)] with τ hτ
  obtain ⟨t1, t2⟩ := abs_lt.mp hτ
  have ga : 0 < Real.Gamma (1 - 2 * τ / 3) := Real.Gamma_pos_of_pos (by linarith)
  have gb : 0 < Real.Gamma (1 - 2 * τ) := Real.Gamma_pos_of_pos (by linarith)
  have hB := B_pos (show τ < 1 by linarith)
  have hA := alg10 τ
  rw [Real.log_mul (by positivity) hB.ne', ← hA]
  ring

lemma F_pos_ev : ∀ᶠ d in nhds (0 : ℝ), 0 < S (1 + d) / (4 * Real.sqrt 3) := by
  filter_upwards [ball_ev (1 / 16) (by norm_num)] with d hd
  obtain ⟨a, b⟩ := abs_lt.mp hd
  have hp : 0 < 1 + d := by linarith
  have hτ : |1 - 1 / (1 + d)| ≤ 1 / 8 := by
    rw [show 1 - 1 / (1 + d) = d / (1 + d) by field_simp; ring, abs_div, abs_of_pos hp,
      div_le_iff₀ hp]
    nlinarith [abs_nonneg d]
  have hfac := factor_id hp hτ
  obtain ⟨u1, u2⟩ := abs_le.mp hτ
  have g1 : 0 < Real.Gamma (1 - 2 * (1 - 1 / (1 + d))) := Real.Gamma_pos_of_pos (by linarith)
  have ga : 0 < Real.Gamma (1 - 2 * (1 - 1 / (1 + d)) / 3) :=
    Real.Gamma_pos_of_pos (by linarith)
  have hB := B_pos (show 1 - 1 / (1 + d) < 1 by linarith)
  have hpos' : 0 < (2 - (1 + d)) * S (1 + d) / (4 * Real.sqrt 3) := by
    rw [hfac]; positivity
  have : (2 - (1 + d)) * S (1 + d) / (4 * Real.sqrt 3) =
      (1 - d) * (S (1 + d) / (4 * Real.sqrt 3)) := by ring
  rw [this] at hpos'
  exact pos_of_mul_pos_right hpos' (by linarith)

end Y05

open Asymptotics in
theorem Y05.proof_10 : Y05.stmt_10 := by
  have hM := mainτ
  constructor
  · -- the `d`-series via T09
    have ht : Tendsto (fun d : ℝ => d / (1 + d)) (nhds 0) (nhds 0) := by
      have h : ContinuousAt (fun d : ℝ => d / (1 + d)) 0 :=
        ContinuousAt.div (f := fun d : ℝ => d) (g := fun d : ℝ => 1 + d) continuousAt_id
          (by fun_prop) (by norm_num)
      simpa using h.tendsto
    have hc := hM.comp_tendsto ht
    have hG : (fun d : ℝ => Real.log ((1 - d) * (S (1 + d) / (4 * Real.sqrt 3)))
        - T09.logSeries (X05.z 3) (d / (1 + d))) =O[nhds 0] fun d : ℝ => (d / (1 + d)) ^ 5 := by
      refine hc.congr' ?_ EventuallyEq.rfl
      filter_upwards [ball_ev (1 / 16) (by norm_num)] with d hd
      obtain ⟨a, b⟩ := abs_lt.mp hd
      have hp : 0 < 1 + d := by linarith
      have e : 1 - 1 / (1 + d) = d / (1 + d) := by field_simp; ring
      have hτ : |1 - 1 / (1 + d)| ≤ 1 / 8 := by
        rw [e, abs_div, abs_of_pos hp, div_le_iff₀ hp]
        nlinarith [abs_nonneg d]
      have hfac := factor_id hp hτ
      rw [e] at hfac
      simp only [Function.comp]
      rw [show (1 - d) * (S (1 + d) / (4 * Real.sqrt 3)) =
        (2 - (1 + d)) * S (1 + d) / (4 * Real.sqrt 3) by ring, hfac]
    have hT := T09.proof_1 (X05.z 3) (fun d => S (1 + d) / (4 * Real.sqrt 3)) F_pos_ev hG
    obtain ⟨c, hc⟩ := hT.bound
    obtain ⟨ε, hε, hεb⟩ := Metric.eventually_nhds_iff.mp hc
    refine ⟨c, ε / 2, by linarith, fun d hd => ?_⟩
    have := hεb (show dist d 0 < ε by rw [Real.dist_eq, sub_zero]; linarith)
    simpa [Real.norm_eq_abs, abs_pow] using this
  · obtain ⟨c, hc⟩ := hM.bound
    obtain ⟨ε, hε, hεb⟩ := Metric.eventually_nhds_iff.mp hc
    refine ⟨c, min (ε / 2) (1 / 8), lt_min (by linarith) (by norm_num), fun p hp hτ => ?_⟩
    have h8 : |1 - 1 / p| ≤ 1 / 8 := hτ.trans (min_le_right _ _)
    have hε2 : |1 - 1 / p| ≤ ε / 2 := hτ.trans (min_le_left _ _)
    rw [factor_id hp h8]
    have := hεb (y := 1 - 1 / p) (by rw [Real.dist_eq, sub_zero]; linarith)
    simpa [Real.norm_eq_abs, abs_pow] using this

/-! ## stmt_11: through `τ⁶` -/

namespace Y05

lemma bern6 : bernoulli' 6 = 1 / 42 := by
  have h5 : bernoulli' 5 = 0 := bernoulli'_eq_zero_of_odd ⟨2, rfl⟩ (by norm_num)
  have c1 : Nat.choose 6 2 = 15 := by decide
  have c2 : Nat.choose 6 3 = 20 := by decide
  have c3 : Nat.choose 6 4 = 15 := by decide
  rw [bernoulli'_def]
  norm_num [Finset.sum_range_succ, Finset.sum_range_zero, c1, c2, c3, h5]

lemma z6_eq : X05.z 6 = Real.pi ^ 6 / 945 := by
  unfold X05.z
  rw [show ((6 : ℕ) : ℂ) = 2 * ((3 : ℕ) : ℂ) by push_cast; ring,
    riemannZeta_two_mul_nat (by norm_num), bernoulli_eq_bernoulli'_of_ne_one (by decide), bern6]
  have : (-1 : ℂ) ^ (3 + 1) * (2 : ℂ) ^ (2 * 3 - 1) * (Real.pi : ℂ) ^ (2 * 3) *
      (((1 / 42 : ℚ)) : ℂ) / ((2 * 3).factorial : ℂ) = ((Real.pi ^ 6 / 945 : ℝ) : ℂ) := by
    push_cast; norm_num [Nat.factorial]; ring
  rw [this, Complex.ofReal_re]

noncomputable def ser (c2 c3 c4 c5 c6 τ : ℝ) : ℝ :=
  c2 * τ ^ 2 + c3 * τ ^ 3 + c4 * τ ^ 4 + c5 * τ ^ 5 + c6 * τ ^ 6

open Asymptotics in
lemma ser_O (c2 c3 c4 c5 c6 : ℝ) :
    (fun τ => ser c2 c3 c4 c5 c6 τ) =O[nhds 0] fun τ : ℝ => τ ^ 2 :=
  O2_of (w := fun τ => c2 + c3 * τ + c4 * τ ^ 2 + c5 * τ ^ 3 + c6 * τ ^ 4) (by fun_prop)
    (fun τ => by unfold ser; ring)

noncomputable def Q3 (P : ℝ) : ℝ := P + P ^ 2 / 2 + P ^ 3 / 6
noncomputable def L3 (P : ℝ) : ℝ := P - P ^ 2 / 2 + P ^ 3 / 3

noncomputable def poly12 (a0 a1 a2 a3 a4 a5 a6 a7 a8 a9 a10 a11 τ : ℝ) : ℝ :=
  a0 + a1 * τ + a2 * τ ^ 2 + a3 * τ ^ 3 + a4 * τ ^ 4 + a5 * τ ^ 5 + a6 * τ ^ 6 + a7 * τ ^ 7
    + a8 * τ ^ 8 + a9 * τ ^ 9 + a10 * τ ^ 10 + a11 * τ ^ 11

lemma poly12_cont (a0 a1 a2 a3 a4 a5 a6 a7 a8 a9 a10 a11 : ℝ) :
    Continuous (fun τ => poly12 a0 a1 a2 a3 a4 a5 a6 a7 a8 a9 a10 a11 τ) := by
  unfold poly12; fun_prop

noncomputable def expR (c2 c3 c4 c5 c6 τ : ℝ) : ℝ :=
  poly12
    ((1/1 : ℝ) * c3 ^ 1 * c4 ^ 1 + (1/1 : ℝ) * c2 ^ 1 * c5 ^ 1 + (1/2 : ℝ) * c2 ^ 2 * c3 ^ 1)
    ((1/2 : ℝ) * c4 ^ 2 + (1/1 : ℝ) * c3 ^ 1 * c5 ^ 1 + (1/1 : ℝ) * c2 ^ 1 * c6 ^ 1 + (1/2 : ℝ) * c2 ^ 1 * c3 ^ 2 + (1/2 : ℝ) * c2 ^ 2 * c4 ^ 1)
    ((1/1 : ℝ) * c4 ^ 1 * c5 ^ 1 + (1/1 : ℝ) * c3 ^ 1 * c6 ^ 1 + (1/6 : ℝ) * c3 ^ 3 + (1/1 : ℝ) * c2 ^ 1 * c3 ^ 1 * c4 ^ 1 + (1/2 : ℝ) * c2 ^ 2 * c5 ^ 1)
    ((1/2 : ℝ) * c5 ^ 2 + (1/1 : ℝ) * c4 ^ 1 * c6 ^ 1 + (1/2 : ℝ) * c3 ^ 2 * c4 ^ 1 + (1/2 : ℝ) * c2 ^ 1 * c4 ^ 2 + (1/1 : ℝ) * c2 ^ 1 * c3 ^ 1 * c5 ^ 1 + (1/2 : ℝ) * c2 ^ 2 * c6 ^ 1)
    ((1/1 : ℝ) * c5 ^ 1 * c6 ^ 1 + (1/2 : ℝ) * c3 ^ 1 * c4 ^ 2 + (1/2 : ℝ) * c3 ^ 2 * c5 ^ 1 + (1/1 : ℝ) * c2 ^ 1 * c4 ^ 1 * c5 ^ 1 + (1/1 : ℝ) * c2 ^ 1 * c3 ^ 1 * c6 ^ 1)
    ((1/2 : ℝ) * c6 ^ 2 + (1/6 : ℝ) * c4 ^ 3 + (1/1 : ℝ) * c3 ^ 1 * c4 ^ 1 * c5 ^ 1 + (1/2 : ℝ) * c3 ^ 2 * c6 ^ 1 + (1/2 : ℝ) * c2 ^ 1 * c5 ^ 2 + (1/1 : ℝ) * c2 ^ 1 * c4 ^ 1 * c6 ^ 1)
    ((1/2 : ℝ) * c4 ^ 2 * c5 ^ 1 + (1/2 : ℝ) * c3 ^ 1 * c5 ^ 2 + (1/1 : ℝ) * c3 ^ 1 * c4 ^ 1 * c6 ^ 1 + (1/1 : ℝ) * c2 ^ 1 * c5 ^ 1 * c6 ^ 1)
    ((1/2 : ℝ) * c4 ^ 1 * c5 ^ 2 + (1/2 : ℝ) * c4 ^ 2 * c6 ^ 1 + (1/1 : ℝ) * c3 ^ 1 * c5 ^ 1 * c6 ^ 1 + (1/2 : ℝ) * c2 ^ 1 * c6 ^ 2)
    ((1/6 : ℝ) * c5 ^ 3 + (1/1 : ℝ) * c4 ^ 1 * c5 ^ 1 * c6 ^ 1 + (1/2 : ℝ) * c3 ^ 1 * c6 ^ 2)
    ((1/2 : ℝ) * c5 ^ 2 * c6 ^ 1 + (1/2 : ℝ) * c4 ^ 1 * c6 ^ 2)
    ((1/2 : ℝ) * c5 ^ 1 * c6 ^ 2)
    ((1/6 : ℝ) * c6 ^ 3)
    τ

noncomputable def logR (c2 c3 c4 c5 c6 τ : ℝ) : ℝ :=
  poly12
    ((-1/1 : ℝ) * c3 ^ 1 * c4 ^ 1 + (-1/1 : ℝ) * c2 ^ 1 * c5 ^ 1 + (1/1 : ℝ) * c2 ^ 2 * c3 ^ 1)
    ((-1/2 : ℝ) * c4 ^ 2 + (-1/1 : ℝ) * c3 ^ 1 * c5 ^ 1 + (-1/1 : ℝ) * c2 ^ 1 * c6 ^ 1 + (1/1 : ℝ) * c2 ^ 1 * c3 ^ 2 + (1/1 : ℝ) * c2 ^ 2 * c4 ^ 1)
    ((-1/1 : ℝ) * c4 ^ 1 * c5 ^ 1 + (-1/1 : ℝ) * c3 ^ 1 * c6 ^ 1 + (1/3 : ℝ) * c3 ^ 3 + (2/1 : ℝ) * c2 ^ 1 * c3 ^ 1 * c4 ^ 1 + (1/1 : ℝ) * c2 ^ 2 * c5 ^ 1)
    ((-1/2 : ℝ) * c5 ^ 2 + (-1/1 : ℝ) * c4 ^ 1 * c6 ^ 1 + (1/1 : ℝ) * c3 ^ 2 * c4 ^ 1 + (1/1 : ℝ) * c2 ^ 1 * c4 ^ 2 + (2/1 : ℝ) * c2 ^ 1 * c3 ^ 1 * c5 ^ 1 + (1/1 : ℝ) * c2 ^ 2 * c6 ^ 1)
    ((-1/1 : ℝ) * c5 ^ 1 * c6 ^ 1 + (1/1 : ℝ) * c3 ^ 1 * c4 ^ 2 + (1/1 : ℝ) * c3 ^ 2 * c5 ^ 1 + (2/1 : ℝ) * c2 ^ 1 * c4 ^ 1 * c5 ^ 1 + (2/1 : ℝ) * c2 ^ 1 * c3 ^ 1 * c6 ^ 1)
    ((-1/2 : ℝ) * c6 ^ 2 + (1/3 : ℝ) * c4 ^ 3 + (2/1 : ℝ) * c3 ^ 1 * c4 ^ 1 * c5 ^ 1 + (1/1 : ℝ) * c3 ^ 2 * c6 ^ 1 + (1/1 : ℝ) * c2 ^ 1 * c5 ^ 2 + (2/1 : ℝ) * c2 ^ 1 * c4 ^ 1 * c6 ^ 1)
    ((1/1 : ℝ) * c4 ^ 2 * c5 ^ 1 + (1/1 : ℝ) * c3 ^ 1 * c5 ^ 2 + (2/1 : ℝ) * c3 ^ 1 * c4 ^ 1 * c6 ^ 1 + (2/1 : ℝ) * c2 ^ 1 * c5 ^ 1 * c6 ^ 1)
    ((1/1 : ℝ) * c4 ^ 1 * c5 ^ 2 + (1/1 : ℝ) * c4 ^ 2 * c6 ^ 1 + (2/1 : ℝ) * c3 ^ 1 * c5 ^ 1 * c6 ^ 1 + (1/1 : ℝ) * c2 ^ 1 * c6 ^ 2)
    ((1/3 : ℝ) * c5 ^ 3 + (2/1 : ℝ) * c4 ^ 1 * c5 ^ 1 * c6 ^ 1 + (1/1 : ℝ) * c3 ^ 1 * c6 ^ 2)
    ((1/1 : ℝ) * c5 ^ 2 * c6 ^ 1 + (1/1 : ℝ) * c4 ^ 1 * c6 ^ 2)
    ((1/1 : ℝ) * c5 ^ 1 * c6 ^ 2)
    ((1/3 : ℝ) * c6 ^ 3)
    τ

set_option maxHeartbeats 1000000 in
open Asymptotics in
lemma expT_O (c2 c3 c4 c5 c6 : ℝ) :
    (fun τ => Q3 (ser c2 c3 c4 c5 c6 τ) - ser c2 c3 (c4 + c2 ^ 2 / 2) (c5 + c2 * c3)
      (c6 + c3 ^ 2 / 2 + c2 * c4 + c2 ^ 3 / 6) τ) =O[nhds 0] fun τ : ℝ => τ ^ 7 := by
  have hc : Continuous (fun τ => expR c2 c3 c4 c5 c6 τ) := poly12_cont _ _ _ _ _ _ _ _ _ _ _ _
  refine (T09.bdd hc.continuousAt (fun τ : ℝ => τ ^ 7)).congr_left (fun τ => ?_)
  unfold Q3 ser expR poly12; ring

set_option maxHeartbeats 1000000 in
open Asymptotics in
lemma logT_O (c2 c3 c4 c5 c6 : ℝ) :
    (fun τ => L3 (ser c2 c3 c4 c5 c6 τ) - ser c2 c3 (c4 - c2 ^ 2 / 2) (c5 - c2 * c3)
      (c6 - c3 ^ 2 / 2 - c2 * c4 + c2 ^ 3 / 3) τ) =O[nhds 0] fun τ : ℝ => τ ^ 7 := by
  have hc : Continuous (fun τ => logR c2 c3 c4 c5 c6 τ) := poly12_cont _ _ _ _ _ _ _ _ _ _ _ _
  refine (T09.bdd hc.continuousAt (fun τ : ℝ => τ ^ 7)).congr_left (fun τ => ?_)
  unfold L3 ser logR poly12; ring

open Asymptotics in
lemma O1_of2 {f : ℝ → ℝ} (h : f =O[nhds 0] fun τ : ℝ => τ ^ 2) : f =O[nhds 0] (1 : ℝ → ℝ) :=
  h.trans (T09.dpow_one 2 (by norm_num))

open Asymptotics in
lemma expO3 {L P : ℝ → ℝ} (hP : P =O[nhds 0] fun τ : ℝ => τ ^ 2)
    (hLP : (fun τ => L τ - P τ) =O[nhds 0] fun τ : ℝ => τ ^ 7) :
    (fun τ => Real.exp (L τ) - (1 + Q3 (P τ))) =O[nhds 0] fun τ : ℝ => τ ^ 7 := by
  have hy : L =O[nhds 0] fun τ : ℝ => τ ^ 2 :=
    (hP.add (hLP.trans (T09.hpow 2 5))).congr_left (fun τ => by ring)
  have hy0 : Tendsto L (nhds 0) (nhds 0) :=
    hy.trans_tendsto (by simpa using ((continuous_pow 2).tendsto (0:ℝ)))
  have hy1 : ∀ᶠ τ in nhds (0:ℝ), |L τ| ≤ 1 := by
    filter_upwards [hy0.eventually (Metric.closedBall_mem_nhds (0:ℝ) one_pos)] with τ hτ
    simpa [Real.norm_eq_abs] using hτ
  have hsum : ∀ x : ℝ, ∑ i ∈ Finset.range 4, x ^ i / (i.factorial : ℝ) =
      1 + x + x ^ 2 / 2 + x ^ 3 / 6 := fun x => by
    simp [Finset.sum_range_succ, Nat.factorial]
  have h4 : (fun τ => Real.exp (L τ) - (1 + Q3 (L τ))) =O[nhds 0] fun τ => L τ ^ 4 := by
    refine IsBigO.of_bound 1 ?_
    filter_upwards [hy1] with τ hτ
    have := Real.exp_bound hτ (n := 4) (by norm_num)
    rw [hsum] at this
    rw [Real.norm_eq_abs, Real.norm_eq_abs, abs_pow]
    have e : (1 + Q3 (L τ)) = 1 + L τ + L τ ^ 2 / 2 + L τ ^ 3 / 6 := by unfold Q3; ring
    rw [e]
    have hc : ((Nat.succ 4 : ℕ) : ℝ) / ((Nat.factorial 4 : ℕ) * (4 : ℕ)) ≤ 1 := by
      norm_num [Nat.factorial]
    calc _ ≤ |L τ| ^ 4 * (((Nat.succ 4 : ℕ) : ℝ) / ((Nat.factorial 4 : ℕ) * (4 : ℕ))) := this
      _ ≤ |L τ| ^ 4 * 1 := mul_le_mul_of_nonneg_left hc (by positivity)
      _ = 1 * |L τ| ^ 4 := by ring
  have h8 : (fun τ => L τ ^ 4) =O[nhds 0] fun τ : ℝ => τ ^ 8 := by
    refine (hy.pow 4).trans (isBigO_of_le _ (fun τ => ?_))
    simp only [Real.norm_eq_abs, ← pow_mul]
    simp
  have hP0 : Tendsto P (nhds 0) (nhds 0) :=
    hP.trans_tendsto (by simpa using ((continuous_pow 2).tendsto (0:ℝ)))
  have hb : (fun τ => 1 + (L τ + P τ) / 2 + (L τ ^ 2 + L τ * P τ + P τ ^ 2) / 6)
      =O[nhds 0] (1 : ℝ → ℝ) := by
    have : Tendsto (fun τ => 1 + (L τ + P τ) / 2 + (L τ ^ 2 + L τ * P τ + P τ ^ 2) / 6) (nhds 0)
        (nhds (1 + (0 + 0) / 2 + (0 ^ 2 + 0 * 0 + 0 ^ 2) / 6)) :=
      ((tendsto_const_nhds.add ((hy0.add hP0).div_const 2)).add
        ((((hy0.pow 2).add (hy0.mul hP0)).add (hP0.pow 2)).div_const 6))
    exact this.isBigO_one ℝ
  have hd := hLP.mul hb
  simp only [Pi.one_apply, mul_one] at hd
  refine ((h4.trans (h8.trans (T09.hpow 7 1))).add hd).congr_left (fun τ => ?_)
  unfold Q3; ring

end Y05

namespace Y05

lemma tail5 {f : ℕ → ℝ} {L A r : ℝ} (hf : HasSum f L) (hr0 : 0 ≤ r) (hr : r ≤ 1 / 2)
    (hb : ∀ k, |f (k + 5)| ≤ A * r ^ k) : |L - (f 0 + f 1 + f 2 + f 3 + f 4)| ≤ 2 * A := by
  have hs : HasSum (fun k => f (k + 5)) (L - (f 0 + f 1 + f 2 + f 3 + f 4)) := by
    have := (hasSum_nat_add_iff' 5).mpr hf
    simpa [Finset.sum_range_succ] using this
  have hA : 0 ≤ A := by
    have := (abs_nonneg _).trans (hb 0); simpa using this
  have hg : HasSum (fun k => A * r ^ k) (A * (1 - r)⁻¹) :=
    (hasSum_geometric_of_lt_one hr0 (by linarith)).mul_left A
  have hup := hasSum_le (fun k => (le_abs_self _).trans (hb k)) hs hg
  have hlo := hasSum_le (fun k => neg_le_of_abs_le (hb k)) hg.neg hs
  have hinv : (1 - r)⁻¹ ≤ 2 := by
    rw [inv_le_comm₀ (by linarith) (by norm_num)]; linarith
  have : A * (1 - r)⁻¹ ≤ 2 * A := by nlinarith
  rw [abs_le]
  constructor <;> linarith

/-- flat-face coefficients -/
noncomputable def fl (k : ℕ) : ℝ :=
  -(2 ^ (k + 2) * (1 - 1 / (3 : ℝ) ^ (k + 1))) / ((k + 2 : ℕ) : ℝ) * X05.z (k + 2)

/-- two-factor log-Gamma coefficients -/
noncomputable def lgc (a b : ℝ) (k : ℕ) : ℝ :=
  (-1 : ℝ) ^ (k + 2) * X05.z (k + 2) / ((k + 2 : ℕ) : ℝ) * (a ^ (k + 2) + b ^ (k + 2) - (a + b) ^ (k + 2))

open Asymptotics in
lemma flat6 : (fun τ : ℝ => Real.log (Real.Gamma (1 - 2 * τ / 3) ^ 3 / Real.Gamma (1 - 2 * τ))
    - ser (fl 0) (fl 1) (fl 2) (fl 3) (fl 4) τ) =O[nhds 0] fun τ => τ ^ 7 := by
  refine IsBigO.of_bound 512 ?_
  filter_upwards [ball_ev (1 / 4) (by norm_num)] with τ hτ
  have hs := Y05.proof_4 τ (by linarith)
  have ht := tail5 hs (by positivity : (0 : ℝ) ≤ 2 * |τ|) (by linarith) (A := 256 * |τ| ^ 7)
    (fun k => by
    have hz := z_le (show 2 ≤ k + 5 + 2 by omega)
    have hz0 := z_nn (show 2 ≤ k + 5 + 2 by omega)
    have h3 : 0 ≤ 1 - 1 / (3 : ℝ) ^ (k + 5 + 1) := by
      have : 1 / (3 : ℝ) ^ (k + 5 + 1) ≤ 1 := by
        rw [div_le_one (by positivity)]; exact one_le_pow₀ (by norm_num)
      linarith
    have h3' : 1 - 1 / (3 : ℝ) ^ (k + 5 + 1) ≤ 1 := by
      have : 0 ≤ 1 / (3 : ℝ) ^ (k + 5 + 1) := by positivity
      linarith
    have hn : (1 : ℝ) ≤ ((k + 5 + 2 : ℕ) : ℝ) := by push_cast; linarith [k.cast_nonneg (α := ℝ)]
    rw [abs_mul, abs_mul, abs_div, abs_neg, abs_of_nonneg (mul_nonneg (by positivity) h3),
      abs_of_pos (by positivity : (0 : ℝ) < ((k + 5 + 2 : ℕ) : ℝ)), abs_of_nonneg hz0, abs_pow]
    calc 2 ^ (k + 5 + 2) * (1 - 1 / 3 ^ (k + 5 + 1)) / ((k + 5 + 2 : ℕ) : ℝ) * X05.z (k + 5 + 2)
          * |τ| ^ (k + 5 + 2)
        ≤ 2 ^ (k + 5 + 2) * 1 / 1 * 2 * |τ| ^ (k + 5 + 2) := by gcongr
      _ = 256 * |τ| ^ 7 * (2 * |τ|) ^ k := by ring)
  rw [Real.norm_eq_abs, Real.norm_eq_abs, abs_pow]
  have e : ser (fl 0) (fl 1) (fl 2) (fl 3) (fl 4) τ =
      fl 0 * τ ^ (0 + 2) + fl 1 * τ ^ (1 + 2) + fl 2 * τ ^ (2 + 2) + fl 3 * τ ^ (3 + 2)
        + fl 4 * τ ^ (4 + 2) := by unfold ser; norm_num
  rw [e]
  unfold fl
  calc _ ≤ 2 * (256 * |τ| ^ 7) := ht
    _ = 512 * |τ| ^ 7 := by ring

lemma lg5 {u v R : ℝ} (hR : R ≤ 1 / 2) (hu : |u| ≤ R) (hv : |v| ≤ R) (huv : |u + v| ≤ R) :
    |Real.log (Real.Gamma (1 + u) * Real.Gamma (1 + v) / Real.Gamma (1 + (u + v)))
      - (∑ k ∈ Finset.range 5, (-1 : ℝ) ^ (k + 2) * X05.z (k + 2) / ((k + 2 : ℕ) : ℝ) *
        (u ^ (k + 2) + v ^ (k + 2) - (u + v) ^ (k + 2)))| ≤ 12 * R ^ 7 := by
  have hR0 : 0 ≤ R := (abs_nonneg u).trans hu
  have hs := logGamma2 (by linarith) hu hv huv
  have pw : ∀ x : ℝ, |x| ≤ R → ∀ n, |x ^ n| ≤ R ^ n := fun x hx n => by
    rw [abs_pow]; exact pow_le_pow_left₀ (abs_nonneg x) hx n
  have ht := tail5 hs hR0 hR (A := 6 * R ^ 7) (fun k => by
    have hz := z_le (show 2 ≤ k + 5 + 2 by omega)
    have hz0 := z_nn (show 2 ≤ k + 5 + 2 by omega)
    have hn : (1 : ℝ) ≤ ((k + 5 + 2 : ℕ) : ℝ) := by push_cast; linarith [k.cast_nonneg (α := ℝ)]
    have h3 : |u ^ (k + 5 + 2) + v ^ (k + 5 + 2) - (u + v) ^ (k + 5 + 2)| ≤
        3 * R ^ (k + 5 + 2) := by
      have a := pw u hu (k + 5 + 2); have b := pw v hv (k + 5 + 2)
      have c := pw _ huv (k + 5 + 2)
      have := abs_sub (u ^ (k + 5 + 2) + v ^ (k + 5 + 2)) ((u + v) ^ (k + 5 + 2))
      have := abs_add_le (u ^ (k + 5 + 2)) (v ^ (k + 5 + 2))
      linarith
    rw [abs_mul, abs_div, abs_mul, abs_pow, abs_neg, abs_one, one_pow, one_mul, abs_of_nonneg hz0,
      abs_of_pos (by positivity : (0 : ℝ) < ((k + 5 + 2 : ℕ) : ℝ))]
    calc X05.z (k + 5 + 2) / ((k + 5 + 2 : ℕ) : ℝ) *
          |u ^ (k + 5 + 2) + v ^ (k + 5 + 2) - (u + v) ^ (k + 5 + 2)|
        ≤ 2 / 1 * (3 * R ^ (k + 5 + 2)) := by gcongr
      _ = 6 * R ^ 7 * R ^ k := by ring)
  simp only [Finset.sum_range_succ, Finset.sum_range_zero, zero_add]
  calc _ ≤ 2 * (6 * R ^ 7) := ht
    _ = 12 * R ^ 7 := by ring

open Asymptotics in
lemma lgO6 (a b : ℝ) : (fun τ : ℝ =>
    Real.log (Real.Gamma (1 + a * τ) * Real.Gamma (1 + b * τ) / Real.Gamma (1 + (a * τ + b * τ)))
      - ser (lgc a b 0) (lgc a b 1) (lgc a b 2) (lgc a b 3) (lgc a b 4) τ)
      =O[nhds 0] fun τ => τ ^ 7 := by
  set M := |a| + |b| + 1 with hM
  have hM0 : 0 < M := by positivity
  refine IsBigO.of_bound (12 * M ^ 7) ?_
  filter_upwards [ball_ev (1 / (2 * M)) (by positivity)] with τ hτ
  have hR : M * |τ| ≤ 1 / 2 := by
    rw [lt_div_iff₀ (by positivity)] at hτ; nlinarith [abs_nonneg τ]
  have hu : |a * τ| ≤ M * |τ| := by
    rw [abs_mul]; exact mul_le_mul_of_nonneg_right (by rw [hM]; linarith [abs_nonneg b])
      (abs_nonneg τ)
  have hv : |b * τ| ≤ M * |τ| := by
    rw [abs_mul]; exact mul_le_mul_of_nonneg_right (by rw [hM]; linarith [abs_nonneg a])
      (abs_nonneg τ)
  have huv : |a * τ + b * τ| ≤ M * |τ| := by
    rw [← add_mul, abs_mul]
    exact mul_le_mul_of_nonneg_right (by rw [hM]; linarith [abs_add_le a b]) (abs_nonneg τ)
  have := lg5 hR hu hv huv
  have e : ser (lgc a b 0) (lgc a b 1) (lgc a b 2) (lgc a b 3) (lgc a b 4) τ =
      ∑ k ∈ Finset.range 5, (-1 : ℝ) ^ (k + 2) * X05.z (k + 2) / ((k + 2 : ℕ) : ℝ) *
        ((a * τ) ^ (k + 2) + (b * τ) ^ (k + 2) - (a * τ + b * τ) ^ (k + 2)) := by
    simp only [Finset.sum_range_succ, Finset.sum_range_zero, zero_add, ser, lgc]
    ring
  rw [Real.norm_eq_abs, Real.norm_eq_abs, abs_pow, e]
  calc _ ≤ 12 * (M * |τ|) ^ 7 := this
    _ = 12 * M ^ 7 * |τ| ^ 7 := by ring

end Y05

namespace Y05

lemma pt9 (τ : ℝ) {g : ℝ × ℝ × ℝ} (hg : g ∈ Q) :
    (g.1 * g.2.1 * g.2.2) ^ (-2 * τ / 3) * Real.exp (-(g.1 + g.2.1 + g.2.2)) * Ymean τ g ^ 3 =
      1 / 27 * (mono (10 * τ / 3) (-(8 * τ / 3)) (-(8 * τ / 3)) g
        + mono (-(8 * τ / 3)) (10 * τ / 3) (-(8 * τ / 3)) g
        + mono (-(8 * τ / 3)) (-(8 * τ / 3)) (10 * τ / 3) g
        + 3 * (mono (4 * τ / 3) (-(2 * τ / 3)) (-(8 * τ / 3)) g
          + mono (4 * τ / 3) (-(8 * τ / 3)) (-(2 * τ / 3)) g
          + mono (-(2 * τ / 3)) (4 * τ / 3) (-(8 * τ / 3)) g
          + mono (-(8 * τ / 3)) (4 * τ / 3) (-(2 * τ / 3)) g
          + mono (-(2 * τ / 3)) (-(8 * τ / 3)) (4 * τ / 3) g
          + mono (-(8 * τ / 3)) (-(2 * τ / 3)) (4 * τ / 3) g)
        + 6 * mono (-(2 * τ / 3)) (-(2 * τ / 3)) (-(2 * τ / 3)) g) := by
  obtain ⟨ha, hb, hc⟩ := hg
  simp only [mem_Ioi] at ha hb hc
  unfold Ymean mono
  set c := 2 * τ / 3 with hcdef
  have e1 : -2 * τ / 3 = -(c * (1 : ℕ)) := by rw [hcdef]; push_cast; ring
  have e2 : 2 * τ = c * (3 : ℕ) := by rw [hcdef]; push_cast; ring
  have e5 : 10 * τ / 3 = c * (5 : ℕ) := by rw [hcdef]; push_cast; ring
  have e4 : -(8 * τ / 3) = -(c * (4 : ℕ)) := by rw [hcdef]; push_cast; ring
  have e3 : 4 * τ / 3 = c * (2 : ℕ) := by rw [hcdef]; push_cast; ring
  have e6 : -c = -(c * (1 : ℕ)) := by push_cast; ring
  rw [Real.mul_rpow (mul_pos ha hb).le hc.le, Real.mul_rpow ha.le hb.le,
    Real.mul_rpow (mul_pos ha hb).le hc.le, Real.mul_rpow ha.le hb.le, e1, e2, e5, e4, e3, e6]
  simp only [rp_pow ha, rp_pow hb, rp_pow hc, rp_neg ha, rp_neg hb, rp_neg hc]
  have qa := Real.rpow_pos_of_pos ha c
  have qb := Real.rpow_pos_of_pos hb c
  have qc := Real.rpow_pos_of_pos hc c
  generalize hA : g.1 ^ c = A at *
  generalize hB : g.2.1 ^ c = B' at *
  generalize hC : g.2.2 ^ c = C at *
  field_simp
  ring

lemma iadd {f h : ℝ × ℝ × ℝ → ℝ} {a b : ℝ} (hf : IntegrableOn f Q ∧ ∫ x in Q, f x = a)
    (hh : IntegrableOn h Q ∧ ∫ x in Q, h x = b) :
    IntegrableOn (fun x => f x + h x) Q ∧ ∫ x in Q, (f x + h x) = a + b :=
  ⟨hf.1.add hh.1, by rw [integral_add hf.1 hh.1, hf.2, hh.2]⟩

lemma ismul {f : ℝ × ℝ × ℝ → ℝ} {a : ℝ} (r : ℝ) (hf : IntegrableOn f Q ∧ ∫ x in Q, f x = a) :
    IntegrableOn (fun x => r * f x) Q ∧ ∫ x in Q, r * f x = r * a :=
  ⟨hf.1.const_mul r, by rw [integral_const_mul, hf.2]⟩

noncomputable def C1val (τ : ℝ) : ℝ :=
  Real.Gamma (1 + 10 * τ / 3) * Real.Gamma (1 - 8 * τ / 3) ^ 2 / Real.Gamma (1 - 2 * τ / 3) ^ 3
noncomputable def C2val (τ : ℝ) : ℝ :=
  Real.Gamma (1 + 4 * τ / 3) * Real.Gamma (1 - 2 * τ / 3) * Real.Gamma (1 - 8 * τ / 3) /
    Real.Gamma (1 - 2 * τ / 3) ^ 3

lemma Y3_pack {τ : ℝ} (hτ : |τ| < 1 / 4) :
    IntegrableOn (fun g => wt τ g * Ymean τ g ^ 3) Q ∧
      gammaE τ (fun g => Ymean τ g ^ 3) = (3 * C1val τ + 18 * C2val τ + 6) / 27 := by
  obtain ⟨t1, t2⟩ := abs_lt.mp hτ
  have a5 : -1 < 10 * τ / 3 := by linarith
  have a4 : -1 < -(8 * τ / 3) := by linarith
  have a2 : -1 < 4 * τ / 3 := by linarith
  have a1 : -1 < -(2 * τ / 3) := by linarith
  have H := ismul (1 / 27) (iadd (iadd (iadd (iadd (mono_int a5 a4 a4) (mono_int a4 a5 a4))
    (mono_int a4 a4 a5)) (ismul 3 (iadd (iadd (iadd (iadd (iadd (mono_int a2 a1 a4)
      (mono_int a2 a4 a1)) (mono_int a1 a2 a4)) (mono_int a4 a2 a1)) (mono_int a1 a4 a2))
      (mono_int a4 a1 a2)))) (ismul 6 (mono_int a1 a1 a1)))
  have hc : ∀ g ∈ Q, wt τ g * Ymean τ g ^ 3 = _ := fun g hg => pt9 τ hg
  refine ⟨IntegrableOn.congr_fun H.1 (fun g hg => (hc g hg).symm) Q_meas, ?_⟩
  rw [gE_wt, setIntegral_congr_fun Q_meas hc, H.2]
  have ga : 0 < Real.Gamma (1 - 2 * τ / 3) := Real.Gamma_pos_of_pos (by linarith)
  unfold C1val C2val
  rw [show 1 + -(8 * τ / 3) = 1 - 8 * τ / 3 by ring, show 1 + -(2 * τ / 3) = 1 - 2 * τ / 3 by ring]
  generalize Real.Gamma (1 - 2 * τ / 3) = G at *
  field_simp
  ring

lemma EX3_eq {τ : ℝ} (hτ : |τ| < 1 / 4) :
    gammaE τ (fun g => (Ymean τ g - 1) ^ 3) =
      (3 * C1val τ + 18 * C2val τ + 6) / 27 - 3 * (1 / 3 * Aval τ + 2 / 3 * Bval τ)
        + 3 * gammaE τ (Ymean τ) - 1 := by
  have hτ2 : |τ| < 1 / 2 := by linarith
  obtain ⟨i0, v0, iY⟩ := wt_pack hτ2
  obtain ⟨iY3, v3⟩ := Y3_pack hτ
  obtain ⟨t1, t2⟩ := abs_lt.mp hτ
  have iY2 : IntegrableOn (fun g => wt τ g * Ymean τ g ^ 2) Q := by
    have h1 : -1 < 2 * τ := by linarith
    have h2 : -1 < -(2 * τ) := by linarith
    have h0 : (-1 : ℝ) < 0 := by norm_num
    have := ((((mono_int h1 h2 h2).1.add (mono_int h2 h1 h2).1).add
      (mono_int h2 h2 h1).1).add ((((mono_int h0 h0 h2).1.add (mono_int h0 h2 h0).1).add
      (mono_int h2 h0 h0).1).const_mul 2)).const_mul (1 / 9)
    exact IntegrableOn.congr_fun this (fun g hg => (pt6 τ hg).symm) Q_meas
  have hG : 0 < Real.Gamma (1 - 2 * τ / 3) := Real.Gamma_pos_of_pos (by linarith)
  have e2 : gammaE τ (fun g => Ymean τ g ^ 2) = 1 / 3 * Aval τ + 2 / 3 * Bval τ := by
    rw [Y05.proof_6 τ hτ2]; unfold Aval Bval
    generalize Real.Gamma (1 - 2 * τ / 3) = G at *
    field_simp; ring
  rw [← e2, ← v3, gE_wt, gE_wt, gE_wt, gE_wt]
  have e : (fun g => wt τ g * (Ymean τ g - 1) ^ 3) =
      fun g => wt τ g * Ymean τ g ^ 3 - 3 * (wt τ g * Ymean τ g ^ 2) + 3 * (wt τ g * Ymean τ g)
        - wt τ g := by
    funext g; ring
  have j2 := iY2.const_mul 3
  have j1 := iY.const_mul 3
  have k1 : IntegrableOn (fun g => wt τ g * Ymean τ g ^ 3 - 3 * (wt τ g * Ymean τ g ^ 2)) Q :=
    Integrable.sub iY3 j2
  have k2 : IntegrableOn (fun g => wt τ g * Ymean τ g ^ 3 - 3 * (wt τ g * Ymean τ g ^ 2)
      + 3 * (wt τ g * Ymean τ g)) Q := Integrable.add k1 j1
  have f1 := integral_sub (μ := volume.restrict Q) (f := fun g => wt τ g * Ymean τ g ^ 3
      - 3 * (wt τ g * Ymean τ g ^ 2) + 3 * (wt τ g * Ymean τ g)) (g := fun g => wt τ g) k2 i0
  have f2 := integral_add (μ := volume.restrict Q)
    (f := fun g => wt τ g * Ymean τ g ^ 3 - 3 * (wt τ g * Ymean τ g ^ 2))
    (g := fun g => 3 * (wt τ g * Ymean τ g)) k1 j1
  have f3 := integral_sub (μ := volume.restrict Q) (f := fun g => wt τ g * Ymean τ g ^ 3)
    (g := fun g => 3 * (wt τ g * Ymean τ g ^ 2)) iY3 j2
  rw [e, f1, f2, f3, integral_const_mul, integral_const_mul, v0]
  generalize Real.Gamma (1 - 2 * τ / 3) = G at *
  have hG' : G ^ 3 ≠ 0 := by positivity
  field_simp

end Y05

namespace Y05

noncomputable def cY (k : ℕ) : ℝ := lgc (2 / 3) (-4 / 3) k - lgc (-2 / 3) (-2 / 3) k
noncomputable def cA (k : ℕ) : ℝ := lgc 2 (-2) k - fl k
noncomputable def cB (k : ℕ) : ℝ := -fl k
noncomputable def cC1 (k : ℕ) : ℝ := lgc (10 / 3) (-8 / 3) k + lgc (2 / 3) (-8 / 3) k - fl k
noncomputable def cC2 (k : ℕ) : ℝ := lgc (4 / 3) (-8 / 3) k - lgc (-2 / 3) (-2 / 3) k

open Asymptotics in
lemma LYO6 : (fun τ => Real.log (gammaE τ (Ymean τ)) - ser (cY 0) (cY 1) (cY 2) (cY 3) (cY 4) τ)
    =O[nhds 0] fun τ : ℝ => τ ^ 7 := by
  refine ((lgO6 (2 / 3) (-4 / 3)).sub (lgO6 (-2 / 3) (-2 / 3))).congr' ?_ EventuallyEq.rfl
  filter_upwards [ball_ev (1 / 8) (by norm_num)] with τ hτ
  obtain ⟨t1, t2⟩ := abs_lt.mp hτ
  show _ = Real.log (gammaE τ (Ymean τ)) - _
  rw [Y05.proof_5 τ (by linarith)]
  rw [show 1 + 2 / 3 * τ = 1 + 2 * τ / 3 by ring, show 1 + -4 / 3 * τ = 1 - 4 * τ / 3 by ring,
    show 1 + (2 / 3 * τ + -4 / 3 * τ) = 1 - 2 * τ / 3 by ring,
    show 1 + -2 / 3 * τ = 1 - 2 * τ / 3 by ring,
    show 1 + (-2 / 3 * τ + -2 / 3 * τ) = 1 - 4 * τ / 3 by ring]
  have ga : 0 < Real.Gamma (1 - 2 * τ / 3) := Real.Gamma_pos_of_pos (by linarith)
  have gb : 0 < Real.Gamma (1 - 4 * τ / 3) := Real.Gamma_pos_of_pos (by linarith)
  have gc : 0 < Real.Gamma (1 + 2 * τ / 3) := Real.Gamma_pos_of_pos (by linarith)
  have hlog : Real.log (Real.Gamma (1 + 2 * τ / 3) * Real.Gamma (1 - 4 * τ / 3) ^ 2 /
      Real.Gamma (1 - 2 * τ / 3) ^ 3) =
      Real.log (Real.Gamma (1 + 2 * τ / 3) * Real.Gamma (1 - 4 * τ / 3) /
        Real.Gamma (1 - 2 * τ / 3)) -
      Real.log (Real.Gamma (1 - 2 * τ / 3) * Real.Gamma (1 - 2 * τ / 3) /
        Real.Gamma (1 - 4 * τ / 3)) := by
    rw [← Real.log_div (by positivity) (by positivity)]
    congr 1
    generalize Real.Gamma (1 - 2 * τ / 3) = a at *
    generalize Real.Gamma (1 - 4 * τ / 3) = b at *
    generalize Real.Gamma (1 + 2 * τ / 3) = c at *
    field_simp
  rw [hlog]
  unfold ser cY
  ring

open Asymptotics in
lemma LAO6 : (fun τ => Real.log (Aval τ) - ser (cA 0) (cA 1) (cA 2) (cA 3) (cA 4) τ)
    =O[nhds 0] fun τ : ℝ => τ ^ 7 := by
  refine ((lgO6 2 (-2)).sub flat6).congr' ?_ EventuallyEq.rfl
  filter_upwards [ball_ev (1 / 8) (by norm_num)] with τ hτ
  obtain ⟨t1, t2⟩ := abs_lt.mp hτ
  show _ = Real.log (Aval τ) - _
  rw [show 1 + -2 * τ = 1 - 2 * τ by ring, show 1 + (2 * τ + -2 * τ) = 1 by ring,
    Real.Gamma_one, div_one]
  have ga : 0 < Real.Gamma (1 - 2 * τ / 3) := Real.Gamma_pos_of_pos (by linarith)
  have gb : 0 < Real.Gamma (1 - 2 * τ) := Real.Gamma_pos_of_pos (by linarith)
  have gc : 0 < Real.Gamma (1 + 2 * τ) := Real.Gamma_pos_of_pos (by linarith)
  have hlog : Real.log (Aval τ) = Real.log (Real.Gamma (1 + 2 * τ) * Real.Gamma (1 - 2 * τ)) -
      Real.log (Real.Gamma (1 - 2 * τ / 3) ^ 3 / Real.Gamma (1 - 2 * τ)) := by
    rw [← Real.log_div (by positivity) (by positivity)]
    unfold Aval
    congr 1
    generalize Real.Gamma (1 - 2 * τ / 3) = a at *
    generalize Real.Gamma (1 - 2 * τ) = b at *
    field_simp
  rw [hlog]
  unfold ser cA
  ring

open Asymptotics in
lemma LBO6 : (fun τ => Real.log (Bval τ) - ser (cB 0) (cB 1) (cB 2) (cB 3) (cB 4) τ)
    =O[nhds 0] fun τ : ℝ => τ ^ 7 := by
  refine flat6.neg_left.congr_left (fun τ => ?_)
  unfold Bval ser cB
  rw [show Real.Gamma (1 - 2 * τ) / Real.Gamma (1 - 2 * τ / 3) ^ 3 =
    (Real.Gamma (1 - 2 * τ / 3) ^ 3 / Real.Gamma (1 - 2 * τ))⁻¹ by rw [inv_div], Real.log_inv]
  ring

open Asymptotics in
lemma LC1O6 : (fun τ => Real.log (C1val τ) - ser (cC1 0) (cC1 1) (cC1 2) (cC1 3) (cC1 4) τ)
    =O[nhds 0] fun τ : ℝ => τ ^ 7 := by
  refine (((lgO6 (10 / 3) (-8 / 3)).add (lgO6 (2 / 3) (-8 / 3))).sub flat6).congr' ?_
    EventuallyEq.rfl
  filter_upwards [ball_ev (1 / 8) (by norm_num)] with τ hτ
  obtain ⟨t1, t2⟩ := abs_lt.mp hτ
  show _ = Real.log (C1val τ) - _
  rw [show 1 + 10 / 3 * τ = 1 + 10 * τ / 3 by ring, show 1 + -8 / 3 * τ = 1 - 8 * τ / 3 by ring,
    show 1 + (10 / 3 * τ + -8 / 3 * τ) = 1 + 2 * τ / 3 by ring,
    show 1 + 2 / 3 * τ = 1 + 2 * τ / 3 by ring,
    show 1 + (2 / 3 * τ + -8 / 3 * τ) = 1 - 2 * τ by ring]
  have ga : 0 < Real.Gamma (1 - 2 * τ / 3) := Real.Gamma_pos_of_pos (by linarith)
  have gb : 0 < Real.Gamma (1 - 8 * τ / 3) := Real.Gamma_pos_of_pos (by linarith)
  have gc : 0 < Real.Gamma (1 + 2 * τ / 3) := Real.Gamma_pos_of_pos (by linarith)
  have gd : 0 < Real.Gamma (1 + 10 * τ / 3) := Real.Gamma_pos_of_pos (by linarith)
  have ge : 0 < Real.Gamma (1 - 2 * τ) := Real.Gamma_pos_of_pos (by linarith)
  have hlog : Real.log (C1val τ) =
      Real.log (Real.Gamma (1 + 10 * τ / 3) * Real.Gamma (1 - 8 * τ / 3) /
        Real.Gamma (1 + 2 * τ / 3)) +
      Real.log (Real.Gamma (1 + 2 * τ / 3) * Real.Gamma (1 - 8 * τ / 3) /
        Real.Gamma (1 - 2 * τ)) -
      Real.log (Real.Gamma (1 - 2 * τ / 3) ^ 3 / Real.Gamma (1 - 2 * τ)) := by
    rw [← Real.log_mul (by positivity) (by positivity), ← Real.log_div (by positivity)
      (by positivity)]
    unfold C1val
    congr 1
    generalize Real.Gamma (1 - 2 * τ / 3) = a at *
    generalize Real.Gamma (1 - 8 * τ / 3) = b at *
    generalize Real.Gamma (1 + 2 * τ / 3) = c at *
    generalize Real.Gamma (1 + 10 * τ / 3) = d at *
    generalize Real.Gamma (1 - 2 * τ) = e at *
    field_simp
  rw [hlog]
  unfold ser cC1
  ring

open Asymptotics in
lemma LC2O6 : (fun τ => Real.log (C2val τ) - ser (cC2 0) (cC2 1) (cC2 2) (cC2 3) (cC2 4) τ)
    =O[nhds 0] fun τ : ℝ => τ ^ 7 := by
  refine ((lgO6 (4 / 3) (-8 / 3)).sub (lgO6 (-2 / 3) (-2 / 3))).congr' ?_ EventuallyEq.rfl
  filter_upwards [ball_ev (1 / 8) (by norm_num)] with τ hτ
  obtain ⟨t1, t2⟩ := abs_lt.mp hτ
  show _ = Real.log (C2val τ) - _
  rw [show 1 + 4 / 3 * τ = 1 + 4 * τ / 3 by ring, show 1 + -8 / 3 * τ = 1 - 8 * τ / 3 by ring,
    show 1 + (4 / 3 * τ + -8 / 3 * τ) = 1 - 4 * τ / 3 by ring,
    show 1 + -2 / 3 * τ = 1 - 2 * τ / 3 by ring,
    show 1 + (-2 / 3 * τ + -2 / 3 * τ) = 1 - 4 * τ / 3 by ring]
  have ga : 0 < Real.Gamma (1 - 2 * τ / 3) := Real.Gamma_pos_of_pos (by linarith)
  have gb : 0 < Real.Gamma (1 - 8 * τ / 3) := Real.Gamma_pos_of_pos (by linarith)
  have gc : 0 < Real.Gamma (1 + 4 * τ / 3) := Real.Gamma_pos_of_pos (by linarith)
  have gd : 0 < Real.Gamma (1 - 4 * τ / 3) := Real.Gamma_pos_of_pos (by linarith)
  have hlog : Real.log (C2val τ) =
      Real.log (Real.Gamma (1 + 4 * τ / 3) * Real.Gamma (1 - 8 * τ / 3) /
        Real.Gamma (1 - 4 * τ / 3)) -
      Real.log (Real.Gamma (1 - 2 * τ / 3) * Real.Gamma (1 - 2 * τ / 3) /
        Real.Gamma (1 - 4 * τ / 3)) := by
    rw [← Real.log_div (by positivity) (by positivity)]
    unfold C2val
    congr 1
    generalize Real.Gamma (1 - 2 * τ / 3) = a at *
    generalize Real.Gamma (1 - 8 * τ / 3) = b at *
    generalize Real.Gamma (1 + 4 * τ / 3) = c at *
    generalize Real.Gamma (1 - 4 * τ / 3) = d at *
    field_simp
  rw [hlog]
  unfold ser cC2
  ring

open Asymptotics in
lemma expS {V : ℝ → ℝ} (c2 c3 c4 c5 c6 : ℝ)
    (hL : (fun τ => Real.log (V τ) - ser c2 c3 c4 c5 c6 τ) =O[nhds 0] fun τ : ℝ => τ ^ 7)
    (hV : ∀ᶠ τ in nhds (0 : ℝ), 0 < V τ) :
    (fun τ => V τ - (1 + ser c2 c3 (c4 + c2 ^ 2 / 2) (c5 + c2 * c3)
      (c6 + c3 ^ 2 / 2 + c2 * c4 + c2 ^ 3 / 6) τ)) =O[nhds 0] fun τ : ℝ => τ ^ 7 := by
  have h1 := expO3 (ser_O c2 c3 c4 c5 c6) hL
  refine (h1.add (expT_O c2 c3 c4 c5 c6)).congr' ?_ EventuallyEq.rfl
  filter_upwards [hV] with τ hτ
  simp only [Real.exp_log hτ]
  ring

lemma ev_pos_EY : ∀ᶠ τ in nhds (0 : ℝ), 0 < gammaE τ (Ymean τ) := by
  filter_upwards [ball_ev (1 / 8) (by norm_num)] with τ hτ
  obtain ⟨t1, t2⟩ := abs_lt.mp hτ
  rw [Y05.proof_5 τ (by linarith)]
  have ga : 0 < Real.Gamma (1 - 2 * τ / 3) := Real.Gamma_pos_of_pos (by linarith)
  have gb : 0 < Real.Gamma (1 - 4 * τ / 3) := Real.Gamma_pos_of_pos (by linarith)
  have gc : 0 < Real.Gamma (1 + 2 * τ / 3) := Real.Gamma_pos_of_pos (by linarith)
  positivity

lemma ev_pos_A : ∀ᶠ τ in nhds (0 : ℝ), 0 < Aval τ := by
  filter_upwards [ball_ev (1 / 8) (by norm_num)] with τ hτ
  obtain ⟨t1, t2⟩ := abs_lt.mp hτ
  unfold Aval
  have ga : 0 < Real.Gamma (1 - 2 * τ / 3) := Real.Gamma_pos_of_pos (by linarith)
  have gb : 0 < Real.Gamma (1 - 2 * τ) := Real.Gamma_pos_of_pos (by linarith)
  have gc : 0 < Real.Gamma (1 + 2 * τ) := Real.Gamma_pos_of_pos (by linarith)
  positivity

lemma ev_pos_B : ∀ᶠ τ in nhds (0 : ℝ), 0 < Bval τ := by
  filter_upwards [ball_ev (1 / 8) (by norm_num)] with τ hτ
  obtain ⟨t1, t2⟩ := abs_lt.mp hτ
  unfold Bval
  have ga : 0 < Real.Gamma (1 - 2 * τ / 3) := Real.Gamma_pos_of_pos (by linarith)
  have gb : 0 < Real.Gamma (1 - 2 * τ) := Real.Gamma_pos_of_pos (by linarith)
  positivity

lemma ev_pos_C1 : ∀ᶠ τ in nhds (0 : ℝ), 0 < C1val τ := by
  filter_upwards [ball_ev (1 / 8) (by norm_num)] with τ hτ
  obtain ⟨t1, t2⟩ := abs_lt.mp hτ
  unfold C1val
  have ga : 0 < Real.Gamma (1 - 2 * τ / 3) := Real.Gamma_pos_of_pos (by linarith)
  have gb : 0 < Real.Gamma (1 - 8 * τ / 3) := Real.Gamma_pos_of_pos (by linarith)
  have gd : 0 < Real.Gamma (1 + 10 * τ / 3) := Real.Gamma_pos_of_pos (by linarith)
  positivity

lemma ev_pos_C2 : ∀ᶠ τ in nhds (0 : ℝ), 0 < C2val τ := by
  filter_upwards [ball_ev (1 / 8) (by norm_num)] with τ hτ
  obtain ⟨t1, t2⟩ := abs_lt.mp hτ
  unfold C2val
  have ga : 0 < Real.Gamma (1 - 2 * τ / 3) := Real.Gamma_pos_of_pos (by linarith)
  have gb : 0 < Real.Gamma (1 - 8 * τ / 3) := Real.Gamma_pos_of_pos (by linarith)
  have gc : 0 < Real.Gamma (1 + 4 * τ / 3) := Real.Gamma_pos_of_pos (by linarith)
  positivity

end Y05

namespace Y05

lemma sqrt_tay3 {x : ℝ} (hx : 0 ≤ x) :
    |Real.sqrt (1 + x) - (1 + x / 2 - x ^ 2 / 8 + x ^ 3 / 16)| ≤ x ^ 4 := by
  have hs2 : Real.sqrt (1 + x) ^ 2 = 1 + x := Real.sq_sqrt (by linarith)
  have hs1 : 1 ≤ Real.sqrt (1 + x) := by
    rw [Real.le_sqrt (by norm_num) (by linarith)]; linarith
  generalize Real.sqrt (1 + x) = s at *
  have hx' : x = s ^ 2 - 1 := by linarith
  rw [hx']
  have h1 : 0 ≤ s - 1 := by linarith
  have e : s - (1 + (s ^ 2 - 1) / 2 - (s ^ 2 - 1) ^ 2 / 8 + (s ^ 2 - 1) ^ 3 / 16) =
      -((s - 1) ^ 4 * ((s ^ 2 + 4 * s + 5) / 16)) := by ring
  rw [e, abs_neg, abs_of_nonneg (mul_nonneg (pow_nonneg h1 4) (by positivity)),
    show (s ^ 2 - 1) ^ 4 = (s - 1) ^ 4 * (s + 1) ^ 4 by ring]
  have h3 : (s ^ 2 + 4 * s + 5) / 16 ≤ (s + 1) ^ 4 := by nlinarith [sq_nonneg (s + 1)]
  exact mul_le_mul_of_nonneg_left h3 (pow_nonneg h1 4)

lemma M8 {M : ℝ} (hM : 0 ≤ M) : M ^ 8 ≤ 16 ^ 8 * Real.exp (M / 2) := by
  have := Real.add_one_le_exp (M / 16)
  have h : M / 16 ≤ Real.exp (M / 16) := by linarith
  have h6 := pow_le_pow_left₀ (by positivity) h 8
  rw [← Real.exp_nat_mul] at h6
  calc M ^ 8 = 16 ^ 8 * (M / 16) ^ 8 := by ring
    _ ≤ 16 ^ 8 * Real.exp (((8 : ℕ) : ℝ) * (M / 16)) := by gcongr
    _ = 16 ^ 8 * Real.exp (M / 2) := by congr 2; push_cast; ring

lemma pt10 {τ : ℝ} (hτ : |τ| ≤ 1 / 90) {g : ℝ × ℝ × ℝ} (hg : g ∈ Q) :
    |wt τ g * (Real.sqrt (Ymean τ g) - (1 + (Ymean τ g - 1) / 2 - (Ymean τ g - 1) ^ 2 / 8
      + (Ymean τ g - 1) ^ 3 / 16))| ≤ 1296 * 16 ^ 8 * τ ^ 8 * (hh g.1 * (hh g.2.1 * hh g.2.2)) := by
  unfold wt
  obtain ⟨ha, hb, hc⟩ := hg
  simp only [mem_Ioi] at ha hb hc
  set c := 2 * τ / 3 with hcdef
  have hY : Ymean τ g = ((Real.exp (c * Real.log g.1)) ^ 3 + (Real.exp (c * Real.log g.2.1)) ^ 3
      + (Real.exp (c * Real.log g.2.2)) ^ 3) / (3 * (Real.exp (c * Real.log g.1) *
        Real.exp (c * Real.log g.2.1) * Real.exp (c * Real.log g.2.2))) := by
    unfold Ymean
    rw [Real.mul_rpow (mul_pos ha hb).le hc.le, Real.mul_rpow ha.le hb.le, ← hcdef,
      rp_q ha, rp_q hb, rp_q hc, rp_q ha, rp_q hb, rp_q hc]
    have k : ∀ x : ℝ, Real.exp (2 * τ * x) = Real.exp (c * x) ^ 3 := fun x => by
      rw [← Real.exp_nat_mul]; congr 1; push_cast; rw [hcdef]; ring
    rw [k, k, k]
  have hw : (g.1 * g.2.1 * g.2.2) ^ (-2 * τ / 3) = (Real.exp (c * Real.log g.1) *
        Real.exp (c * Real.log g.2.1) * Real.exp (c * Real.log g.2.2))⁻¹ := by
    rw [Real.mul_rpow (mul_pos ha hb).le hc.le, Real.mul_rpow ha.le hb.le,
      rp_q ha, rp_q hb, rp_q hc, mul_inv, mul_inv, ← Real.exp_neg, ← Real.exp_neg,
      ← Real.exp_neg, hcdef]
    congr 3 <;> ring
  rw [hY, hw]
  set l1 := Real.log g.1
  set l2 := Real.log g.2.1
  set l3 := Real.log g.2.2
  set M := |l1| + |l2| + |l3| with hM
  have hM0 : 0 ≤ M := by positivity
  set E := Real.exp (|c| * M) with hE
  have hc0 : 0 ≤ |c| := abs_nonneg c
  have up : ∀ l : ℝ, |l| ≤ M → Real.exp (c * l) ≤ E := fun l hl => by
    apply Real.exp_le_exp.mpr
    calc c * l ≤ |c * l| := le_abs_self _
      _ = |c| * |l| := abs_mul c l
      _ ≤ |c| * M := mul_le_mul_of_nonneg_left hl hc0
  have upi : ∀ l : ℝ, |l| ≤ M → (Real.exp (c * l))⁻¹ ≤ E := fun l hl => by
    rw [← Real.exp_neg, ← mul_neg]; apply up; rwa [abs_neg]
  have b1 : |l1| ≤ M := by have := abs_nonneg l2; have := abs_nonneg l3; linarith
  have b2 : |l2| ≤ M := by have := abs_nonneg l1; have := abs_nonneg l3; linarith
  have b3 : |l3| ≤ M := by have := abs_nonneg l2; have := abs_nonneg l1; linarith
  have dd : ∀ l l' : ℝ, |l| ≤ M → |l'| ≤ M → |l| + |l'| ≤ M →
      |Real.exp (c * l) - Real.exp (c * l')| ≤ 2 * |c| * M * E := fun l l' hl hl' hll => by
    calc |Real.exp (c * l) - Real.exp (c * l')|
        ≤ |c * l - c * l'| * (Real.exp (c * l) + Real.exp (c * l')) := exp_sub_le _ _
      _ ≤ (|c| * M) * (E + E) := by
          apply mul_le_mul _ (add_le_add (up l hl) (up l' hl')) (by positivity) (by positivity)
          rw [← mul_sub, abs_mul]
          exact mul_le_mul_of_nonneg_left ((abs_sub _ _).trans hll) hc0
      _ = 2 * |c| * M * E := by ring
  obtain ⟨X0, Xb, hP⟩ := algX (Real.exp_pos (c * l1)) (Real.exp_pos (c * l2))
    (Real.exp_pos (c * l3))
    (up l1 b1) (up l2 b2) (up l3 b3) (upi l1 b1) (upi l2 b2) (upi l3 b3)
    (dd l1 l2 b1 b2 (by have := abs_nonneg l3; linarith))
    (dd l2 l3 b2 b3 (by have := abs_nonneg l1; linarith))
    (dd l1 l3 b1 b3 (by have := abs_nonneg l2; linarith))
  have hW0 : 0 ≤ (Real.exp (c * l1) * Real.exp (c * l2) * Real.exp (c * l3))⁻¹ := by positivity
  have hprod : Real.exp (-(g.1 + g.2.1 + g.2.2)) * Real.exp (7 / 10 * M) ≤
      hh g.1 * (hh g.2.1 * hh g.2.2) := by
    have e1 : Real.exp (-(g.1 + g.2.1 + g.2.2)) * Real.exp (7 / 10 * M) =
        (Real.exp (-g.1) * Real.exp (7 / 10 * |l1|)) * ((Real.exp (-g.2.1) *
          Real.exp (7 / 10 * |l2|)) * (Real.exp (-g.2.2) * Real.exp (7 / 10 * |l3|))) := by
      simp only [← Real.exp_add]; congr 1; rw [hM]; ring
    rw [e1]
    unfold hh
    simp only [← mul_add]
    have k1 := ekap (κ := 7 / 10) ha
    have k2 := ekap (κ := 7 / 10) hb
    have k3 := ekap (κ := 7 / 10) hc
    gcongr
  have hc8 : |c| ^ 8 ≤ τ ^ 8 := by
    rw [hcdef, abs_div, abs_mul, show |(2 : ℝ)| = 2 by norm_num, show |(3 : ℝ)| = 3 by norm_num,
      div_pow, mul_pow, show |τ| ^ 8 = τ ^ 8 from Even.pow_abs ⟨4, rfl⟩ τ]
    have := pow_nonneg (abs_nonneg τ) 8
    rw [show |τ| ^ 8 = τ ^ 8 from Even.pow_abs ⟨4, rfl⟩ τ] at this
    linarith
  have h27 : 27 * |c| * M ≤ M / 5 := by
    have : |c| = 2 * |τ| / 3 := by
      rw [hcdef, abs_div, abs_mul]; norm_num
    rw [this]; nlinarith [abs_nonneg τ]
  have hEpow : E ^ 3 * ((3 / 2) * (2 * |c| * M * E) ^ 2 * E ^ 4) ^ 4 =
      1296 * |c| ^ 8 * M ^ 8 * Real.exp (27 * |c| * M) := by
    rw [hE, show 27 * |c| * M = ((27 : ℕ) : ℝ) * (|c| * M) by push_cast; ring, Real.exp_nat_mul]
    ring
  have hm8 := M8 hM0
  have hτ8 : 0 ≤ τ ^ 8 := by positivity
  set e := Real.exp (-(g.1 + g.2.1 + g.2.2))
  have he : 0 < e := Real.exp_pos _
  have key := sqrt_tay3 X0
  rw [show 1 + ((Real.exp (c * l1) ^ 3 + Real.exp (c * l2) ^ 3 + Real.exp (c * l3) ^ 3) /
    (3 * (Real.exp (c * l1) * Real.exp (c * l2) * Real.exp (c * l3))) - 1) =
    (Real.exp (c * l1) ^ 3 + Real.exp (c * l2) ^ 3 + Real.exp (c * l3) ^ 3) /
    (3 * (Real.exp (c * l1) * Real.exp (c * l2) * Real.exp (c * l3))) by ring] at key
  generalize (Real.exp (c * l1) ^ 3 + Real.exp (c * l2) ^ 3 + Real.exp (c * l3) ^ 3) /
    (3 * (Real.exp (c * l1) * Real.exp (c * l2) * Real.exp (c * l3))) = Y at *
  generalize (Real.exp (c * l1) * Real.exp (c * l2) * Real.exp (c * l3))⁻¹ = W at *
  rw [abs_mul, abs_of_nonneg (mul_nonneg hW0 he.le)]
  calc W * e * |Real.sqrt Y - (1 + (Y - 1) / 2 - (Y - 1) ^ 2 / 8 + (Y - 1) ^ 3 / 16)|
      ≤ W * e * (Y - 1) ^ 4 := mul_le_mul_of_nonneg_left key (mul_nonneg hW0 he.le)
    _ ≤ E ^ 3 * e * ((3 / 2) * (2 * |c| * M * E) ^ 2 * E ^ 4) ^ 4 :=
        mul_le_mul (mul_le_mul_of_nonneg_right hP he.le) (pow_le_pow_left₀ X0 Xb 4)
          (by positivity) (by positivity)
    _ = e * (1296 * |c| ^ 8 * M ^ 8 * Real.exp (27 * |c| * M)) := by rw [← hEpow]; ring
    _ ≤ e * (1296 * τ ^ 8 * (16 ^ 8 * Real.exp (M / 2)) * Real.exp (M / 5)) := by
        gcongr
    _ = 1296 * 16 ^ 8 * τ ^ 8 * (e * Real.exp (7 / 10 * M)) := by
        rw [show 7 / 10 * M = M / 2 + M / 5 by ring, Real.exp_add]; ring
    _ ≤ 1296 * 16 ^ 8 * τ ^ 8 * (hh g.1 * (hh g.2.1 * hh g.2.2)) := by gcongr


lemma B3_bd : ∃ C δ : ℝ, 0 < δ ∧ ∀ τ : ℝ, |τ| ≤ δ →
    |B τ - (1 + gammaE τ (fun g => Ymean τ g - 1) / 2
      - gammaE τ (fun g => (Ymean τ g - 1) ^ 2) / 8
      + gammaE τ (fun g => (Ymean τ g - 1) ^ 3) / 16)| ≤ C * τ ^ 8 := by
  obtain ⟨ε, hε, hGl⟩ := gam_lb
  set K := ∫ g in Q, hh g.1 * (hh g.2.1 * hh g.2.2)
  refine ⟨8 * (1296 * 16 ^ 8) * K, min (1 / 90) ε, lt_min (by norm_num) hε, fun τ hτ => ?_⟩
  have hτ1 : |τ| ≤ 1 / 90 := hτ.trans (min_le_left _ _)
  have hG := hGl τ (hτ.trans (min_le_right _ _))
  have hτ4 : |τ| < 1 / 4 := by linarith
  have hτ2 : |τ| < 1 / 2 := by linarith
  obtain ⟨i0, v0, iY⟩ := wt_pack hτ2
  obtain ⟨iY3, -⟩ := Y3_pack hτ4
  obtain ⟨t1, t2⟩ := abs_lt.mp hτ4
  have iY2 : IntegrableOn (fun g => wt τ g * Ymean τ g ^ 2) Q := by
    have h1 : -1 < 2 * τ := by linarith
    have h2 : -1 < -(2 * τ) := by linarith
    have h0 : (-1 : ℝ) < 0 := by norm_num
    have := ((((mono_int h1 h2 h2).1.add (mono_int h2 h1 h2).1).add
      (mono_int h2 h2 h1).1).add ((((mono_int h0 h0 h2).1.add (mono_int h0 h2 h0).1).add
      (mono_int h2 h0 h0).1).const_mul 2)).const_mul (1 / 9)
    exact IntegrableOn.congr_fun this (fun g hg => (pt6 τ hg).symm) Q_meas
  have i1 : IntegrableOn (fun g => wt τ g * Real.sqrt (Ymean τ g)) Q := Fb_int (by linarith)
  have i2 : IntegrableOn (fun g => wt τ g * (Ymean τ g - 1)) Q :=
    (iY.sub i0).congr_fun (fun g _ => by simp only [Pi.sub_apply]; ring) Q_meas
  have i3 : IntegrableOn (fun g => wt τ g * (Ymean τ g - 1) ^ 2) Q :=
    ((iY2.sub (iY.const_mul 2)).add i0).congr_fun
      (fun g _ => by simp only [Pi.sub_apply, Pi.add_apply]; ring) Q_meas
  have i4 : IntegrableOn (fun g => wt τ g * (Ymean τ g - 1) ^ 3) Q :=
    (((iY3.sub (iY2.const_mul 3)).add (iY.const_mul 3)).sub i0).congr_fun
      (fun g _ => by simp only [Pi.sub_apply, Pi.add_apply]; ring) Q_meas
  set f : ℝ × ℝ × ℝ → ℝ := fun g => wt τ g * (Real.sqrt (Ymean τ g) -
    (1 + (Ymean τ g - 1) / 2 - (Ymean τ g - 1) ^ 2 / 8 + (Ymean τ g - 1) ^ 3 / 16)) with hfdef
  have hf : f = fun g => wt τ g * Real.sqrt (Ymean τ g) - wt τ g
      - (wt τ g * (Ymean τ g - 1)) / 2 + (wt τ g * (Ymean τ g - 1) ^ 2) / 8
      - (wt τ g * (Ymean τ g - 1) ^ 3) / 16 := by
    funext g; rw [hfdef]; ring
  have hint : ∫ g in Q, f g = (∫ g in Q, wt τ g * Real.sqrt (Ymean τ g)) - (∫ g in Q, wt τ g)
      - (∫ g in Q, wt τ g * (Ymean τ g - 1)) / 2 + (∫ g in Q, wt τ g * (Ymean τ g - 1) ^ 2) / 8
      - (∫ g in Q, wt τ g * (Ymean τ g - 1) ^ 3) / 16 := by
    have k2 : IntegrableOn (fun g => wt τ g * Real.sqrt (Ymean τ g) - wt τ g) Q :=
      Integrable.sub i1 i0
    have k3 : IntegrableOn (fun g => wt τ g * Real.sqrt (Ymean τ g) - wt τ g
        - wt τ g * (Ymean τ g - 1) / 2) Q := Integrable.sub k2 (i2.div_const 2)
    have k4 : IntegrableOn (fun g => wt τ g * Real.sqrt (Ymean τ g) - wt τ g
        - wt τ g * (Ymean τ g - 1) / 2 + wt τ g * (Ymean τ g - 1) ^ 2 / 8) Q :=
      Integrable.add k3 (i3.div_const 8)
    have e0 := integral_sub (μ := volume.restrict Q)
      (f := fun g => wt τ g * Real.sqrt (Ymean τ g) - wt τ g - wt τ g * (Ymean τ g - 1) / 2
        + wt τ g * (Ymean τ g - 1) ^ 2 / 8)
      (g := fun g => wt τ g * (Ymean τ g - 1) ^ 3 / 16) k4 (i4.div_const 16)
    have e1 := integral_add (μ := volume.restrict Q)
      (f := fun g => wt τ g * Real.sqrt (Ymean τ g) - wt τ g - wt τ g * (Ymean τ g - 1) / 2)
      (g := fun g => wt τ g * (Ymean τ g - 1) ^ 2 / 8) k3 (i3.div_const 8)
    have e2 := integral_sub (μ := volume.restrict Q)
      (f := fun g => wt τ g * Real.sqrt (Ymean τ g) - wt τ g)
      (g := fun g => wt τ g * (Ymean τ g - 1) / 2) k2 (i2.div_const 2)
    have e3 := integral_sub (μ := volume.restrict Q)
      (f := fun g => wt τ g * Real.sqrt (Ymean τ g)) (g := fun g => wt τ g) i1 i0
    rw [hf]
    simp only at e0 e1 e2 e3 ⊢
    rw [e0, e1, e2, e3, integral_div, integral_div, integral_div]
  have hGp : 0 < Real.Gamma (1 - 2 * τ / 3) := by linarith
  have hrep : B τ - (1 + gammaE τ (fun g => Ymean τ g - 1) / 2
      - gammaE τ (fun g => (Ymean τ g - 1) ^ 2) / 8
      + gammaE τ (fun g => (Ymean τ g - 1) ^ 3) / 16) =
      (Real.Gamma (1 - 2 * τ / 3) ^ 3)⁻¹ * ∫ g in Q, f g := by
    rw [hint, B, gE_wt, gE_wt, gE_wt, gE_wt, v0]
    have hGn : Real.Gamma (1 - 2 * τ / 3) ≠ 0 := hGp.ne'
    generalize Real.Gamma (1 - 2 * τ / 3) = G at *
    field_simp
    ring
  rw [hrep, abs_mul, abs_of_pos (by positivity)]
  have hbd : |∫ g in Q, f g| ≤
      ∫ g in Q, 1296 * 16 ^ 8 * τ ^ 8 * (hh g.1 * (hh g.2.1 * hh g.2.2)) := by
    rw [← Real.norm_eq_abs]
    refine norm_integral_le_of_norm_le (hh3_int.const_mul _) ?_
    refine (ae_restrict_iff' Q_meas).mpr (Eventually.of_forall fun g hg => ?_)
    rw [Real.norm_eq_abs]
    exact pt10 hτ1 hg
  rw [integral_const_mul] at hbd
  have hinv : (Real.Gamma (1 - 2 * τ / 3) ^ 3)⁻¹ ≤ 8 := by
    rw [inv_le_comm₀ (by positivity) (by norm_num)]
    calc 8⁻¹ = (1 / 2 : ℝ) ^ 3 := by norm_num
      _ ≤ _ := pow_le_pow_left₀ (by norm_num) hG 3
  calc (Real.Gamma (1 - 2 * τ / 3) ^ 3)⁻¹ * |∫ g in Q, f g|
      ≤ 8 * (1296 * 16 ^ 8 * τ ^ 8 * K) := mul_le_mul hinv hbd (abs_nonneg _) (by norm_num)
    _ = 8 * (1296 * 16 ^ 8) * K * τ ^ 8 := by ring

end Y05

namespace Y05

noncomputable def ex (c : ℕ → ℝ) : ℕ → ℝ
  | 0 => c 0
  | 1 => c 1
  | 2 => c 2 + c 0 ^ 2 / 2
  | 3 => c 3 + c 0 * c 1
  | _ => c 4 + c 1 ^ 2 / 2 + c 0 * c 2 + c 0 ^ 3 / 6

noncomputable def lx (c : ℕ → ℝ) : ℕ → ℝ
  | 0 => c 0
  | 1 => c 1
  | 2 => c 2 - c 0 ^ 2 / 2
  | 3 => c 3 - c 0 * c 1
  | _ => c 4 - c 1 ^ 2 / 2 - c 0 * c 2 + c 0 ^ 3 / 3

noncomputable def cy (k : ℕ) : ℝ :=
  ex cY k / 2 - (ex cA k / 3 + 2 * ex cB k / 3 - 2 * ex cY k) / 8
    + ((3 * ex cC1 k + 18 * ex cC2 k) / 27 - (ex cA k + 2 * ex cB k) + 3 * ex cY k) / 16

noncomputable def tot (k : ℕ) : ℝ := fl k + lx cy k

set_option maxHeartbeats 4000000 in
lemma tot_eq (τ : ℝ) : ser (tot 0) (tot 1) (tot 2) (tot 3) (tot 4) τ =
    T09.logSeries (X05.z 3) τ - 160 * X05.z 5 / 27 * τ ^ 5
      + (16 / 81 * X05.z 3 ^ 2 - 572 / 76545 * Real.pi ^ 6) * τ ^ 6 := by
  simp only [ser, tot, lx, cy, ex, cY, cA, cB, cC1, cC2, lgc, fl, T09.logSeries]
  norm_num
  rw [z2_eq, z4_eq, z6_eq]
  ring

end Y05

namespace Y05

open Asymptotics in
lemma h8O : (fun τ => B τ - (1 + gammaE τ (fun g => Ymean τ g - 1) / 2
    - gammaE τ (fun g => (Ymean τ g - 1) ^ 2) / 8
    + gammaE τ (fun g => (Ymean τ g - 1) ^ 3) / 16)) =O[nhds 0] fun τ : ℝ => τ ^ 7 := by
  obtain ⟨C, δ, hδ, h⟩ := B3_bd
  refine IsBigO.of_bound |C| ?_
  filter_upwards [ball_ev (min δ 1) (lt_min hδ one_pos)] with τ hτ
  have h1 : |τ| ≤ δ := (le_of_lt hτ).trans (min_le_left _ _)
  have h2 : |τ| ≤ 1 := (le_of_lt hτ).trans (min_le_right _ _)
  rw [Real.norm_eq_abs, Real.norm_eq_abs, abs_pow]
  calc _ ≤ C * τ ^ 8 := h τ h1
    _ ≤ |C| * |τ| ^ 7 := by
      rw [show τ ^ 8 = |τ| ^ 8 from (Even.pow_abs ⟨4, rfl⟩ τ).symm]
      exact mul_le_mul (le_abs_self _) (pow_le_pow_of_le_one (abs_nonneg _) h2 (by norm_num))
        (by positivity) (abs_nonneg _)

open Asymptotics in
lemma yO6 : (fun τ => (B τ - 1) - ser (cy 0) (cy 1) (cy 2) (cy 3) (cy 4) τ)
    =O[nhds 0] fun τ : ℝ => τ ^ 7 := by
  have hEY : (fun τ => gammaE τ (Ymean τ) - (1 + ser (ex cY 0) (ex cY 1) (ex cY 2) (ex cY 3)
      (ex cY 4) τ)) =O[nhds 0] fun τ : ℝ => τ ^ 7 := expS _ _ _ _ _ LYO6 ev_pos_EY
  have hA : (fun τ => Aval τ - (1 + ser (ex cA 0) (ex cA 1) (ex cA 2) (ex cA 3)
      (ex cA 4) τ)) =O[nhds 0] fun τ : ℝ => τ ^ 7 := expS _ _ _ _ _ LAO6 ev_pos_A
  have hB : (fun τ => Bval τ - (1 + ser (ex cB 0) (ex cB 1) (ex cB 2) (ex cB 3)
      (ex cB 4) τ)) =O[nhds 0] fun τ : ℝ => τ ^ 7 := expS _ _ _ _ _ LBO6 ev_pos_B
  have hC1 : (fun τ => C1val τ - (1 + ser (ex cC1 0) (ex cC1 1) (ex cC1 2) (ex cC1 3)
      (ex cC1 4) τ)) =O[nhds 0] fun τ : ℝ => τ ^ 7 := expS _ _ _ _ _ LC1O6 ev_pos_C1
  have hC2 : (fun τ => C2val τ - (1 + ser (ex cC2 0) (ex cC2 1) (ex cC2 2) (ex cC2 3)
      (ex cC2 4) τ)) =O[nhds 0] fun τ : ℝ => τ ^ 7 := expS _ _ _ _ _ LC2O6 ev_pos_C2
  have := ((h8O.add (hEY.const_mul_left (1 / 2))).sub
    ((((hA.const_mul_left (1 / 3)).add (hB.const_mul_left (2 / 3))).sub
      (hEY.const_mul_left 2)).const_mul_left (1 / 8))).add
    (((((hC1.const_mul_left (3 / 27)).add (hC2.const_mul_left (18 / 27))).sub
      (hA.add (hB.const_mul_left 2))).add (hEY.const_mul_left 3)).const_mul_left (1 / 16))
  refine this.congr' ?_ EventuallyEq.rfl
  filter_upwards [ball_ev (1 / 8) (by norm_num)] with τ hτ
  have hτ' : |τ| < 1 / 2 := by linarith
  have hτ4 : |τ| < 1 / 4 := by linarith
  rw [EX_eq hτ', EX2_eq hτ', EX3_eq hτ4]
  simp only [ser, cy]
  ring

open Asymptotics in
lemma yO6_2 : (fun τ => B τ - 1) =O[nhds 0] fun τ : ℝ => τ ^ 2 :=
  ((yO6.trans (T09.hpow 2 5)).add (ser_O _ _ _ _ _)).congr_left (fun τ => by ring)

lemma log_cub {y : ℝ} (hy : |y| ≤ 1 / 2) :
    |Real.log (1 + y) - L3 y| ≤ 2 * |y| ^ 4 := by
  have := Real.abs_log_sub_add_sum_range_le (x := -y) (by rw [abs_neg]; linarith) 3
  have e : (∑ i ∈ Finset.range 3, (-y) ^ (i + 1) / ((i : ℝ) + 1)) + Real.log (1 - -y) =
      Real.log (1 + y) - L3 y := by
    simp [Finset.sum_range_succ, L3]; ring
  rw [e, abs_neg] at this
  refine this.trans ?_
  rw [show |y| ^ (3 + 1) = |y| ^ 4 from rfl, div_le_iff₀ (by linarith)]
  have := pow_nonneg (abs_nonneg y) 4
  have := mul_le_mul_of_nonneg_left (show 1 ≤ 2 * (1 - |y|) by linarith) this
  linarith

open Asymptotics in
lemma logBO6 : (fun τ => Real.log (B τ) - ser (lx cy 0) (lx cy 1) (lx cy 2) (lx cy 3) (lx cy 4) τ)
    =O[nhds 0] fun τ : ℝ => τ ^ 7 := by
  have ht : Tendsto (fun τ => B τ - 1) (nhds 0) (nhds 0) :=
    yO6_2.trans_tendsto (by simpa using ((continuous_pow 2).tendsto (0:ℝ)))
  have hev : ∀ᶠ τ in nhds (0 : ℝ), |B τ - 1| ≤ 1 / 2 := by
    filter_upwards [ht.eventually (Metric.closedBall_mem_nhds (0:ℝ) (by norm_num : (0:ℝ) < 1 / 2))]
      with τ hτ
    simpa [Real.norm_eq_abs] using hτ
  have h4 : (fun τ => Real.log (B τ) - L3 (B τ - 1)) =O[nhds 0] fun τ => (B τ - 1) ^ 4 := by
    refine IsBigO.of_bound 2 ?_
    filter_upwards [hev] with τ hτ
    have := log_cub hτ
    rw [show 1 + (B τ - 1) = B τ by ring] at this
    rw [Real.norm_eq_abs, Real.norm_eq_abs, abs_pow]
    exact this
  have h8 : (fun τ => (B τ - 1) ^ 4) =O[nhds 0] fun τ : ℝ => τ ^ 8 := by
    refine (yO6_2.pow 4).trans (isBigO_of_le _ (fun τ => ?_))
    simp only [Real.norm_eq_abs, ← pow_mul]
    simp
  set s := fun τ => ser (cy 0) (cy 1) (cy 2) (cy 3) (cy 4) τ with hs
  have hsO : s =O[nhds 0] fun τ : ℝ => τ ^ 2 := ser_O _ _ _ _ _
  have hs0 : Tendsto s (nhds 0) (nhds 0) :=
    hsO.trans_tendsto (by simpa using ((continuous_pow 2).tendsto (0:ℝ)))
  have hb : (fun τ => 1 - ((B τ - 1) + s τ) / 2 + ((B τ - 1) ^ 2 + (B τ - 1) * s τ + s τ ^ 2) / 3)
      =O[nhds 0] (1 : ℝ → ℝ) := by
    have : Tendsto (fun τ => 1 - ((B τ - 1) + s τ) / 2
        + ((B τ - 1) ^ 2 + (B τ - 1) * s τ + s τ ^ 2) / 3) (nhds 0)
        (nhds (1 - (0 + 0) / 2 + (0 ^ 2 + 0 * 0 + 0 ^ 2) / 3)) :=
      ((tendsto_const_nhds.sub ((ht.add hs0).div_const 2)).add
        ((((ht.pow 2).add (ht.mul hs0)).add (hs0.pow 2)).div_const 3))
    exact this.isBigO_one ℝ
  have hd := yO6.mul hb
  simp only [Pi.one_apply, mul_one] at hd
  have hT := logT_O (cy 0) (cy 1) (cy 2) (cy 3) (cy 4)
  refine (((h4.trans (h8.trans (T09.hpow 7 1))).add hd).add hT).congr_left (fun τ => ?_)
  simp only [hs, lx]
  unfold L3
  ring

open Asymptotics in
lemma mainτ6 : (fun τ => Real.log (Real.Gamma (1 - 2 * τ / 3) ^ 3 / Real.Gamma (1 - 2 * τ) * B τ)
    - (T09.logSeries (X05.z 3) τ - 160 * X05.z 5 / 27 * τ ^ 5
      + (16 / 81 * X05.z 3 ^ 2 - 572 / 76545 * Real.pi ^ 6) * τ ^ 6))
    =O[nhds 0] fun τ : ℝ => τ ^ 7 := by
  refine (flat6.add logBO6).congr' ?_ EventuallyEq.rfl
  filter_upwards [ball_ev (1 / 8) (by norm_num)] with τ hτ
  obtain ⟨t1, t2⟩ := abs_lt.mp hτ
  have ga : 0 < Real.Gamma (1 - 2 * τ / 3) := Real.Gamma_pos_of_pos (by linarith)
  have gb : 0 < Real.Gamma (1 - 2 * τ) := Real.Gamma_pos_of_pos (by linarith)
  have hB := B_pos (show τ < 1 by linarith)
  rw [Real.log_mul (by positivity) hB.ne', ← tot_eq τ]
  simp only [ser, tot]
  ring

end Y05

open Asymptotics in
theorem Y05.proof_11 : Y05.stmt_11 := by
  refine ⟨Y05.proof_10.1, ?_⟩
  obtain ⟨c, hc⟩ := mainτ6.bound
  obtain ⟨ε, hε, hεb⟩ := Metric.eventually_nhds_iff.mp hc
  refine ⟨c, min (ε / 2) (1 / 8), lt_min (by linarith) (by norm_num), fun p hp hτ => ?_⟩
  have h8 : |1 - 1 / p| ≤ 1 / 8 := hτ.trans (min_le_right _ _)
  have hε2 : |1 - 1 / p| ≤ ε / 2 := hτ.trans (min_le_left _ _)
  rw [factor_id hp h8]
  have := hεb (y := 1 - 1 / p) (by rw [Real.dist_eq, sub_zero]; linarith)
  simpa [Real.norm_eq_abs, abs_pow] using this
