import Mathlib
import Zeta3.Wave3.W02_Holo
import Zeta3.Wave5.Y05_Octahedron
/-!
# Z04 — the area of the superellipsoid (paper label `thm:S`)

Target id: Z04.
Paper labels: `thm:S` (paper/ieee/main.tex lines 1921–1945; setting line 1909–1912), with the
detailed proof paper/ieee/proofs/sec_g3_S.tex (Lemmas `g3:lem:area`, `g3:lem:LD`,
`g3:lem:exact`, `g3:lem:smooth`, `g3:lem:coeffs`).

Quoted verbatim.

main.tex 1909:

  The three-dimensional squircle is the superellipsoid $|x|^p+|y|^p+|z|^p=1$
  (Fig.~\ref{fig:superellipsoid}). Its surface area $S(p)$ goes from $4\sqrt3$ (the octahedron,
  $p=1$) through $4\pi$ (the sphere) to $24$ (the cube).

main.tex 1921–1945 (`thm:S`, the parts formalized here):

  Let $\varepsilon=1/p$ and $G(\varepsilon)=\Gamma(1+\varepsilon)^2/\Gamma(1+2\varepsilon)$.
  For $p\ge1$,
  \[ S(p)=6\,G(\varepsilon)\bigl(P(p)-4\bigr)-8\varepsilon^2Q(\varepsilon), \]
  where $P(p)$ is the squircle perimeter and
  \begin{gather*}
  Q(\varepsilon)=-\frac{1}{\sqrt\pi\,\Gamma(1+2\varepsilon)}\int_0^\infty h_\varepsilon(y)^3\dd y,\\
  h_\varepsilon(y)=\int_0^\infty\bigl(1-e^{-x^{2-2\varepsilon}/y^2}\bigr)x^{\varepsilon-1}e^{-x}\dd x .
  \end{gather*}
  Moreover $S(p)=\sum_{k\le K}s_kp^{-k}+O(p^{-K-1})$ for every $K$, with no $\log p$ terms,
  $s_0=24$ and
  \begin{gather*}
  s_1=48M_0,\qquad s_2=48M_1-4\pi^2-8Q_0,\\
  s_3=48M_2-8\pi^2M_0+48\zeta(3)-8Q_1,
  \end{gather*}

sec_g3_S.tex, Lemma `g3:lem:area`:

  $S(p)=8\varepsilon^2\int_{\Delta}(\sigma_1\sigma_2\sigma_3)^{\varepsilon-1}\|X\|\dd\sigma$ for
  $d=3$, and $P(p)=4\varepsilon\int_0^1(\sigma_1\sigma_2)^{\varepsilon-1}\|X\|\dd\sigma_1$ for
  $d=2$ ($\sigma_2=1-\sigma_1$).

sec_g3_S.tex, before Lemma `g3:lem:LD`:

  \mathrm h(t)=\int_0^\infty\bigl(1-e^{-tx^{2-2\varepsilon}}\bigr)x^{\varepsilon-1}e^{-x}\dd x,
  \qquad H_m(\varepsilon)=\frac1{2\sqrt\pi}\int_0^\infty t^{-3/2}\mathrm h(t)^m\dd t .
  ... so $H_m=\frac1{\sqrt\pi}\int_0^\infty h_\varepsilon(y)^m\dd y$. By Tonelli and
  $\frac1{2\sqrt\pi}\int_0^\infty(1-e^{-tr})t^{-3/2}\dd t=\sqrt r$ ... we get
  $H_1=\int_0^\infty x^{1-\varepsilon}x^{\varepsilon-1}e^{-x}\dd x=1$.

sec_g3_S.tex, Lemma `g3:lem:exact`:

  For $0<\varepsilon\le1$: $P(p)=8-4\varepsilon\frac{H_2(\varepsilon)}{\Gamma(1+\varepsilon)}$ and
  $S(p)=8\,\frac{3\Gamma(1+\varepsilon)^2-3\varepsilon\Gamma(1+\varepsilon)H_2+\varepsilon^2H_3}
  {\Gamma(1+2\varepsilon)}$. Consequently $S(p)=6G(\varepsilon)(P(p)-4)-8\varepsilon^2Q(\varepsilon)$
  with $Q=-H_3/\Gamma(1+2\varepsilon)$ as in the theorem.

sec_g3_S.tex, proof of Lemma `g3:lem:smooth` (case `j = 0`):

  |\partial_\varepsilon^j\mathrm h|\le C_j\,(1+|\log y|)^{j+1}\ (y\le1),\qquad
  |\partial_\varepsilon^j\mathrm h|\le C_jy^{-2}\ (y\ge1),
  uniformly in $\varepsilon\in[0,\frac14]$.

Formalization choices.
* `S = Y05.S` (imported). `Y05.S p` is DEFINED as `8 (1/p)² ∫_Δ (t₁t₂t₃)^{-τ}
  (t₁^{2τ}+t₂^{2τ}+t₃^{2τ})^{1/2}`, `τ = 1 - 1/p`, which is literally the simplex integral of
  Lemma `g3:lem:area` (`(σ₁σ₂σ₃)^{ε-1} = (σ₁σ₂σ₃)^{-τ}`, `‖X‖² = Σ σᵢ^{2-2ε} = Σ σᵢ^{2τ}`).
  Y05 proves (Y05.stmt_1, Y05.stmt_2) that its integrand is integrable and that `S` is `8` times
  the graph-area integral of the first-octant piece `z = (1-x^p-y^p)^{1/p}`.
* `P = W02.Pint` (imported): `P(p) = (4/p) ∫_0^1 √(t^{2/p-2} + (1-t)^{2/p-2}) dt` (eq:Pint);
  W04 proves it is the true length (`eVariationOn`) of the squircle for `p ≥ 1`.
* `h ε y = ∫_{(0,∞)} (1 - e^{-x^{2-2ε}/y²}) x^{ε-1} e^{-x} dx` (Bochner integral over `Ioi 0`,
  so `x > 0` and real `rpow` is unambiguous; only `y > 0` is ever used). Every integral whose
  value is asserted comes with an explicit `IntegrableOn` claim, so no statement holds through
  the "non-integrable ↦ 0" convention.
* `H m ε = (1/√π) ∫_{(0,∞)} h_ε(y)^m dy` (the supplement's `y`-form of `H_m`).
  `Q ε` is defined literally as in the theorem, NOT via `H`.
* `G ε = Γ(1+ε)²/Γ(1+2ε)` with `Real.Gamma`.
* `M_0` is `(W02.Mk 0).re`, the real part of `M_0 = -F(0)/2` of eq:expansion (`W02.Mk` is
  complex; `W02.stmt_4` shows `F` is real on real `ε < 1`, so `M_0` is real).
* `s_0 = 24` is stated as `S(p) → 24` (`p → ∞`); `s_1 = 48 M_0` as the case `K = 1` of the
  expansion: `|S(p) - 24 - 48 M_0/p| ≤ C/p²` for `p ≥ N`.
* Numerical check before freezing (mpmath, 12–15 digits; face integral of `Y05.S` computed by
  2D quadrature, independently of the identity): `p = 1`: both sides `6.92820323028 = 4√3`;
  `p = 1.5`: both `10.1175144429`; `p = 2`: `12.56637048` vs `12.56637044` (`4π = 12.5663706`);
  `p = 3`: `15.72798812` vs `15.72798696` (agreement limited by the 2D quadrature).

Statements.
* stmt_1: the Gaussian linearization `r = (1/(2√π)) ∫_0^∞ (1-e^{-tr²}) t^{-3/2} dt`, `r ≥ 0`,
  with integrability.
* stmt_2: Lemma `g3:lem:smooth`, case `j = 0`, for `h` itself: uniformly in `ε ∈ [0,1/4]`,
  the `x`-integrand is integrable, `|h_ε(y)| ≤ C(1+|log y|)` for `y > 0` and
  `|h_ε(y)| ≤ C y^{-2}` for `y ≥ 1`.
* stmt_3: `H_1 = 1` for `0 < ε ≤ 1` (with integrability of `h_ε` in `y`).
* stmt_4: Lemma `g3:lem:exact`, perimeter part: `P(p) = 8 - 4ε H_2/Γ(1+ε)`, `p ≥ 1`.
* stmt_5: Lemma `g3:lem:exact`, area part:
  `S(p) = 8 (3Γ(1+ε)² - 3εΓ(1+ε)H_2 + ε²H_3)/Γ(1+2ε)`, `p ≥ 1`.
* stmt_6: `thm:S`, first claim (THE TARGET): for `p ≥ 1`, all integrals involved are integrable
  (the face integrand of `S`, the perimeter integrand, the `x`-integrand of `h_ε(y)` for every
  `y > 0`, and `h_ε³` on `(0,∞)`) and `S(p) = 6G(ε)(P(p)-4) - 8ε²Q(ε)`.
* stmt_7: `thm:S`, "Moreover", `s_0 = 24`: `S(p) → 24` as `p → ∞`.
* stmt_8: `thm:S`, "Moreover", `K = 1`, `s_0 = 24`, `s_1 = 48 M_0`:
  `|S(p) - 24 - 48M_0/p| ≤ C/p²` for `p ≥ N`.

The order follows the theorem: the identity first, then the expansion of its "Moreover"
sentence; stmt_7, stmt_8 are further claims of `thm:S`, not steps towards stmt_6.

NOT COVERED:
* That `Y05.S p` is the surface area (2-dimensional Hausdorff measure) of
  `{|x|^p+|y|^p+|z|^p=1}`: the 8-fold symmetry and the identification of the graph-area integral
  with Hausdorff measure (Y05 covers only the graph-area formula of the first-octant piece).
* Lemma `g3:lem:LD` (Laplace–Dirichlet convolution) as a separate general statement; it is used
  inside stmt_4, stmt_5. Lemma `g3:lem:area` for `P` (a rewriting of eq:Pint) is not stated.
* The expansion for general `K`, "no `log p` terms", the smoothness of `Q` on `[0,1/4]`
  (derivatives `j ≥ 1`), `s_2`, `s_3`, `Q_1`, and the closed form of `Q_0`
  (Section `g3:sec:S:Q0`); the numerical values of `s_1, s_2, s_3, Q_1`; the values
  `S(1) = 4√3`, `S(2) = 4π`; the `n`-dimensional `S_n(p)` and Proposition `prop:cuberadius`.
-/
-- BEGIN STATEMENTS
namespace Z04

open MeasureTheory Real Filter Topology

/-- `G(ε) = Γ(1+ε)²/Γ(1+2ε)`. -/
noncomputable def G (ε : ℝ) : ℝ := Real.Gamma (1 + ε) ^ 2 / Real.Gamma (1 + 2 * ε)

/-- The integrand of `h_ε(y)`: `(1 - e^{-x^{2-2ε}/y²}) x^{ε-1} e^{-x}`. -/
noncomputable def hIntegrand (ε y x : ℝ) : ℝ :=
  (1 - Real.exp (-(x ^ (2 - 2 * ε) / y ^ 2))) * x ^ (ε - 1) * Real.exp (-x)

/-- `h_ε(y) = ∫_0^∞ (1 - e^{-x^{2-2ε}/y²}) x^{ε-1} e^{-x} dx`. -/
noncomputable def h (ε y : ℝ) : ℝ := ∫ x in Set.Ioi (0 : ℝ), hIntegrand ε y x

/-- `H_m(ε) = (1/√π) ∫_0^∞ h_ε(y)^m dy`. -/
noncomputable def H (m : ℕ) (ε : ℝ) : ℝ :=
  1 / Real.sqrt π * ∫ y in Set.Ioi (0 : ℝ), h ε y ^ m

/-- `Q(ε) = -(1/(√π Γ(1+2ε))) ∫_0^∞ h_ε(y)³ dy`, literally as in `thm:S`. -/
noncomputable def Q (ε : ℝ) : ℝ :=
  -(1 / (Real.sqrt π * Real.Gamma (1 + 2 * ε))) * ∫ y in Set.Ioi (0 : ℝ), h ε y ^ 3

/-- sec_g3_S.tex, before Lemma `g3:lem:LD` and in the proof of Lemma `g3:lem:exact`:
"For $r\ge0$, $r=\frac1{2\sqrt\pi}\int_0^\infty(1-e^{-tr^2})t^{-3/2}\dd t$", with integrability. -/
def stmt_1 : Prop :=
  ∀ r : ℝ, 0 ≤ r →
    IntegrableOn (fun t : ℝ => (1 - Real.exp (-(t * r ^ 2))) * t ^ (-(3 / 2 : ℝ)))
        (Set.Ioi 0) ∧
      1 / (2 * Real.sqrt π) *
          ∫ t in Set.Ioi (0 : ℝ), (1 - Real.exp (-(t * r ^ 2))) * t ^ (-(3 / 2 : ℝ)) = r

/-- sec_g3_S.tex, proof of Lemma `g3:lem:smooth`, case `j = 0`: uniformly in `ε ∈ [0,1/4]`,
the integral defining `h_ε(y)` converges, `|h_ε(y)| ≤ C (1+|log y|)` and, for `y ≥ 1`,
`|h_ε(y)| ≤ C y^{-2}`. -/
def stmt_2 : Prop :=
  ∃ C : ℝ, ∀ ε : ℝ, 0 ≤ ε → ε ≤ 1 / 4 → ∀ y : ℝ, 0 < y →
    IntegrableOn (hIntegrand ε y) (Set.Ioi 0) ∧
      |h ε y| ≤ C * (1 + |Real.log y|) ∧ (1 ≤ y → |h ε y| ≤ C / y ^ 2)

/-- sec_g3_S.tex, before Lemma `g3:lem:LD`: "we get $H_1=\int_0^\infty
x^{1-\varepsilon}x^{\varepsilon-1}e^{-x}\dd x=1$", for `0 < ε ≤ 1`, with the integrability of
`h_ε(y)` in `x` and in `y`. -/
def stmt_3 : Prop :=
  ∀ ε : ℝ, 0 < ε → ε ≤ 1 →
    (∀ y : ℝ, 0 < y → IntegrableOn (hIntegrand ε y) (Set.Ioi 0)) ∧
      IntegrableOn (h ε) (Set.Ioi 0) ∧ H 1 ε = 1

/-- sec_g3_S.tex, Lemma `g3:lem:exact`, perimeter part: for `p ≥ 1`, `ε = 1/p`,
"$P(p)=8-4\varepsilon\frac{H_2(\varepsilon)}{\Gamma(1+\varepsilon)}$", with the integrability of
`h_ε²` on `(0,∞)`. -/
def stmt_4 : Prop :=
  ∀ p : ℝ, 1 ≤ p →
    IntegrableOn (fun y => h (1 / p) y ^ 2) (Set.Ioi 0) ∧
      W02.Pint p = 8 - 4 * (1 / p) * H 2 (1 / p) / Real.Gamma (1 + 1 / p)

/-- sec_g3_S.tex, Lemma `g3:lem:exact`, area part: for `p ≥ 1`, `ε = 1/p`,
"$S(p)=8\,\frac{3\Gamma(1+\varepsilon)^2-3\varepsilon\Gamma(1+\varepsilon)H_2+\varepsilon^2H_3}
{\Gamma(1+2\varepsilon)}$", with the integrability of `h_ε²` and `h_ε³` on `(0,∞)`. -/
def stmt_5 : Prop :=
  ∀ p : ℝ, 1 ≤ p →
    IntegrableOn (fun y => h (1 / p) y ^ 2) (Set.Ioi 0) ∧
      IntegrableOn (fun y => h (1 / p) y ^ 3) (Set.Ioi 0) ∧
      Y05.S p = 8 * (3 * Real.Gamma (1 + 1 / p) ^ 2
          - 3 * (1 / p) * Real.Gamma (1 + 1 / p) * H 2 (1 / p)
          + (1 / p) ^ 2 * H 3 (1 / p)) / Real.Gamma (1 + 2 * (1 / p))

/-- `thm:S`, first claim (the target): for `p ≥ 1`, `ε = 1/p`,
"$S(p)=6\,G(\varepsilon)\bigl(P(p)-4\bigr)-8\varepsilon^2Q(\varepsilon)$", together with the
integrability of every integral involved: the face integrand of `S` on `Δ`, the perimeter
integrand on `[0,1]`, the integrand of `h_ε(y)` on `(0,∞)` for each `y > 0`, and `h_ε³` on
`(0,∞)`. -/
def stmt_6 : Prop :=
  ∀ p : ℝ, 1 ≤ p →
    IntegrableOn (Y05.faceIntegrand (1 - 1 / p)) Y05.Δ ∧
      IntervalIntegrable (W02.pIntegrand p) volume 0 1 ∧
      (∀ y : ℝ, 0 < y → IntegrableOn (hIntegrand (1 / p) y) (Set.Ioi 0)) ∧
      IntegrableOn (fun y => h (1 / p) y ^ 3) (Set.Ioi 0) ∧
      Y05.S p = 6 * G (1 / p) * (W02.Pint p - 4) - 8 * (1 / p) ^ 2 * Q (1 / p)

/-- `thm:S`, "Moreover", `s_0 = 24`: `S(p) → 24` as `p → ∞` (the cube). -/
def stmt_7 : Prop :=
  Tendsto Y05.S atTop (𝓝 24)

/-- `thm:S`, "Moreover", case `K = 1` with `s_0 = 24`, `s_1 = 48 M_0`
(`M_0 = (W02.Mk 0).re`): there are `C`, `N` with `|S(p) - 24 - 48 M_0/p| ≤ C/p²` for `p ≥ N`. -/
def stmt_8 : Prop :=
  ∃ C N : ℝ, ∀ p : ℝ, N ≤ p →
    |Y05.S p - 24 - 48 * (W02.Mk 0).re / p| ≤ C / p ^ 2

end Z04
-- END STATEMENTS

open MeasureTheory Real Set Filter Topology

namespace Z04

lemma one_sub_exp_neg_le {x : ℝ} : 1 - Real.exp (-x) ≤ x := by
  have := Real.add_one_le_exp (-x); linarith

lemma one_sub_exp_neg_nn {x : ℝ} (hx : 0 ≤ x) : 0 ≤ 1 - Real.exp (-x) := by
  have : Real.exp (-x) ≤ 1 := Real.exp_le_one_iff.mpr (by linarith)
  linarith

lemma one_sub_exp_neg_le_one (x : ℝ) : 1 - Real.exp (-x) ≤ 1 := by
  have := Real.exp_pos (-x); linarith

lemma gt_int {c : ℝ} (hc : 0 ≤ c) :
    IntegrableOn (fun t : ℝ => (1 - Real.exp (-(t * c))) * t ^ (-(3 / 2 : ℝ))) (Ioi 0) := by
  have hm : Measurable (fun t : ℝ => (1 - Real.exp (-(t * c))) * t ^ (-(3 / 2 : ℝ))) := by
    fun_prop
  have hsplit : Ioi (0 : ℝ) = Ioc 0 1 ∪ Ioi 1 := (Ioc_union_Ioi_eq_Ioi zero_le_one).symm
  rw [hsplit]
  refine IntegrableOn.union ?_ ?_
  · have hi : IntegrableOn (fun t : ℝ => c * t ^ (-(1 / 2 : ℝ))) (Ioc 0 1) := by
      have := (intervalIntegral.intervalIntegrable_rpow' (r := -(1 / 2 : ℝ)) (by norm_num)
        (a := 0) (b := 1))
      rw [intervalIntegrable_iff_integrableOn_Ioc_of_le zero_le_one] at this
      exact this.const_mul c
    refine hi.mono' hm.aestronglyMeasurable ?_
    refine (ae_restrict_iff' measurableSet_Ioc).mpr (Eventually.of_forall fun t ht => ?_)
    obtain ⟨ht0, ht1⟩ := ht
    have h1 := one_sub_exp_neg_nn (mul_nonneg ht0.le hc)
    have h2 := one_sub_exp_neg_le (x := t * c)
    have hp : 0 < t ^ (-(3 / 2 : ℝ)) := Real.rpow_pos_of_pos ht0 _
    rw [Real.norm_eq_abs, abs_of_nonneg (mul_nonneg h1 hp.le)]
    have e : t ^ (-(1 / 2 : ℝ)) = t * t ^ (-(3 / 2 : ℝ)) := by
      rw [show (-(1 / 2 : ℝ)) = 1 + (-(3 / 2)) by norm_num, Real.rpow_add ht0, Real.rpow_one]
    rw [e]
    calc (1 - Real.exp (-(t * c))) * t ^ (-(3 / 2 : ℝ)) ≤ (t * c) * t ^ (-(3 / 2 : ℝ)) :=
          mul_le_mul_of_nonneg_right h2 hp.le
      _ = _ := by ring
  · have hi : IntegrableOn (fun t : ℝ => t ^ (-(3 / 2 : ℝ))) (Ioi 1) :=
      integrableOn_Ioi_rpow_of_lt (by norm_num) zero_lt_one
    refine hi.mono' hm.aestronglyMeasurable ?_
    refine (ae_restrict_iff' measurableSet_Ioi).mpr (Eventually.of_forall fun t ht => ?_)
    have ht0 : (0 : ℝ) < t := lt_trans zero_lt_one ht
    have h1 := one_sub_exp_neg_nn (mul_nonneg ht0.le hc)
    have h2 := one_sub_exp_neg_le_one (t * c)
    have hp : 0 < t ^ (-(3 / 2 : ℝ)) := Real.rpow_pos_of_pos ht0 _
    rw [Real.norm_eq_abs, abs_of_nonneg (mul_nonneg h1 hp.le)]
    calc (1 - Real.exp (-(t * c))) * t ^ (-(3 / 2 : ℝ)) ≤ 1 * t ^ (-(3 / 2 : ℝ)) :=
          mul_le_mul_of_nonneg_right h2 hp.le
      _ = _ := by ring

lemma gt_val {c : ℝ} (hc : 0 < c) :
    ∫ t in Ioi (0 : ℝ), (1 - Real.exp (-(t * c))) * t ^ (-(3 / 2 : ℝ)) =
      2 * Real.sqrt π * Real.sqrt c := by
  set u : ℝ → ℝ := fun t => 1 - Real.exp (-(t * c))
  set u' : ℝ → ℝ := fun t => c * Real.exp (-(t * c))
  set v : ℝ → ℝ := fun t => -2 * t ^ (-(1 / 2 : ℝ))
  set v' : ℝ → ℝ := fun t => t ^ (-(3 / 2 : ℝ))
  have hu : ∀ x ∈ Ioi (0 : ℝ), HasDerivAt u (u' x) x := by
    intro x _
    have h1 : HasDerivAt (fun t : ℝ => -(t * c)) (-c) x := by
      have := ((hasDerivAt_id x).mul_const c).neg
      simp only [id, one_mul] at this
      exact this
    have := h1.exp.const_sub 1
    convert this using 1
    simp only [u']; ring
  have hv : ∀ x ∈ Ioi (0 : ℝ), HasDerivAt v (v' x) x := by
    intro x hx
    have := (Real.hasDerivAt_rpow_const (p := -(1 / 2 : ℝ)) (Or.inl (ne_of_gt hx))).const_mul (-2)
    convert this using 1
    simp only [v']
    rw [show (-(1 / 2 : ℝ)) - 1 = -(3 / 2) by norm_num]; ring
  have huv' : IntegrableOn (u * v') (Ioi 0) := gt_int hc.le
  have hg : IntegrableOn (fun x : ℝ => x ^ (-(1 / 2 : ℝ)) * Real.exp (-c * x ^ (1 : ℝ))) (Ioi 0) :=
    integrableOn_rpow_mul_exp_neg_mul_rpow (by norm_num) zero_lt_one hc
  have hu'v : IntegrableOn (u' * v) (Ioi 0) := by
    refine IntegrableOn.congr_fun (hg.const_mul (-2 * c)) (fun x _ => ?_) measurableSet_Ioi
    simp only [Pi.mul_apply, u', v, Real.rpow_one]
    rw [show -c * x = -(x * c) by ring]; ring
  have h0 : Tendsto (u * v) (𝓝[>] 0) (𝓝 0) := by
    have hb : Tendsto (fun t : ℝ => 2 * c * t ^ (1 / 2 : ℝ)) (𝓝[>] 0) (𝓝 0) := by
      have : Tendsto (fun t : ℝ => t ^ (1 / 2 : ℝ)) (𝓝[>] 0) (𝓝 0) := by
        have := (Real.continuousAt_rpow_const 0 (1 / 2 : ℝ) (Or.inr (by norm_num))).tendsto
        rw [Real.zero_rpow (by norm_num)] at this
        exact this.mono_left nhdsWithin_le_nhds
      simpa using this.const_mul (2 * c)
    refine squeeze_zero_norm' ?_ hb
    filter_upwards [self_mem_nhdsWithin] with t ht
    have ht0 : (0 : ℝ) < t := ht
    simp only [Pi.mul_apply, u, v, Real.norm_eq_abs]
    have h1 := one_sub_exp_neg_nn (mul_nonneg ht0.le hc.le)
    have h2 := one_sub_exp_neg_le (x := t * c)
    have hp : 0 < t ^ (-(1 / 2 : ℝ)) := Real.rpow_pos_of_pos ht0 _
    have e : t ^ (1 / 2 : ℝ) = t * t ^ (-(1 / 2 : ℝ)) := by
      conv_lhs => rw [show (1 / 2 : ℝ) = 1 + (-(1 / 2)) by norm_num]
      rw [Real.rpow_add ht0, Real.rpow_one]
    rw [e, abs_mul, abs_of_nonneg h1, abs_mul, abs_of_pos hp]
    norm_num
    nlinarith
  have hinf : Tendsto (u * v) atTop (𝓝 0) := by
    have hb : Tendsto (fun t : ℝ => 2 * t ^ (-(1 / 2 : ℝ))) atTop (𝓝 0) := by
      simpa using (tendsto_rpow_neg_atTop (y := (1 / 2 : ℝ)) (by norm_num)).const_mul 2
    refine squeeze_zero_norm' ?_ hb
    filter_upwards [eventually_gt_atTop 0] with t ht0
    simp only [Pi.mul_apply, u, v, Real.norm_eq_abs]
    have h1 := one_sub_exp_neg_nn (mul_nonneg ht0.le hc.le)
    have h2 := one_sub_exp_neg_le_one (t * c)
    have hp : 0 < t ^ (-(1 / 2 : ℝ)) := Real.rpow_pos_of_pos ht0 _
    rw [abs_mul, abs_of_nonneg h1, abs_mul, abs_of_pos hp]
    norm_num
    nlinarith
  have key := integral_Ioi_mul_deriv_eq_deriv_mul hu hv huv' hu'v h0 hinf
  simp only [u, v', u', v] at key
  rw [key]
  have e2 : ∫ x in Ioi (0 : ℝ), c * Real.exp (-(x * c)) * (-2 * x ^ (-(1 / 2 : ℝ))) =
      -2 * c * ∫ x in Ioi (0 : ℝ), x ^ ((1 / 2 : ℝ) - 1) * Real.exp (-(c * x)) := by
    rw [← integral_const_mul]
    refine setIntegral_congr_fun measurableSet_Ioi (fun x _ => ?_)
    rw [show (1 / 2 : ℝ) - 1 = -(1 / 2) by norm_num, mul_comm c x]; ring
  rw [e2, integral_rpow_mul_exp_neg_mul_Ioi (by norm_num) hc, Real.Gamma_one_half_eq,
    Real.div_rpow zero_le_one hc.le, Real.one_rpow, ← Real.sqrt_eq_rpow]
  have hs : 0 < Real.sqrt c := Real.sqrt_pos.mpr hc
  have hsq : Real.sqrt c ^ 2 = c := Real.sq_sqrt hc.le
  field_simp
  rw [hsq]; ring

lemma gy_pt (c : ℝ) {x : ℝ} (hx : 0 < x) :
    (|(-2 : ℝ)| * x ^ ((-2 : ℝ) - 1)) •
        ((fun t : ℝ => (1 - Real.exp (-(t * c))) * t ^ (-(3 / 2 : ℝ))) (x ^ (-2 : ℝ))) =
      2 * (1 - Real.exp (-(c / x ^ 2))) := by
  simp only [smul_eq_mul]
  have h1 : x ^ (-2 : ℝ) * c = c / x ^ 2 := by
    rw [Real.rpow_neg hx.le, Real.rpow_two]; field_simp
  have h2 : (x ^ (-2 : ℝ)) ^ (-(3 / 2 : ℝ)) = x ^ (3 : ℝ) := by
    rw [← Real.rpow_mul hx.le]; norm_num
  have h3 : x ^ ((-2 : ℝ) - 1) * x ^ (3 : ℝ) = 1 := by
    rw [← Real.rpow_add hx]; norm_num
  rw [h1, h2, abs_of_neg (by norm_num : (-2 : ℝ) < 0)]
  linear_combination (2 * (1 - Real.exp (-(c / x ^ 2)))) * h3

lemma gy {c : ℝ} (hc : 0 ≤ c) :
    IntegrableOn (fun y : ℝ => 1 - Real.exp (-(c / y ^ 2))) (Ioi 0) ∧
      ∫ y in Ioi (0 : ℝ), (1 - Real.exp (-(c / y ^ 2))) = Real.sqrt π * Real.sqrt c := by
  set g : ℝ → ℝ := fun t => (1 - Real.exp (-(t * c))) * t ^ (-(3 / 2 : ℝ)) with hg
  have hI := (integrableOn_Ioi_comp_rpow_iff g (p := -2) (by norm_num)).mpr (gt_int hc)
  have hV := integral_comp_rpow_Ioi g (p := -2) (by norm_num)
  have hI2 : IntegrableOn (fun y : ℝ => 2 * (1 - Real.exp (-(c / y ^ 2)))) (Ioi 0) :=
    hI.congr_fun (fun x hx => gy_pt c hx) measurableSet_Ioi
  rw [setIntegral_congr_fun measurableSet_Ioi (fun x hx => gy_pt c hx)] at hV
  refine ⟨?_, ?_⟩
  · have := hI2.const_mul (1 / 2)
    simp only [← mul_assoc, one_div, inv_mul_cancel₀ (two_ne_zero' ℝ), one_mul] at this
    exact this
  · rw [integral_const_mul] at hV
    rcases hc.lt_or_eq with hc | hc
    · rw [hg, gt_val hc] at hV; linarith
    · subst hc; simp

lemma hI_meas (ε y : ℝ) : Measurable (hIntegrand ε y) := by
  unfold hIntegrand; fun_prop

lemma hI_nn (ε y : ℝ) {x : ℝ} (hx : 0 < x) : 0 ≤ hIntegrand ε y x := by
  unfold hIntegrand
  have := one_sub_exp_neg_nn (x := x ^ (2 - 2 * ε) / y ^ 2) (by positivity)
  positivity

lemma hI_le (ε : ℝ) {y x : ℝ} (hy : 0 < y) (hx : 0 < x) :
    hIntegrand ε y x ≤ 1 / y ^ 2 * (Real.exp (-x) * x ^ ((2 - ε) - 1)) := by
  unfold hIntegrand
  have h1 := one_sub_exp_neg_le (x := x ^ (2 - 2 * ε) / y ^ 2)
  have hB : 0 < x ^ (ε - 1) := Real.rpow_pos_of_pos hx _
  have hAB : x ^ (2 - 2 * ε) * x ^ (ε - 1) = x ^ ((2 - ε) - 1) := by
    rw [← Real.rpow_add hx]; ring_nf
  calc (1 - Real.exp (-(x ^ (2 - 2 * ε) / y ^ 2))) * x ^ (ε - 1) * Real.exp (-x)
      ≤ (x ^ (2 - 2 * ε) / y ^ 2) * x ^ (ε - 1) * Real.exp (-x) := by gcongr
    _ = _ := by rw [← hAB]; field_simp

lemma hI_int {ε : ℝ} (hε : ε < 2) {y : ℝ} (hy : 0 < y) :
    IntegrableOn (hIntegrand ε y) (Ioi 0) := by
  have hg := (Real.GammaIntegral_convergent (s := 2 - ε) (by linarith)).const_mul (1 / y ^ 2)
  refine Integrable.mono' hg (hI_meas ε y).aestronglyMeasurable ?_
  refine (ae_restrict_iff' measurableSet_Ioi).mpr (Eventually.of_forall fun x hx => ?_)
  rw [Real.norm_eq_abs, abs_of_nonneg (hI_nn ε y hx)]
  exact hI_le ε hy hx

lemma h_nn (ε y : ℝ) : 0 ≤ h ε y :=
  setIntegral_nonneg measurableSet_Ioi (fun x hx => hI_nn ε y hx)

/-- `h_ε(y) ≤ 2/y²` for `0 ≤ ε ≤ 1`. -/
lemma h_le_inv_sq {ε : ℝ} (h0 : 0 ≤ ε) (h1 : ε ≤ 1) {y : ℝ} (hy : 0 < y) :
    h ε y ≤ 2 / y ^ 2 := by
  obtain ⟨i0, v0⟩ := Y05.gamma1 (x := 0) (by norm_num)
  obtain ⟨i1, v1⟩ := Y05.gamma1 (x := 1) (by norm_num)
  have hi := (i0.add i1).const_mul (1 / y ^ 2)
  have hv : ∫ x in Ioi (0 : ℝ), 1 / y ^ 2 * (Real.exp (-x) * x ^ (0 : ℝ) +
      Real.exp (-x) * x ^ (1 : ℝ)) = 2 / y ^ 2 := by
    rw [integral_const_mul, integral_add i0 i1, v0, v1]
    norm_num [Real.Gamma_two]; ring
  rw [← hv]
  refine setIntegral_mono_on (hI_int (by linarith) hy) hi measurableSet_Ioi (fun x hx => ?_)
  have hx0 : (0 : ℝ) < x := hx
  refine (hI_le ε hy hx0).trans ?_
  have hp : x ^ ((2 - ε) - 1) ≤ x ^ (0 : ℝ) + x ^ (1 : ℝ) := by
    rcases le_total x 1 with hx1 | hx1
    · have ha := Real.rpow_le_one hx0.le hx1 (show 0 ≤ (2 - ε) - 1 by linarith)
      have hb := Real.rpow_nonneg hx0.le (1 : ℝ)
      rw [Real.rpow_zero]; exact le_add_of_le_of_nonneg ha hb
    · have ha := Real.rpow_le_rpow_of_exponent_le hx1 (show (2 - ε) - 1 ≤ 1 by linarith)
      have hb := Real.rpow_nonneg hx0.le (0 : ℝ)
      exact le_add_of_nonneg_of_le hb ha
  have he := Real.exp_pos (-x)
  have : Real.exp (-x) * x ^ ((2 - ε) - 1) ≤ Real.exp (-x) * x ^ (0 : ℝ) +
      Real.exp (-x) * x ^ (1 : ℝ) := by nlinarith
  exact mul_le_mul_of_nonneg_left this (by positivity)

/-- the pointwise bound on `(0,1]`. -/
lemma hI_le_k {ε : ℝ} (h0 : 0 ≤ ε) {y x : ℝ} (hy : 0 < y) (hx : 0 < x) (hx1 : x ≤ 1) :
    hIntegrand ε y x ≤ 2 * x ^ (1 - 2 * ε) / (y ^ 2 + x ^ (2 - 2 * ε)) := by
  unfold hIntegrand
  have hA : 0 ≤ x ^ (2 - 2 * ε) := Real.rpow_nonneg hx.le _
  have hB : 0 < x ^ (ε - 1) := Real.rpow_pos_of_pos hx _
  have hAB : x ^ (2 - 2 * ε) * x ^ (ε - 1) = x ^ (1 - ε) := by
    rw [← Real.rpow_add hx]; ring_nf
  have hCD : x ^ (1 - ε) ≤ x ^ (1 - 2 * ε) :=
    Real.rpow_le_rpow_of_exponent_ge hx hx1 (by linarith)
  have hy2 : 0 < y ^ 2 := by positivity
  have hE1 := one_sub_exp_neg_le (x := x ^ (2 - 2 * ε) / y ^ 2)
  have hE2 := one_sub_exp_neg_le_one (x ^ (2 - 2 * ε) / y ^ 2)
  have hE0 := one_sub_exp_neg_nn (x := x ^ (2 - 2 * ε) / y ^ 2) (by positivity)
  have hex : Real.exp (-x) ≤ 1 := Real.exp_le_one_iff.mpr (by linarith)
  generalize 1 - Real.exp (-(x ^ (2 - 2 * ε) / y ^ 2)) = E at *
  generalize x ^ (2 - 2 * ε) = A at *
  generalize x ^ (ε - 1) = B at *
  generalize x ^ (1 - ε) = C at *
  generalize x ^ (1 - 2 * ε) = D at *
  have hEy : E * y ^ 2 ≤ A := by
    have := mul_le_mul_of_nonneg_right hE1 hy2.le
    rwa [div_mul_cancel₀ _ hy2.ne'] at this
  rw [le_div_iff₀ (by positivity)]
  have hEB : E * B * Real.exp (-x) ≤ E * B := by
    have : 0 ≤ E * B := by positivity
    nlinarith [Real.exp_pos (-x)]
  have key : E * B * (y ^ 2 + A) ≤ 2 * C := by
    have : E * (y ^ 2 + A) ≤ 2 * A := by nlinarith
    calc E * B * (y ^ 2 + A) = E * (y ^ 2 + A) * B := by ring
      _ ≤ 2 * A * B := mul_le_mul_of_nonneg_right this hB.le
      _ = 2 * C := by rw [← hAB]; ring
  have : 0 ≤ y ^ 2 + A := by positivity
  nlinarith

lemma k_int {ε : ℝ} (h1 : ε ≤ 1 / 2) {y : ℝ} (hy : 0 < y) :
    IntervalIntegrable (fun x : ℝ => 2 * x ^ (1 - 2 * ε) / (y ^ 2 + x ^ (2 - 2 * ε))) volume 0 1 ∧
    ∫ x in (0 : ℝ)..1, 2 * x ^ (1 - 2 * ε) / (y ^ 2 + x ^ (2 - 2 * ε)) =
      (Real.log (y ^ 2 + 1) - Real.log (y ^ 2)) / (1 - ε) := by
  have hy2 : 0 < y ^ 2 := by positivity
  have hc1 : Continuous (fun x : ℝ => x ^ (1 - 2 * ε)) := Real.continuous_rpow_const (by linarith)
  have hc2 : Continuous (fun x : ℝ => x ^ (2 - 2 * ε)) := Real.continuous_rpow_const (by linarith)
  have hden : ∀ x ∈ Icc (0 : ℝ) 1, y ^ 2 + x ^ (2 - 2 * ε) ≠ 0 := fun x hx => by
    have := Real.rpow_nonneg hx.1 (2 - 2 * ε); positivity
  have hk : ContinuousOn (fun x : ℝ => 2 * x ^ (1 - 2 * ε) / (y ^ 2 + x ^ (2 - 2 * ε)))
      (Icc 0 1) :=
    ((continuous_const.mul hc1).continuousOn).div ((continuous_const.add hc2).continuousOn) hden
  have hint : IntervalIntegrable (fun x : ℝ => 2 * x ^ (1 - 2 * ε) / (y ^ 2 + x ^ (2 - 2 * ε)))
      volume 0 1 := (hk.mono (by rw [uIcc_of_le zero_le_one])).intervalIntegrable
  refine ⟨hint, ?_⟩
  have hF : ContinuousOn (fun x : ℝ => Real.log (y ^ 2 + x ^ (2 - 2 * ε)) / (1 - ε)) (Icc 0 1) :=
    (((continuous_const.add hc2).continuousOn).log hden).div_const _
  have hd : ∀ x ∈ Ioo (0 : ℝ) 1, HasDerivAt (fun x : ℝ => Real.log (y ^ 2 + x ^ (2 - 2 * ε)) /
      (1 - ε)) (2 * x ^ (1 - 2 * ε) / (y ^ 2 + x ^ (2 - 2 * ε))) x := by
    intro x hx
    have hx0 : 0 < x := hx.1
    have h := (((Real.hasDerivAt_rpow_const (p := 2 - 2 * ε) (Or.inl hx0.ne')).const_add
      (y ^ 2)).log (hden x ⟨hx0.le, hx.2.le⟩)).div_const (1 - ε)
    convert h using 1
    rw [show 2 - 2 * ε - 1 = 1 - 2 * ε by ring]
    have : (1 - ε) ≠ 0 := by intro h; linarith
    have := hden x ⟨hx0.le, hx.2.le⟩
    field_simp
  rw [intervalIntegral.integral_eq_sub_of_hasDerivAt_of_le zero_le_one hF hd hint]
  simp only [Real.one_rpow]
  rw [Real.zero_rpow (by linarith), add_zero, sub_div]

lemma log_diff_le {y : ℝ} (hy : 0 < y) :
    Real.log (y ^ 2 + 1) - Real.log (y ^ 2) ≤ Real.log 2 + 2 * |Real.log y| := by
  have hy2 : 0 < y ^ 2 := by positivity
  have hl : Real.log (y ^ 2) = 2 * Real.log y := by
    rw [Real.log_pow]; norm_num
  rcases le_total 1 y with h1 | h1
  · have : Real.log (y ^ 2 + 1) ≤ Real.log (2 * y ^ 2) :=
      Real.log_le_log (by positivity) (by nlinarith)
    rw [Real.log_mul (by norm_num) hy2.ne'] at this
    have := abs_nonneg (Real.log y)
    linarith
  · have : Real.log (y ^ 2 + 1) ≤ Real.log 2 :=
      Real.log_le_log (by positivity) (by nlinarith)
    have hn : Real.log y ≤ 0 := Real.log_nonpos hy.le h1
    rw [abs_of_nonpos hn]
    linarith

lemma h_le_log {ε : ℝ} (h0 : 0 ≤ ε) (h1 : ε ≤ 1 / 2) {y : ℝ} (hy : 0 < y) :
    h ε y ≤ 4 * (1 + |Real.log y|) := by
  have hint := hI_int (ε := ε) (by linarith) hy
  have hsplit : Ioi (0 : ℝ) = Ioc 0 1 ∪ Ioi 1 := (Ioc_union_Ioi_eq_Ioi zero_le_one).symm
  have hU : h ε y = (∫ x in Ioc (0 : ℝ) 1, hIntegrand ε y x) + ∫ x in Ioi (1 : ℝ), hIntegrand ε y x := by
    unfold h
    rw [hsplit] at hint ⊢
    exact setIntegral_union (Ioc_disjoint_Ioi le_rfl) measurableSet_Ioi
      (hint.mono_set subset_union_left) (hint.mono_set subset_union_right)
  obtain ⟨ki, kv⟩ := k_int h1 hy
  have hA : ∫ x in Ioc (0 : ℝ) 1, hIntegrand ε y x ≤
      (Real.log (y ^ 2 + 1) - Real.log (y ^ 2)) / (1 - ε) := by
    rw [← kv, intervalIntegral.integral_of_le zero_le_one]
    refine setIntegral_mono_on (hint.mono_set (by rw [hsplit]; exact subset_union_left))
      ((intervalIntegrable_iff_integrableOn_Ioc_of_le zero_le_one).mp ki) measurableSet_Ioc
      (fun x hx => hI_le_k h0 hy hx.1 hx.2)
  obtain ⟨g0, gv⟩ := Y05.gamma1 (x := 0) (by norm_num)
  have hB : ∫ x in Ioi (1 : ℝ), hIntegrand ε y x ≤ 1 := by
    calc ∫ x in Ioi (1 : ℝ), hIntegrand ε y x ≤ ∫ x in Ioi (1 : ℝ), Real.exp (-x) * x ^ (0 : ℝ) := by
          refine setIntegral_mono_on (hint.mono_set (Ioi_subset_Ioi zero_le_one))
            (g0.mono_set (Ioi_subset_Ioi zero_le_one)) measurableSet_Ioi (fun x hx => ?_)
          have hx1 : (1 : ℝ) < x := hx
          unfold hIntegrand
          rw [Real.rpow_zero, mul_one]
          have e1 := one_sub_exp_neg_le_one (x ^ (2 - 2 * ε) / y ^ 2)
          have e0 := one_sub_exp_neg_nn (x := x ^ (2 - 2 * ε) / y ^ 2)
            (by have := Real.rpow_nonneg (show (0:ℝ) ≤ x by linarith) (2 - 2 * ε); positivity)
          have hb := Real.rpow_le_one_of_one_le_of_nonpos hx1.le (show ε - 1 ≤ 0 by linarith)
          have hb0 := Real.rpow_nonneg (show (0:ℝ) ≤ x by linarith) (ε - 1)
          have := Real.exp_pos (-x)
          calc (1 - Real.exp (-(x ^ (2 - 2 * ε) / y ^ 2))) * x ^ (ε - 1) * Real.exp (-x)
              ≤ 1 * 1 * Real.exp (-x) := by gcongr
            _ = _ := by ring
      _ ≤ ∫ x in Ioi (0 : ℝ), Real.exp (-x) * x ^ (0 : ℝ) :=
          setIntegral_mono_set g0 (Eventually.of_forall fun x => by
            simp only [Pi.zero_apply, Real.rpow_zero, mul_one]; exact (Real.exp_pos _).le)
            (Eventually.of_forall (Ioi_subset_Ioi zero_le_one))
      _ = 1 := by rw [gv]; simp
  have hL := log_diff_le hy
  have hl2 := Real.log_two_lt_d9
  have hpos : 0 ≤ Real.log (y ^ 2 + 1) - Real.log (y ^ 2) := by
    have := Real.log_le_log (by positivity : 0 < y ^ 2) (by linarith : y ^ 2 ≤ y ^ 2 + 1)
    linarith
  have hdiv : (Real.log (y ^ 2 + 1) - Real.log (y ^ 2)) / (1 - ε) ≤
      2 * (Real.log (y ^ 2 + 1) - Real.log (y ^ 2)) := by
    rw [div_le_iff₀ (by linarith)]; nlinarith
  have := abs_nonneg (Real.log y)
  rw [hU]; nlinarith

end Z04

theorem Z04.proof_1 : Z04.stmt_1 := by
  intro r hr
  refine ⟨Z04.gt_int (sq_nonneg r), ?_⟩
  rcases hr.lt_or_eq with hr | hr
  · rw [Z04.gt_val (by positivity), Real.sqrt_sq hr.le]
    have : 0 < Real.sqrt π := Real.sqrt_pos.mpr Real.pi_pos
    field_simp
  · subst hr
    simp

theorem Z04.proof_2 : Z04.stmt_2 := by
  refine ⟨4, fun ε h0 h1 y hy => ⟨Z04.hI_int (by linarith) hy, ?_, fun hy1 => ?_⟩⟩
  · rw [abs_of_nonneg (Z04.h_nn ε y)]; exact Z04.h_le_log h0 (by linarith) hy
  · rw [abs_of_nonneg (Z04.h_nn ε y)]
    refine (Z04.h_le_inv_sq h0 (by linarith) hy).trans ?_
    have : 0 < y ^ 2 := by positivity
    gcongr; norm_num

namespace Z04

/-- the `y`-integral of the `h`-integrand at fixed `x`. -/
lemma hI_y {ε : ℝ} {x : ℝ} (hx : 0 < x) :
    IntegrableOn (fun y => hIntegrand ε y x) (Ioi 0) ∧
      ∫ y in Ioi (0 : ℝ), hIntegrand ε y x = Real.sqrt π * Real.exp (-x) := by
  obtain ⟨gi, gv⟩ := gy (c := x ^ (2 - 2 * ε)) (Real.rpow_nonneg hx.le _)
  have e : (fun y => hIntegrand ε y x) =
      fun y => x ^ (ε - 1) * Real.exp (-x) * (1 - Real.exp (-(x ^ (2 - 2 * ε) / y ^ 2))) := by
    funext y; unfold hIntegrand; ring
  rw [e]
  refine ⟨gi.const_mul _, ?_⟩
  rw [integral_const_mul, gv, Real.sqrt_eq_rpow (x ^ (2 - 2 * ε)), ← Real.rpow_mul hx.le]
  have : x ^ (ε - 1) * x ^ ((2 - 2 * ε) * (1 / 2)) = 1 := by
    rw [← Real.rpow_add hx]; ring_nf; simp
  linear_combination (Real.exp (-x) * Real.sqrt π) * this

lemma hF_int {ε : ℝ} (hε : ε < 2) :
    Integrable (fun z : ℝ × ℝ => hIntegrand ε z.2 z.1)
      ((volume.restrict (Ioi (0 : ℝ))).prod (volume.restrict (Ioi (0 : ℝ)))) := by
  have hm : AEStronglyMeasurable (fun z : ℝ × ℝ => hIntegrand ε z.2 z.1)
      ((volume.restrict (Ioi (0 : ℝ))).prod (volume.restrict (Ioi (0 : ℝ)))) := by
    unfold hIntegrand; fun_prop
  rw [integrable_prod_iff hm]
  refine ⟨(ae_restrict_iff' measurableSet_Ioi).mpr (Eventually.of_forall fun x hx => (hI_y hx).1), ?_⟩
  have hg := (Real.GammaIntegral_convergent (s := 1) one_pos).const_mul (Real.sqrt π)
  refine IntegrableOn.congr_fun hg (fun x hx => ?_) measurableSet_Ioi
  have hx0 : (0 : ℝ) < x := hx
  simp only
  rw [setIntegral_congr_fun measurableSet_Ioi (fun y hy => by
      rw [Real.norm_eq_abs, abs_of_nonneg (hI_nn ε y hx0)]), (hI_y hx0).2]
  simp

lemma h_int {ε : ℝ} (hε : ε < 2) : IntegrableOn (h ε) (Ioi 0) :=
  (hF_int hε).integral_prod_right

lemma h_integral {ε : ℝ} (hε : ε < 2) :
    ∫ y in Ioi (0 : ℝ), h ε y = Real.sqrt π := by
  unfold h
  have := integral_integral_swap (f := fun x y => hIntegrand ε y x) (hF_int hε)
  rw [← this, setIntegral_congr_fun measurableSet_Ioi (fun x hx => (hI_y hx).2),
    integral_const_mul]
  have := (Y05.gamma1 (x := 0) (by norm_num)).2
  simp only [Real.rpow_zero, mul_one, add_zero, Real.Gamma_one] at this
  rw [this, mul_one]

end Z04

theorem Z04.proof_3 : Z04.stmt_3 := by
  intro ε h0 h1
  refine ⟨fun y hy => Z04.hI_int (by linarith) hy, Z04.h_int (by linarith), ?_⟩
  unfold Z04.H
  simp only [pow_one]
  rw [Z04.h_integral (by linarith)]
  have : 0 < Real.sqrt π := Real.sqrt_pos.mpr Real.pi_pos
  field_simp

namespace Z04

/-- `u₀(x) = e^{-x} x^{-τ}`. -/
noncomputable def u0 (τ x : ℝ) : ℝ := Real.exp (-x) * x ^ (-τ)

/-- `v(x) = u₀(x) - (h-integrand)(x) = u₀(x) e^{-x^{2τ}/y²}`. -/
noncomputable def vv (ε y x : ℝ) : ℝ := u0 (1 - ε) x - hIntegrand ε y x

lemma u0_int {ε : ℝ} (hε : 0 < ε) :
    IntegrableOn (u0 (1 - ε)) (Ioi 0) ∧ ∫ x in Ioi (0 : ℝ), u0 (1 - ε) x = Real.Gamma ε := by
  obtain ⟨i, v⟩ := Y05.gamma1 (x := -(1 - ε)) (by linarith)
  exact ⟨i, by rw [show ∫ x in Ioi (0 : ℝ), u0 (1 - ε) x = _ from v]; congr 1; ring⟩

lemma vv_int {ε : ℝ} (hε : 0 < ε) (hε1 : ε ≤ 1) {y : ℝ} (hy : 0 < y) :
    IntegrableOn (vv ε y) (Ioi 0) ∧
      ∫ x in Ioi (0 : ℝ), vv ε y x = Real.Gamma ε - h ε y := by
  obtain ⟨i0, v0⟩ := u0_int hε
  have i1 := hI_int (ε := ε) (by linarith) hy
  refine ⟨i0.sub i1, ?_⟩
  unfold vv; rw [integral_sub i0 i1, v0]; rfl

lemma h_le_gamma {ε : ℝ} (hε : 0 < ε) (hε1 : ε ≤ 1) {y : ℝ} (hy : 0 < y) : h ε y ≤ Real.Gamma ε := by
  obtain ⟨i0, v0⟩ := u0_int hε
  rw [← v0]
  refine setIntegral_mono_on (hI_int (by linarith) hy) i0 measurableSet_Ioi (fun x hx => ?_)
  have hx0 : (0 : ℝ) < x := hx
  unfold hIntegrand u0
  have e1 := one_sub_exp_neg_le_one (x ^ (2 - 2 * ε) / y ^ 2)
  have hb : 0 < x ^ (ε - 1) := Real.rpow_pos_of_pos hx0 _
  rw [show -(1 - ε) = ε - 1 by ring]
  have := Real.exp_pos (-x)
  calc (1 - Real.exp (-(x ^ (2 - 2 * ε) / y ^ 2))) * x ^ (ε - 1) * Real.exp (-x)
      ≤ 1 * x ^ (ε - 1) * Real.exp (-x) := by gcongr
    _ = _ := by ring

lemma prod3 {f1 f2 f3 : ℝ → ℝ} (i1 : IntegrableOn f1 (Ioi 0)) (i2 : IntegrableOn f2 (Ioi 0))
    (i3 : IntegrableOn f3 (Ioi 0)) :
    IntegrableOn (fun g : ℝ × ℝ × ℝ => f1 g.1 * f2 g.2.1 * f3 g.2.2) Y05.Q ∧
      ∫ g in Y05.Q, f1 g.1 * f2 g.2.1 * f3 g.2.2 =
        (∫ x in Ioi (0 : ℝ), f1 x) * (∫ x in Ioi (0 : ℝ), f2 x) * ∫ x in Ioi (0 : ℝ), f3 x := by
  have e : (fun g : ℝ × ℝ × ℝ => f1 g.1 * f2 g.2.1 * f3 g.2.2) =
      fun g => f1 g.1 * ((fun k : ℝ × ℝ => f2 k.1 * f3 k.2) g.2) := by
    funext g; simp only; ring
  have i23 : IntegrableOn (fun k : ℝ × ℝ => f2 k.1 * f3 k.2) (Ioi (0 : ℝ) ×ˢ Ioi (0 : ℝ))
      (volume.prod volume) := by
    rw [IntegrableOn, ← Measure.prod_restrict]
    exact Integrable.mul_prod i2 i3
  rw [e]
  constructor
  · show IntegrableOn _ _ (volume.prod volume)
    rw [IntegrableOn, ← Measure.prod_restrict]
    exact Integrable.mul_prod i1 i23
  · refine Eq.trans (setIntegral_prod_mul (μ := (volume : Measure ℝ))
      (ν := (volume : Measure (ℝ × ℝ))) f1 (fun k : ℝ × ℝ => f2 k.1 * f3 k.2)
      (Ioi (0 : ℝ)) (Ioi (0 : ℝ) ×ˢ Ioi (0 : ℝ))) ?_
    refine Eq.trans (congrArg _ (setIntegral_prod_mul (μ := (volume : Measure ℝ))
      (ν := (volume : Measure ℝ)) f2 f3 (Ioi (0 : ℝ)) (Ioi (0 : ℝ)))) ?_
    ring

/-- the pointwise expansion `w(g)(1 - e^{-|g^τ|²/y²}) = u₀u₀u₀ - vvv`. -/
lemma pt3 {ε y : ℝ} {g : ℝ × ℝ × ℝ} (hg : g ∈ Y05.Q) :
    Real.exp (-(g.1 + g.2.1 + g.2.2)) * (g.1 * g.2.1 * g.2.2) ^ (-(1 - ε)) *
        (1 - Real.exp (-((g.1 ^ (2 * (1 - ε)) + g.2.1 ^ (2 * (1 - ε)) + g.2.2 ^ (2 * (1 - ε)))
          / y ^ 2))) =
      u0 (1 - ε) g.1 * u0 (1 - ε) g.2.1 * u0 (1 - ε) g.2.2 -
        vv ε y g.1 * vv ε y g.2.1 * vv ε y g.2.2 := by
  obtain ⟨ha, hb, hc⟩ := hg
  simp only [mem_Ioi] at ha hb hc
  unfold vv u0 hIntegrand
  rw [Real.mul_rpow (by positivity) hc.le, Real.mul_rpow ha.le hb.le,
    show (2 - 2 * ε) = 2 * (1 - ε) by ring, show ε - 1 = -(1 - ε) by ring,
    add_div, add_div, neg_add, neg_add, Real.exp_add, Real.exp_add, neg_add, neg_add,
    Real.exp_add, Real.exp_add]
  ring

end Z04

namespace Z04

/-- `w(g) |g^τ|`, the Laplace-side integrand of `S` (`= √3 · Y05.Fb`). -/
noncomputable def W3 (ε : ℝ) (g : ℝ × ℝ × ℝ) : ℝ :=
  Real.exp (-(g.1 + g.2.1 + g.2.2)) * (g.1 * g.2.1 * g.2.2) ^ (-(1 - ε)) *
    Real.sqrt (g.1 ^ (2 * (1 - ε)) + g.2.1 ^ (2 * (1 - ε)) + g.2.2 ^ (2 * (1 - ε)))

/-- `w(g) (1 - e^{-|g^τ|²/y²})`. -/
noncomputable def Phi3 (ε : ℝ) (z : (ℝ × ℝ × ℝ) × ℝ) : ℝ :=
  Real.exp (-(z.1.1 + z.1.2.1 + z.1.2.2)) * (z.1.1 * z.1.2.1 * z.1.2.2) ^ (-(1 - ε)) *
    (1 - Real.exp (-((z.1.1 ^ (2 * (1 - ε)) + z.1.2.1 ^ (2 * (1 - ε)) + z.1.2.2 ^ (2 * (1 - ε)))
      / z.2 ^ 2)))

lemma W3_eq {ε : ℝ} {g : ℝ × ℝ × ℝ} (hg : g ∈ Y05.Q) :
    W3 ε g = Real.sqrt 3 * Y05.Fb (1 - ε) g := by
  rw [Y05.Fb_form (1 - ε) hg, W3]
  have : 0 < Real.sqrt 3 := by positivity
  field_simp

lemma W3_int {ε : ℝ} (hε : 0 < ε) : IntegrableOn (W3 ε) Y05.Q := by
  have := (Y05.Fb_int (τ := 1 - ε) (by linarith)).const_mul (Real.sqrt 3)
  exact IntegrableOn.congr_fun this (fun g hg => (W3_eq hg).symm) Y05.Q_meas

lemma Phi3_y {ε : ℝ} {g : ℝ × ℝ × ℝ} (hg : g ∈ Y05.Q) :
    IntegrableOn (fun y => Phi3 ε (g, y)) (Ioi 0) ∧
      ∫ y in Ioi (0 : ℝ), Phi3 ε (g, y) = Real.sqrt π * W3 ε g := by
  obtain ⟨ha, hb, hc⟩ := hg
  simp only [mem_Ioi] at ha hb hc
  obtain ⟨gi, gv⟩ := gy (c := g.1 ^ (2 * (1 - ε)) + g.2.1 ^ (2 * (1 - ε)) + g.2.2 ^ (2 * (1 - ε)))
    (by positivity)
  unfold Phi3 W3
  simp only
  refine ⟨gi.const_mul _, ?_⟩
  rw [integral_const_mul, gv]; ring

lemma Phi3_int {ε : ℝ} (hε : 0 < ε) :
    Integrable (Phi3 ε) ((volume.restrict Y05.Q).prod (volume.restrict (Ioi (0 : ℝ)))) := by
  have hm : AEStronglyMeasurable (Phi3 ε)
      ((volume.restrict Y05.Q).prod (volume.restrict (Ioi (0 : ℝ)))) := by
    unfold Phi3; fun_prop
  rw [integrable_prod_iff hm]
  refine ⟨(ae_restrict_iff' Y05.Q_meas).mpr (Eventually.of_forall fun g hg => (Phi3_y hg).1), ?_⟩
  have := (W3_int hε).const_mul (Real.sqrt π)
  refine IntegrableOn.congr_fun this (fun g hg => ?_) Y05.Q_meas
  have hg' := hg
  obtain ⟨ha, hb, hc⟩ := hg
  simp only [mem_Ioi] at ha hb hc
  rw [setIntegral_congr_fun measurableSet_Ioi (fun y hy => by
      rw [Real.norm_eq_abs, abs_of_nonneg (by
        unfold Phi3
        have := one_sub_exp_neg_nn (x := (g.1 ^ (2 * (1 - ε)) + g.2.1 ^ (2 * (1 - ε)) +
          g.2.2 ^ (2 * (1 - ε))) / y ^ 2) (by positivity)
        positivity)]), (Phi3_y hg').2]

lemma Phi3_g {ε : ℝ} (hε : 0 < ε) (hε1 : ε ≤ 1) {y : ℝ} (hy : 0 < y) :
    ∫ g in Y05.Q, Phi3 ε (g, y) = Real.Gamma ε ^ 3 - (Real.Gamma ε - h ε y) ^ 3 := by
  obtain ⟨i0, v0⟩ := u0_int hε
  obtain ⟨i1, v1⟩ := vv_int hε hε1 hy
  obtain ⟨j0, w0⟩ := prod3 i0 i0 i0
  obtain ⟨j1, w1⟩ := prod3 i1 i1 i1
  rw [setIntegral_congr_fun Y05.Q_meas (fun g hg => by unfold Phi3; exact pt3 hg),
    integral_sub j0 j1, w0, w1, v0, v1]
  ring

lemma hpow_int {ε : ℝ} (hε : 0 < ε) (hε1 : ε ≤ 1) (m : ℕ) :
    IntegrableOn (fun y => h ε y ^ (m + 1)) (Ioi 0) := by
  have hi := (h_int (ε := ε) (by linarith)).const_mul (Real.Gamma ε ^ m)
  refine Integrable.mono' hi ((h_int (ε := ε) (by linarith)).aestronglyMeasurable.pow _) ?_
  refine (ae_restrict_iff' measurableSet_Ioi).mpr (Eventually.of_forall fun y hy => ?_)
  have h0 := h_nn ε y
  have h1 := h_le_gamma hε hε1 hy
  rw [Real.norm_eq_abs, abs_of_nonneg (by positivity), pow_succ]
  exact mul_le_mul_of_nonneg_right (pow_le_pow_left₀ h0 h1 m) h0

/-- the Laplace side of `S`: `∫_Q w |g^τ| = 3Γ(ε)² - 3Γ(ε)H₂ + H₃`. -/
lemma W3_integral {ε : ℝ} (hε : 0 < ε) (hε1 : ε ≤ 1) :
    ∫ g in Y05.Q, W3 ε g =
      3 * Real.Gamma ε ^ 2 - 3 * Real.Gamma ε * H 2 ε + H 3 ε := by
  have hsw := integral_integral_swap (f := fun g y => Phi3 ε (g, y)) (Phi3_int hε)
  rw [setIntegral_congr_fun Y05.Q_meas (fun g hg => (Phi3_y hg).2), integral_const_mul,
    setIntegral_congr_fun measurableSet_Ioi (fun y hy => Phi3_g hε hε1 hy)] at hsw
  have i1 := h_int (ε := ε) (by linarith)
  have i2 : IntegrableOn (fun y => h ε y ^ 2) (Ioi 0) := hpow_int hε hε1 1
  have i3 : IntegrableOn (fun y => h ε y ^ 3) (Ioi 0) := hpow_int hε hε1 2
  have e : ∀ y, Real.Gamma ε ^ 3 - (Real.Gamma ε - h ε y) ^ 3 =
      3 * Real.Gamma ε ^ 2 * h ε y - 3 * Real.Gamma ε * h ε y ^ 2 + h ε y ^ 3 := fun y => by ring
  simp only [e] at hsw
  have j1 : IntegrableOn (fun y => 3 * Real.Gamma ε ^ 2 * h ε y) (Ioi 0) := i1.const_mul _
  have j2 : IntegrableOn (fun y => 3 * Real.Gamma ε * h ε y ^ 2) (Ioi 0) := i2.const_mul _
  have j12 : IntegrableOn (fun y => 3 * Real.Gamma ε ^ 2 * h ε y - 3 * Real.Gamma ε * h ε y ^ 2)
      (Ioi 0) := j1.sub j2
  rw [integral_add j12 i3,
    integral_sub j1 j2, integral_const_mul, integral_const_mul,
    h_integral (by linarith)] at hsw
  unfold H
  have hs : 0 < Real.sqrt π := Real.sqrt_pos.mpr Real.pi_pos
  field_simp
  linarith

end Z04

namespace Z04

lemma face_eq {ε : ℝ} (hε : 0 < ε) (hε1 : ε ≤ 1) :
    ∫ t in Y05.Δ, Y05.faceIntegrand (1 - ε) t =
      (3 * Real.Gamma ε ^ 2 - 3 * Real.Gamma ε * H 2 ε + H 3 ε) / Real.Gamma (1 + 2 * ε) := by
  have hd := Y05.dirichlet (τ := 1 - ε) (by linarith)
  have hw : ∫ g in Y05.Q, W3 ε g = Real.sqrt 3 * ∫ g in Y05.Q, Y05.Fb (1 - ε) g := by
    rw [← integral_const_mul]
    exact setIntegral_congr_fun Y05.Q_meas (fun g hg => W3_eq hg)
  rw [W3_integral hε hε1, hd, show 3 - 2 * (1 - ε) = 1 + 2 * ε by ring] at hw
  have hg : 0 < Real.Gamma (1 + 2 * ε) := Real.Gamma_pos_of_pos (by linarith)
  have h3 : 0 < Real.sqrt 3 := by positivity
  rw [hw]
  field_simp

end Z04

theorem Z04.proof_5 : Z04.stmt_5 := by
  intro p hp
  have hε : 0 < 1 / p := by positivity
  have hε1 : 1 / p ≤ 1 := by rw [div_le_one (by linarith)]; exact hp
  refine ⟨Z04.hpow_int hε hε1 1, Z04.hpow_int hε hε1 2, ?_⟩
  unfold Y05.S
  rw [Z04.face_eq hε hε1, add_comm (1 : ℝ) (1 / p), Real.Gamma_add_one hε.ne']
  have hg : 0 < Real.Gamma (1 + 2 * (1 / p)) := Real.Gamma_pos_of_pos (by linarith)
  field_simp

namespace Z04

/-- the positive quadrant -/
abbrev Q2 : Set (ℝ × ℝ) := Ioi (0 : ℝ) ×ˢ Ioi (0 : ℝ)

lemma Q2_meas : MeasurableSet Q2 := measurableSet_Ioi.prod measurableSet_Ioi

lemma prod2 {f1 f2 : ℝ → ℝ} (i1 : IntegrableOn f1 (Ioi 0)) (i2 : IntegrableOn f2 (Ioi 0)) :
    IntegrableOn (fun k : ℝ × ℝ => f1 k.1 * f2 k.2) Q2 ∧
      ∫ k in Q2, f1 k.1 * f2 k.2 = (∫ x in Ioi (0 : ℝ), f1 x) * ∫ x in Ioi (0 : ℝ), f2 x := by
  constructor
  · show IntegrableOn _ _ (volume.prod volume)
    rw [IntegrableOn, ← Measure.prod_restrict]
    exact Integrable.mul_prod i1 i2
  · exact setIntegral_prod_mul (μ := (volume : Measure ℝ)) (ν := (volume : Measure ℝ)) f1 f2
      (Ioi (0 : ℝ)) (Ioi (0 : ℝ))

/-- `w(k) |k^τ|` in two dimensions. -/
noncomputable def W2 (ε : ℝ) (k : ℝ × ℝ) : ℝ :=
  Real.exp (-(k.1 + k.2)) * (k.1 * k.2) ^ (-(1 - ε)) *
    Real.sqrt (k.1 ^ (2 * (1 - ε)) + k.2 ^ (2 * (1 - ε)))

/-- `w(k) (1 - e^{-|k^τ|²/y²})` in two dimensions. -/
noncomputable def Phi2 (ε : ℝ) (z : (ℝ × ℝ) × ℝ) : ℝ :=
  Real.exp (-(z.1.1 + z.1.2)) * (z.1.1 * z.1.2) ^ (-(1 - ε)) *
    (1 - Real.exp (-((z.1.1 ^ (2 * (1 - ε)) + z.1.2 ^ (2 * (1 - ε))) / z.2 ^ 2)))

lemma W2_meas (ε : ℝ) : Measurable (W2 ε) := by unfold W2; fun_prop

lemma W2_nn (ε : ℝ) {k : ℝ × ℝ} (hk : k ∈ Q2) : 0 ≤ W2 ε k := by
  obtain ⟨ha, hb⟩ := hk
  simp only [mem_Ioi] at ha hb
  unfold W2; positivity

lemma W2_bd {ε : ℝ} {k : ℝ × ℝ} (hk : k ∈ Q2) :
    W2 ε k ≤ (Real.exp (-k.1) * k.1 ^ (0 : ℝ)) * u0 (1 - ε) k.2 +
      u0 (1 - ε) k.1 * (Real.exp (-k.2) * k.2 ^ (0 : ℝ)) := by
  obtain ⟨ha, hb⟩ := hk
  simp only [mem_Ioi] at ha hb
  unfold W2 u0
  set τ := 1 - ε
  rw [Real.mul_rpow ha.le hb.le, Real.rpow_zero, Real.rpow_zero]
  have e : ∀ x : ℝ, 0 < x → x ^ (2 * τ) = (x ^ τ) ^ 2 := fun x hx => by
    rw [← Real.rpow_natCast, ← Real.rpow_mul hx.le]; congr 1; push_cast; ring
  have ek : ∀ x : ℝ, 0 < x → x ^ (-τ) * x ^ τ = 1 := fun x hx => by
    rw [← Real.rpow_add hx]; simp
  rw [e _ ha, e _ hb]
  have pa := Real.rpow_pos_of_pos ha τ; have pb := Real.rpow_pos_of_pos hb τ
  have na := Real.rpow_pos_of_pos ha (-τ); have nb := Real.rpow_pos_of_pos hb (-τ)
  have ka := ek _ ha; have kb := ek _ hb
  have hs : Real.sqrt ((k.1 ^ τ) ^ 2 + (k.2 ^ τ) ^ 2) ≤ k.1 ^ τ + k.2 ^ τ := by
    rw [Real.sqrt_le_iff]; constructor
    · positivity
    · nlinarith [mul_pos pa pb]
  rw [neg_add, Real.exp_add]
  have hE1 := Real.exp_pos (-k.1); have hE2 := Real.exp_pos (-k.2)
  generalize Real.exp (-k.1) = E1 at *
  generalize Real.exp (-k.2) = E2 at *
  generalize k.1 ^ τ = A at *
  generalize k.2 ^ τ = B at *
  generalize k.1 ^ (-τ) = A' at *
  generalize k.2 ^ (-τ) = B' at *
  have : E1 * E2 * (A' * B') * Real.sqrt (A ^ 2 + B ^ 2) ≤ E1 * E2 * (A' * B') * (A + B) :=
    mul_le_mul_of_nonneg_left hs (by positivity)
  calc E1 * E2 * (A' * B') * Real.sqrt (A ^ 2 + B ^ 2) ≤ E1 * E2 * (A' * B') * (A + B) := this
    _ = E1 * E2 * B' * (A' * A) + E1 * E2 * A' * (B' * B) := by ring
    _ = _ := by rw [ka, kb]; ring

lemma W2_int {ε : ℝ} (hε : 0 < ε) : IntegrableOn (W2 ε) Q2 := by
  obtain ⟨g0, _⟩ := Y05.gamma1 (x := 0) (by norm_num)
  obtain ⟨i0, _⟩ := u0_int hε
  have i := (prod2 g0 i0).1.add (prod2 i0 g0).1
  refine Integrable.mono' i (W2_meas ε).aestronglyMeasurable ?_
  refine (ae_restrict_iff' Q2_meas).mpr (Eventually.of_forall fun k hk => ?_)
  rw [Real.norm_eq_abs, abs_of_nonneg (W2_nn ε hk)]
  exact W2_bd hk

lemma pt2 {ε y : ℝ} {k : ℝ × ℝ} (hk : k ∈ Q2) :
    Real.exp (-(k.1 + k.2)) * (k.1 * k.2) ^ (-(1 - ε)) *
        (1 - Real.exp (-((k.1 ^ (2 * (1 - ε)) + k.2 ^ (2 * (1 - ε))) / y ^ 2))) =
      u0 (1 - ε) k.1 * u0 (1 - ε) k.2 - vv ε y k.1 * vv ε y k.2 := by
  obtain ⟨ha, hb⟩ := hk
  simp only [mem_Ioi] at ha hb
  unfold vv u0 hIntegrand
  rw [Real.mul_rpow ha.le hb.le,
    show (2 - 2 * ε) = 2 * (1 - ε) by ring, show ε - 1 = -(1 - ε) by ring,
    add_div, neg_add, Real.exp_add, neg_add, Real.exp_add]
  ring

lemma Phi2_y {ε : ℝ} {k : ℝ × ℝ} (hk : k ∈ Q2) :
    IntegrableOn (fun y => Phi2 ε (k, y)) (Ioi 0) ∧
      ∫ y in Ioi (0 : ℝ), Phi2 ε (k, y) = Real.sqrt π * W2 ε k := by
  obtain ⟨ha, hb⟩ := hk
  simp only [mem_Ioi] at ha hb
  obtain ⟨gi, gv⟩ := gy (c := k.1 ^ (2 * (1 - ε)) + k.2 ^ (2 * (1 - ε))) (by positivity)
  unfold Phi2 W2
  simp only
  refine ⟨gi.const_mul _, ?_⟩
  rw [integral_const_mul, gv]; ring

lemma Phi2_int {ε : ℝ} (hε : 0 < ε) :
    Integrable (Phi2 ε) ((volume.restrict Q2).prod (volume.restrict (Ioi (0 : ℝ)))) := by
  have hm : AEStronglyMeasurable (Phi2 ε)
      ((volume.restrict Q2).prod (volume.restrict (Ioi (0 : ℝ)))) := by
    unfold Phi2; fun_prop
  rw [integrable_prod_iff hm]
  refine ⟨(ae_restrict_iff' Q2_meas).mpr (Eventually.of_forall fun k hk => (Phi2_y hk).1), ?_⟩
  have := (W2_int hε).const_mul (Real.sqrt π)
  refine IntegrableOn.congr_fun this (fun k hk => ?_) Q2_meas
  have hk' := hk
  obtain ⟨ha, hb⟩ := hk
  simp only [mem_Ioi] at ha hb
  rw [setIntegral_congr_fun measurableSet_Ioi (fun y hy => by
      rw [Real.norm_eq_abs, abs_of_nonneg (by
        unfold Phi2
        have := one_sub_exp_neg_nn (x := (k.1 ^ (2 * (1 - ε)) + k.2 ^ (2 * (1 - ε))) / y ^ 2)
          (by positivity)
        positivity)]), (Phi2_y hk').2]

lemma Phi2_g {ε : ℝ} (hε : 0 < ε) (hε1 : ε ≤ 1) {y : ℝ} (hy : 0 < y) :
    ∫ k in Q2, Phi2 ε (k, y) = Real.Gamma ε ^ 2 - (Real.Gamma ε - h ε y) ^ 2 := by
  obtain ⟨i0, v0⟩ := u0_int hε
  obtain ⟨i1, v1⟩ := vv_int hε hε1 hy
  obtain ⟨j0, w0⟩ := prod2 i0 i0
  obtain ⟨j1, w1⟩ := prod2 i1 i1
  rw [setIntegral_congr_fun Q2_meas (fun k hk => by unfold Phi2; exact pt2 hk),
    integral_sub j0 j1, w0, w1, v0, v1]
  ring

/-- the Laplace side of `P`: `∫_{Q₂} w |k^τ| = 2Γ(ε) - H₂`. -/
lemma W2_integral {ε : ℝ} (hε : 0 < ε) (hε1 : ε ≤ 1) :
    ∫ k in Q2, W2 ε k = 2 * Real.Gamma ε - H 2 ε := by
  have hsw := integral_integral_swap (f := fun k y => Phi2 ε (k, y)) (Phi2_int hε)
  rw [setIntegral_congr_fun Q2_meas (fun k hk => (Phi2_y hk).2), integral_const_mul,
    setIntegral_congr_fun measurableSet_Ioi (fun y hy => Phi2_g hε hε1 hy)] at hsw
  have i1 := h_int (ε := ε) (by linarith)
  have i2 : IntegrableOn (fun y => h ε y ^ 2) (Ioi 0) := hpow_int hε hε1 1
  have e : ∀ y, Real.Gamma ε ^ 2 - (Real.Gamma ε - h ε y) ^ 2 =
      2 * Real.Gamma ε * h ε y - h ε y ^ 2 := fun y => by ring
  simp only [e] at hsw
  have j1 : IntegrableOn (fun y => 2 * Real.Gamma ε * h ε y) (Ioi 0) := i1.const_mul _
  rw [integral_sub j1 i2, integral_const_mul, h_integral (by linarith)] at hsw
  unfold H
  have hs : 0 < Real.sqrt π := Real.sqrt_pos.mpr Real.pi_pos
  field_simp
  linarith

end Z04

namespace Z04

lemma W2_scale {p : ℝ} (hp : 0 < p) {r s : ℝ} (hr : 0 < r) (hs : s ∈ Ioo (0 : ℝ) 1) :
    W2 (1 / p) (r * s, r - r * s) =
      Real.exp (-r) * r ^ (-(1 - 1 / p)) * W02.pIntegrand p s := by
  obtain ⟨hs0, hs1⟩ := hs
  have h1s : 0 < 1 - s := by linarith
  set τ := 1 - 1 / p with hτ
  unfold W2 W02.pIntegrand
  rw [show r - r * s = r * (1 - s) by ring, show r * s + r * (1 - s) = r by ring,
    show r * s * (r * (1 - s)) = r * s * r * (1 - s) by ring,
    Real.mul_rpow (by positivity) h1s.le, Real.mul_rpow (by positivity) hr.le,
    Real.mul_rpow hr.le hs0.le, Real.mul_rpow hr.le hs0.le, Real.mul_rpow hr.le h1s.le,
    show 2 / p - 2 = -2 * τ by rw [hτ]; ring, ← hτ]
  have e2 : ∀ x : ℝ, 0 < x → x ^ (2 * τ) = (x ^ τ) ^ 2 := fun x hx => by
    rw [← Real.rpow_natCast, ← Real.rpow_mul hx.le]; congr 1; push_cast; ring
  have em2 : ∀ x : ℝ, 0 < x → x ^ (-2 * τ) = ((x ^ τ) ^ 2)⁻¹ := fun x hx => by
    rw [show -2 * τ = -(2 * τ) by ring, Real.rpow_neg hx.le, e2 x hx]
  have em : ∀ x : ℝ, 0 < x → x ^ (-τ) = (x ^ τ)⁻¹ := fun x hx => Real.rpow_neg hx.le τ
  rw [e2 _ hr, e2 _ hs0, e2 _ h1s, em2 _ hs0, em2 _ h1s, em _ hr, em _ hs0, em _ h1s]
  have pa := Real.rpow_pos_of_pos hs0 τ
  have pb := Real.rpow_pos_of_pos h1s τ
  have pR := Real.rpow_pos_of_pos hr τ
  generalize s ^ τ = a at *
  generalize (1 - s) ^ τ = b at *
  generalize r ^ τ = R at *
  have hsq : Real.sqrt (((a ^ 2)⁻¹ + (b ^ 2)⁻¹)) = Real.sqrt (a ^ 2 + b ^ 2) / (a * b) := by
    rw [show (a ^ 2)⁻¹ + (b ^ 2)⁻¹ = (a ^ 2 + b ^ 2) / (a * b) ^ 2 by field_simp; ring,
      Real.sqrt_div' _ (by positivity), Real.sqrt_sq (by positivity)]
  rw [hsq, show R ^ 2 * a ^ 2 + R ^ 2 * b ^ 2 = R ^ 2 * (a ^ 2 + b ^ 2) by ring,
    Real.sqrt_mul (by positivity), Real.sqrt_sq pR.le]
  field_simp

lemma dir2 {p : ℝ} (hp : 1 ≤ p) :
    ∫ k in Q2, W2 (1 / p) k = Real.Gamma (1 + 1 / p) * ∫ s in (0 : ℝ)..1, W02.pIntegrand p s := by
  have hp0 : 0 < p := by linarith
  have hε : 0 < 1 / p := by positivity
  set F := Q2.indicator (W2 (1 / p)) with hF
  have hFm : Measurable F := (W2_meas _).indicator Q2_meas
  have hFi : Integrable F := (integrable_indicator_iff Q2_meas).mpr (W2_int hε)
  rw [← integral_indicator Q2_meas]
  have hFi' : Integrable F (volume.prod volume) := hFi
  have hS : MeasurePreserving (fun w : ℝ × ℝ => (w.1, w.2 - w.1))
      (volume.prod volume) (volume.prod volume) := by
    refine MeasurePreserving.skew_product (MeasurePreserving.id _)
      (g := fun (u : ℝ) (r : ℝ) => r - u) (by fun_prop) (Eventually.of_forall fun u => ?_)
    exact map_sub_right_eq_self (volume : Measure ℝ) u
  have e2 : ∫ v, F v ∂(volume.prod volume) =
      ∫ w, F (w.1, w.2 - w.1) ∂(volume.prod volume) := by
    conv_lhs => rw [← hS.map_eq]
    rw [integral_map hS.measurable.aemeasurable hFm.aestronglyMeasurable]
  have hHi : Integrable (fun w : ℝ × ℝ => F (w.1, w.2 - w.1)) (volume.prod volume) :=
    hS.integrable_comp_of_integrable hFi'
  show ∫ v, F v ∂(volume.prod volume) = _
  rw [e2, integral_prod_symm _ hHi]
  have e3 : ∀ r : ℝ, ∫ u : ℝ, F (u, r - u) =
      (Ioi (0 : ℝ)).indicator (fun r => Real.exp (-r) * r ^ ((1 + 1 / p) - 1) *
        ∫ s in (0 : ℝ)..1, W02.pIntegrand p s) r := by
    intro r
    by_cases hr : 0 < r
    · rw [indicator_of_mem (show r ∈ Ioi 0 from hr)]
      have hs := Measure.integral_comp_mul_left (fun u : ℝ => F (u, r - u)) r
      have hpt : (fun x : ℝ => F (r * x, r - r * x)) = (Ioo (0 : ℝ) 1).indicator
          (fun s => Real.exp (-r) * r ^ (-(1 - 1 / p)) * W02.pIntegrand p s) := by
        funext x
        by_cases hx : x ∈ Ioo (0 : ℝ) 1
        · rw [indicator_of_mem hx, hF, indicator_of_mem, W2_scale hp0 hr hx]
          obtain ⟨h0, h1⟩ := hx
          refine ⟨?_, ?_⟩ <;> simp only [mem_Ioi]
          · positivity
          · nlinarith
        · rw [indicator_of_notMem hx, hF, indicator_of_notMem]
          rintro ⟨ha, hb⟩
          simp only [mem_Ioi] at ha hb
          apply hx
          refine ⟨pos_of_mul_pos_right ha hr.le, ?_⟩
          by_contra hc
          push Not at hc
          nlinarith
      rw [hpt, integral_indicator measurableSet_Ioo, integral_const_mul,
        abs_of_pos (inv_pos.mpr hr), smul_eq_mul, ← integral_Ioc_eq_integral_Ioo,
        ← intervalIntegral.integral_of_le zero_le_one] at hs
      rw [show (1 + 1 / p) - 1 = 1 + (-(1 - 1 / p)) by ring, Real.rpow_add hr, Real.rpow_one]
      have h' : ∫ u : ℝ, F (u, r - u) = r * (Real.exp (-r) * r ^ (-(1 - 1 / p)) *
          ∫ s in (0 : ℝ)..1, W02.pIntegrand p s) := by
        rw [hs]; field_simp
      rw [h']; ring
    · rw [indicator_of_notMem (show r ∉ Ioi 0 from hr)]
      have : ∀ u : ℝ, F (u, r - u) = 0 := fun u => by
        rw [hF, indicator_of_notMem]
        rintro ⟨ha, hb⟩
        simp only [mem_Ioi] at ha hb
        exact hr (by linarith)
      simp only [this, integral_zero]
  simp_rw [e3]
  rw [integral_indicator measurableSet_Ioi, integral_mul_const,
    ← Real.Gamma_eq_integral (by linarith)]

end Z04

theorem Z04.proof_4 : Z04.stmt_4 := by
  intro p hp
  have hp0 : 0 < p := by linarith
  have hε : 0 < 1 / p := by positivity
  have hε1 : 1 / p ≤ 1 := by rw [div_le_one hp0]; exact hp
  refine ⟨Z04.hpow_int hε hε1 1, ?_⟩
  have h1 := Z04.dir2 hp
  rw [Z04.W2_integral hε hε1] at h1
  have hg : 0 < Real.Gamma (1 / p) := Real.Gamma_pos_of_pos hε
  have hG : Real.Gamma (1 + 1 / p) = 1 / p * Real.Gamma (1 / p) := by
    rw [add_comm, Real.Gamma_add_one hε.ne']
  unfold W02.Pint
  rw [hG] at h1 ⊢
  have hK : ∫ s in (0 : ℝ)..1, W02.pIntegrand p s =
      (2 * Real.Gamma (1 / p) - Z04.H 2 (1 / p)) / (1 / p * Real.Gamma (1 / p)) := by
    rw [h1]; field_simp
  rw [hK]
  field_simp
  ring

theorem Z04.proof_6 : Z04.stmt_6 := by
  intro p hp
  have hp0 : 0 < p := by linarith
  have hε : 0 < 1 / p := by positivity
  have hε1 : 1 / p ≤ 1 := by rw [div_le_one hp0]; exact hp
  refine ⟨Y05.proof_1 p hp0, ?_, fun y hy => Z04.hI_int (by linarith) hy,
    Z04.hpow_int hε hε1 2, ?_⟩
  · rcases hp.lt_or_eq with hp1 | hp1
    · exact (W02.proof_2 p hp1).1
    · subst hp1
      have : W02.pIntegrand 1 = fun _ => Real.sqrt 2 := by
        funext t; unfold W02.pIntegrand; norm_num
      rw [this]; exact intervalIntegrable_const
  · obtain ⟨_, hP⟩ := Z04.proof_4 p hp
    obtain ⟨_, _, hS⟩ := Z04.proof_5 p hp
    rw [hS, hP]
    unfold Z04.G Z04.Q Z04.H
    have g1 : 0 < Real.Gamma (1 + 1 / p) := Real.Gamma_pos_of_pos (by linarith)
    have g2 : 0 < Real.Gamma (1 + 2 * (1 / p)) := Real.Gamma_pos_of_pos (by linarith)
    have hs : 0 < Real.sqrt π := Real.sqrt_pos.mpr Real.pi_pos
    field_simp
    ring

namespace Z04

lemma G_le_one {ε : ℝ} (hε : 0 ≤ ε) : G ε ≤ 1 := by
  have hc := Real.convexOn_log_Gamma.2 (show (1 : ℝ) ∈ Ioi 0 by simp)
    (show 1 + 2 * ε ∈ Ioi (0 : ℝ) by simp only [mem_Ioi]; linarith)
    (show (0 : ℝ) ≤ 1 / 2 by norm_num) (show (0 : ℝ) ≤ 1 / 2 by norm_num)
    (by norm_num : (1 / 2 : ℝ) + 1 / 2 = 1)
  simp only [Function.comp, smul_eq_mul, Real.Gamma_one, Real.log_one, mul_zero, zero_add] at hc
  rw [show 1 / 2 * 1 + 1 / 2 * (1 + 2 * ε) = 1 + ε by ring] at hc
  have g1 : 0 < Real.Gamma (1 + ε) := Real.Gamma_pos_of_pos (by linarith)
  have g2 : 0 < Real.Gamma (1 + 2 * ε) := Real.Gamma_pos_of_pos (by linarith)
  have : Real.log (Real.Gamma (1 + ε) ^ 2) ≤ Real.log (Real.Gamma (1 + 2 * ε)) := by
    rw [Real.log_pow]; push_cast; linarith
  have := (Real.log_le_log_iff (by positivity) g2).mp this
  unfold G
  rw [div_le_one g2]; exact this

lemma G_ge {ε : ℝ} (hε : 0 < ε) : 1 - 4 * ε ^ 2 ≤ G ε := by
  have hb := X05.beta_real (a := 1 + ε) (b := 1 + ε) (by linarith) (by linarith)
  rw [show 1 + ε - 1 = ε by ring, show 1 + ε + (1 + ε) = (1 + 2 * ε) + 1 by ring,
    Real.Gamma_add_one (by linarith)] at hb
  have g2 : 0 < Real.Gamma (1 + 2 * ε) := Real.Gamma_pos_of_pos (by linarith)
  set B := ∫ t in (0 : ℝ)..1, t ^ ε * (1 - t) ^ ε
  have hG : G ε = (1 + 2 * ε) * B := by
    unfold G; rw [sq, ← hb]; field_simp
  have hlog : IntervalIntegrable Real.log volume 0 1 := intervalIntegral.intervalIntegrable_log'
  have hlog' : IntervalIntegrable (fun t : ℝ => Real.log (1 - t)) volume 0 1 := by
    have := hlog.comp_sub_left 1
    simpa using this.symm
  have hl1 : ∫ t in (0 : ℝ)..1, Real.log t = -1 := by
    rw [integral_log_from_zero]; simp
  have hl2 : ∫ t in (0 : ℝ)..1, Real.log (1 - t) = -1 := by
    rw [intervalIntegral.integral_comp_sub_left (fun t => Real.log t) 1]; simpa using hl1
  have hlow : IntervalIntegrable (fun t : ℝ => 1 + ε * (Real.log t + Real.log (1 - t)))
      volume 0 1 := intervalIntegrable_const.add ((hlog.add hlog').const_mul ε)
  have hc : Continuous (fun t : ℝ => t ^ ε * (1 - t) ^ ε) :=
    (Real.continuous_rpow_const hε.le).mul
      ((Real.continuous_rpow_const hε.le).comp (continuous_const.sub continuous_id))
  have hB : 1 - 2 * ε ≤ B := by
    have hm := intervalIntegral.integral_mono_on_of_le_Ioo zero_le_one hlow
      (hc.intervalIntegrable 0 1) (fun t ht => by
        obtain ⟨h0, h1⟩ := ht
        have h1t : 0 < 1 - t := by linarith
        rw [← Real.mul_rpow h0.le h1t.le, Real.rpow_def_of_pos (by positivity),
          Real.log_mul h0.ne' h1t.ne']
        have := Real.add_one_le_exp (Real.log t + Real.log (1 - t)) |>.trans_eq' rfl
        have := Real.add_one_le_exp ((Real.log t + Real.log (1 - t)) * ε)
        linarith)
    rw [intervalIntegral.integral_add intervalIntegrable_const ((hlog.add hlog').const_mul ε),
      intervalIntegral.integral_const_mul, intervalIntegral.integral_add hlog hlog', hl1, hl2]
      at hm
    simp at hm
    linarith
  rw [hG]; nlinarith

end Z04

namespace Z04

lemma FR_lip : ∃ C δ : ℝ, 0 < δ ∧ ∀ ε : ℝ, |ε| < δ → |W02.FR ε - W02.FR 0| ≤ C * |ε| := by
  have hopen : IsOpen {z : ℂ | z.re < 1} := isOpen_lt Complex.continuous_re continuous_const
  have hd : DifferentiableAt ℂ W02.FC 0 :=
    W02.proof_5.2.differentiableAt (hopen.mem_nhds (by simp))
  obtain ⟨c, hc⟩ := hd.isBigO_sub.bound
  obtain ⟨δ, hδ, hδP⟩ := Metric.eventually_nhds_iff.mp hc
  refine ⟨c, min δ 1, lt_min hδ one_pos, fun ε hε => ?_⟩
  have hε1 : |ε| < 1 := lt_of_lt_of_le hε (min_le_right _ _)
  have hεδ : |ε| < δ := lt_of_lt_of_le hε (min_le_left _ _)
  have h1 := hδP (y := (ε : ℂ)) (by simpa [dist_eq_norm] using hεδ)
  have e1 := (W02.proof_4 ε (lt_of_le_of_lt (le_abs_self ε) hε1)).2
  have e0 := (W02.proof_4 0 (by norm_num)).2
  simp only [Complex.ofReal_zero] at e0
  rw [e1, e0, sub_zero, ← Complex.ofReal_sub, Complex.norm_real, Complex.norm_real,
    Real.norm_eq_abs, Real.norm_eq_abs] at h1
  exact h1

lemma M0_eq : (W02.Mk 0).re = -(W02.FR 0) / 2 := by
  have e0 := (W02.proof_4 0 (by norm_num)).2
  simp only [Complex.ofReal_zero] at e0
  unfold W02.Mk
  rw [iteratedDeriv_zero, e0]
  simp

lemma P_exp : ∃ C δ : ℝ, 0 < δ ∧ ∀ p : ℝ, 1 < p → 1 / p < δ →
    |W02.Pint p - 8 - 8 * (W02.Mk 0).re / p| ≤ C / p ^ 2 := by
  obtain ⟨C, δ, hδ, hL⟩ := FR_lip
  refine ⟨4 * C, δ, hδ, fun p hp hpδ => ?_⟩
  have hp0 : 0 < p := by linarith
  have hε : 0 < 1 / p := by positivity
  have hP := (W02.proof_2 p hp).2.2
  have hL' := hL (1 / p) (by rw [abs_of_pos hε]; exact hpδ)
  rw [abs_of_pos hε] at hL'
  have e : W02.Pint p - 8 - 8 * (W02.Mk 0).re / p =
      -(4 * (1 / p)) * (W02.FR (1 / p) - W02.FR 0) := by
    rw [M0_eq]
    have : W02.Pint p = 8 - 4 * (1 / p) * W02.FR (1 / p) := by linarith
    rw [this]; field_simp; ring
  rw [e, abs_mul, abs_neg, abs_of_pos (by positivity : 0 < 4 * (1 / p))]
  calc 4 * (1 / p) * |W02.FR (1 / p) - W02.FR 0| ≤ 4 * (1 / p) * (C * (1 / p)) :=
        mul_le_mul_of_nonneg_left hL' (by positivity)
    _ = 4 * C / p ^ 2 := by field_simp

lemma gamma_lower : ∃ m : ℝ, 0 < m ∧ ∀ x ∈ Icc (1 : ℝ) (3 / 2), m ≤ Real.Gamma x := by
  have hc : ContinuousOn Real.Gamma (Icc (1 : ℝ) (3 / 2)) := fun x hx =>
    (Real.differentiableAt_Gamma (fun m => by
      have := hx.1; have : (0 : ℝ) ≤ m := m.cast_nonneg; intro h; linarith)).continuousAt
      |>.continuousWithinAt
  obtain ⟨x0, hx0, hmin⟩ := isCompact_Icc.exists_isMinOn ⟨1, by norm_num⟩ hc
  refine ⟨Real.Gamma x0, Real.Gamma_pos_of_pos (by linarith [hx0.1]), fun x hx => hmin hx⟩

lemma h_le_small {ε : ℝ} (h0 : 0 ≤ ε) (h1 : ε ≤ 1 / 4) {y : ℝ} (hy : 0 < y) (hy1 : y ≤ 1) :
    h ε y ≤ 20 * y ^ (-(1 / 4 : ℝ)) := by
  have hb := h_le_log h0 (by linarith) hy
  have hlog : |Real.log y| ≤ 4 * y ^ (-(1 / 4 : ℝ)) := by
    have hn : Real.log y ≤ 0 := Real.log_nonpos hy.le hy1
    rw [abs_of_nonpos hn, ← Real.log_inv]
    have := Real.log_le_rpow_div (x := y⁻¹) (by positivity) (ε := 1 / 4) (by norm_num)
    rw [Real.inv_rpow hy.le, ← Real.rpow_neg hy.le] at this
    linarith
  have hone : 1 ≤ y ^ (-(1 / 4 : ℝ)) :=
    Real.one_le_rpow_of_pos_of_le_one_of_nonpos hy hy1 (by norm_num)
  linarith

lemma int_h3_le : ∃ B : ℝ, ∀ ε : ℝ, 0 < ε → ε ≤ 1 / 4 →
    ∫ y in Ioi (0 : ℝ), h ε y ^ 3 ≤ B := by
  have i1 : IntegrableOn (fun y : ℝ => 8000 * y ^ (-(3 / 4 : ℝ))) (Ioc 0 1) := by
    have := intervalIntegral.intervalIntegrable_rpow' (r := -(3 / 4 : ℝ)) (by norm_num)
      (a := 0) (b := 1)
    rw [intervalIntegrable_iff_integrableOn_Ioc_of_le zero_le_one] at this
    exact this.const_mul _
  have i2 : IntegrableOn (fun y : ℝ => 8 * y ^ (-(6 : ℝ))) (Ioi 1) :=
    (integrableOn_Ioi_rpow_of_lt (by norm_num) zero_lt_one).const_mul _
  refine ⟨(∫ y in Ioc (0 : ℝ) 1, 8000 * y ^ (-(3 / 4 : ℝ))) +
    ∫ y in Ioi (1 : ℝ), 8 * y ^ (-(6 : ℝ)), fun ε h0 h1 => ?_⟩
  have hint : IntegrableOn (fun y => h ε y ^ 3) (Ioi 0) := hpow_int h0 (by linarith) 2
  have hsplit : Ioi (0 : ℝ) = Ioc 0 1 ∪ Ioi 1 := (Ioc_union_Ioi_eq_Ioi zero_le_one).symm
  rw [hsplit] at hint ⊢
  rw [setIntegral_union (Ioc_disjoint_Ioi le_rfl) measurableSet_Ioi
      (hint.mono_set subset_union_left) (hint.mono_set subset_union_right)]
  gcongr ?_ + ?_
  · refine setIntegral_mono_on (hint.mono_set subset_union_left) i1 measurableSet_Ioc
      (fun y hy => ?_)
    have hb := h_le_small h0.le h1 hy.1 hy.2
    have hn := h_nn ε y
    have e : (20 * y ^ (-(1 / 4 : ℝ))) ^ 3 = 8000 * y ^ (-(3 / 4 : ℝ)) := by
      rw [mul_pow, ← Real.rpow_natCast (y ^ (-(1 / 4 : ℝ))) 3, ← Real.rpow_mul hy.1.le]; norm_num
    rw [← e]; exact pow_le_pow_left₀ hn hb 3
  · refine setIntegral_mono_on (hint.mono_set subset_union_right) i2 measurableSet_Ioi
      (fun y hy => ?_)
    have hy0 : (0 : ℝ) < y := lt_trans zero_lt_one hy
    have hb := h_le_inv_sq h0.le (by linarith) hy0
    have hn := h_nn ε y
    have e : (2 / y ^ 2) ^ 3 = 8 * y ^ (-(6 : ℝ)) := by
      rw [Real.rpow_neg hy0.le, show (6 : ℝ) = ((6 : ℕ) : ℝ) by norm_num, Real.rpow_natCast]
      field_simp; ring
    rw [← e]; exact pow_le_pow_left₀ hn hb 3

lemma Q_bd : ∃ C : ℝ, ∀ ε : ℝ, 0 < ε → ε ≤ 1 / 4 → |Q ε| ≤ C := by
  obtain ⟨B, hB⟩ := int_h3_le
  obtain ⟨m, hm, hmin⟩ := gamma_lower
  have hs : 0 < Real.sqrt π := Real.sqrt_pos.mpr Real.pi_pos
  refine ⟨B / (Real.sqrt π * m), fun ε h0 h1 => ?_⟩
  have hg := hmin (1 + 2 * ε) ⟨by linarith, by linarith⟩
  have hg0 : 0 < Real.Gamma (1 + 2 * ε) := lt_of_lt_of_le hm hg
  have hI0 : 0 ≤ ∫ y in Ioi (0 : ℝ), h ε y ^ 3 :=
    setIntegral_nonneg measurableSet_Ioi (fun y _ => pow_nonneg (h_nn ε y) 3)
  have hIB := hB ε h0 h1
  unfold Q
  rw [abs_mul, abs_neg, abs_of_pos (by positivity), abs_of_nonneg hI0]
  rw [one_div_mul_eq_div, div_le_div_iff₀ (by positivity) (by positivity)]
  have : Real.sqrt π * m ≤ Real.sqrt π * Real.Gamma (1 + 2 * ε) :=
    mul_le_mul_of_nonneg_left hg hs.le
  have hB0 : 0 ≤ B := le_trans hI0 hIB
  have := mul_le_mul hIB this (by positivity) hB0
  linarith

end Z04

namespace Z04

/-- the algebra of the `K = 1` estimate. -/
lemma alg8 {e g P q M C1 C2 : ℝ} (he0 : 0 < e) (he4 : e ≤ 1 / 4) (hg1 : g ≤ 1) (hg2 : 1 - 4 * e ^ 2 ≤ g)
    (hP : |P - 8 - 8 * M * e| ≤ C1 * e ^ 2) (hq : |q| ≤ C2) :
    |6 * g * (P - 4) - 8 * e ^ 2 * q - 24 - 48 * M * e| ≤
      (96 + 192 * |M| + 6 * C1 + 8 * C2) * e ^ 2 := by
  have he1 : e ≤ 1 := by linarith
  set D := P - 8 - 8 * M * e with hD
  have hg0 : 0 ≤ g := by nlinarith
  have hgm : |g - 1| ≤ 4 * e ^ 2 := by rw [abs_le]; constructor <;> nlinarith
  have e1 : 6 * g * (P - 4) - 8 * e ^ 2 * q - 24 - 48 * M * e =
      24 * (g - 1) + 48 * M * e * (g - 1) + 6 * g * D - 8 * e ^ 2 * q := by rw [hD]; ring
  rw [e1]
  have t1 : |24 * (g - 1)| ≤ 96 * e ^ 2 := by
    rw [abs_mul, abs_of_pos (by norm_num : (0:ℝ) < 24)]; linarith
  have t2 : |48 * M * e * (g - 1)| ≤ 192 * |M| * e ^ 2 := by
    rw [abs_mul, abs_mul, abs_mul, abs_of_pos (by norm_num : (0:ℝ) < 48), abs_of_pos he0]
    have hM := abs_nonneg M
    have : e * |g - 1| ≤ 4 * e ^ 2 := by
      calc e * |g - 1| ≤ 1 * (4 * e ^ 2) :=
            mul_le_mul he1 hgm (abs_nonneg _) zero_le_one
        _ = 4 * e ^ 2 := one_mul _
    nlinarith
  have t3 : |6 * g * D| ≤ 6 * C1 * e ^ 2 := by
    rw [abs_mul, abs_mul, abs_of_pos (by norm_num : (0:ℝ) < 6), abs_of_nonneg hg0]
    have := abs_nonneg D
    have : g * |D| ≤ 1 * (C1 * e ^ 2) := mul_le_mul hg1 hP (abs_nonneg _) zero_le_one
    nlinarith
  have t4 : |8 * e ^ 2 * q| ≤ 8 * C2 * e ^ 2 := by
    rw [abs_mul, abs_of_pos (by positivity : (0:ℝ) < 8 * e ^ 2)]
    nlinarith [sq_nonneg e]
  calc |24 * (g - 1) + 48 * M * e * (g - 1) + 6 * g * D - 8 * e ^ 2 * q|
      ≤ |24 * (g - 1)| + |48 * M * e * (g - 1)| + |6 * g * D| + |8 * e ^ 2 * q| := by
        have a1 := abs_sub (24 * (g - 1) + 48 * M * e * (g - 1) + 6 * g * D) (8 * e ^ 2 * q)
        have a2 := abs_add_le (24 * (g - 1) + 48 * M * e * (g - 1)) (6 * g * D)
        have a3 := abs_add_le (24 * (g - 1)) (48 * M * e * (g - 1))
        linarith
    _ ≤ _ := by linarith

end Z04

theorem Z04.proof_8 : Z04.stmt_8 := by
  obtain ⟨C1, δ, hδ, hPe⟩ := Z04.P_exp
  obtain ⟨C2, hQ⟩ := Z04.Q_bd
  set M := (W02.Mk 0).re
  refine ⟨96 + 192 * |M| + 6 * C1 + 8 * C2, max 4 (2 / δ), fun p hp => ?_⟩
  have hp4 : 4 ≤ p := le_trans (le_max_left _ _) hp
  have hpd : 2 / δ ≤ p := le_trans (le_max_right _ _) hp
  have hp0 : 0 < p := by linarith
  have he0 : 0 < 1 / p := by positivity
  have he4 : 1 / p ≤ 1 / 4 := one_div_le_one_div_of_le (by norm_num) hp4
  have heδ : 1 / p < δ := by
    rw [div_lt_iff₀ hp0]
    have := (div_le_iff₀ hδ).mp hpd
    nlinarith
  obtain ⟨_, _, _, _, hS⟩ := Z04.proof_6 p (by linarith)
  have hP := hPe p (by linarith) heδ
  have hP' : |W02.Pint p - 8 - 8 * M * (1 / p)| ≤ C1 * (1 / p) ^ 2 := by
    rw [show 8 * M * (1 / p) = 8 * M / p by ring, show C1 * (1 / p) ^ 2 = C1 / p ^ 2 by ring]
    exact hP
  have key := Z04.alg8 he0 he4 (Z04.G_le_one he0.le) (Z04.G_ge he0) hP'
    (hQ (1 / p) he0 he4)
  rw [← hS] at key
  rw [show 48 * M / p = 48 * M * (1 / p) by ring,
    show (96 + 192 * |M| + 6 * C1 + 8 * C2) / p ^ 2 =
      (96 + 192 * |M| + 6 * C1 + 8 * C2) * (1 / p) ^ 2 by ring]
  exact key

theorem Z04.proof_7 : Z04.stmt_7 := by
  obtain ⟨C, N, hN⟩ := Z04.proof_8
  set M := (W02.Mk 0).re
  show Tendsto Y05.S atTop (𝓝 24)
  have hlim' : Tendsto (fun p : ℝ => |48 * M / p| + |C| / p ^ 2) atTop (𝓝 0) := by
    have h1 : Tendsto (fun p : ℝ => |48 * M / p|) atTop (𝓝 0) := by
      simpa using ((tendsto_const_nhds (x := 48 * M)).div_atTop tendsto_id).abs
    have h2 : Tendsto (fun p : ℝ => |C| / p ^ 2) atTop (𝓝 0) :=
      tendsto_const_nhds.div_atTop (tendsto_pow_atTop (by norm_num))
    simpa using h1.add h2
  rw [tendsto_iff_norm_sub_tendsto_zero]
  refine squeeze_zero' (Eventually.of_forall fun p => norm_nonneg _) ?_ hlim'
  filter_upwards [eventually_ge_atTop N, eventually_gt_atTop 0] with p hp hp0
  have := hN p hp
  rw [Real.norm_eq_abs]
  have hC : C / p ^ 2 ≤ |C| / p ^ 2 := div_le_div_of_nonneg_right (le_abs_self C) (by positivity)
  calc |Y05.S p - 24| = |(Y05.S p - 24 - 48 * M / p) + 48 * M / p| := by ring_nf
    _ ≤ |Y05.S p - 24 - 48 * M / p| + |48 * M / p| := abs_add_le _ _
    _ ≤ _ := by linarith
