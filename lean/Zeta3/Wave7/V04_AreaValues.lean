import Mathlib
import Zeta3.Wave3.W02_Holo
import Zeta3.Wave5.Y05_Octahedron
/-!
# V04 — exact values of the area and perimeter functions (sanity anchors)

Target id: V04.
Paper labels: the text around `thm:S` (paper/ieee/main.tex line 1909, theorem lines 1921–1945),
the landmarks of the perimeter (main.tex lines 802–806) and the checks after `eq:Pint`
(main.tex line 829), `eq:pw:Sn` (paper/ieee/proofs/sec_pureweight.tex lines 163–166), the
normalisation checks of the proof of `thm:S` (paper/ieee/proofs/sec_g3_S.tex line 87) and the
octahedron value in paper/ieee/proofs/sec_g7_octa.tex line 6.

Quoted verbatim.

main.tex 1909:

  The three-dimensional squircle is the superellipsoid $|x|^p+|y|^p+|z|^p=1$
  (Fig.~\ref{fig:superellipsoid}). Its surface area $S(p)$ goes from $4\sqrt3$ (the octahedron,
  $p=1$) through $4\pi$ (the sphere) to $24$ (the cube).

main.tex 802–806:

  The perimeter $P(p)$ has no such closed form. We only know the landmarks
  \[
  P(1)=4\sqrt2,\qquad P(2)=2\pi,\qquad \lim_{p\to\infty}P(p)=8 ,
  \]

main.tex 829 (after `eq:Pint`):

  \textbf{Case 2:} $p=2$. Here $t^{-1}+(1-t)^{-1}=1/(t(1-t))$, so
  $P(2)=2\int_0^1dt/\sqrt{t(1-t)}=2B(\frac12,\frac12)=2\pi$. \textbf{Case 3:} $p\to\infty$ gives
  $P\to8$; exact hypergeometric series are in~\cite{MRL,RM}.

sec_pureweight.tex 163–166 (`eq:pw:Sn`):

  S_n(p)=2^n\varepsilon^{n-1}\int_\Delta\Bigl(\prod_{i=1}^nt_i\Bigr)^{-\tau}
  \Bigl(\sum_{i=1}^nt_i^{2\tau}\Bigr)^{1/2}dt_1\cdots dt_{n-1},
  $\Delta=\{t_i\ge0,\ \sum_{i<n}t_i\le1\}$. The integral converges since $\tau<1$.

sec_g3_S.tex 87:

  As checks of the normalisation, at $p=1$ and $p=2$ the formula gives $4\sqrt3$ and $4\pi$, and
  $P(1)=4\sqrt2$, $P(2)=2\pi$ (Checks).

sec_g7_octa.tex 6:

  ... $S_3(1)=2^3\sqrt3/2!=4\sqrt3$, the area of the octahedron.

Imported objects (nothing is redefined).
* `S = Y05.S`: `Y05.S p = 8 (1/p)² ∫_Δ Y05.faceIntegrand (1 - 1/p)`, the paper's `eq:pw:Sn`
  for `n = 3` (`ε = 1/p`, `τ = 1 - 1/p`), with `Y05.faceIntegrand τ t =
  (t₁t₂t₃)^{-τ} (t₁^{2τ}+t₂^{2τ}+t₃^{2τ})^{1/2}`, `t₃ = 1 - t₁ - t₂`, and
  `Y05.Δ = {t₁, t₂ ≥ 0, t₁ + t₂ ≤ 1}` (Lebesgue measure on `ℝ × ℝ`).
* `P = W02.Pint`: `W02.Pint p = (4/p) ∫_0^1 W02.pIntegrand p`,
  `W02.pIntegrand p t = √(t^{2/p-2} + (1-t)^{2/p-2})` (`eq:Pint`).

Conventions. Every integral whose value is asserted comes with an explicit integrability claim,
so no statement holds through the "non-integrable ↦ 0" convention. At `p = 1` the exponents are
`τ = 0`, and `x ^ (0 : ℝ) = 1` for every real `x` (also `x = 0`), which is the paper's reading
(the integrand is `√3`, the octahedron's faces are flat). At `p = 2` the integrands are unbounded
near the boundary of `Δ` resp. near `t = 0, 1`; the boundary values (where `0 ^ (-1/2) = 0` and
`1/√0 = 0` in Lean) form a null set and do not matter. No `ζ` value appears in this target.

Statements.
* stmt_1: at `p = 1` (`τ = 0`) the face integrand is integrable on `Δ` and its integral is
  `√3 · area(Δ) = √3/2`.
* stmt_2: `S(1) = 4√3` (the octahedron), with the integrability of the face integrand.
* stmt_3: the Dirichlet integral used at `p = 2` (`τ = 1/2`):
  `∫_Δ (t₁t₂t₃)^{-1/2} = Γ(1/2)³/Γ(3/2) = 2π`, with integrability.
* stmt_4: `S(2) = 4π` (the sphere), with the integrability of the face integrand.
* stmt_5: `∫_0^1 dt/√(t(1-t)) = B(1/2,1/2) = π`, with integrability (main.tex 829).
* stmt_6: `P(2) = 2π` (the circle), with the integrability of the perimeter integrand.
* stmt_7: `lim_{p→∞} P(p) = 8` (main.tex 804 and 829, "Case 3").

NOT COVERED:
* `P(1) = 4√2`: it belongs to target Z01 (as the task assigns), not restated here.
* `S(p) → 24` (the cube) is Z04.stmt_7, not restated here.
* That `Y05.S` / `W02.Pint` are the true surface area / length (Hausdorff measure, resp.
  `eVariationOn`; the latter is W04) is not part of these value checks.
* The closed form `A(p) = 4Γ(1+1/p)²/Γ(1+2/p)` of the enclosed area and its values
  `A(1) = 2`, `A(2) = π`, `A(p) → 4` (main.tex 799–801): the enclosed area is not an object of
  the imported library.
* The general octahedron value `S_n(1) = 2^n √n/(n-1)!` for `n ≠ 3`.
-/
-- BEGIN STATEMENTS
namespace V04

open MeasureTheory Real Filter Topology

/-- Proof of `S(1) = 4√3`: at `p = 1`, `τ = 1 - 1/1 = 0`, the face integrand of `eq:pw:Sn` is the
constant `√3` on `Δ`; it is integrable and `∫_Δ = √3 · (1/2)`. -/
def stmt_1 : Prop :=
  IntegrableOn (Y05.faceIntegrand (1 - 1 / 1)) Y05.Δ ∧
    ∫ t in Y05.Δ, Y05.faceIntegrand (1 - 1 / 1) t = Real.sqrt 3 / 2

/-- main.tex 1909 "from $4\sqrt3$ (the octahedron, $p=1$)", sec_g7_octa.tex 6
"$S_3(1)=2^3\sqrt3/2!=4\sqrt3$": `S(1) = 4√3`, with the integrability of the face integrand. -/
def stmt_2 : Prop :=
  IntegrableOn (Y05.faceIntegrand (1 - 1 / 1)) Y05.Δ ∧ Y05.S 1 = 4 * Real.sqrt 3

/-- Proof of `S(2) = 4π` (`τ = 1/2`, where `Σ tᵢ^{2τ} = t₁+t₂+t₃ = 1`): the Dirichlet integral
`∫_Δ (t₁t₂t₃)^{-1/2} dt = Γ(1/2)³/Γ(3/2) = 2π`, with integrability. -/
def stmt_3 : Prop :=
  IntegrableOn (fun t : ℝ × ℝ => (t.1 * t.2 * (1 - t.1 - t.2)) ^ (-(1 / 2 : ℝ))) Y05.Δ ∧
    ∫ t in Y05.Δ, (t.1 * t.2 * (1 - t.1 - t.2)) ^ (-(1 / 2 : ℝ)) = 2 * π

/-- main.tex 1909 "through $4\pi$ (the sphere)", sec_g3_S.tex 87 "at $p=1$ and $p=2$ the formula
gives $4\sqrt3$ and $4\pi$": `S(2) = 4π`, with the integrability of the face integrand. -/
def stmt_4 : Prop :=
  IntegrableOn (Y05.faceIntegrand (1 - 1 / 2)) Y05.Δ ∧ Y05.S 2 = 4 * π

/-- main.tex 829: "$\int_0^1dt/\sqrt{t(1-t)}=B(\frac12,\frac12)$" `= π`, with integrability. -/
def stmt_5 : Prop :=
  IntervalIntegrable (fun t : ℝ => 1 / Real.sqrt (t * (1 - t))) volume 0 1 ∧
    ∫ t in (0 : ℝ)..1, 1 / Real.sqrt (t * (1 - t)) = π

/-- main.tex 804 and 829 "Case 2": `P(2) = 2π` (the circle), with the integrability of the
integrand of `eq:Pint` at `p = 2`. -/
def stmt_6 : Prop :=
  IntervalIntegrable (W02.pIntegrand 2) volume 0 1 ∧ W02.Pint 2 = 2 * π

/-- main.tex 804 "$\lim_{p\to\infty}P(p)=8$" and 829 "Case 3: $p\to\infty$ gives $P\to8$". -/
def stmt_7 : Prop :=
  Tendsto W02.Pint atTop (𝓝 8)

end V04
-- END STATEMENTS

open MeasureTheory Real Set Filter Topology

namespace V04

lemma arc_deriv {c b : ℝ} (hc : 0 < c) (hb0 : 0 < b) (hb1 : b < c) :
    HasDerivAt (fun x => Real.arcsin (2 * x / c - 1)) (1 / Real.sqrt (b * (c - b))) b := by
  have h1 : 2 * b / c - 1 ≠ -1 := by
    intro h; have : 2 * b / c = 0 := by linarith
    rw [div_eq_zero_iff] at this; rcases this with h | h <;> linarith
  have h2 : 2 * b / c - 1 ≠ 1 := by
    intro h; have : 2 * b / c = 2 := by linarith
    rw [div_eq_iff hc.ne'] at this; linarith
  have hd : HasDerivAt (fun x => 2 * x / c - 1) (2 / c) b := by
    have := ((hasDerivAt_id b).const_mul 2).div_const c |>.sub_const 1
    simpa using this
  have := (Real.hasDerivAt_arcsin h1 h2).comp b hd
  refine this.congr_deriv ?_
  have e : 1 - (2 * b / c - 1) ^ 2 = (2 / c) ^ 2 * (b * (c - b)) := by field_simp; ring
  rw [e, Real.sqrt_mul (by positivity), Real.sqrt_sq (by positivity)]
  have : 0 < Real.sqrt (b * (c - b)) := Real.sqrt_pos.mpr (mul_pos hb0 (by linarith))
  field_simp

lemma arc_int {c : ℝ} (hc : 0 < c) :
    IntervalIntegrable (fun b : ℝ => 1 / Real.sqrt (b * (c - b))) volume 0 c ∧
      ∫ b in (0 : ℝ)..c, 1 / Real.sqrt (b * (c - b)) = π := by
  have hcont : ContinuousOn (fun x => Real.arcsin (2 * x / c - 1)) (Icc 0 c) := by
    fun_prop
  have hder : ∀ x ∈ Ioo 0 c, HasDerivAt (fun x => Real.arcsin (2 * x / c - 1))
      (1 / Real.sqrt (x * (c - x))) x := fun x hx => arc_deriv hc hx.1 hx.2
  have hint : IntervalIntegrable (fun b : ℝ => 1 / Real.sqrt (b * (c - b))) volume 0 c := by
    rw [intervalIntegrable_iff_integrableOn_Ioc_of_le hc.le]
    exact intervalIntegral.integrableOn_deriv_of_nonneg hcont hder (fun x _ => by positivity)
  refine ⟨hint, ?_⟩
  rw [intervalIntegral.integral_eq_sub_of_hasDerivAt_of_le hc.le hcont hder hint]
  simp only [mul_zero, zero_div, zero_sub]
  rw [show 2 * c / c - 1 = 1 by field_simp; ring, Real.arcsin_one, Real.arcsin_neg_one]
  ring

/-- the Dirichlet integrand at `τ = 1/2` -/
noncomputable def fD (t : ℝ × ℝ) : ℝ := (t.1 * t.2 * (1 - t.1 - t.2)) ^ (-(1 / 2 : ℝ))

lemma face_half {t : ℝ × ℝ} (ht : t ∈ Y05.Δ) : Y05.faceIntegrand (1 - 1 / 2) t = fD t := by
  obtain ⟨ha, hb, hab⟩ := ht
  have hc : 0 ≤ 1 - t.1 - t.2 := by linarith
  unfold Y05.faceIntegrand fD
  rw [show (1 : ℝ) - 1 / 2 = 1 / 2 by norm_num, show (2 : ℝ) * (1 / 2) = 1 by norm_num,
    Real.rpow_one, Real.rpow_one, Real.rpow_one,
    show t.1 + t.2 + (1 - t.1 - t.2) = 1 by ring, Real.sqrt_one, mul_one]

lemma fD_int : IntegrableOn fD Y05.Δ :=
  (Y05.proof_1 2 two_pos).congr_fun (fun _ ht => face_half ht) Y05.Δ_meas

lemma rpow_half {x : ℝ} (hx : 0 ≤ x) : x ^ (-(1 / 2 : ℝ)) = 1 / Real.sqrt x := by
  rw [Real.rpow_neg hx, Real.sqrt_eq_rpow, inv_eq_one_div]

lemma slice_zero {a : ℝ} (ha : a ≤ 0 ∨ 1 ≤ a) (b : ℝ) : Y05.Δ.indicator fD (a, b) = 0 := by
  by_cases h : (a, b) ∈ Y05.Δ
  · rw [indicator_of_mem h]
    obtain ⟨h1, h2, h3⟩ := h
    simp only at h1 h2 h3
    have : a * b * (1 - a - b) = 0 := by
      rcases ha with ha | ha
      · rw [show a = 0 by linarith]; ring
      · rw [show b = 0 by linarith]; ring
    unfold fD; simp only
    rw [this, Real.zero_rpow (by norm_num)]
  · exact indicator_of_notMem h _

lemma slice_val {a : ℝ} (ha0 : 0 < a) (ha1 : a < 1) :
    ∫ b, Y05.Δ.indicator fD (a, b) = π * a ^ (-(1 / 2 : ℝ)) := by
  have heq : (fun b => Y05.Δ.indicator fD (a, b)) = (Icc 0 (1 - a)).indicator
      (fun b => a ^ (-(1 / 2 : ℝ)) * (1 / Real.sqrt (b * ((1 - a) - b)))) := by
    funext b
    by_cases hb : b ∈ Icc 0 (1 - a)
    · have hm : (a, b) ∈ Y05.Δ := ⟨ha0.le, hb.1, by linarith [hb.2]⟩
      rw [indicator_of_mem hm, indicator_of_mem hb]
      unfold fD; simp only
      have hc : 0 ≤ b * (1 - a - b) := mul_nonneg hb.1 (by linarith [hb.2])
      rw [show a * b * (1 - a - b) = a * (b * (1 - a - b)) by ring,
        Real.mul_rpow ha0.le hc, rpow_half hc]
    · rw [indicator_of_notMem hb]
      apply indicator_of_notMem
      rintro ⟨h1, h2, h3⟩
      exact hb ⟨h2, by simp only at h3; linarith⟩
  rw [heq, integral_indicator measurableSet_Icc, integral_Icc_eq_integral_Ioc,
    ← intervalIntegral.integral_of_le (by linarith), intervalIntegral.integral_const_mul,
    (arc_int (by linarith : (0:ℝ) < 1 - a)).2]
  ring

theorem dirichlet_half : IntegrableOn fD Y05.Δ ∧ ∫ t in Y05.Δ, fD t = 2 * π := by
  refine ⟨fD_int, ?_⟩
  have hI : Integrable (Y05.Δ.indicator fD) (volume.prod volume) := by
    rw [← Measure.volume_eq_prod]
    exact (integrable_indicator_iff Y05.Δ_meas).mpr fD_int
  rw [← integral_indicator Y05.Δ_meas, Measure.volume_eq_prod, integral_prod _ hI]
  have hg : ∀ a : ℝ, ∫ b, Y05.Δ.indicator fD (a, b) =
      (Ioo 0 1).indicator (fun a => π * a ^ (-(1 / 2 : ℝ))) a := by
    intro a
    by_cases ha : a ∈ Ioo (0 : ℝ) 1
    · rw [indicator_of_mem ha, slice_val ha.1 ha.2]
    · rw [indicator_of_notMem ha]
      have : a ≤ 0 ∨ 1 ≤ a := by
        rw [mem_Ioo, not_and_or, not_lt, not_lt] at ha; exact ha
      simp only [slice_zero this, integral_zero]
  simp_rw [hg]
  rw [integral_indicator measurableSet_Ioo, ← integral_Ioc_eq_integral_Ioo,
    ← intervalIntegral.integral_of_le zero_le_one, intervalIntegral.integral_const_mul,
    integral_rpow (Or.inl (by norm_num))]
  norm_num
  ring

end V04

lemma V04.S_one : Y05.S 1 = 4 * Real.sqrt 3 := by
  rw [Y05.proof_3 1 one_pos]
  have hB : Y05.B (1 - 1 / 1) = 1 := by
    have h5 := Y05.proof_5 0 (by norm_num)
    have : (fun g => Real.sqrt (Y05.Ymean (1 - 1 / 1) g)) = Y05.Ymean 0 := by
      funext g; simp [Y05.Ymean]; norm_num
    rw [Y05.B, this, show (1:ℝ) - 1/1 = 0 by norm_num, h5]; simp
  rw [hB]
  have g3 : Real.Gamma 3 = 2 := by
    rw [show (3 : ℝ) = 2 + 1 by norm_num, Real.Gamma_add_one (by norm_num), Real.Gamma_two]; norm_num
  norm_num [g3]
  ring

theorem V04.proof_1 : V04.stmt_1 := by
  refine ⟨Y05.proof_1 1 one_pos, ?_⟩
  have h := V04.S_one
  unfold Y05.S at h
  norm_num at h ⊢
  linarith

theorem V04.proof_2 : V04.stmt_2 := ⟨Y05.proof_1 1 one_pos, V04.S_one⟩

theorem V04.proof_3 : V04.stmt_3 := V04.dirichlet_half

theorem V04.proof_4 : V04.stmt_4 := by
  refine ⟨Y05.proof_1 2 two_pos, ?_⟩
  unfold Y05.S
  rw [setIntegral_congr_fun Y05.Δ_meas (fun t ht => V04.face_half ht), V04.dirichlet_half.2]
  ring

theorem V04.proof_5 : V04.stmt_5 := V04.arc_int one_pos

theorem V04.proof_6 : V04.stmt_6 := by
  have hEq : EqOn (fun t : ℝ => 1 / Real.sqrt (t * (1 - t))) (W02.pIntegrand 2) (uIoo 0 1) := by
    intro t ht
    rw [uIoo_of_le zero_le_one] at ht
    obtain ⟨h0, h1⟩ := ht
    have h1' : 0 < 1 - t := by linarith
    simp only [W02.pIntegrand]
    rw [show (2 : ℝ) / 2 - 2 = -1 by norm_num, Real.rpow_neg_one, Real.rpow_neg_one,
      show t⁻¹ + (1 - t)⁻¹ = (t * (1 - t))⁻¹ by field_simp; ring, Real.sqrt_inv, one_div]
  obtain ⟨hi, hv⟩ := V04.arc_int one_pos
  refine ⟨hi.congr_uIoo hEq, ?_⟩
  unfold W02.Pint
  rw [← intervalIntegral.integral_congr_uIoo hEq, hv]
  ring

theorem V04.proof_7 : V04.stmt_7 := by
  have hc : ContinuousAt W02.FC 0 :=
    (W02.proof_5.2.differentiableAt (by
      have : IsOpen {ε : ℂ | ε.re < 1} := isOpen_lt Complex.continuous_re continuous_const
      exact this.mem_nhds (by simp))).continuousAt
  have h0 : Tendsto (fun p : ℝ => 1 / p) atTop (𝓝 0) := by
    simpa using tendsto_inv_atTop_zero
  have h0c : Tendsto (fun p : ℝ => ((1 / p : ℝ) : ℂ)) atTop (𝓝 0) := by
    have := (Complex.continuous_ofReal.tendsto 0).comp h0
    rwa [Complex.ofReal_zero] at this
  have hF : Tendsto (fun p : ℝ => (W02.FC ((1 / p : ℝ) : ℂ)).re) atTop (𝓝 (W02.FC 0).re) :=
    (Complex.continuous_re.tendsto _).comp (hc.tendsto.comp h0c)
  have hlim : Tendsto (fun p : ℝ => 4 * (2 - (1 / p) * (W02.FC ((1 / p : ℝ) : ℂ)).re)) atTop
      (𝓝 (4 * (2 - 0 * (W02.FC 0).re))) :=
    (tendsto_const_nhds.sub (h0.mul hF)).const_mul 4
  rw [show (4 : ℝ) * (2 - 0 * (W02.FC 0).re) = 8 by ring] at hlim
  refine hlim.congr' ?_
  filter_upwards [eventually_gt_atTop 1] with p hp
  obtain ⟨-, -, h2⟩ := W02.proof_2 p hp
  have h4 := (W02.proof_4 (1 / p) (by rw [div_lt_one (by linarith)]; exact hp)).2
  rw [h4, Complex.ofReal_re]
  linarith
